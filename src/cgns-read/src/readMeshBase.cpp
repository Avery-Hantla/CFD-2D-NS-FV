#include "../include/readMeshBase.hpp"

ReadMeshBase::ReadMeshBase(std::filesystem::path cgnsFilePath) : _cgnsFilePath(cgnsFilePath) {
  if (cg_open(_cgnsFilePath.string().c_str(), CG_MODE_READ, &_fileIndex)) cg_error_exit();
  readBases();
  createBoundaryNameLookupTable();
}

ReadMeshBase::~ReadMeshBase() {
  cg_close(_fileIndex);
}

void ReadMeshBase::readBases() {
  // check how many bases there are. Only a single base per mesh is supported, which should always be true. 
  int numberOfBases;
  if (cg_nbases(_fileIndex, &numberOfBases)) cg_error_exit();
  assert(numberOfBases == 1 && "Mesh reading only supported for a mesh with a single base");
}

void ReadMeshBase::readZones(ZoneType_t zoneTypeRequired) {
  char zonename[64];
  cgsize_t sizes[3][3];

  // check how many zones there are in the CGNS file
  if (cg_nzones(_fileIndex, 1, &_numberOfZones)) cg_error_exit();
  _zoneSizesMax.resize(_numberOfZones);
  _zoneSizesMin.resize(_numberOfZones);

  for (int zone = 0; zone < _numberOfZones; ++zone) {
    // check that the zone within the mesh file contains an structured/unstructured mesh, depending on function argument
    ZoneType_t zoneType;
    if (cg_zone_type(_fileIndex, 1, zone + 1, &zoneType)) cg_error_exit();
    assert(zoneType == zoneTypeRequired && "Wrong mesh type read, expected either structured or unstrucutred");

    // ensure that we only have 1 grid per zone. This should always be the case but we check to avoid nasty surprises
    int numberOfGrids;
    if (cg_ngrids(_fileIndex, 1, zone + 1, &numberOfGrids)) cg_error_exit();
    assert(numberOfGrids == 1 && "Mesh reading only supported for a mesh with a single grid per zone");

    // get the name of the zone and store it in a temporary map. We'll need that for reading boundary conditions later
    if (cg_zone_read(_fileIndex, 1, zone + 1, zonename, sizes[0])) cg_error_exit();
    _zoneNamesMap.emplace(std::string(zonename), zone);

    // we only support 2D meshes for now, make sure this is the case
    int numberOfCoordinatesToRead;
    if (cg_ncoords(_fileIndex, 1, zone + 1, &numberOfCoordinatesToRead)) cg_error_exit();
    assert(numberOfCoordinatesToRead == 2 && "Mesh reading only supported for 2D meshes");

    // store the number of vertices (start and end) for each coordiante (e.g. x and y). Start is always 1
    _zoneSizesMin[zone].resize(numberOfCoordinatesToRead);
    _zoneSizesMax[zone].resize(numberOfCoordinatesToRead);
    for (int coordinate = 0; coordinate < numberOfCoordinatesToRead; ++coordinate) {
      _zoneSizesMin[zone][coordinate] = 1;
      _zoneSizesMax[zone][coordinate] = sizes[0][coordinate];
    }
  }
}

void ReadMeshBase::createBoundaryNameLookupTable() {
  // In CGNS, we can read boundary conditions in two ways; either by placing them directly onto the grid, or by
  // inserting an additional layer (called a family) which connects the grid with a boundary condition. The
  // advantage of using a family node is that the boundary condition remains intact, even if the mesh is changed.
  // The disadvantage is that we have to do a bit more work to read the boundary condition information. First, we
  // need to read the family information and then get the boundary condition name and type from the family. Later,
  // we need to to map that information to the correct boundary condition. In the following, we check if the
  // boundary conditions are stored under a family (i.e. the bc type is FamilySpecified) and then we read the family
  // node with all of its boundary information, which we'll store in the boundaryConditionMap for later lookup.
  int numberOfFamilies = 0;
  if (cg_nfamilies(_fileIndex, 1, &numberOfFamilies)) cg_error_exit();

  for (int family = 0; family < numberOfFamilies; ++family) {
    char familyName[128];
    int isBoundaryCondition, nGeo;
    if (cg_family_read(_fileIndex, 1, family + 1, familyName, &isBoundaryCondition, &nGeo)) cg_error_exit();

    char familyBoundaryName[128];
    BCType_t familyBoundaryType;

    // only if isBoundaryCondition is equal to 1 are boundary conditions stored. For a value other than 1, other
    // family information may be stored e.g. different volume conditions (useful for separating zones/volumes to
    // apply MRF or sliding grids, for example). We are only interested in reading boundary conditions, though.
    if (isBoundaryCondition == 1) {
      if (cg_fambc_read(_fileIndex, 1, family + 1, 1, familyBoundaryName, &familyBoundaryType)) cg_error_exit();
      _boundaryConditionMap.emplace(familyName, familyBoundaryType);
    }
  }
}

unsigned ReadMeshBase::getBoundaryConditionID(BCType_t bcType) {
  if (bcType == 20 || bcType == 21 || bcType == 22 || bcType == 23 || bcType == 24)
    return BC::WALL;
  else if (bcType == 9 || bcType == 10 || bcType == 11 || bcType == 18)
    return BC::INLET;
  else if (bcType == 13 || bcType == 14 || bcType == 15 || bcType == 19)
    return BC::OUTLET;
  else if (bcType == 16 || bcType == 17)
    return BC::SYMMETRY;
  else {
    throw std::runtime_error("Unsupported boundary condition type specified");
    exit(-1);
  }
}
