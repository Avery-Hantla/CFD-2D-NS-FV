#include "../include/readStructuredMesh.hpp"

ReadStructuredMesh::ReadStructuredMesh(std::filesystem::path cgnsFilePath) : ReadMeshBase(cgnsFilePath) { }

void ReadStructuredMesh::readMesh() {
  readZones(CGNS_ENUMV(Structured));
  readCoorinates();
  readInterfaceConnectivity();
  readBoundaries();
}

void ReadStructuredMesh::readCoorinates() {
  DataType_t type; char coordname[128];
  int numberOfDimensions = _zoneSizesMax[0].size();

  // resize the coordinates array, we store one per zone for both x and y
  _coordinates.resize(_numberOfZones);
  
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    // set the number of max vertices for x and y. Coordinate reading will require this data structure later.
    cgsize_t maxVertices[2] = { _zoneSizesMax[zone][COORDINATE::X], _zoneSizesMax[zone][COORDINATE::Y] };
    cgsize_t minVertices[2] = { _zoneSizesMin[zone][COORDINATE::X], _zoneSizesMin[zone][COORDINATE::Y] };

    // based on the max vertices, resize the coordinates array for both x and y
    _coordinates[zone].resize(maxVertices[COORDINATE::X]);
    for (int i = 0; i < maxVertices[0]; ++i) {
      _coordinates[zone][i].resize(maxVertices[COORDINATE::Y]);
      for (int j = 0; j < maxVertices[1]; ++j) {
        _coordinates[zone][i][j].resize(numberOfDimensions);
      }
    }

    // loop over each direction, e.g. x and y (numberOfDimensions should be 2)
    for (int dimension = 0; dimension < numberOfDimensions; ++dimension) {
      // get the coordinates information, mainly the type and name. The type is either single or double preciosion and
      // the coordname is a SIDS-compliant identification, e.g. 'CoordinateX' and 'CoordinateY'      
      if (cg_coord_info(_fileIndex, 1, zone + 1, dimension + 1, &type, coordname)) cg_error_exit();

      // temporary (1D) coordinate array to store coordinates in. 2D std::vectors do no guarantee contiguous memory 
      // layout, thus we need to first read in all coordinates into a 1D array and then map that to a 2D array 
      std::vector<double> temp(maxVertices[COORDINATE::X] * maxVertices[COORDINATE::Y]);      
      if (cg_coord_read(_fileIndex, 1, zone + 1, coordname, type, minVertices, maxVertices, &temp[0])) cg_error_exit();

      // now map temporary (1D) coordinate array to 2D coordinate array
      for (int vertex = 0; vertex < maxVertices[COORDINATE::X] * maxVertices[COORDINATE::Y]; ++vertex) {
        unsigned i = vertex % maxVertices[0];
        unsigned j = vertex / maxVertices[1];
        _coordinates[zone][i][j][dimension] = temp[vertex];
      }
    }
  }
}

void ReadStructuredMesh::readInterfaceConnectivity() {
  int numberOfInterfaces;

  // read the number of interfaces in the mesh
  if (cg_n1to1_global(_fileIndex, 1, &numberOfInterfaces)) cg_error_exit();
  
  // allocate temporary memory. This is required to ensure memory layout is contiguous
  char **interfaceName = new char*[numberOfInterfaces];
  char **zoneName = new char*[numberOfInterfaces];
  char **donorName = new char*[numberOfInterfaces];
  cgsize_t **owner = new cgsize_t*[numberOfInterfaces];
  cgsize_t **neighbour = new cgsize_t*[numberOfInterfaces];
  int **transform = new int*[numberOfInterfaces];

  int numberOfDimensions = _zoneSizesMax[0].size();
  for (int interface = 0; interface < numberOfInterfaces; ++interface) {
    interfaceName[interface] = new char[128];
    zoneName[interface] = new char[128];
    donorName[interface] = new char[128];
    owner[interface] = new cgsize_t[2 * numberOfDimensions];
    neighbour[interface] = new cgsize_t[2 * numberOfDimensions];
    transform[interface] = new int[numberOfDimensions];
  }
  
  // read all connectivity (interface) information all at once
  if (cg_1to1_read_global(_fileIndex, 1, interfaceName, zoneName, donorName, owner, neighbour, transform))
    cg_error_exit();

  _interfaceConnectivity.resize(numberOfInterfaces);
  for (int interface = 0; interface < numberOfInterfaces; ++interface) {
    // store zone index of the interface owner and neighbour 
    _interfaceConnectivity[interface].zones[INTERFACE::OWNER] = _zoneNamesMap[std::string(zoneName[0])];
    _interfaceConnectivity[interface].zones[INTERFACE::NEIGHBOUR] = _zoneNamesMap[std::string(donorName[0])];

    // construct (temporary) 2x2 transformation matrix based on information read from the interface. Details about the
    // transformation matrix can be found here: https://cgns.github.io/CGNS_docs_current/sids/cnct.html#Transform
    std::vector<std::vector<int>> tMatrix(numberOfDimensions, std::vector<int>(numberOfDimensions));
    for (int dimension = 0; dimension < numberOfDimensions; ++dimension) {
      auto value = transform[interface][dimension];
      tMatrix[std::abs(value) - 1][dimension] = (value > 0) ? 1 : ((value < 0) ? -1 : 0);
    }
    
    // number of vertices in the interface alogn the x and y direction
    unsigned numberOfIndicesX = std::abs(owner[interface][2] - owner[interface][0]);
    unsigned numberOfIndicesY = std::abs(owner[interface][3] - owner[interface][1]);
    unsigned numberOfIndices  = std::max(numberOfIndicesX + 1, numberOfIndicesY + 1);

    // a map of indices, mapping i,j indices from the owning zone to the i,j indices of the neighbouring zone
    _interfaceConnectivity[interface].ownerIndex.resize(numberOfIndices);
    _interfaceConnectivity[interface].neighbourIndex.resize(numberOfIndices);

    for (unsigned index = 0; index < numberOfIndices; ++index) {
      const unsigned ownerStartX = owner[interface][COORDINATE::X] - 1;
      const unsigned ownerStartY = owner[interface][COORDINATE::Y] - 1;
      const unsigned neighbourStartX = neighbour[interface][COORDINATE::X] - 1;
      const unsigned neighbourStartY = neighbour[interface][COORDINATE::Y] - 1;

      const unsigned indexX = numberOfIndicesX == 0 ? 0 : index;
      const unsigned indexY = numberOfIndicesY == 0 ? 0 : index;

      unsigned ownerIndexX = ownerStartX + indexX;
      unsigned ownerIndexY = ownerStartY + indexY;

      // calculate the i, j indices of the neighbour zone with the transformation matrix and owner indices. See
      // https://cgns.github.io/CGNS_docs_current/sids/cnct.html#Transform for more information about this equation
      unsigned neighbourIndexX = tMatrix[0][0] * indexX + tMatrix[0][1] * indexY + neighbourStartX;
      unsigned neighbourIndexY = tMatrix[1][0] * indexX + tMatrix[1][1] * indexY + neighbourStartY;

      // write connectivty information to interface connectivity structure
      _interfaceConnectivity[interface].ownerIndex[index][COORDINATE::X] = ownerIndexX;
      _interfaceConnectivity[interface].ownerIndex[index][COORDINATE::Y] = ownerIndexY;
      _interfaceConnectivity[interface].neighbourIndex[index][COORDINATE::X] = neighbourIndexX;
      _interfaceConnectivity[interface].neighbourIndex[index][COORDINATE::Y] = neighbourIndexY;
    }
  }

  // deallocate memory to avoid memory leakage
  for (int interface = 0; interface < numberOfInterfaces; ++interface) {
    delete [] interfaceName[interface];
    delete [] zoneName[interface];
    delete [] donorName[interface];
    delete [] owner[interface];
    delete [] neighbour[interface];
    delete [] transform[interface];
  }
  delete [] interfaceName;
  delete [] zoneName;
  delete [] donorName;
  delete [] owner;
  delete [] neighbour;
  delete [] transform;
}

void ReadStructuredMesh::readBoundaries() {
  _boundaryConditions.resize(_numberOfZones);
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    // get the number of boundary conditions for the current zone
    int numberOfBoundaries;
    if (cg_nbocos(_fileIndex, 1, zone + 1, &numberOfBoundaries)) cg_error_exit();
    _boundaryConditions[zone].resize(numberOfBoundaries);

    for (int boundary = 0; boundary < numberOfBoundaries; ++boundary) {
      // temporary boundary condition information;
      char boundaryName[128]; BCType_t boundaryType; PointSetType_t pointSetType; cgsize_t numberOfPointsInBoundary;
      int normalIndex, numberOfDataset, normalList; cgsize_t normalListSize; DataType_t normalDataType;

      // read boundary condition information
      if (cg_boco_info(_fileIndex, 1, zone + 1, boundary + 1, boundaryName, &boundaryType, &pointSetType,
        &numberOfPointsInBoundary, &normalIndex, &normalListSize, &normalDataType, &numberOfDataset)) cg_error_exit();

      assert(pointSetType == PointRange && "Expected to read a start and end location for boundary conditions only");

      // start and end location for boundary conditions in i, j index
      cgsize_t boundaryPoints[2][2]; 

      // read boundary condition start and end location
      if (cg_boco_read(_fileIndex, 1, zone + 1, boundary + 1, boundaryPoints[0], &normalList)) cg_error_exit();

      // number of vertices in the interface along the x and y direction
      unsigned numberOfIndicesX = std::abs(boundaryPoints[1][0] - boundaryPoints[0][0]);
      unsigned numberOfIndicesY = std::abs(boundaryPoints[1][1] - boundaryPoints[0][1]);
      unsigned numberOfIndices  = std::max(numberOfIndicesX + 1, numberOfIndicesY + 1);
      _boundaryConditions[zone][boundary].indices.resize(numberOfIndices);

      // get all indices along current boundary condition and store them in boundary structure. Reconstructed based on
      // the start and end location of the current boundary.
      for (unsigned index = 0; index < numberOfIndices; ++index) {
        const unsigned indexX = numberOfIndicesX == 0 ? 0 : index;
        const unsigned indexY = numberOfIndicesY == 0 ? 0 : index;
        _boundaryConditions[zone][boundary].indices[index][COORDINATE::X] = boundaryPoints[0][COORDINATE::X] + indexX - 1;
        _boundaryConditions[zone][boundary].indices[index][COORDINATE::Y] = boundaryPoints[0][COORDINATE::Y] + indexY - 1;
      }
      
      // Writing the boundary name and type to the boundary structure. We need to check if we need to get this
      // information from a family or if we can read this directly.
      if (boundaryType == CGNS_ENUMV(FamilySpecified)) {
        // go to the current boundary node within the ZoneBC_t node, which must be present for each zone
        if (cg_goto(_fileIndex, 1, "Zone_t", zone + 1, "ZoneBC_t", 1, "BC_t", boundary + 1,"end")) cg_error_exit();
      
        // once we have navigated to the boundary node, read the family name, this will be checked against the family
        // names we read previously. From the family name, we can get the boundary condition type through the map we set
        // up earlier.
        char familyName[128];
        if (cg_famname_read(familyName)) cg_error_exit();

        // store boundary condition name and type in the boundary condition structure
        _boundaryConditions[zone][boundary].boundaryName = familyName;
        _boundaryConditions[zone][boundary].boundaryType = getBoundaryConditionID(_boundaryConditionMap[familyName]);
      } else {
        // store boundary condition name and type in the boundary condition structure
        _boundaryConditions[zone][boundary].boundaryName = boundaryName;
        _boundaryConditions[zone][boundary].boundaryType = getBoundaryConditionID(boundaryType);
      }
    }
  }
}