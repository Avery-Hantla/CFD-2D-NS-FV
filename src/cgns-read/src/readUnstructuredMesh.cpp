#include "../include/readUnstructuredMesh.hpp"

ReadUnstructuredMesh::ReadUnstructuredMesh(std::filesystem::path cgnsFilePath) : ReadMeshBase(cgnsFilePath) { }

void ReadUnstructuredMesh::readMesh() {
  readZones(CGNS_ENUMV(Unstructured));
  readCoorinates();
  readCellConnectivity();
  readInterfaceConnectivity();
  readBoundaries();
}

void ReadUnstructuredMesh::readCoorinates() {
  DataType_t type; char coordname[128];
  int numberOfDimensions = _zoneSizesMax[0].size();

  // resize the coordinates array, we store one per zone for both x and y
  _coordinates.resize(_numberOfZones);
  
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    // set the number of max vertices for x and y. Coordinate reading will require this data structure later.
    cgsize_t maxVertices = _zoneSizesMax[zone][0];
    cgsize_t minVertices = _zoneSizesMin[zone][0];

    // based on the max vertices, resize the coordinates array.
    _coordinates[zone].resize(maxVertices);
    for (int i = 0; i < maxVertices; ++i) {
      _coordinates[zone][i].resize(numberOfDimensions);
    }

    // loop over each direction, e.g. x and y (numberOfDimensions should be 2)
    for (int dimension = 0; dimension < numberOfDimensions; ++dimension) {
      // get the coordinates information, mainly the type and name. The type is either single or double preciosion and
      // the coordname is a SIDS-compliant identification, e.g. 'CoordinateX' and 'CoordinateY'      
      if (cg_coord_info(_fileIndex, 1, zone + 1, dimension + 1, &type, coordname)) cg_error_exit();

      // temporary (1D) coordinate array to store coordinates in
      std::vector<double> temp(maxVertices);      
      if (cg_coord_read(_fileIndex, 1, zone + 1, coordname, type, &minVertices, &maxVertices, &temp[0]))
        cg_error_exit();

      // store temporary coordinates in final coordinates array
      for (int vertex = 0; vertex < maxVertices; ++vertex) {
        _coordinates[zone][vertex][dimension] = temp[vertex];
      }
    }
  }
}

void ReadUnstructuredMesh::readInterfaceConnectivity() {
  char interfaceName[128], donorName[128]; GridLocation_t location; GridConnectivityType_t connectType;
  PointSetType_t pointSet, donorPointSet; cgsize_t numberOfPoints, numberOfDonorCells;
  ZoneType_t donorZoneType; DataType_t donorDataType;

  // create map that stores, for each interface (identified by name), the zones that it connects
  std::unordered_map<std::string, std::vector<std::array<int, 2>>> interfaceNameToZone;
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    int numberOfInterfaces;
    if (cg_nconns(_fileIndex, 1, zone + 1, &numberOfInterfaces)) cg_error_exit();
    for (int interface = 0; interface < numberOfInterfaces; ++interface) {
      if (cg_conn_info(_fileIndex, 1, zone + 1, interface + 1, interfaceName, &location, &connectType, &pointSet,
        &numberOfPoints, donorName, &donorZoneType, &donorPointSet, &donorDataType, &numberOfDonorCells))
        cg_error_exit();
      
      if (interfaceNameToZone.find(interfaceName) == interfaceNameToZone.end())
        interfaceNameToZone.emplace(interfaceName, std::vector<std::array<int, 2>>{{zone, interface}});
      else
        interfaceNameToZone[interfaceName].push_back({zone, interface});
    }
  }

  // read interface connectivity for each zone
  _interfaceConnectivity.resize(_numberOfZones);
  for (int zone = 0; zone < _numberOfZones; ++zone) {    
    int numberOfInterfaces;
    if (cg_nconns(_fileIndex, 1, zone + 1, &numberOfInterfaces)) cg_error_exit();
    assert(numberOfInterfaces == _interfaceElementConnectivity[zone].size() && "The number of interfaces should be the"
      " same as the number of interfaces for which elements have been read");

    // for each zone, go through all interfaces stored in this zone and store the connectivity information
    _interfaceConnectivity[zone].resize(numberOfInterfaces);
    for (int interface = 0; interface < numberOfInterfaces; ++interface) {
      if (cg_conn_info(_fileIndex, 1, zone + 1, interface + 1, interfaceName, &location, &connectType, &pointSet,
        &numberOfPoints, donorName, &donorZoneType, &donorPointSet, &donorDataType, &numberOfDonorCells))
        cg_error_exit();
        
      assert(donorZoneType == CGNS_ENUMV(Unstructured) && "Expected unstructured interface");

      auto zones = interfaceNameToZone[interfaceName];
      assert(zones.size() == 2 && "Expected 2 zones per interface");

      const auto &zoneOwner = zones[0][0];
      const auto &interfaceOwner = zones[0][1];
      const auto &zoneNeighbour = zones[1][0];
      const auto &interfaceNeighbour = zones[1][1];

      _interfaceConnectivity[zone][interface].zones = {zoneOwner, zoneNeighbour};
      _interfaceConnectivity[zone][interface].ownerIndex = _interfaceElementConnectivity[zoneOwner][interfaceOwner];
      _interfaceConnectivity[zone][interface].neighbourIndex =
        _interfaceElementConnectivity[zoneNeighbour][interfaceNeighbour];
    }
  }
}

void ReadUnstructuredMesh::readBoundaries() {
  // go through all zones and check which boundary conditions are stored under each zone
  _boundaryConditions.resize(_numberOfZones);
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    int numberOfBoundaries;
    if (cg_nbocos(_fileIndex, 1, zone + 1, &numberOfBoundaries)) cg_error_exit();
    _boundaryConditions[zone].resize(numberOfBoundaries);

    // for each zone, go through all boundaries stored in this zone and store the boundary information
    for (int boundary = 0; boundary < numberOfBoundaries; ++boundary) {
      char boundaryName[128]; int normalIndex, numberOfDataset; PointSetType_t pointSetType;
      cgsize_t numberOfPointsInBoundary, normalListSize; DataType_t normalDataType; BCType_t boundaryType;

      if (cg_boco_info(_fileIndex, 1, zone + 1, boundary + 1, boundaryName, &boundaryType, &pointSetType,
        &numberOfPointsInBoundary, &normalIndex, &normalListSize, &normalDataType, &numberOfDataset)) cg_error_exit();

      assert(pointSetType == CGNS_ENUMV(PointRange) &&
        "Expected to read a start and end location for boundary conditions only");

      _boundaryConditions[zone][boundary].indices = _boundaryElementConnectivity[zone][boundary];

      // Writing the boundary name and type to the boundary structure. We need to check if we need to get this
      // information from a family node or if we can read this directly from the boundary node.
      if (boundaryType == CGNS_ENUMV(FamilySpecified)) {
        // go to the current boundary node within the ZoneBC_t node, which must be present for each zone
        if (cg_goto(_fileIndex, 1, "Zone_t", zone + 1, "ZoneBC_t", 1, "BC_t", boundary + 1, "end")) cg_error_exit();
      
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

void ReadUnstructuredMesh::readCellConnectivity() {
  _internalElementConnectivity.resize(_numberOfZones);
  _interfaceElementConnectivity.resize(_numberOfZones);
  _boundaryElementConnectivity.resize(_numberOfZones);

  // the element connectivity contains the internal, interface, and boundary cells. First, we need to find out what the
  // starting location is for each cell type so that we no which part of the connectivy array to start reading from to
  // get internal, interface, or boudnary cell connectivity.
  auto [internalElementStart, interfaceElementStart, boundaryElementStart] = getStartLocationForElementConnectivity();

  // go through all zones and read the connectivity
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    int numberOfSections;
    if (cg_nsections(_fileIndex, 1, zone + 1, &numberOfSections)) cg_error_exit();

    _interfaceElementConnectivity[zone].resize(interfaceElementStart[zone].size());
    _boundaryElementConnectivity[zone].resize(boundaryElementStart[zone].size());
    
    unsigned internalElementIndex = 0;
    unsigned interfaceElementIndex = 0;
    unsigned boundaryElementIndex = 0;

    // for each zone, go through all sections and read the element connectivity information
    for (int section = 0; section < numberOfSections; ++section) {
      char sectionName[128]; ElementType_t elementType; cgsize_t start, end, parentData, elementDataSize;
      int boundaryFlag, parentFlag, numberOfVerticesPerElement;

      // this routine provides us with the variable 'start', which tells us where to start reading element connectivity
      // from. We use this value later together with the above received start locations for internal, interface, and
      // boundary start locations to figure out what element connectivity we have just read
      if (cg_section_read(_fileIndex, 1, zone + 1, section + 1, sectionName, &elementType, &start, &end, &boundaryFlag,
        &parentFlag)) cg_error_exit();

      // the numbe rof elements in the element connectivity array
      if (cg_ElementDataSize(_fileIndex, 1, zone + 1, section + 1, &elementDataSize)) cg_error_exit();

      // this temporary array will hold all element connectivity information in this section
      std::vector<cgsize_t> elements1D(elementDataSize);
      if (cg_elements_read(_fileIndex, 1, zone + 1, section + 1, &elements1D[0], &parentData)) cg_error_exit();
      if (cg_npe(elementType, &numberOfVerticesPerElement)) cg_error_exit();

      assert(elementDataSize % numberOfVerticesPerElement == 0 &&
        "The number of elements in the section is not a multiple of the number of vertices per element");
      
      unsigned numberOfElements = elementDataSize / numberOfVerticesPerElement;

      // the element connectivity is received from CGNS in a 1D array. It has a total size of
      // size = [numberOfElements * numberOfVerticesPerElement]. What we want, for easier indexing, is a 2D array where
      // each dimension will have a size of [numberOfElements][numberOfVerticesPerElement]. This lambda expression does
      // that for us. We use a lambda here because we want to use this routine 3 time (interal, interface, and boundary)
      // cells. 
      auto map1DElementsto2DVector = [numberOfVerticesPerElement, numberOfElements, &elements1D]
        (std::vector<ReadUnstructuredMesh::IndexType> &target) {
        unsigned counter = 0;
        unsigned elementStartLocation = target.size();
        target.resize(elementStartLocation + numberOfElements);
        for (int element = 0; element < numberOfElements; ++element) {
          unsigned elementIndex = elementStartLocation + element;
          target[elementIndex].resize(numberOfVerticesPerElement);
          for (int vertex = 0; vertex < numberOfVerticesPerElement; ++vertex) {
            target[elementIndex][vertex] = elements1D[counter++] - 1;
          }
        }
      };

      auto &internal = internalElementStart[zone];
      auto &interface = interfaceElementStart[zone];
      auto &boundary = boundaryElementStart[zone];

      // now we make use of the variable 'start' and the starting location for the internal, interface, and boundary
      // cells. We go through all of them and check if they contain the value for 'start'. If so, then we know that
      // we have read the connectivity for the internal, interface, or boundary cells. We then map the 1D connectivity
      // array to the 2D element connectivity array using our above defined lambda expression.
      if (std::find(internal.begin(), internal.end(), start) != internal.end()) {
        map1DElementsto2DVector(_internalElementConnectivity[zone]);
        internalElementIndex++;
      } else if (std::find(interface.begin(), interface.end(), start) != interface.end()) {
        map1DElementsto2DVector(_interfaceElementConnectivity[zone][interfaceElementIndex]);
        interfaceElementIndex++;
      } else if (std::find(boundary.begin(), boundary.end(), start) != boundary.end()) {
        map1DElementsto2DVector(_boundaryElementConnectivity[zone][boundaryElementIndex]);
        boundaryElementIndex++;
      } else {
        throw std::runtime_error("Could not itendify which element connectivity was read");
      }
    }
    int totalSectionsProcessed = internalElementIndex + interfaceElementIndex + boundaryElementIndex;
    assert(numberOfSections == totalSectionsProcessed && "Did not process all sections, some information is lost");
  }
}

std::tuple<std::vector<std::vector<cgsize_t>>, std::vector<std::vector<cgsize_t>>, std::vector<std::vector<cgsize_t>>>
  ReadUnstructuredMesh::getStartLocationForElementConnectivity() const {

  std::vector<std::vector<cgsize_t>> internalElementStart(_numberOfZones);
  std::vector<std::vector<cgsize_t>> interfaceElementStart(_numberOfZones);
  std::vector<std::vector<cgsize_t>> boundaryElementStart(_numberOfZones);

  // for each zone, read the starting location of the element connectivity for internal, interface, and boundary cells
  unsigned totalNumberOfInterfaces = 0;
  for (int zone = 0; zone < _numberOfZones; ++zone) {
    // read interface connectivity 
    int numberOfInterfaces = 0;
    if (cg_nconns(_fileIndex, 1, zone + 1, &numberOfInterfaces)) cg_error_exit();
    for (int interface = 0; interface < numberOfInterfaces; ++interface) {
      char interfaceName[128], donorName[128]; GridLocation_t location; GridConnectivityType_t connectType;
      cgsize_t donorData, numberOfPoints, numberOfDonorCells, start; DataType_t donorDataType; ZoneType_t donorZoneType;
      PointSetType_t pointSet, donorPointSet; 

      if (cg_conn_info(_fileIndex, 1, zone + 1, interface + 1, interfaceName, &location, &connectType, &pointSet,
        &numberOfPoints, donorName, &donorZoneType, &donorPointSet, &donorDataType, &numberOfDonorCells))
        cg_error_exit();

      if (cg_conn_read(_fileIndex, 1, zone + 1, interface + 1, &start, donorDataType, &donorData)) cg_error_exit();
      interfaceElementStart[zone].push_back(start);
      totalNumberOfInterfaces++;
    }

    // read boundary connectivity 
    int numberOfBoundaries = 0;
    if (cg_nbocos(_fileIndex, 1, zone + 1, &numberOfBoundaries)) cg_error_exit();
    for (int boundary = 0; boundary < numberOfBoundaries; ++boundary) {
      int normalList; cgsize_t start;
      if (cg_boco_read(_fileIndex, 1, zone + 1, boundary + 1, &start, &normalList)) cg_error_exit();
      boundaryElementStart[zone].push_back(start);
    }

    // read all connectivity. we can't directly read only connectivity of internal cells but instead have to read all of
    // them, which will include the interface and boundary connectivity, as they are all stored under the same node (we 
    // just happen to have separate functions to read boundary and interface connectivity but not for internal cells).
    // Once we have all connectivity read, we simple remove the interface and boundary connectivity from it and are left
    // with internal connectivity information only.
    int numberOfSections = 0;
    if (cg_nsections(_fileIndex, 1, zone + 1, &numberOfSections)) cg_error_exit();
    for (int section = 0; section < numberOfSections; ++section) {
      char sectionName[128]; ElementType_t elementType; cgsize_t start, end; int nbndry, parentFlag;
      if (cg_section_read(_fileIndex, 1, zone + 1, section + 1, sectionName, &elementType, &start, &end, &nbndry,
        &parentFlag)) cg_error_exit();
      internalElementStart[zone].push_back(start);
    }
    // erase all boundary connectivity start locations from internal elements
    for (const auto &boundary : boundaryElementStart[zone]) {
      internalElementStart[zone].erase(std::remove(internalElementStart[zone].begin(), internalElementStart[zone].end(),
        boundary), internalElementStart[zone].end());
    }
    // erase all interface connectivity start locations from internal elements
    for (const auto &interface : interfaceElementStart[zone]) {
      internalElementStart[zone].erase(std::remove(internalElementStart[zone].begin(), internalElementStart[zone].end(),
        interface), internalElementStart[zone].end());
    }
  }
  assert(totalNumberOfInterfaces % 2 == 0 && "Expected an even number of interface pairs");

  return {internalElementStart, interfaceElementStart, boundaryElementStart};
}
