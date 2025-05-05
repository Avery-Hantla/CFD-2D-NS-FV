/////////////////////////////////////////////////////////
//               Function To Read CGNS File
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include <numeric>
//#include <algorithm>

#include "../support/class_mesh.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_BC.hpp"
#include "../support/struct_inputs.hpp"

#include "../cgns-read/meshReader.hpp"

void readCGNS() {
    // unstructured mesh reading test, with family-based boundary condition
    auto unstructuredMeshPathFamily = std::filesystem::path("../../grids/unstructured2D.cgns");
    ReadUnstructuredMesh unstructuredMeshFamily(unstructuredMeshPathFamily);
    unstructuredMeshFamily.readMesh();

    // unstructured mesh reading test, without family-based boundary condition
    auto unstructuredMeshPathNoFamily = std::filesystem::path("../grids/unstructured2DNoFamily.cgns");
    ReadUnstructuredMesh unstructuredMeshNoFamily(unstructuredMeshPathNoFamily);
    unstructuredMeshNoFamily.readMesh();


    // Check number of zones
    int num_zones = coordinates.size();
    std::cout << "Read " num_zones << " Zones" << std::endl;

    // Read coordinates 
    for (int idx = 0; idx < num_zones; idx ++) {
        int num_verts = coordinates[idx].size();
        for (int jdx = 0; jdx < num_verts; jdx ++) {
            coordinates[idx][jdx][COORDINATE::X];
            coordinates[idx][jdx][COORDINATE::Y];
        }
    }

    assert(coordinates[0].size() == 13 && "Expected 13 vertices in zone 0");
    assert(coordinates[1].size() == 9 && "Expected 9 vertices in zone 1");
}