#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cgnslib.h>

#include "class_mesh.hpp"
#include "struct_size.hpp"

void write_cg(class_mesh* mesh, struct_size* size) {
    int fileId, baseId = 1, zoneId = 1;
    cgsize_t sizes[2];// one = 1;

    cg_open("test", CG_MODE_WRITE, &fileId);
    cg_base_write(fileId, "mesh", 2, 2, &baseId);
    cg_zone_write(fileId, baseId, "Zone 1", sizes, Unstructured, &zoneId);

    // Make X and Y coords vectors
    std::vector<double> xCoords(sizes[0]), yCoords(sizes[0]);
    for (int idx = 0; idx < size->num_points; idx++) {
        xCoords[idx] = mesh->x[idx];
        yCoords[idx] = mesh->y[idx];
        // output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << std::endl;
    }
    for (int idx = 0; idx < size->num_cells; idx ++) {
        xCoords[size->num_points + idx] = mesh->cell_centerx[idx]; 
        yCoords[size->num_points + idx] = mesh->cell_centery[idx];
        // output << Qbar->rho[idx] << " " << Qbar->u[idx] << " " << Qbar->v[idx] << " " << Qbar->P[idx] << std::endl;
    }

    std::vector<cgsize_t> localElems;
    // Make conectivity matrix
    for (int idx = 0; idx < size->num_cells; idx ++) {
        for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx ++) {
          int face_num = mesh->find_cell_face(idx, jdx);
    
          int point1 = mesh->face_point1[face_num];
          int point2 = mesh->face_point2[face_num];
    
          localElems.push_back(point1+1);
          localElems.push_back(point2+1);
          localElems.push_back(size->num_points + idx + 1);
        }
      }

    int xid;
    cg_coord_write(fileId, baseId, zoneId, RealDouble, "CoordinateX", xCoords.data(), &xid);
    cg_coord_write(fileId, baseId, zoneId, RealDouble, "CoordinateY", yCoords.data(), &xid);

    cgsize_t start = 1, end = (size->num_cells+size->num_faces)/3;

    // end = start + localElems.size()/3 - 1;

    int secId;

    char buff[32];
    // sprintf(buff, "%d", 1);
    std::string secName = std::string("Proc") + buff;

    cg_section_write(fileId, baseId, zoneId, secName.c_str(), TRI_3, start, end, 0, localElems.data(), &secId);

    cg_close(fileId);
}