/////////////////////////////////////////////////////////
//           Function To Compute Residual
/////////////////////////////////////////////////////////
#include <vector> 

#include "class_q.hpp"
#include "class_mesh.hpp"
#include "class_f.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"
#include "class_flow.hpp"
#include "struct_BC.hpp"
#include "struct_residual.hpp"

#include "reconstruction.cpp"
#include "rusanov.cpp"

void res(struct_residual* residual, class_mesh* mesh, struct_size* size, struct_inputs* inputs, class_Q* Qbar, class_Q* Qface_c1, class_Q* Qface_c2, class_flow* freestream, struct_BC* BC) {
    class_F F_reimann;
    F_reimann.init(size->num_cells*6, freestream->gamma);

    // Reconstuct Cell Faces
    reconstruction(mesh, Qbar, Qface_c1, Qface_c2, freestream, size, inputs, BC);
    // Find Reimann Flux
    rusanov(&F_reimann, Qface_c1, Qface_c2, mesh, freestream, size, BC); 

    // Compute Residual
    int face_num;
    for (int idx = 0; idx < size->num_cells; idx++) {
        double sum_temp1 = 0, sum_temp2 = 0, sum_temp3 = 0, sum_temp4 = 0;
        for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx ++) {
            face_num = mesh->find_cell_face(idx,jdx);

            sum_temp1 += mesh->get_cell_face_rx(idx,jdx)*(F_reimann.p1[face_num]*mesh->face_area[face_num]);
            sum_temp2 += mesh->get_cell_face_rx(idx,jdx)*(F_reimann.p2[face_num]*mesh->face_area[face_num]);
            sum_temp3 += mesh->get_cell_face_rx(idx,jdx)*(F_reimann.p3[face_num]*mesh->face_area[face_num]);
            sum_temp4 += mesh->get_cell_face_rx(idx,jdx)*(F_reimann.p4[face_num]*mesh->face_area[face_num]);
        }
        residual->p1[idx] = (-1/mesh->cell_vol[idx]) * sum_temp1;
        residual->p2[idx] = (-1/mesh->cell_vol[idx]) * sum_temp2;
        residual->p3[idx] = (-1/mesh->cell_vol[idx]) * sum_temp3;
        residual->p4[idx] = (-1/mesh->cell_vol[idx]) * sum_temp4;
    }
}