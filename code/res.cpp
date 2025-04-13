/////////////////////////////////////////////////////////
//           Function To Compute Residual
/////////////////////////////////////////////////////////
#include <vector> 

#include "class_q.hpp"
#include "class_mesh.hpp"
#include "class_f.hpp"
#include "class_residual.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"
#include "class_flow.hpp"
#include "struct_BC.hpp"

#include "reconstruction.cpp"
#include "rusanov.cpp"
#include "roe_flux.cpp"
#include "squeeze.cpp"

void res(class_residual* residual, class_mesh* mesh, struct_size* size, struct_inputs* inputs, class_Q* Qbar, class_Q* Qface, class_Q* Qface_c1, class_Q* Qface_c2, class_flow* freestream, struct_BC* BC) {
    class_F F_reimann;
    F_reimann.init(size->num_faces, freestream->gamma);

    // Reconstuct Cell Faces
    switch (inputs->order) {
        case 1:
            reconstruction(mesh, Qbar, Qface, Qface_c1, Qface_c2, freestream, size, inputs, BC);
            break;
        case 2:
            compute_gradient(Qbar, mesh, freestream, size, inputs);
            reconstruction(mesh, Qbar, Qface, Qface_c1, Qface_c2, freestream, size, inputs, BC);

            switch (inputs->limiter) {
                case 1: // Minmod

                    break;
                case 2: // Squeeze
                    squeeze(mesh, Qbar, Qface_c1, Qface_c2, size);
                    break;
            }

            switch (inputs->eqn) {
                case 2: // NS
                F_reimann.update_vis(mesh, Qface);
            }

            break;
        default: 
            std::cout << "ERROR: ORDER NOT IMPLEMENTED";
            break;
    }

    // Find Reimann Flux
    switch (inputs->flux_solver) {
        case 1: // Rusanov Flux
            rusanov(&F_reimann, Qface_c1, Qface_c2, mesh, freestream, size, BC); 
            break;
        case 2: // Roe Flux
            break;
    }

    // Compute Residual
    int c1_num, c2_num;
    residual->reset();
    for (int idx = 0; idx < size->num_faces; idx++) {
        double sum_temp1 = 0, sum_temp2 = 0, sum_temp3 = 0, sum_temp4 = 0;
        double sum_temp_vis1 = 0, sum_temp_vis2 = 0, sum_temp_vis3 = 0, sum_temp_vis4 = 0;

        c1_num = mesh->face_cell1[idx];
        c2_num = mesh->face_cell2[idx];

        sum_temp1 = (F_reimann.p1[idx]*mesh->face_area[idx]);
        sum_temp2 = (F_reimann.p2[idx]*mesh->face_area[idx]);
        sum_temp3 = (F_reimann.p3[idx]*mesh->face_area[idx]);
        sum_temp4 = (F_reimann.p4[idx]*mesh->face_area[idx]);

        sum_temp_vis1 = (F_reimann.v1[idx]*mesh->face_area[idx]);
        sum_temp_vis2 = (F_reimann.v2[idx]*mesh->face_area[idx]);
        sum_temp_vis3 = (F_reimann.v3[idx]*mesh->face_area[idx]);
        sum_temp_vis4 = (F_reimann.v4[idx]*mesh->face_area[idx]);

        residual->p1[c1_num] += sum_temp_vis1 - sum_temp1;
        residual->p2[c1_num] += sum_temp_vis2 - sum_temp2;
        residual->p3[c1_num] += sum_temp_vis3 - sum_temp3;
        residual->p4[c1_num] += sum_temp_vis4 - sum_temp4;

        if (c2_num >= 0) {
            residual->p1[c2_num] -= sum_temp_vis1 - sum_temp1;
            residual->p2[c2_num] -= sum_temp_vis2 - sum_temp2;
            residual->p3[c2_num] -= sum_temp_vis3 - sum_temp3;
            residual->p4[c2_num] -= sum_temp_vis4 - sum_temp4;
        }
    }
    
    for (int idx = 0; idx < size->num_cells; idx ++) {
        residual->p1[idx] = residual->p1[idx] / mesh->cell_vol[idx];
        residual->p2[idx] = residual->p2[idx] / mesh->cell_vol[idx];
        residual->p3[idx] = residual->p3[idx] / mesh->cell_vol[idx];
        residual->p4[idx] = residual->p4[idx] / mesh->cell_vol[idx];
    }
}