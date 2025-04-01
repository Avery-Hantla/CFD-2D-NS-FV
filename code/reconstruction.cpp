/////////////////////////////////////////////////////////
//       Function To Reconstruct Solution in Cell
/////////////////////////////////////////////////////////
#include <iostream>

#include "class_mesh.hpp"
#include "class_q.hpp"
#include "class_flow.hpp"
#include "struct_inputs.hpp"
#include "struct_size.hpp"
#include "struct_BC.hpp"

void reconstruction(class_mesh* mesh, class_Q* Qbar, class_Q* Qfaces_c1, class_Q* Qfaces_c2, class_flow* freestream, struct_size* size, struct_inputs* inputs, struct_BC* BC) {
    switch (inputs->order) {
        case 1: // First Order FV Scheme
            int cell_1, cell_2;
            double Vn;
            class_flow wall_BC;
            wall_BC.gamma = freestream->gamma;
            for (int idx = 0; idx < size->num_faces; idx ++) {
                cell_1 = mesh->face_cell1[idx];

                Qfaces_c1->p1[idx] = Qbar->p1[cell_1];
                Qfaces_c1->p2[idx] = Qbar->p2[cell_1];
                Qfaces_c1->p3[idx] = Qbar->p3[cell_1];
                Qfaces_c1->p4[idx] = Qbar->p4[cell_1];

                cell_2 = mesh->face_cell2[idx];
                // Apply BC
                if (cell_2 == -1) { // Freestream BC
                    Qfaces_c2->p1[idx] = freestream->Q1;
                    Qfaces_c2->p2[idx] = freestream->Q2;
                    Qfaces_c2->p3[idx] = freestream->Q3;
                    Qfaces_c2->p4[idx] = freestream->Q4;
                } else if (cell_2 == -2) { // Wall BC
                    wall_BC.rho = Qbar->rho[cell_2];
                    wall_BC.P = Qbar->P[cell_2];
                    Vn = Qbar->u[cell_2]*mesh->face_nx[idx] + Qbar->v[cell_2]*mesh->face_ny[idx];
                    wall_BC.u = Qbar->u[cell_2] - 2*Vn*mesh->face_nx[idx];
                    wall_BC.v = Qbar->v[cell_2] - 2*Vn*mesh->face_ny[idx];
                    wall_BC.updateQ();

                    Qfaces_c2->p1[idx] = wall_BC.Q1;
                    Qfaces_c2->p2[idx] = wall_BC.Q2;
                    Qfaces_c2->p3[idx] = wall_BC.Q3;
                    Qfaces_c2->p4[idx] = wall_BC.Q4;
                } else if (cell_2 == -3) { // Extrapolation BC
                    Qfaces_c2->p1[idx] = Qbar->p1[cell_1];
                    Qfaces_c2->p2[idx] = Qbar->p2[cell_1];
                    Qfaces_c2->p3[idx] = Qbar->p3[cell_1];
                    Qfaces_c2->p4[idx] = Qbar->p4[cell_1];
                } else {
                    Qfaces_c2->p1[idx] = Qbar->p1[cell_2];
                    Qfaces_c2->p2[idx] = Qbar->p2[cell_2];
                    Qfaces_c2->p3[idx] = Qbar->p3[cell_2];
                    Qfaces_c2->p4[idx] = Qbar->p4[cell_2];
                }
            }
            Qfaces_c1->updateflow();
            Qfaces_c2->updateflow();
            break;
        case 2: // Second Order FV Scheme
            
        default:
            std::cout << "ERROR: ORDER NOT IMPLEMENTED";
    };
}