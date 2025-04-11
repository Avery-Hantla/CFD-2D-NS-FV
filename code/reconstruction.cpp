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

#include "cell_gradient.cpp"

void reconstruction(class_mesh* mesh, class_Q* Qbar, class_Q* Qfaces_c1, class_Q* Qfaces_c2, class_flow* freestream, struct_size* size, struct_inputs* inputs, struct_BC* BC) {
    int cell_1, cell_2;
    double Vn;
    class_flow wall_BC;
    wall_BC.gamma = freestream->gamma;

    // Reconstruct Faces
    for (int idx = 0; idx < size->num_faces; idx ++) {
        cell_1 = mesh->face_cell1[idx];
        cell_2 = mesh->face_cell2[idx];

        Qfaces_c1->p1[idx] = Qbar->p1[cell_1] + Qbar->Qxp1[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp1[cell_1]*mesh->face_cell1_dy[idx];
        Qfaces_c1->p2[idx] = Qbar->p2[cell_1] + Qbar->Qxp2[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp2[cell_1]*mesh->face_cell1_dy[idx];
        Qfaces_c1->p3[idx] = Qbar->p3[cell_1] + Qbar->Qxp3[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp3[cell_1]*mesh->face_cell1_dy[idx];
        Qfaces_c1->p4[idx] = Qbar->p4[cell_1] + Qbar->Qxp4[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp4[cell_1]*mesh->face_cell1_dy[idx];

        Qfaces_c1->updateflow_int(idx);
        
        // Apply BC
        switch (cell_2) { 
            case -1: // Freestream BC
                Qfaces_c2->p1[idx] = freestream->Q1;
                Qfaces_c2->p2[idx] = freestream->Q2;
                Qfaces_c2->p3[idx] = freestream->Q3;
                Qfaces_c2->p4[idx] = freestream->Q4;
                break;
            case -2:  // Wall BC
                wall_BC.rho = Qfaces_c1->rho[idx];
                wall_BC.P = Qfaces_c1->P[idx];
                Vn = Qfaces_c1->u[idx]*mesh->face_nx[idx] + Qfaces_c1->v[idx]*mesh->face_ny[idx];
                wall_BC.u = Qfaces_c1->u[idx] - 2.0*Vn*mesh->face_nx[idx];
                wall_BC.v = Qfaces_c1->v[idx] - 2.0*Vn*mesh->face_ny[idx];
                wall_BC.updateQ();

                Qfaces_c2->p1[idx] = wall_BC.Q1;
                Qfaces_c2->p2[idx] = wall_BC.Q2;
                Qfaces_c2->p3[idx] = wall_BC.Q3;
                Qfaces_c2->p4[idx] = wall_BC.Q4;
                break;
            case -3: // Extrapolation BC
                Qfaces_c2->p1[idx] = Qfaces_c1->p1[idx];
                Qfaces_c2->p2[idx] = Qfaces_c1->p2[idx];
                Qfaces_c2->p3[idx] = Qfaces_c1->p3[idx];
                Qfaces_c2->p4[idx] = Qfaces_c1->p4[idx];
                break;
            default:
                Qfaces_c2->p1[idx] = Qbar->p1[cell_2] + Qbar->Qxp1[cell_2]*mesh->face_cell2_dx[idx] + Qbar->Qyp1[cell_2]*mesh->face_cell2_dy[idx];
                Qfaces_c2->p2[idx] = Qbar->p2[cell_2] + Qbar->Qxp2[cell_2]*mesh->face_cell2_dx[idx] + Qbar->Qyp2[cell_2]*mesh->face_cell2_dy[idx];
                Qfaces_c2->p3[idx] = Qbar->p3[cell_2] + Qbar->Qxp3[cell_2]*mesh->face_cell2_dx[idx] + Qbar->Qyp3[cell_2]*mesh->face_cell2_dy[idx];
                Qfaces_c2->p4[idx] = Qbar->p4[cell_2] + Qbar->Qxp4[cell_2]*mesh->face_cell2_dx[idx] + Qbar->Qyp4[cell_2]*mesh->face_cell2_dy[idx];
                break;
        }
    }
    Qfaces_c2->updateflow();

    // Find Viscous Fluxes for NS
    switch (inputs->eqn) {
        case 2: // If eqn == NS

            break;
    }
}