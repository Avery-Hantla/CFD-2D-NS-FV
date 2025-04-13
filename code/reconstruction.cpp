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

void reconstruction(class_mesh* mesh, class_Q* Qbar, class_Q* Qface, class_Q* Qfaces_c1, class_Q* Qfaces_c2, class_flow* freestream, struct_size* size, struct_inputs* inputs, struct_BC* BC) {
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
            case -4: // NO SLIP WALL BC
                wall_BC.rho = Qfaces_c1->rho[idx];
                wall_BC.P = Qfaces_c1->P[idx];
                wall_BC.u = -Qfaces_c1->u[idx];
                wall_BC.v = -Qfaces_c1->v[idx];
                wall_BC.updateQ();

                Qfaces_c2->p1[idx] = wall_BC.Q1;
                Qfaces_c2->p2[idx] = wall_BC.Q2;
                Qfaces_c2->p3[idx] = wall_BC.Q3;
                Qfaces_c2->p4[idx] = wall_BC.Q4;
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

    // Reconstruct Viscous Terms for NS
    switch (inputs->eqn) {
        case 2: // If eqn == NS
            // Average the face values
            for (int idx = 0; idx < size->num_faces; idx++) {
                Qface->p1[idx] = 0.5*(Qfaces_c1->p1[idx] + Qfaces_c2->p1[idx]);
                Qface->p2[idx] = 0.5*(Qfaces_c1->p2[idx] + Qfaces_c2->p2[idx]);
                Qface->p3[idx] = 0.5*(Qfaces_c1->p3[idx] + Qfaces_c2->p3[idx]);
                Qface->p4[idx] = 0.5*(Qfaces_c1->p4[idx] + Qfaces_c2->p4[idx]);
            }
            Qface->updateflow();

            // Find dQdl and dQdm
            for (int idx = 0; idx < size->num_faces; idx ++) {
                int cell1 = mesh->face_cell1[idx];
                int cell2 = mesh->face_cell2[idx];

                switch (cell2) {
                    case -1: // Freestream
                        Qface->comp_vis(idx, 0, 0, 0, 0, 0, 0, 0, 0);
                        break;
                    case -2: // SLIP_WALL
                        std::cout << "ERROR: SLIP_WALL CANNOT BE USED WITH NS \n";
                        break; 
                    case -3: // EXTRAPOLATION
                        Qface->comp_vis(idx, 0, 0, 0, 0, 0, 0, 0, 0);
                        break;
                    case -4: // NO_SLIP_WALL
                        
                        // Asumming Adiabatic
                        Qface->dTdx[idx] = 0;
                        Qface->dTdy[idx] = 0;
                        break;
                    default:
                        // Find dQdl
                        double dQdl_p1 = (Qbar->p1[cell2] - Qbar->p1[cell1]) / std::abs(mesh->face_dist[idx]);
                        double dQdl_p2 = (Qbar->p2[cell2] - Qbar->p2[cell1]) / std::abs(mesh->face_dist[idx]);
                        double dQdl_p3 = (Qbar->p3[cell2] - Qbar->p3[cell1]) / std::abs(mesh->face_dist[idx]);
                        double dQdl_p4 = (Qbar->p4[cell2] - Qbar->p4[cell1]) / std::abs(mesh->face_dist[idx]);

                        // Find dQdm
                        double dQdm_p1 = 0.5*( Qbar->Qxp1[cell1]*mesh->face_mx[idx]+Qbar->Qyp1[cell1]*mesh->face_my[idx] + Qbar->Qxp1[cell2]*mesh->face_mx[idx]+Qbar->Qyp1[cell2]*mesh->face_my[idx] );
                        double dQdm_p2 = 0.5*( Qbar->Qxp2[cell1]*mesh->face_mx[idx]+Qbar->Qyp2[cell1]*mesh->face_my[idx] + Qbar->Qxp2[cell2]*mesh->face_mx[idx]+Qbar->Qyp2[cell2]*mesh->face_my[idx] );
                        double dQdm_p3 = 0.5*( Qbar->Qxp3[cell1]*mesh->face_mx[idx]+Qbar->Qyp3[cell1]*mesh->face_my[idx] + Qbar->Qxp3[cell2]*mesh->face_mx[idx]+Qbar->Qyp3[cell2]*mesh->face_my[idx] );
                        double dQdm_p4 = 0.5*( Qbar->Qxp4[cell1]*mesh->face_mx[idx]+Qbar->Qyp4[cell1]*mesh->face_my[idx] + Qbar->Qxp4[cell2]*mesh->face_mx[idx]+Qbar->Qyp4[cell2]*mesh->face_my[idx] );

                        // Find dQdx and dQdy
                        double dQdy_p1 = dQdl_p1*mesh->face_ly[idx] - (dQdm_p1*mesh->face_lx[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdy_p2 = dQdl_p2*mesh->face_ly[idx] - (dQdm_p2*mesh->face_lx[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdy_p3 = dQdl_p3*mesh->face_ly[idx] - (dQdm_p3*mesh->face_lx[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdy_p4 = dQdl_p4*mesh->face_ly[idx] - (dQdm_p4*mesh->face_lx[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];

                        double dQdx_p1 = dQdl_p1*mesh->face_lx[idx] + (dQdm_p1*mesh->face_ly[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdx_p2 = dQdl_p2*mesh->face_lx[idx] + (dQdm_p2*mesh->face_ly[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdx_p3 = dQdl_p3*mesh->face_lx[idx] + (dQdm_p3*mesh->face_ly[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];
                        double dQdx_p4 = dQdl_p4*mesh->face_lx[idx] + (dQdm_p4*mesh->face_ly[idx]*mesh->face_ly[idx]) / mesh->face_mx[idx];

                        // Find ux, uy, vx, vy, dTdx, dTdy
                        double drhodx = dQdx_p1;
                        double dudx = dQdx_p2 - drhodx * Qface->u[idx] / Qface->rho[idx];
                        double dvdx = dQdx_p3 - drhodx * Qface->v[idx] / Qface->rho[idx];

                        double drhody = dQdy_p1;
                        double dudy = dQdy_p2 - drhody * Qface->u[idx] / Qface->rho[idx];
                        double dvdy = dQdy_p3 - drhody * Qface->u[idx] / Qface->rho[idx];

                        double dPdx = (freestream->gamma-1)*(dQdx_p4 - 0.5*drhodx*(Qface->u[idx]*Qface->u[idx]+Qface->v[idx]*Qface->v[idx]) - Qface->rho[idx]*(Qface->u[idx]*dudx + Qface->v[idx]*dvdx) );
                        double dPdy = (freestream->gamma-1)*(dQdy_p4 - 0.5*drhody*(Qface->u[idx]*Qface->u[idx]+Qface->v[idx]*Qface->v[idx]) - Qface->rho[idx]*(Qface->u[idx]*dudy + Qface->v[idx]*dvdy) );

                        // Compute tau xx, tau xy, tau yy, dTdx, dTdy
                        Qface->comp_vis(idx, dudx, dudy, dvdx, dvdy, dPdx, dPdy, drhodx, drhody);
                        break;
                }
            }   

            break;
    }
}