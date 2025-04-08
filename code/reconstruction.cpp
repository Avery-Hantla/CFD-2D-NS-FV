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
    std::vector<double> sum_temp1x_p1(size->num_faces), sum_temp1x_p2(size->num_faces), sum_temp1x_p3(size->num_faces), sum_temp1x_p4(size->num_faces);
    std::vector<double> sum_temp1y_p1(size->num_faces), sum_temp1y_p2(size->num_faces), sum_temp1y_p3(size->num_faces), sum_temp1y_p4(size->num_faces);

    int cell_1, cell_2;
    double Vn;
    class_flow wall_BC;
    wall_BC.gamma = freestream->gamma;
    switch (inputs->order) {
        case 1: // First Order FV Scheme
            for (int idx = 0; idx < size->num_faces; idx ++) {
                cell_1 = mesh->face_cell1[idx];

                Qfaces_c1->p1[idx] = Qbar->p1[cell_1];
                Qfaces_c1->p2[idx] = Qbar->p2[cell_1];
                Qfaces_c1->p3[idx] = Qbar->p3[cell_1];
                Qfaces_c1->p4[idx] = Qbar->p4[cell_1];

                cell_2 = mesh->face_cell2[idx];
                // Apply BC
                switch (cell_2) { 
                    case -1: // Freestream BC
                        Qfaces_c2->p1[idx] = freestream->Q1;
                        Qfaces_c2->p2[idx] = freestream->Q2;
                        Qfaces_c2->p3[idx] = freestream->Q3;
                        Qfaces_c2->p4[idx] = freestream->Q4;
                        break;
                    case -2:  // Wall BC
                        wall_BC.rho = Qbar->rho[cell_1];
                        wall_BC.P = Qbar->P[cell_1];
                        Vn = Qbar->u[cell_1]*mesh->face_nx[idx] + Qbar->v[cell_1]*mesh->face_ny[idx];
                        wall_BC.u = Qbar->u[cell_1] - 2.0*Vn*mesh->face_nx[idx];
                        wall_BC.v = Qbar->v[cell_1] - 2.0*Vn*mesh->face_ny[idx];
                        wall_BC.updateQ();

                        Qfaces_c2->p1[idx] = wall_BC.Q1;
                        Qfaces_c2->p2[idx] = wall_BC.Q2;
                        Qfaces_c2->p3[idx] = wall_BC.Q3;
                        Qfaces_c2->p4[idx] = wall_BC.Q4;
                        break;
                    case -3: // Extrapolation BC
                        Qfaces_c2->p1[idx] = Qbar->p1[cell_1];
                        Qfaces_c2->p2[idx] = Qbar->p2[cell_1];
                        Qfaces_c2->p3[idx] = Qbar->p3[cell_1];
                        Qfaces_c2->p4[idx] = Qbar->p4[cell_1];
                        break;
                    default:
                        Qfaces_c2->p1[idx] = Qbar->p1[cell_2];
                        Qfaces_c2->p2[idx] = Qbar->p2[cell_2];
                        Qfaces_c2->p3[idx] = Qbar->p3[cell_2];
                        Qfaces_c2->p4[idx] = Qbar->p4[cell_2];
                        break;
                }
            }
            Qfaces_c1->updateflow();
            Qfaces_c2->updateflow();
            break;

        case 2: // Second Order FV Scheme
            std::fill(sum_temp1x_p1.begin(), sum_temp1x_p1.end(), 0);
            std::fill(sum_temp1x_p2.begin(), sum_temp1x_p2.end(), 0);
            std::fill(sum_temp1x_p3.begin(), sum_temp1x_p3.end(), 0);
            std::fill(sum_temp1x_p4.begin(), sum_temp1x_p4.end(), 0);

            std::fill(sum_temp1y_p1.begin(), sum_temp1y_p1.end(), 0);
            std::fill(sum_temp1y_p2.begin(), sum_temp1y_p2.end(), 0);
            std::fill(sum_temp1y_p3.begin(), sum_temp1y_p3.end(), 0);
            std::fill(sum_temp1y_p4.begin(), sum_temp1y_p4.end(), 0);
            // Find dQdx and dQdy
            // double dQdx_p1[size->num_faces], dQdx_p2[size->num_faces], dQdx_p3[size->num_faces], dQdx_p4[size->num_faces];
            // double dQdy_p1[size->num_faces], dQdy_p2[size->num_faces], dQdy_p3[size->num_faces], dQdy_p4[size->num_faces];
            // double sum_temp1x_p1[size->num_faces], sum_temp1x_p2[size->num_faces], sum_temp1x_p3[size->num_faces], sum_temp1x_p4[size->num_faces];
            // double sum_temp1y_p1[size->num_faces], sum_temp1y_p2[size->num_faces], sum_temp1y_p3[size->num_faces], sum_temp1y_p4[size->num_faces];

            for (int idx = 0; idx < size->num_faces; idx ++) {
                cell_1 = mesh->face_cell1[idx];
                cell_2 = mesh->face_cell2[idx];

                double Qjmi_p1, Qjmi_p2, Qjmi_p3, Qjmi_p4;
                
                switch (cell_2) {
                    case -1: // Freestream
                        Qjmi_p1 = freestream->Q1 - Qbar->p1[cell_1];
                        Qjmi_p2 = freestream->Q2 - Qbar->p2[cell_1];
                        Qjmi_p3 = freestream->Q3 - Qbar->p3[cell_1];
                        Qjmi_p4 = freestream->Q4 - Qbar->p4[cell_1];
                        break;
                    case -2: // Wall BC
                        wall_BC.rho = Qbar->rho[cell_1];
                        wall_BC.P = Qbar->P[cell_1];
                        Vn = Qbar->u[cell_1]*mesh->face_nx[idx] + Qbar->v[cell_1]*mesh->face_ny[idx];
                        wall_BC.u = Qbar->u[cell_1] - 2.0*Vn*mesh->face_nx[idx];
                        wall_BC.v = Qbar->v[cell_1] - 2.0*Vn*mesh->face_ny[idx];
                        wall_BC.updateQ();

                        Qjmi_p1 = wall_BC.Q1 - Qbar->p1[cell_1];
                        Qjmi_p2 = wall_BC.Q2 - Qbar->p2[cell_1];
                        Qjmi_p3 = wall_BC.Q3 - Qbar->p3[cell_1];
                        Qjmi_p4 = wall_BC.Q4 - Qbar->p4[cell_1];
                        break;
                    case -3: // Extrapolation BC
                        Qjmi_p1 = 0;
                        Qjmi_p2 = 0;
                        Qjmi_p3 = 0;
                        Qjmi_p4 = 0;
                        break;
                    default:
                        Qjmi_p1 = Qbar->p1[cell_2] - Qbar->p1[cell_1];
                        Qjmi_p2 = Qbar->p2[cell_2] - Qbar->p2[cell_1];
                        Qjmi_p3 = Qbar->p3[cell_2] - Qbar->p3[cell_1];
                        Qjmi_p4 = Qbar->p4[cell_2] - Qbar->p4[cell_1];
                        break;
                } 

                sum_temp1x_p1[cell_1] += Qjmi_p1*mesh->face_dxj[idx];
                sum_temp1x_p2[cell_1] += Qjmi_p2*mesh->face_dxj[idx];
                sum_temp1x_p3[cell_1] += Qjmi_p3*mesh->face_dxj[idx];
                sum_temp1x_p4[cell_1] += Qjmi_p4*mesh->face_dxj[idx];

                sum_temp1y_p1[cell_1] += Qjmi_p1*mesh->face_dyj[idx];
                sum_temp1y_p2[cell_1] += Qjmi_p2*mesh->face_dyj[idx];
                sum_temp1y_p3[cell_1] += Qjmi_p3*mesh->face_dyj[idx];
                sum_temp1y_p4[cell_1] += Qjmi_p4*mesh->face_dyj[idx];

                if (cell_2 >= 0) {
                    sum_temp1x_p1[cell_2] += Qjmi_p1*mesh->face_dxj[idx];
                    sum_temp1x_p2[cell_2] += Qjmi_p2*mesh->face_dxj[idx];
                    sum_temp1x_p3[cell_2] += Qjmi_p3*mesh->face_dxj[idx];
                    sum_temp1x_p4[cell_2] += Qjmi_p4*mesh->face_dxj[idx];

                    sum_temp1y_p1[cell_2] += Qjmi_p1*mesh->face_dyj[idx];
                    sum_temp1y_p2[cell_2] += Qjmi_p2*mesh->face_dyj[idx];
                    sum_temp1y_p3[cell_2] += Qjmi_p3*mesh->face_dyj[idx];
                    sum_temp1y_p4[cell_2] += Qjmi_p4*mesh->face_dyj[idx];
                }
            }

            for (int idx = 0; idx < size->num_cells; idx ++) {
                Qbar->Qxp1[idx] = (sum_temp1x_p1[idx]*mesh->Iyy[idx] - sum_temp1y_p1[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qxp2[idx] = (sum_temp1x_p2[idx]*mesh->Iyy[idx] - sum_temp1y_p2[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qxp3[idx] = (sum_temp1x_p3[idx]*mesh->Iyy[idx] - sum_temp1y_p3[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qxp4[idx] = (sum_temp1x_p4[idx]*mesh->Iyy[idx] - sum_temp1y_p4[idx]*mesh->Ixy[idx]);//mesh->delta[idx];

                Qbar->Qyp1[idx] = (sum_temp1y_p1[idx]*mesh->Ixx[idx] - sum_temp1x_p1[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qyp2[idx] = (sum_temp1y_p2[idx]*mesh->Ixx[idx] - sum_temp1x_p2[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qyp3[idx] = (sum_temp1y_p3[idx]*mesh->Ixx[idx] - sum_temp1x_p3[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
                Qbar->Qyp4[idx] = (sum_temp1y_p4[idx]*mesh->Ixx[idx] - sum_temp1x_p4[idx]*mesh->Ixy[idx]);//mesh->delta[idx];
            }

            // Reconstruct Faces
            for (int idx = 0; idx < size->num_faces; idx ++) {
                cell_1 = mesh->face_cell1[idx];
                cell_2 = mesh->face_cell2[idx];

                Qfaces_c1->p1[idx] = Qbar->p1[cell_1] + Qbar->Qxp1[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp1[cell_1]*mesh->face_cell1_dy[idx];
                Qfaces_c1->p2[idx] = Qbar->p2[cell_1] + Qbar->Qxp2[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp2[cell_1]*mesh->face_cell1_dy[idx];
                Qfaces_c1->p3[idx] = Qbar->p3[cell_1] + Qbar->Qxp3[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp3[cell_1]*mesh->face_cell1_dy[idx];
                Qfaces_c1->p4[idx] = Qbar->p4[cell_1] + Qbar->Qxp4[cell_1]*mesh->face_cell1_dx[idx] + Qbar->Qyp4[cell_1]*mesh->face_cell1_dy[idx];

                Qfaces_c1->updateflow_int(idx);

                cell_2 = mesh->face_cell2[idx];
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
            break;

        default:
            std::cout << "ERROR: ORDER NOT IMPLEMENTED \n";
    };
}