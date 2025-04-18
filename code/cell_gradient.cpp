#include <vector>

#include "class_mesh.hpp"
#include "class_q.hpp"
#include "class_flow.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"

#include "minmod.cpp"

void compute_gradient(class_Q* Qbar, class_mesh* mesh, class_flow* freestream, struct_size* size, struct_inputs* inputs) {
    int cell_1, cell_2;
    double Vn;
    class_flow wall_BC;
    wall_BC.gamma = freestream->gamma;
    std::vector<double> sum_temp1x_p1(size->num_faces), sum_temp1x_p2(size->num_faces), sum_temp1x_p3(size->num_faces), sum_temp1x_p4(size->num_faces);
    std::vector<double> sum_temp1y_p1(size->num_faces), sum_temp1y_p2(size->num_faces), sum_temp1y_p3(size->num_faces), sum_temp1y_p4(size->num_faces);

    std::vector<double> face_slopex_p1(size->num_faces), face_slopex_p2(size->num_faces), face_slopex_p3(size->num_faces), face_slopex_p4(size->num_faces);
    std::vector<double> face_slopey_p1(size->num_faces), face_slopey_p2(size->num_faces), face_slopey_p3(size->num_faces), face_slopey_p4(size->num_faces);

    std::fill(sum_temp1x_p1.begin(), sum_temp1x_p1.end(), 0);
    std::fill(sum_temp1x_p2.begin(), sum_temp1x_p2.end(), 0);
    std::fill(sum_temp1x_p3.begin(), sum_temp1x_p3.end(), 0);
    std::fill(sum_temp1x_p4.begin(), sum_temp1x_p4.end(), 0);

    std::fill(sum_temp1y_p1.begin(), sum_temp1y_p1.end(), 0);
    std::fill(sum_temp1y_p2.begin(), sum_temp1y_p2.end(), 0);
    std::fill(sum_temp1y_p3.begin(), sum_temp1y_p3.end(), 0);
    std::fill(sum_temp1y_p4.begin(), sum_temp1y_p4.end(), 0);

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
            case -4: // NO SLIP WALL BC
                wall_BC.rho = Qbar->rho[idx];
                wall_BC.P = Qbar->P[idx];
                wall_BC.u = -Qbar->u[idx];
                wall_BC.v = -Qbar->v[idx];
                wall_BC.updateQ();

                Qjmi_p1 = wall_BC.Q1 - Qbar->p1[cell_1];
                Qjmi_p2 = wall_BC.Q2 - Qbar->p2[cell_1];
                Qjmi_p3 = wall_BC.Q3 - Qbar->p3[cell_1];
                Qjmi_p4 = wall_BC.Q4 - Qbar->p4[cell_1];
                break;
            default:
                Qjmi_p1 = Qbar->p1[cell_2] - Qbar->p1[cell_1];
                Qjmi_p2 = Qbar->p2[cell_2] - Qbar->p2[cell_1];
                Qjmi_p3 = Qbar->p3[cell_2] - Qbar->p3[cell_1];
                Qjmi_p4 = Qbar->p4[cell_2] - Qbar->p4[cell_1];
                break;
        } 

        face_slopex_p1[idx] = Qjmi_p1*mesh->face_dxj[idx];
        face_slopex_p2[idx] = Qjmi_p2*mesh->face_dxj[idx];
        face_slopex_p3[idx] = Qjmi_p3*mesh->face_dxj[idx];
        face_slopex_p4[idx] = Qjmi_p4*mesh->face_dxj[idx];

        face_slopey_p1[idx] = Qjmi_p1*mesh->face_dyj[idx];
        face_slopey_p2[idx] = Qjmi_p2*mesh->face_dyj[idx];
        face_slopey_p3[idx] = Qjmi_p3*mesh->face_dyj[idx];
        face_slopey_p4[idx] = Qjmi_p4*mesh->face_dyj[idx];
    }

    for (int idx = 0; idx < size->num_faces; idx ++) {
        cell_1 = mesh->face_cell1[idx];
        cell_2 = mesh->face_cell2[idx];

        sum_temp1x_p1[cell_1] += face_slopex_p1[idx];
        sum_temp1x_p2[cell_1] += face_slopex_p2[idx];
        sum_temp1x_p3[cell_1] += face_slopex_p3[idx];
        sum_temp1x_p4[cell_1] += face_slopex_p4[idx];

        sum_temp1y_p1[cell_1] += face_slopey_p1[idx];
        sum_temp1y_p2[cell_1] += face_slopey_p2[idx];
        sum_temp1y_p3[cell_1] += face_slopey_p3[idx];
        sum_temp1y_p4[cell_1] += face_slopey_p4[idx];

        if (cell_2 >= 0) {
            sum_temp1x_p1[cell_2] += face_slopex_p1[idx];
            sum_temp1x_p2[cell_2] += face_slopex_p2[idx];
            sum_temp1x_p3[cell_2] += face_slopex_p3[idx];
            sum_temp1x_p4[cell_2] += face_slopex_p4[idx];

            sum_temp1y_p1[cell_2] += face_slopey_p1[idx];
            sum_temp1y_p2[cell_2] += face_slopey_p2[idx];
            sum_temp1y_p3[cell_2] += face_slopey_p3[idx];
            sum_temp1y_p4[cell_2] += face_slopey_p4[idx];
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
}