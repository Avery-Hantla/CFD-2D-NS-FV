/////////////////////////////////////////////////////////
//           Function To Compute Roe Flux
/////////////////////////////////////////////////////////

#include <Eigen>

#include "../support/class_q.hpp"
#include "../support/class_mesh.hpp"
#include "../support/class_f.hpp"
#include "../support/class_flow.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_BC.hpp"

void roe(class_F* F_roe, class_Q* Qface_c1, class_Q* Qface_c2, class_mesh* mesh, class_flow* freestream, struct_size* size, struct_BC* BC) {
    class_F F_c1, F_c2;
    F_c1.init(size->num_faces, freestream->gamma);
    F_c2.init(size->num_faces, freestream->gamma);
    
    F_c1.update(mesh, Qface_c1, size); 
    F_c2.update(mesh, Qface_c2, size);


    Eigen::MatrixXf GAMMA(3,3);
    Eigen::MatrixXf R(3,3);
    Eigen::MatrixXf R_inv(3,3);
    Eigen::MatrixXf A(3,3);

    for (int idx = 0; idx < size->num_faces; idx++) {
        // Find Left Values
        double rhoL = Qface_c1->rho[idx];
        double uL = Qface_c1->u[idx];
        double PL = Qface_c1->P[idx];
        double EL = Qface_c1->E[idx];
        double HL = (EL + PL)/rhoL;

        // Find Right Values
        double rhoR = Qface_c2->rho[idx];
        double uR = Qface_c2->u[idx];
        double PR = Qface_c2->P[idx];
        double ER = Qface_c2->E[idx];
        double HR = (ER + PR)/rhoR;

        // Find roes average u and c
        double rho_bar = std::sqrt(rhoL*rhoR);
        double Vn_bar = (std::sqrt(rhoL)*uL+std::sqrt(rhoR)*uR)/(std::sqrt(rhoL)+sqrt(rhoR));
        double H_bar = (std::sqrt(rhoL)*HL+std::sqrt(rhoR)*HR)/(std::sqrt(rhoL)+sqrt(rhoR));
        double c_bar = std::sqrt((freestream->gamma-1)*(H_bar-(Vn_bar*Vn_bar)/2));
        double P_bar = (c_bar*c_bar)*rho_bar/freestream->gamma;

        // Eigen::MatrixXf GAMMA(3,3);
        // Eigen::MatrixXf R(3,3);
        // Eigen::MatrixXf R_inv(3,3);
        GAMMA << Vn_bar, 0, 0, 
                0, Vn_bar + c_bar, 0, 
                0, 0, Vn_bar - c_bar;      
        R << 2/(Vn_bar*Vn_bar), (2*rho_bar*(35*P_bar + 5*rho_bar*Vn_bar*Vn_bar - 2*std::sqrt(35)*Vn_bar*std::sqrt(P_bar*rho_bar)))/(245*P_bar*P_bar + 42*P_bar*rho_bar*Vn_bar*Vn_bar + 5*rho_bar*rho_bar*std::pow(Vn_bar,4)),  (2*rho_bar*(35*P_bar + 5*rho_bar*Vn_bar*Vn_bar + 2*std::sqrt(35)*Vn_bar*std::sqrt(P_bar*rho_bar)))/(245*P_bar*P_bar + 42*P_bar*rho_bar*Vn_bar*Vn_bar + 5*rho_bar*rho_bar*std::pow(Vn_bar,4)),
            2/Vn_bar,  (10*rho_bar*rho_bar*std::pow(Vn_bar,3) + 14*std::sqrt(7)*P_bar*std::sqrt(5*P_bar*rho_bar) + 42*P_bar*rho_bar*Vn_bar - 2*std::sqrt(7)*rho_bar*Vn_bar*Vn_bar*std::sqrt(5*P_bar*rho_bar))/(245*P_bar*P_bar + 42*P_bar*rho_bar*Vn_bar*Vn_bar + 5*rho_bar*rho_bar*std::pow(Vn_bar,4)),   (10*rho_bar*rho_bar*std::pow(Vn_bar,3) - 14*std::sqrt(7)*P_bar*std::sqrt(5*P_bar*rho_bar) + 42*P_bar*rho_bar*Vn_bar + 2*std::sqrt(7)*rho_bar*Vn_bar*Vn_bar*std::sqrt(5*P_bar*rho_bar))/(245*P_bar*P_bar + 42*P_bar*rho_bar*Vn_bar*Vn_bar + 5*rho_bar*rho_bar*std::pow(Vn_bar,4)),
            1, 1, 1;
        R_inv = R.inverse();

        A = R*GAMMA*R_inv;

        double dis_temp_1 = A(0,1)*

        std::cout << "Here is mat*mat:\n" << A(1,1) << std::endl;
        // double GAMMA[3][3] = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
        // double R[3][3] = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
        // double A[3][3];

        // for (int idx = 0; idx < 3; idx++)
        //     for (int jdx = 0; jdx < 3; jdx++) {
        //         for (int udx = 0; udx < 3; udx++)
        //             results[idx][jdx] += GAMMA[idx][udx] * matrix2[udx][jdx];
        //     }
        // }

        // Find Vn and c averag 
        double Vn_avg = (F_c1.Vn[idx] + F_c2.Vn[idx])/2;
        double c_avg = (Qface_c1->c[idx] + Qface_c2->c[idx])/2;

        Qface_c1->Vn_avg[idx] = Vn_avg;
        Qface_c1->c_avg[idx] = c_avg;
        Qface_c2->Vn_avg[idx] = Vn_avg;
        Qface_c2->c_avg[idx] = c_avg;

        // Compute Roe Flux
        F_roe->p1[idx] = (F_c1.p1[idx]+F_c2.p1[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p1[idx] - Qface_c1->p1[idx]);
        F_roe->p2[idx] = (F_c1.p2[idx]+F_c2.p2[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p2[idx] - Qface_c1->p2[idx]);
        F_roe->p3[idx] = (F_c1.p3[idx]+F_c2.p3[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p3[idx] - Qface_c1->p3[idx]);
        F_roe->p4[idx] = (F_c1.p4[idx]+F_c2.p4[idx])/2 - (std::abs(Vn_avg) + c_avg)/2 * (Qface_c2->p4[idx] - Qface_c1->p4[idx]);
    }
}