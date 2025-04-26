/////////////////////////////////////////////////////////
//         Function Monitor and Output Solution
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <filesystem> 

#include "../support/class_mesh.hpp"
#include "../support/class_q.hpp"
#include "../support/class_residual.hpp"
#include "../support/class_f.hpp"
#include "../support/struct_inputs.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_report.hpp"

void save(class_Q* Qbar, class_mesh* mesh, class_residual* residual, class_flow* freestream, struct_inputs* inputs, struct_size* size, struct_report* report, int ndx) {
    // Monitor 
    if ((ndx % inputs->monitor_step) == 0) {

        /////////////////////// Reconstruct to nodes ///////////////////////
        double face_reconst_p1, face_reconst_p2, face_reconst_p3, face_reconst_p4;
        double rho, E, u, v, P, force, Fy = 0, Fx = 0, L, D, CL, CD, tau_w;
        for (int idx = 0; idx < size->numWALL; idx++) {
            int face = mesh->BC_faces[idx];

            int cell1 = mesh->face_cell1[face];

            int cellx = mesh->cell_centerx[cell1];
            int celly = mesh->cell_centery[cell1];

            int facex = mesh->face_centerx[face];
            int facey = mesh->face_centery[face];

            face_reconst_p1 = Qbar->p1[cell1] + ((facex-cellx) * Qbar->Qxp1[cell1] + (facey-celly) * Qbar->Qyp1[cell1]);
            face_reconst_p2 = Qbar->p2[cell1] + ((facex-cellx) * Qbar->Qxp2[cell1] + (facey-celly) * Qbar->Qyp2[cell1]);
            face_reconst_p3 = Qbar->p3[cell1] + ((facex-cellx) * Qbar->Qxp3[cell1] + (facey-celly) * Qbar->Qyp3[cell1]);
            face_reconst_p4 = Qbar->p4[cell1] + ((facex-cellx) * Qbar->Qxp4[cell1] + (facey-celly) * Qbar->Qyp4[cell1]);

            // Find flow parameters
            rho = face_reconst_p1;
            u = face_reconst_p2/rho;
            v = face_reconst_p3/rho;
            E = face_reconst_p4;
            P = (E-(0.5*rho*(u*u+v*v)))*(Qbar->gamma-1);

            if (inputs->eqn == 2) {
                double dn = std::abs(std::sqrt(mesh->face_cell1_dx[face] * mesh->face_cell1_dx[face] + mesh->face_cell1_dy[face] * mesh->face_cell1_dy[face]));
                double dudn = -Qbar->u[cell1]/dn;
                double dvdn = -Qbar->v[cell1]/dn;
                tau_w = Qbar->mu * std::sqrt(std::abs(dudn*dudn + dvdn*dvdn));
                if (mesh->face_ny[face] < 0) {
                    Fx += tau_w*mesh->face_area[face]*(-mesh->face_ny[face]);
                    Fy += tau_w*mesh->face_area[face]*(mesh->face_nx[face]);
                } else {
                    Fx += tau_w*mesh->face_area[face]*(mesh->face_ny[face]);
                    Fy += tau_w*mesh->face_area[face]*(-mesh->face_nx[face]);
                }
            }

            // Find Fx and Fy
            force = P*mesh->face_area[face];
            Fy += force*mesh->face_ny[face];
            Fx += force*mesh->face_nx[face];
        }

        L = Fx*(-freestream->flowy) + Fy*freestream->flowx;
        D = Fy*freestream->flowy + Fx*freestream->flowx;
        CL = L/(0.5*report->rho*std::sqrt(report->u * report->u + report->v * report->v)*std::sqrt(report->u * report->u + report->v * report->v)*report->length);
        CD = D/(0.5*report->rho*std::sqrt(report->u * report->u + report->v * report->v)*std::sqrt(report->u * report->u + report->v * report->v)*report->length);

        //////////////// Output ////////////////
        // Find L2 Norm of residuals
        double res1 = 0, res2 = 0, res3 = 0, res4 = 0;
        for (int idx = 0; idx < size->num_cells; idx++) {
            res1 += residual->p1[idx] * residual->p1[idx];
            res2 += residual->p2[idx] * residual->p2[idx];
            res3 += residual->p3[idx] * residual->p3[idx];
            res4 += residual->p4[idx] * residual->p4[idx];
        }
        res1 = std::sqrt(res1);
        res2 = std::sqrt(res2);
        res3 = std::sqrt(res3);
        res4 = std::sqrt(res4);

        // Print to terminal
        std::cout << "iter: " << ndx << ", CL: " << CL << ", CD: " << CD << ", res: " << res1 << ", " << res2 << ", " << res3 << ", " << res4 << std::endl;

        // Save to output file
        std::ofstream output("sol/output.dat", std::ios_base::app | std::ios_base::out);
        output << ndx << " " << CL << " " << CD << std::endl;
        output.close();

        // Save to residual file
        std::ofstream res_history("sol/res_history.dat", std::ios_base::app | std::ios_base::out);
        res_history << ndx << " " << res1 << " " << res2 << " " << res3 << " " << res4 << std::endl;
        res_history.close();
    }

    // Output
    if ((ndx % inputs->output_step) == 0) {
        std::cout <<  "Saving Solution \n";

        /////////////////////// Reconstruct to nodes ///////////////////////
        std::vector<int> pointcounter(size->num_points,0);
        std::vector<double> point_avg_p1(size->num_points,0.0), point_avg_p2(size->num_points,0.0), point_avg_p3(size->num_points,0.0), point_avg_p4(size->num_points,0.0);
        std::vector<double> rho(size->num_points,0.0), E(size->num_points,0.0), u(size->num_points,0.0), v(size->num_points,0.0), P(size->num_points,0.0);
        for (int idx = 0; idx < size->num_faces; idx++) {
            int cell1 = mesh->face_cell1[idx];
            int cell2 = mesh->face_cell2[idx];

            int point1 = mesh->face_point1[idx];
            int point2 = mesh->face_point2[idx];

            pointcounter[point1]++;
            pointcounter[point2]++;

            point_avg_p1[point1] += Qbar->p1[cell1] + ((mesh->x[point1]-mesh->cell_centerx[cell1]) * Qbar->Qxp1[cell1] + (mesh->y[point1]-mesh->cell_centery[cell1]) * Qbar->Qyp1[cell1]);
            point_avg_p2[point1] += Qbar->p2[cell1] + ((mesh->x[point1]-mesh->cell_centerx[cell1]) * Qbar->Qxp2[cell1] + (mesh->y[point1]-mesh->cell_centery[cell1]) * Qbar->Qyp2[cell1]);
            point_avg_p3[point1] += Qbar->p3[cell1] + ((mesh->x[point1]-mesh->cell_centerx[cell1]) * Qbar->Qxp3[cell1] + (mesh->y[point1]-mesh->cell_centery[cell1]) * Qbar->Qyp3[cell1]);
            point_avg_p4[point1] += Qbar->p4[cell1] + ((mesh->x[point1]-mesh->cell_centerx[cell1]) * Qbar->Qxp4[cell1] + (mesh->y[point1]-mesh->cell_centery[cell1]) * Qbar->Qyp4[cell1]);

            point_avg_p1[point2] += Qbar->p1[cell1] + ((mesh->x[point2]-mesh->cell_centerx[cell1]) * Qbar->Qxp1[cell1] + (mesh->y[point2]-mesh->cell_centery[cell1]) * Qbar->Qyp1[cell1]);
            point_avg_p2[point2] += Qbar->p2[cell1] + ((mesh->x[point2]-mesh->cell_centerx[cell1]) * Qbar->Qxp2[cell1] + (mesh->y[point2]-mesh->cell_centery[cell1]) * Qbar->Qyp2[cell1]);
            point_avg_p3[point2] += Qbar->p3[cell1] + ((mesh->x[point2]-mesh->cell_centerx[cell1]) * Qbar->Qxp3[cell1] + (mesh->y[point2]-mesh->cell_centery[cell1]) * Qbar->Qyp3[cell1]);
            point_avg_p4[point2] += Qbar->p4[cell1] + ((mesh->x[point2]-mesh->cell_centerx[cell1]) * Qbar->Qxp4[cell1] + (mesh->y[point2]-mesh->cell_centery[cell1]) * Qbar->Qyp4[cell1]);

            if (cell2 >= 0) {
                point_avg_p1[point1] += Qbar->p1[cell2] + ((mesh->x[point1]-mesh->cell_centerx[cell2]) * Qbar->Qxp1[cell2] + (mesh->y[point1]-mesh->cell_centery[cell2]) * Qbar->Qyp1[cell2]);
                point_avg_p2[point1] += Qbar->p2[cell2] + ((mesh->x[point1]-mesh->cell_centerx[cell2]) * Qbar->Qxp2[cell2] + (mesh->y[point1]-mesh->cell_centery[cell2]) * Qbar->Qyp2[cell2]);
                point_avg_p3[point1] += Qbar->p3[cell2] + ((mesh->x[point1]-mesh->cell_centerx[cell2]) * Qbar->Qxp3[cell2] + (mesh->y[point1]-mesh->cell_centery[cell2]) * Qbar->Qyp3[cell2]);
                point_avg_p4[point1] += Qbar->p4[cell2] + ((mesh->x[point1]-mesh->cell_centerx[cell2]) * Qbar->Qxp4[cell2] + (mesh->y[point1]-mesh->cell_centery[cell2]) * Qbar->Qyp4[cell2]);

                point_avg_p1[point2] += Qbar->p1[cell2] + ((mesh->x[point2]-mesh->cell_centerx[cell2]) * Qbar->Qxp1[cell2] + (mesh->y[point2]-mesh->cell_centery[cell2]) * Qbar->Qyp1[cell2]);
                point_avg_p2[point2] += Qbar->p2[cell2] + ((mesh->x[point2]-mesh->cell_centerx[cell2]) * Qbar->Qxp2[cell2] + (mesh->y[point2]-mesh->cell_centery[cell2]) * Qbar->Qyp2[cell2]);
                point_avg_p3[point2] += Qbar->p3[cell2] + ((mesh->x[point2]-mesh->cell_centerx[cell2]) * Qbar->Qxp3[cell2] + (mesh->y[point2]-mesh->cell_centery[cell2]) * Qbar->Qyp3[cell2]);
                point_avg_p4[point2] += Qbar->p4[cell2] + ((mesh->x[point2]-mesh->cell_centerx[cell2]) * Qbar->Qxp4[cell2] + (mesh->y[point2]-mesh->cell_centery[cell2]) * Qbar->Qyp4[cell2]);
                
                pointcounter[point1]++;
                pointcounter[point2]++;
            }
        }

        for (int idx = 0; idx < size->num_points; idx ++) {
            point_avg_p1[idx] = point_avg_p1[idx]/(pointcounter[idx]);
            point_avg_p2[idx] = point_avg_p2[idx]/(pointcounter[idx]);
            point_avg_p3[idx] = point_avg_p3[idx]/(pointcounter[idx]);
            point_avg_p4[idx] = point_avg_p4[idx]/(pointcounter[idx]);
        }
        // Find flow parameters
        for (int idx = 0; idx < size->num_points; idx++) {
            rho[idx] = point_avg_p1[idx];
            E[idx] = point_avg_p4[idx];
            u[idx] = point_avg_p2[idx]/rho[idx];
            v[idx] = point_avg_p3[idx]/rho[idx];
            P[idx] = (E[idx]-(0.5*rho[idx]*(u[idx]*u[idx]+v[idx]*v[idx])))*(Qbar->gamma-1);
        }

        std::vector<double> tau_w(size->numWALL,0);
        if (inputs->eqn == 2) {
            for (int idx = 0; idx < size->numWALL; idx++) {
                int face = mesh->BC_faces[idx];
                int cell1 = mesh->face_cell1[face];

                double dn = std::abs(std::sqrt(mesh->face_cell1_dx[face] * mesh->face_cell1_dx[face] + mesh->face_cell1_dy[face] * mesh->face_cell1_dy[face]));
                double dudn = -Qbar->u[cell1]/dn;
                double dvdn = -Qbar->v[cell1]/dn;
                tau_w[idx] = Qbar->mu * std::sqrt(std::abs(dudn*dudn + dvdn*dvdn));
            }
        }
        
        /////////////////////// Save Tecplot Solution  ///////////////////////
        // Save volume data
        std::ofstream output;
        if (inputs->nmax == ndx) {
            output.open ("sol/tec_sol.dat");
        } else {
            output.open ("sol/tec_sol_" + std::to_string(ndx) + ".dat");
        }
        if (output.is_open()) {
            // Output Header
            output << "VARIABLES = \"X\", \"Y\", \"density\", \"u\", \"v\", \"P\" \n";
            output << "Zone STRANDID=1, ZoneType=FETRIANGLE, SOLUTIONTIME=" << ndx << std::endl;
            output << "n=" << size->num_faces << " e=" << size->num_faces*2-mesh->num_of_BC << " DATAPACKING=POINT" <<std::endl;
            output << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )" << std::endl;

            // Output variables
            for (int idx = 0; idx < size->num_points; idx++) {
                output << mesh->x[idx] << " " << mesh->y[idx] << " ";
                output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << std::endl;
            }
            for (int idx = 0; idx < size->num_cells; idx ++) {
                output << mesh->cell_centerx[idx] << " " << mesh->cell_centery[idx] << " ";
                output << Qbar->rho[idx] << " " << Qbar->u[idx] << " " << Qbar->v[idx] << " " << Qbar->P[idx] << std::endl;
            }
            
            // Output connectivity matrix
            for (int idx = 0; idx < mesh->connect_out.size()/3; idx ++) {
                output << mesh->connect_out[3*idx + 0] << " " << mesh->connect_out[3*idx + 1] << " " << mesh->connect_out[3*idx + 2] << std::endl;
            }

        } else {
            std::cout << "ERROR: Cannot Save Tecplot File \n";
        }
        output.close();
        
        // Save surface data
        if (inputs->nmax == ndx) {
            output.open ("sol/tec_surf.dat");
        } else {
            output.open ("sol/tec_surf_" + std::to_string(ndx) + ".dat");
        }
        if (output.is_open()) {
            output << "TITLE     = \"Translation of CGNS file merged_sol_aver.cgns\"" << std::endl;
            output << "VARIABLES = \"CoordinateX\" \n \"CoordinateY\" \n \"rho\" \n \"u\" \n \"v\" \n \"P\" \n \"tau_w\"" << std::endl;
            output << "ZONE T=\"ZONE1\"\nSTRANDID=0, SOLUTIONTIME=" << ndx << "\nI=" << size->numWALL+1 << ", J=1, K=1, ZONETYPE=Ordered\n";
            output << "DATAPACKING=POINT\nDT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n";

            for (int idx = 0; idx < size->numWALL+1; idx++) {
                output << mesh->x[idx] << " " << mesh->y[idx] << " ";
                output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << " " << tau_w[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save Tecplot Surface File \n";
        }
        output.close();

        /////////////////////// Save MATLAB Solution  ///////////////////////
        // std::ofstream output;
        // Save volume data
        if (inputs->nmax == ndx) {
            output.open ("sol/mat_sol.dat");
        } else {
            output.open ("sol/mat_sol_" + std::to_string(ndx) + ".dat");
        }
        if (output.is_open()) {
            for (int idx = 0; idx < size->num_cells; idx++) {
                output << Qbar->rho[idx] << " " << Qbar->u[idx] << " " << Qbar->v[idx] << " " << Qbar->P[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save MATLAB File \n";
        }
        output.close();

        // Save surface data
        if (inputs->nmax == ndx) {
            output.open ("sol/mat_surf.dat");
        } else {
            output.open ("sol/mat_surf_" + std::to_string(ndx) + ".dat");
        }
        if (output.is_open()) {
            for (int idx = 0; idx < size->numWALL+1; idx++) {
                output << mesh->x[idx] << " " << mesh->y[idx] << " ";
                output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << " " << tau_w[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save MATLAB Surface File \n";
        }
        output.close();
    }
}