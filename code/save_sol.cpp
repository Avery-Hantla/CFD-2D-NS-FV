/////////////////////////////////////////////////////////
//         Function Monitor and Output Solution
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <filesystem> 

#include "class_mesh.hpp"
#include "class_q.hpp"
#include "class_residual.hpp"
#include "class_f.hpp"
#include "struct_inputs.hpp"
#include "struct_size.hpp"
#include "struct_report.hpp"

void save(class_Q* Qbar, class_mesh* mesh, class_residual* residual, class_flow* freestream, struct_inputs* inputs, struct_size* size, struct_report* report, int ndx) {
    // Monitor 
    if ((ndx % inputs->monitor_step) == 0) {

        /////////////////////// Reconstruct to nodes ///////////////////////
        double face_reconst_p1, face_reconst_p2, face_reconst_p3, face_reconst_p4;
        double rho, E, u, v, P, force, L = 0, D = 0, CL, CD;
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
            E = face_reconst_p2;
            u = face_reconst_p3/rho;
            v = face_reconst_p4/rho;
            P = (E-(0.5*rho*(u*u+v*v)))*(Qbar->gamma-1);

            // Find CL and CD
            force = P*mesh->face_area[face];
            L -= force*mesh->face_ny[face];
            D -= force*mesh->face_nx[face];
        }

        CL = L/(0.5*report->rho*std::sqrt(report->u * report->u + report->v * report->v)*std::sqrt(report->u * report->u + report->v * report->v)*report->area);
        CD = D/(0.5*report->rho*std::sqrt(report->u * report->u + report->v * report->v)*std::sqrt(report->u * report->u + report->v * report->v)*report->area);

        //////////////// Output ////////////////
        double res1 = *std::max_element(residual->p1.begin(),residual->p1.end());
        double res2 = *std::max_element(residual->p2.begin(),residual->p2.end());
        double res3 = *std::max_element(residual->p3.begin(),residual->p3.end());
        double res4 = *std::max_element(residual->p4.begin(),residual->p4.end());

        // Print to terminal
        std::cout << "iter: " << ndx << ", CL: " << CL << ", CD: " << CD << ", res: " << res1 << ", " << res2 << ", " << res3 << ", " << res4 << std::endl;

        // Save to output file
        std::ofstream output("../sol/output.dat", std::ios_base::app | std::ios_base::out);
        output << ndx << " " << CL << " " << CD << std::endl;
        output.close();

        // Save to residual file
        std::ofstream res_history("../sol/res_history.dat", std::ios_base::app | std::ios_base::out);
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
        
        /////////////////////// Save Tecplot Solution  ///////////////////////
        // Save volume data
        std::ofstream output;
        output.open ("../sol/tec_sol_" + std::to_string(ndx) + ".dat");
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
        output.open ("../sol/tec_surf_" + std::to_string(ndx) + ".dat");
        if (output.is_open()) {
            output << "TITLE     = \"Translation of CGNS file merged_sol_aver.cgns\"" << std::endl;
            output << "VARIABLES = \"CoordinateX\" \n \"CoordinateY\" \n \"rho\" \n \"u\" \n \"v\" \n \"P\" " << std::endl;
            output << "ZONE T=\"ZONE1\"\nSTRANDID=0, SOLUTIONTIME=" << ndx << "\nI=" << size->numWALL+1 << ", J=1, K=1, ZONETYPE=Ordered\n";
            output << "DATAPACKING=POINT\nDT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )\n";

            for (int idx = 0; idx < size->numWALL+1; idx++) {
                output << mesh->x[idx] << " " << mesh->y[idx] << " ";
                output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save Tecplot  File \n";
        }
        output.close();

        /////////////////////// Save MATLAB Solution  ///////////////////////
        // std::ofstream output;
        // Save volume data
        output.open ("../sol/mat_sol_" + std::to_string(ndx) + ".dat");
        if (output.is_open()) {
            for (int idx = 0; idx < size->num_cells; idx++) {
                output << Qbar->rho[idx] << " " << Qbar->u[idx] << " " << Qbar->v[idx] << " " << Qbar->P[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save MATLAB File \n";
        }
        output.close();

        // Save surface data
        output.open ("../sol/mat_surf_" + std::to_string(ndx) + ".dat");
        if (output.is_open()) {
            for (int idx = 0; idx < size->numWALL+1; idx++) {
                output << mesh->x[idx] << " " << mesh->y[idx] << " ";
                output << rho[idx] << " " << u[idx] << " " << v[idx] << " " << P[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save MATLAB Surface File \n";
        }
        output.close();
    }
}