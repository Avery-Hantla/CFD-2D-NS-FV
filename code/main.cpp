/* ////////////////////////////////////////////////////////// 
              FV 2D Euler Code for Mixed Meshes
                        Avery Hantla
                         March 2025
*/ //////////////////////////////////////////////////////////
#include <iostream>
// #include <cgnslib.h>

// Load classes, structs, etc.
#include "class_mesh.hpp"
#include "class_q.hpp"
#include "class_f.hpp"
#include "class_flow.hpp"

#include "struct_BC.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"
#include "struct_time.hpp"
#include "struct_report.hpp"

// Load Functions 
#include "read_inputs.cpp"
#include "readmesh.cpp"
// #include "ssp_rk2.cpp"
#include "ssp_rk3.cpp"
#include "save_sol.cpp"
#include "post_process.cpp"
// #include "writecgns.cpp"

int main() { 
    //////////// Initilize Variables, Objects, Structs /////////////
    class_mesh mesh; 
    class_Q Qbar;
    class_Q Qface_c1;
    class_Q Qface_c2;
    class_flow freestream;
    class_residual residual;

    struct_size size;
    struct_inputs inputs; 
    struct_BC BC;
    struct_time time;
    struct_report report;
    
    ////////// Load Variables, Read Mesh, Set up Outputs ///////////
    read_inputs(&freestream, &inputs, &BC, &time, &report); // Load Input variables
    readmesh(&mesh,inputs.grid_file,&size, &BC, &inputs); // Import Mesh

    if (inputs.restart == 0) {
        std::ofstream res_history;
        res_history.open ("../sol/res_history.dat");
        if (res_history.is_open()) {
            res_history << "#iter, res1, res2, res3, res4" << std::endl;
        }
        res_history.close();

        std::ofstream output;
        output.open ("../sol/output.dat");
        if (output.is_open()) {
            output << "#iter, CL, CD" << std::endl;
        }
        output.close();
    }

    ////////////////////// Initilize Domain ////////////////////////
    int ndx;
    if (inputs.restart == 1) {
        std::ifstream restart;
        restart.open ("../sol/restart_sol.dat");
        if (restart.is_open()) {
            int num_cells_restart;
            double Q1, Q2, Q3, Q4;

            restart >> ndx >> num_cells_restart;
            if (num_cells_restart != size.num_cells) {
                std::cout << "ERROR: Restart File Does Not Match Grid Input";
            }

            Qbar.init(size.num_cells, freestream.gamma);
            for (int idx = 0; idx < num_cells_restart; idx ++) {
                restart >> Q1 >> Q2 >> Q3 >> Q4;
                Qbar.p1[idx] = Q1;
                Qbar.p2[idx] = Q2;
                Qbar.p3[idx] = Q3;
                Qbar.p4[idx] = Q4;
            }
        }
        restart.close();

        Qbar.updateflow();
        Qface_c1.init(size.num_faces, &freestream);
        Qface_c2.init(size.num_faces, &freestream);
    } else {
        Qbar.init(size.num_cells, &freestream);
        Qface_c1.init(size.num_faces, &freestream);
        Qface_c2.init(size.num_faces, &freestream);

        Qbar.updateQ();

        ndx = 1;
    }
    freestream.updateQ();

    residual.init(size.num_cells);

    //////////////////////// Run Simulation ////////////////////////
    while (ndx <= inputs.nmax) {
        switch (time.scheme) {
            case 1:
                // ssp_rk2(&mesh, &Qbar, &Qface_c1, &Qface_c2, &residual, &size, &inputs, &freestream, &BC, &time);
            case 2:
                ssp_rk3(&mesh, &Qbar, &Qface_c1, &Qface_c2, &residual, &size, &inputs, &freestream, &BC, &time);
        }

        // Monitor and Output Solution
        save(&Qbar, &mesh, &residual, &freestream, &inputs, &size, &report, ndx);

        ndx ++;
    }

    //////////////////////// Post Process ////////////////////////
    post(&Qbar, &residual, &inputs, &size, inputs.nmax);
    // write_cg(&mesh, &size);

    return 0;
}