/* ////////////////////////////////////////////////////////// 
              FV 2D Euler Code for Mixed Meshes
                        Avery Hantla
                         March 2025
*/ //////////////////////////////////////////////////////////
#include <iostream>

// Load classes, structs, etc.
#include "class_mesh.hpp"
#include "class_q.hpp"
#include "class_f.hpp"
#include "class_flow.hpp"

#include "struct_BC.hpp"
#include "struct_size.hpp"
#include "struct_inputs.hpp"

// Load Functions 
#include "read_inputs.cpp"
#include "readmesh.cpp"
#include "ssp_rk2.cpp"
#include "save_sol.cpp"

int main() { 
    //////////// Initilize Variables, Objects, Structs /////////////
    class_mesh mesh; 
    class_Q Qbar;
    class_Q Qface_c1;
    class_Q Qface_c2;
    class_F Fface;
    class_flow freestream;

    struct_size size;
    struct_inputs inputs; 
    struct_BC BC;
    
    ////////// Load Variables, Read Mesh, Set up Outputs ///////////
    read_inputs(&freestream, &inputs, &BC); // Load Input variables
    readmesh(&mesh,inputs.grid_file,&size, &BC); // Import Mesh

    std::ofstream output;
    output.open ("../sol/output.dat");
    if (output.is_open()) {
        output << "#iter, res1, res2, res3, res4" << std::endl;
    }

    ////////////////////// Initilize Domain ////////////////////////
    Qbar.init(size.num_cells, &freestream);
    Qbar.updateQ();
    Qface_c1.init(size.num_cells*6);
    Qface_c2.init(size.num_cells*6);

    freestream.updateQ();

    struct_residual residual;
    residual.p1.assign(size.num_cells, -101);
    residual.p2.assign(size.num_cells, -101);
    residual.p3.assign(size.num_cells, -101);
    residual.p4.assign(size.num_cells, -101);

    //////////////////////// Run Simulation ////////////////////////
    for (int ndx = 1; ndx <= inputs.nmax; ndx++) {
        ssp_rk2(&mesh, &Qbar, &Qface_c1, &Qface_c2, &residual, &size, &inputs, &freestream, &BC);

        // Monitor and Output Solution
        save(&Qbar, &residual, &inputs, &size, ndx);
    }

    return 0;
}