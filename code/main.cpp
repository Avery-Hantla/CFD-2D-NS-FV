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
#include "ssp_rk3.cpp"
#include "save_sol.cpp"

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
    
    ////////// Load Variables, Read Mesh, Set up Outputs ///////////
    read_inputs(&freestream, &inputs, &BC); // Load Input variables
    readmesh(&mesh,inputs.grid_file,&size, &BC); // Import Mesh

    std::ofstream output;
    output.open ("../sol/output.dat");
    if (output.is_open()) {
        output << "#iter, res1, res2, res3, res4" << std::endl;
    }

    ////////////////////// Initilize Domain ////////////////////////
    if (inputs.restart == 1) {
        Qbar.init(size.num_cells, &freestream);
        Qface_c1.init(size.num_faces, &freestream);
        Qface_c2.init(size.num_faces, &freestream);
    } else {
        Qbar.init(size.num_cells, &freestream);
        Qface_c1.init(size.num_faces, &freestream);
        Qface_c2.init(size.num_faces, &freestream);
    }

    Qbar.updateQ();
    freestream.updateQ();

    residual.init(size.num_cells);

    //////////////////////// Run Simulation ////////////////////////
    for (int ndx = 1; ndx <= inputs.nmax; ndx++) {
        ssp_rk3(&mesh, &Qbar, &Qface_c1, &Qface_c2, &residual, &size, &inputs, &freestream, &BC);

        // Monitor and Output Solution
        save(&Qbar, &residual, &inputs, &size, ndx);
    }

    return 0;
}