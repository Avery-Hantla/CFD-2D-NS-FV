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
    
    ///////////////// Load Variables and Read Mesh /////////////////
    read_inputs(&freestream, &inputs, &BC); // Load Input variables
    readmesh(&mesh,inputs.grid_file,&size, &BC); // Import Mesh

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

        // Print Information to terminal 
        if ((ndx % inputs.monitor_step) == 0) {
            //double max1 = * std::max_element(residual.p1.begin(),residual.p1.end());
            std::cout << "iter: " << ndx << " res: " << *std::max_element(residual.p1.begin(),residual.p1.end()) << " ";
            std::cout << ", " << *std::max_element(residual.p2.begin(),residual.p2.end()) << " ";
            std::cout << ", " << *std::max_element(residual.p3.begin(),residual.p3.end()) << " ";
            std::cout << ", " << *std::max_element(residual.p4.begin(),residual.p4.end()) << std::endl;
        }
    }
    // std::ofstream myfile;
    // myfile.open ("centroids.txt");
    // for (int idx = 0; idx < size.num_cells; idx++) {
    //     myfile << mesh.cell_centerx[idx] << " " << mesh.cell_centery[idx] << std::endl;
    // }
    // myfile.close();

    return 0;
}