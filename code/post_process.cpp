/////////////////////////////////////////////////////////
//                 Function Post Process
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <filesystem> 

#include "class_q.hpp"
#include "class_residual.hpp"
#include "struct_inputs.hpp"
#include "struct_size.hpp"

void post(class_Q* Qbar, class_residual* residual, struct_inputs* inputs, struct_size* size, int ndx) {
    // Save Restart File
    std::cout <<  "Writing Restart File \n";

    std::ofstream restart;
    restart.open ("sol/restart_sol.dat");
    if (restart.is_open()) {
        restart << ndx << " " << size->num_cells <<std::endl;
        for (int idx = 0; idx < size->num_cells; idx++) {
            restart << Qbar->p1[idx] << " " << Qbar->p2[idx] << " " << Qbar->p3[idx] << " " << Qbar->p4[idx] << std::endl;
        }
    } else {
        std::cout << "ERROR: Cannot Write Restart File \n";
    }
    restart.close();

    // Ou
}