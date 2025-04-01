/////////////////////////////////////////////////////////
//         Function Monitor and Output Solution
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <filesystem> 

#include "class_q.hpp"
#include "struct_residual.hpp"
#include "struct_inputs.hpp"
#include "struct_size.hpp"

void save(class_Q* Qbar, struct_residual* residual, struct_inputs* inputs, struct_size* size, int ndx) {
    // Monitor 
    if ((ndx % inputs->monitor_step) == 0) {
        std::string iter = std::to_string(ndx);
        std::string res1 = std::to_string(*std::max_element(residual->p1.begin(),residual->p1.end()));
        std::string res2 = std::to_string(*std::max_element(residual->p2.begin(),residual->p2.end()));
        std::string res3 = std::to_string(*std::max_element(residual->p3.begin(),residual->p3.end()));
        std::string res4 = std::to_string(*std::max_element(residual->p4.begin(),residual->p4.end()));

        // Print to terminal
        std::cout << "iter: " << iter << " res: " << res1 << ", " + res2 << ", " << res3 << ", " << res4 << std::endl;

        // Save to output file
        std::ofstream output("../sol/output.dat", std::ios_base::app | std::ios_base::out);
        output << iter << " " << res1 << " " + res2 << " " << res3 << " " << res4 << std::endl;
    }

    // Output
    if ((ndx % inputs->output_step) == 0) {
        std::cout <<  "Saving Solution \n";

        std::ofstream output;
        output.open ("../sol/sol_" + std::to_string(ndx) + ".dat");
        if (output.is_open()) {
            for (int idx = 0; idx < size->num_cells; idx++) {
                output << Qbar->p1[idx] << " " << Qbar->p2[idx] << " " << Qbar->p3[idx] << " " << Qbar->p4[idx] << std::endl;
            }
        } else {
            std::cout << "ERROR: Cannot Save File \n";
        }
        output.close();
    }
}