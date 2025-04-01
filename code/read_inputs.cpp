#include <iostream>
#include <fstream>
#include "class_flow.hpp"
#include "struct_BC.hpp"
#include "struct_inputs.hpp"

void read_inputs(class_flow* freestream, struct_inputs* inputs, struct_BC* BC) {
    std::ifstream input_file;
        input_file.open ("../input.in");
        if (input_file.is_open()) {
            std::cout << "Reading input.in \n";
            std::string in_string1, in_string2, throw_away;
            for (int idx = 0; idx < 3; idx++) { // Look for inputs and flow
                input_file >> in_string1;
                if (in_string1 == "[inputs]") { // Input inputs secction
                    for (int jdx = 0; jdx < 3*6; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "CFL") ? inputs->CFL = stod(in_string2): false;
                        (in_string1 == "order") ? inputs->order = stoi(in_string2): false;
                        (in_string1 == "islimiteron") ? (in_string1 == "true") ? inputs->islimiteron = true: inputs->islimiteron = false : false;
                        (in_string1 == "nmax") ? inputs->nmax = stoi(in_string2): false;
                        (in_string1 == "output_step") ? inputs->output_step = stoi(in_string2): false;
                        (in_string1 == "grid_file") ? inputs->grid_file = in_string2: "false";
                    }
                    std::cout << "[inputs] \n";
                    std::cout << "  CFL = " << inputs->CFL << std::endl;
                    std::cout << "  order = " << inputs->order << std::endl;
                    std::cout << "  islimiteron = " << inputs->islimiteron << std::endl;
                    std::cout << "  nmax = " << inputs->nmax << std::endl;
                    std::cout << "  output_step = " << inputs->output_step << std::endl;
                    std::cout << "  grid_file = " << inputs->grid_file<< std::endl;
                }
                if (in_string1 == "[flow]") { // flow Flow section
                    for (int jdx = 0; jdx < 3*5; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "gamma") ? freestream->gamma = stod(in_string2): false;
                        (in_string1 == "P") ? freestream->P = stod(in_string2): false;
                        (in_string1 == "rho") ? freestream->rho = stod(in_string2): false;
                        (in_string1 == "u") ? freestream->u = stod(in_string2): false;
                        (in_string1 == "v") ? freestream->v = stod(in_string2): false;
                    }
                    std::cout << "[flow] \n";
                    std::cout << "  gamma = " << freestream->gamma << std::endl;
                    std::cout << "  P = " << freestream->P << std::endl;
                    std::cout << "  rho = " << freestream->rho << std::endl;
                    std::cout << "  u = " << freestream->u << std::endl;
                    std::cout << "  v = " << freestream->v << std::endl;
                }
                if (in_string1 == "[BC]") { // BC Flow section
                    for (int jdx = 0; jdx < 3*2; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "freestream_patch") ? BC->freestream_patch = stoi(in_string2)-1: false;
                        (in_string1 == "wall_patch") ? BC->wall_patch = stoi(in_string2)-1: false;
                    }
                    std::cout << "[BC] \n";
                    std::cout << "  freestream_patch = " << BC->freestream_patch+1 << std::endl;
                    std::cout << "  wall_patch = " << BC->wall_patch+1 << std::endl;
                }
            }
        } else {
            std::cout << "ERROR: Cannot Open Input File\n";
        }
    input_file.close();
}