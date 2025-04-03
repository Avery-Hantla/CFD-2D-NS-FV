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
                    for (int jdx = 0; jdx < 3*8; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "CFL") ? inputs->CFL = stod(in_string2): false;
                        (in_string1 == "order") ? inputs->order = stoi(in_string2): false;
                        (in_string1 == "islimiteron") ? (in_string1 == "true") ? inputs->islimiteron = true: inputs->islimiteron = false : false;
                        (in_string1 == "nmax") ? inputs->nmax = stoi(in_string2): false;
                        (in_string1 == "monitor_step") ? inputs->monitor_step = stoi(in_string2): false;
                        (in_string1 == "output_step") ? inputs->output_step = stoi(in_string2): false;
                        (in_string1 == "grid_file") ? inputs->grid_file = in_string2: "false";
                        (in_string1 == "restart") ? inputs->restart = stoi(in_string2): false;
                    }
                    std::cout << "[inputs] \n";
                    std::cout << "  CFL = " << inputs->CFL << std::endl;
                    std::cout << "  restart = " << inputs->restart << std::endl;
                    std::cout << "  order = " << inputs->order << std::endl;
                    std::cout << "  islimiteron = " << inputs->islimiteron << std::endl;
                    std::cout << "  nmax = " << inputs->nmax << std::endl;
                    std::cout << "  monitor_step = " << inputs->monitor_step << std::endl;
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
                    std::cout << "[BC] \n";
                    for (int jdx = 0; jdx < 3*2; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        if (in_string1 == "BC1") {
                            if (in_string2 == "FREESTREAM") {
                                std::cout << "  BC1 = FREESTREAM \n";
                                BC->BC1 = -1;
                            }
                            if (in_string2 == "SLIP_WALL") {
                                std::cout << "  BC1 = SLIP_WALL \n";
                                BC->BC1 = -2;
                            }
                            if (in_string2 == "EXTRAPOLATION") {
                                std::cout << "  BC1 = EXTRAPOLATION \n";
                                BC->BC1 = -3;
                            }
                        }
                        if (in_string1 == "BC2") {
                            if (in_string2 == "FREESTREAM") {
                                std::cout << "  BC2 = FREESTREAM \n";
                                BC->BC2 = -1;
                            }
                            if (in_string2 == "SLIP_WALL") {
                                std::cout << "  BC2 = SLIP_WALL \n";
                                BC->BC2 = -2;
                            }
                            if (in_string2 == "EXTRAPOLATION") {
                                std::cout << "  BC2 = EXTRAPOLATION \n";
                                BC->BC2 = -3;
                            }
                        }
                    }
                }
            }
        } else {
            std::cout << "ERROR: Cannot Open Input File\n";
        }
    input_file.close();
}