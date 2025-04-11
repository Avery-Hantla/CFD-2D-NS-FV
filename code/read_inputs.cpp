#include <iostream>
#include <fstream>
#include "class_flow.hpp"
#include "struct_BC.hpp"
#include "struct_inputs.hpp"
#include "struct_time.hpp"
#include "struct_report.hpp"

void read_inputs(class_flow* freestream, struct_inputs* inputs, struct_BC* BC, struct_time* time, struct_report* report, int argc, char** argv) {
    std::ifstream input_file;
        if (argc == 2) {
            input_file.open (argv[1]);
        } else {
            input_file.open ("input.in");
        }
        if (input_file.is_open()) {
            std::cout << "\nReading input.in \n";
            std::string in_string1, in_string2, throw_away;
            for (int idx = 0; idx < 8; idx++) { // Look for inputs and flow
                input_file >> in_string1;
                /////////////////////////////////// Case ///////////////////////////////////
                if (in_string1 == "[case]") { // Input inputs secction
                    for (int jdx = 0; jdx < 3*5; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "solution_order") ? inputs->order = stoi(in_string2): false;
                        (in_string1 == "nmax") ? inputs->nmax = stoi(in_string2): false;
                        (in_string1 == "monitor_step") ? inputs->monitor_step = stoi(in_string2): false;
                        (in_string1 == "output_step") ? inputs->output_step = stoi(in_string2): false;
                        (in_string1 == "restart") ? inputs->restart = stoi(in_string2): false;
                    }
                    std::cout << "[case] \n";
                    std::cout << "  nmax = " << inputs->nmax << std::endl;
                    std::cout << "  monitor_step = " << inputs->monitor_step << std::endl;
                    std::cout << "  output_step = " << inputs->output_step << std::endl;
                    std::cout << "  solution_order = " << inputs->order << std::endl;
                    std::cout << "  restart = " << inputs->restart << std::endl;
                }
                /////////////////////////////////// Limiter ///////////////////////////////////
                if (in_string1 == "[limiter]") { // Input inputs secction
                    std::cout << "[limiter] \n";
                    for (int jdx = 0; jdx < 3*1; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        if (in_string1 == "limiter") {
                            if (in_string2 == "NONE") {
                                std::cout << "  limiter = NONE \n";
                                inputs->limiter = 0;
                            }
                            if (in_string2 == "MINMOD") {
                                std::cout << "  limiter = MINMOD \n";
                                inputs->limiter = 1;
                            }
                            if (in_string2 == "SQUEEZE") {
                                std::cout << "  limiter = SQUEEZE \n";
                                inputs->limiter = 2;
                            }
                        }
                    }
                }
                /////////////////////////////////// Grid ///////////////////////////////////
                if (in_string1 == "[grid]") { // Input inputs secction
                    for (int jdx = 0; jdx < 3*1; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "grid_file") ? inputs->grid_file = in_string2: "false";
                    }
                    std::cout << "[grid] \n";
                    std::cout << "  grid_file = " << inputs->grid_file<< std::endl;
                }
                /////////////////////////////////// Equations ///////////////////////////////////
                if (in_string1 == "[equations]") { // Input inputs secction
                    std::cout << "[equations] \n";
                    for (int jdx = 0; jdx < 3*2; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        if (in_string1 == "eqn") {
                            if (in_string2 == "EULER") {
                                std::cout << "  eqn = EULER \n";
                                inputs->eqn = 1;
                            }
                            if (in_string2 == "NS") {
                                std::cout << "  eqn = NS \n";
                                inputs->eqn = 2;
                            }
                        }
                        if (in_string1 == "flux_solver") {
                            if (in_string2 == "RUSANOV") {
                                std::cout << "  flux_solver = RUSANOV \n";
                                inputs->flux_solver = 1;
                            }
                            if (in_string2 == "ROE") {
                                std::cout << "  flux_solver = ROE \n";
                                inputs->flux_solver = 2;
                            }
                        }
                    }
                }
                /////////////////////////////////// Flow ///////////////////////////////////
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
                /////////////////////////////////// Time ///////////////////////////////////
                if (in_string1 == "[time]") { // Input inputs secction
                    std::cout << "[time] \n";
                    for (int jdx = 0; jdx < 3*2; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "CFL") ? time->CFL = stod(in_string2): false;
                        if (in_string1 == "time_scheme") {
                            if (in_string2 == "SSP_RK2") {
                                std::cout << "  scheme = SSP_RK2 \n";
                                time->scheme = 1;
                            }
                            if (in_string2 == "SSP_RK3") {
                                std::cout << "  scheme = SSP_RK3 \n";
                                time->scheme = 2;
                            }
                        }

                    }
                    std::cout << "  CFL = " << time->CFL << std::endl;
                }
                /////////////////////////////////// report ///////////////////////////////////
                if (in_string1 == "[report]") { // flow Flow section
                    for (int jdx = 0; jdx < 3*6; jdx+=3) { // Less than 3*(num of inputs)
                        input_file >> in_string1 >> throw_away >> in_string2;
                        (in_string1 == "u") ? report->u = stod(in_string2): false;
                        (in_string1 == "v") ? report->v = stod(in_string2): false;
                        (in_string1 == "P") ? report->P = stod(in_string2): false;
                        (in_string1 == "rho") ? report->rho = stod(in_string2): false;
                        (in_string1 == "area") ? report->area = stod(in_string2): false;
                        (in_string1 == "length") ? report->length = stod(in_string2): false;
                    }
                    std::cout << "[report] \n";
                    std::cout << "  u = " << report->u << std::endl;
                    std::cout << "  v = " << report->v << std::endl;
                    std::cout << "  P = " << report->P << std::endl;
                    std::cout << "  rho = " << report->rho << std::endl;
                    std::cout << "  area = " << report->area << std::endl;
                    std::cout << "  length = " << report->length << std::endl;
                }
                /////////////////////////////////// BC ///////////////////////////////////
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