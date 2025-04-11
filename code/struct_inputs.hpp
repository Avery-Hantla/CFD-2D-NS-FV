#ifndef STRUCT_INPUTS_HPP
    #define STRUCT_INPUTS_HPP
    #include <iostream>
    struct struct_inputs{
        int order; 
        int nmax; 
        int monitor_step;
        int output_step;
        int restart;

        // Grid
        std::string grid_file; 

        // Equation
        int eqn;
        int flux_solver;

        // Limiter
        int limiter; // 0 is no limiter, 1 is minmod, 2 is squeeze, 
    };
#endif