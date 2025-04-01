#ifndef STRUCT_INPUTS_HPP
    #define STRUCT_INPUTS_HPP
    #include <iostream>
    struct struct_inputs{
        double CFL;
        int order; 
        bool islimiteron;
        int nmax; 
        int output_step;
        std::string grid_file; 
    };
#endif