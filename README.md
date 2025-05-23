# CFD 2D Finite Volume Euler and Navier-Stokes Code with Arbitrary Mesh 

This code is capable of handling an arbitrary mesh with cells containing n number of sides. There is both first and second-order reconstruction implemented along with the squeeze limiter. 

Note that some features such as CGNS, GPU, and ROE Flux are still under development  

## Requirements  
 - CMake  

## Installation  
To install and build the program, clone the repo then make a build directory and use CMake to build  

```bash
git clone https://github.com/Avery-Hantla/CFD-2D-NS-FV.git
cd CFD-2D-NS-FV
mkdir build
cd build
cmake ../
make
```

The executable will be located in   

CFD-2D-NS-FV/build/bin  

## Inputs
An input file can be directly passed to the binary and can be ran using.

```bash
./fv2d input.in
```

If no input file is specified, then the binary will search for an "input.in" file in the working directory. A sample input file is provided in the repository, an input file is also shown below with all working options.  

[case]  
    nmax = Number of iterations to run   
    monitor_step = Number of iterations between each job output monitor  
    output_step = Number of iterations between each solution file dump  
    solution_poly_order = 0, 1  
    restart = 0, 1  

[grid]  
    grid_file = PATH_TO_GRID_FILE  

[equations]   
    eqn = EULER, NS  
    flux_solver = RUSANOV  

[fluid]  
    gamma = Ratio of Specific Heats  
    R = Gas Constant  
    mu = Dynamic Viscosity   
    Pr = Prandtl Number   

[init]  
    P = Initial Pressure  
    rho = Initial Density  
    u = Initial X velocity  
    v = Initial Y Velocity  

[time]  
    time_scheme = SSP_RK3  
    time_step_type = GLOBAL, LOCAL  
    CFL = CFL_NUM  

[limiter]  
    limiter = NONE, SQUEEZE  

[report]  
    length = REF_LENGTH  

[BC]  
    BC1 = SLIP_WALL, NO_SLIP_WALL, FREESTREAM, EXTRAPOLATION  
    BC2 = SLIP_WALL, NO_SLIP_WALL, FREESTREAM, EXTRAPOLATION  
