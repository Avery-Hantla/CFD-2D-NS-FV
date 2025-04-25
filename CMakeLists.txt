cmake_minimum_required(VERSION 3.10.0)
project(fv2d VERSION 0.1.0 LANGUAGES C CXX)

# Find Header Files
find_path(CGNS_HEAD_DIR cgnslib.h)
find_path(EIGEN_HEAD_DIR Eigen)

# Define preprocessor macro
option(USE_GPU OFF)
if(USE_GPU)
    add_definitions(-DUSE_GPU)

    # Build executable
    add_executable(fv2d_gpu src/main.cpp)

    # Link header files
    target_include_directories(fv2d_gpu PUBLIC "${CGNS_HEAD_DIR}")
    target_include_directories(fv2d_gpu PUBLIC "${EIGEN_HEAD_DIR}")
else()
    # Build executable
    add_executable(fv2d src/main.cpp)

    # Link header files
    target_include_directories(fv2d PUBLIC "${CGNS_HEAD_DIR}")
    target_include_directories(fv2d PUBLIC "${EIGEN_HEAD_DIR}")
endif()

include(CTest)
enable_testing()