cmake_minimum_required(VERSION 3.22)
project(fv2d VERSION 0.1.0 LANGUAGES CXX CUDA)

# Add subdirectory
# add_subdirectory(src/cuda)

# Link static libarys
# find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})

# Define preprocessor macro
option(USE_GPU OFF)
if(USE_GPU)
    find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    add_subdirectory(src/cuda)
    add_definitions(-DUSE_GPU)
else()
    add_subdirectory(src/solver)
    # find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()

# Find Header Files
find_path(CGNS_HEAD_DIR cgnslib.h)
find_path(EIGEN_HEAD_DIR Eigen)

# Build executable
# set(executable_output_path "../build/bin/")
add_executable(fv2d src/base/main.cpp)

# Link header files
target_include_directories(fv2d PUBLIC "${CGNS_HEAD_DIR}")
target_include_directories(fv2d PUBLIC "${EIGEN_HEAD_DIR}")

# Link static librarys
if(USE_GPU)
    target_link_libraries(fv2d cuda_lib ${CUDART_LIBRARY})
else()
    target_link_libraries(fv2d solverlib)
endif()

include(CTest)
enable_testing()