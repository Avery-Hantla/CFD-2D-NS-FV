cmake_minimum_required(VERSION 3.22)
project(fv2d VERSION 0.1.0 LANGUAGES CXX CUDA)

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
# find_path(CGNS_HEAD_DIR cgnslib.h)

# Build executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/bin/)
add_executable(fv2d src/base/main.cpp)

# Link config
target_include_directories(fv2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config/cgns/include/)
target_include_directories(fv2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config/eigen/)

# Link static librarys
if(USE_GPU)
    target_link_libraries(fv2d cuda_lib ${CUDART_LIBRARY})
else()
    target_link_libraries(fv2d solverlib)
endif()