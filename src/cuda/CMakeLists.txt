cmake_minimum_required(VERSION 3.10.0)
# set(CMAKE_C_COMPILER "nvcc")
# set(CMAKE_CXX_COMPILER "nvcc")
project(fv2d VERSION 0.1.0 LANGUAGES CXX CUDA)

# if(NOT DEFINED CMAKE_CUDA_STANDARD)
#     set(CMAKE_CUDA_STANDARD 11)
#     set(CMAKE_CUDA_STANDARD_REQUIRED ON)
# endif()

# Add static libarys
ADD_LIBRARY(_cuda STATIC rusanov.cu)

include(CTest)
enable_testing()