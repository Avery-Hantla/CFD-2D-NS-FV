cmake_minimum_required(VERSION 3.22)
project(CUDA VERSION 0.1.0 LANGUAGES CXX CUDA)

ADD_LIBRARY(cuda_lib STATIC rusanov.cu)

set_target_properties(cuda_lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(cuda_lib PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# set(CUDA_LIB ${CMAKE_CURRENT_BINARY_DIR}/cuda.lib PARENT_SCOPE)