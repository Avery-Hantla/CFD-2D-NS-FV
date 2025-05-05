cmake_minimum_required(VERSION 3.22)
project(solverlib)

ADD_LIBRARY(solverlib STATIC res.cpp rusanov.cpp reconstruction.cpp squeeze.cpp ssp_rk3.cpp cell_gradient.cpp)

# set_target_properties(solverlib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# set_target_properties(solverlib PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# set(SOLVER_LIB ${CMAKE_CURRENT_BINARY_DIR}/solverlib.lib PARENT_SCOPE)