cmake_minimum_required(VERSION 3.10.0)
# set(CMAKE_C_COMPILER "nvcc")
# set(CMAKE_CXX_COMPILER "nvcc")
project(fv2d VERSION 0.1.0 LANGUAGES C CXX)

# Find Header Files
find_path(CGNS_HEAD_DIR cgnslib.h)
find_path(EIGEN_HEAD_DIR Eigen)

# Add static libarys
# ADD_LIBRARY(_rusanov STATIC src/rusanov.cpp)

# Define preprocessor macro
option(USE_GPU OFF)
if(USE_GPU)
    add_subdirectory(src/cuda)
    add_definitions(-DUSE_GPU)
endif()

# Build executable
set(executable_output_path "../build/bin/")
add_executable(fv2d src/base/main.cpp)

# Link header files
target_include_directories(fv2d PUBLIC "${CGNS_HEAD_DIR}")
target_include_directories(fv2d PUBLIC "${EIGEN_HEAD_DIR}")

# Link static librarys
if(USE_GPU)
    TARGET_LINK_LIBRARIES(fv2d _cuda)
endif()

include(CTest)
enable_testing()