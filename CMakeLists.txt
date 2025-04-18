cmake_minimum_required(VERSION 3.10.0)
project(fv2d VERSION 0.1.0 LANGUAGES C CXX)

# Append PREFIX_PATH
# list(APPEND CMAKE_SYSTEM_INCLUDE_PATH "/home/a568h695/work/846/libraries/include")
# set(CMAKE_PREFIX_PATH "/home/a568h695/work/846/libraries/include" FORCE)
# string(REPLACE ":" ";" CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")
# list(APPEND CMAKE_PREFIX_PATH "/home/a568h695/work/846/libraries/include")

# Find Header Files
find_path(CGNS_HEAD_DIR cgnslib.h)
find_path(EIGEN_HEAD_DIR Eigen)

add_executable(fv2d code/main.cpp)

# Link header files
target_include_directories(fv2d PUBLIC "${CGNS_HEAD_DIR}")
target_include_directories(fv2d PUBLIC "${EIGEN_HEAD_DIR}")

include(CTest)
enable_testing()