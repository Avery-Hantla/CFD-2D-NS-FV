cmake_minimum_required(VERSION 3.22)
project(fv2d)

# Add subdirectorys
add_subdirectory(src/cgns-read/src)
# find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})



# Contrib
include_directories(${CMAKE_SOURCE_DIR}/contrib/cgnslib/include)
include_directories(${CMAKE_SOURCE_DIR}/contrib/eigen)
include_directories(${CMAKE_SOURCE_DIR}/contrib/hdf5/include)


# Define preprocessor macro
option(USE_GPU OFF)
if(USE_GPU) 
    # find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    add_subdirectory(src/cuda)
    add_definitions(-DUSE_GPU)
else()
    add_subdirectory(src/solver)
    # find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()



# Find librarys
add_library( libcgns SHARED IMPORTED GLOBAL)
set_target_properties( libcgns PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/contrib/cgnslib/lib/libcgns.a )

add_library( libhdf5 SHARED IMPORTED GLOBAL)
set_target_properties( libhdf5 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/contrib/hdf5/lib/libhdf5.a )

add_library( libzlib SHARED IMPORTED GLOBAL)
set_target_properties( libzlib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/contrib/zlib/lib/libz.a )

# find_path(CGNS_LIBRARY_DIR
#     NAMES libcgns.a
#     HINTS $ENV{CGNS_DIR}/lib ${CMAKE_SOURCE_DIR}/contrib/cgnslib/lib/
#     REQUIRED
# )
# if (CGNS_LIBRARY)
# cmake_path(GET CGNS_LIBRARY PARENT_PATH CGNS_LIB_DIR)
# cmake_path(GET CGNS_LIB_DIR PARENT_PATH CGNS_DIR)
# add_library(cgns::cgns STATIC IMPORTED)
# set_target_properties(cgns::cgns PROPERTIES
# INTERFACE_INCLUDE_DIRECTORIES "${CGNS_DIR}/include"
# INTERFACE_LINK_LIBRARIES ${CGNS_LIBRARY}
# IMPORTED_LOCATION ${CGNS_LIBRARY})
# # set(CGNS_FOUND TRUE)
# endif()

# set(HDF5_USE_STATIC_LIBRARIES true)

# if (DEFINED ENV{HDF5_ROOT})
#     set (HDF5_ROOT $ENV{HDF5_ROOT})
#     # we want to bypass the builtin FindHDF5.cmake as it may not match our external one. Instead rely on CMake installation
#     list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakefile/modules/")
#     find_package(HDF5 COMPONENTS C CXX REQUIRED)
#     message ("Found HDF5 = ${HDF5_FOUND} , ${HDF5_CXX_LIBRARIES} , ${HDF5_C_LIBRARIES}, ${HDF5_LIBRARIES}, ${HDF5_VERSION}, ${HDF5_INCLUDE_DIRS}")
#     # add the modules back after we've found cmake
#     list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmakefile/modules/")
# endif()

# if (NOT DEFINED HDF5_ROOT)
#     set (HDF5_ROOT ${CMAKE_SOURCE_DIR}/contrib/hdf5)
#     set (HDF5_C_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/contrib/hdf5/include)
#     set (HDF5_CXX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/contrib/hdf5/include/cpp)
#     find_package(HDF5 COMPONENTS C CXX REQUIRED)
# endif()

# set(ZLIB_USE_STATIC_LIBS true)
# find_package(ZLIB)



# Build executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/bin/)
add_executable(fv2d src/base/main.cpp)



# Link config
# target_include_directories(fv2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config/cgnslib/include/)
# target_include_directories(fv2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config/eigen/)



# Link static librarys
if(USE_GPU)
    target_link_libraries(fv2d cuda_lib ${CUDART_LIBRARY})
else()
    target_link_libraries(fv2d solverlib)
endif()
# target_link_libraries(fv2d libcgns libhdf5)
target_link_libraries(fv2d readCGNS)