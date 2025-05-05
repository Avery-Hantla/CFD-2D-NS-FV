cmake_minimum_required(VERSION 3.22)
project(CGNSLIB VERSION 0.1.0 LANGUAGES CXX CUDA)


include_directories(${CMAKE_SOURCE_DIR}/contrib/cgnslib/include)


# find_path(CGNS_LIBRARY
#     NAMES libcgns.a
#     HINTS $ENV{CGNS_DIR}/lib ${CMAKE_SOURCE_DIR}/contrib/cgnslib/lib/libcgns.a
#     REQUIRED
# )
# if (CGNS_LIBRARY)
#     cmake_path(GET CGNS_LIBRARY PARENT_PATH CGNS_LIB_DIR)
#     cmake_path(GET CGNS_LIB_DIR PARENT_PATH CGNS_DIR)
#     add_library(cgns::cgns STATIC IMPORTED)
#     set_target_properties(cgns::cgns PROPERTIES
#     INTERFACE_INCLUDE_DIRECTORIES "${CGNS_DIR}/include"
#     INTERFACE_LINK_LIBRARIES ${CGNS_LIBRARY}
#     IMPORTED_LOCATION ${CGNS_LIBRARY})
#     set(CGNS_FOUND TRUE)
# endif()

# if (CGNS_FOUND)
# set(EXTRA_LIBS ${EXTRA_LIBS} cgns::cgns)
# else()
# set(EXTRA_LIBS ${EXTRA_LIBS} libcgns.a)
# endif()

ADD_LIBRARY(readCGNS_func STATIC readMeshBase.cpp readStructuredMesh.cpp readUnstructuredMesh.cpp)

# add_library( libcgns SHARED IMPORTED GLOBAL)
# set_target_properties( libcgns PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/contrib/cgnslib/lib/libcgns.a )

# add_library( libhdf5 SHARED IMPORTED GLOBAL)
# set_target_properties( libhdf5 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/contrib/hdf5/lib/libhdf5.a )

# find_library(CGNS_LIBRARY libcgns.a ${CMAKE_SOURCE_DIR}/contrib/cgnslib/lib/)
add_library(readCGNS INTERFACE)
target_link_libraries(readCGNS INTERFACE readCGNS_func libcgns libhdf5 libzlib)

# target_link_libraries(readCGNS PRIVATE CGNS_LIBRARY)