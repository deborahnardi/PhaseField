cmake_minimum_required (VERSION 3.11.0)
project (s) # Name of the project here

list (APPEND CMAKE_MODULE_PATH "${s_SOURCE_DIR}/cmake-modules")

find_package (PETSc REQUIRED)   
find_package (HDF5 REQUIRED)
find_package (METIS REQUIRED)
find_package (MPI REQUIRED)

include_directories(include ${MPI_INCLUDE_PATH} ${PETSC_INCLUDES})

# Essential: include our directories first otherwise we can get internal headers from some installed path
include_directories ("${s_SOURCE_DIR}/include" "${s_SOURCE_DIR}/include" ${PETSC_INCLUDES} ${PETSC_KSP_LIB})
add_definitions(${PETSC_DEFINITIONS})
set (s_DEPENDENT_LIBRARIES "${PETSC_LIBRARIES}")
set (s_LIBRARIES s ${s_DEPENDENT_LIBRARIES})
include_directories (${HDF5_INCLUDE_DIRS} ${HDF5_DEFINITIONS} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
include_directories(${MPI_INCLUDE_PATH})

#hide warnings
add_definitions(-w) 

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#[[
If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
to Debug prior to calling PROJECT(). This assumes you're only build one type of solution
a time. If you're working with multiple solutions like Visual Studio Professional can do,
you need to comment this out.
]]
if(DEFINED CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: Debug, Release.")
else() # CMAKE_BUILD_TYPE
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug, Release.")
endif() # CMAKE_BUILD_TYPE
 
set(CMAKE_CXX_STANDARD 17)

# Eigen
include(cmake-modules/EnableEigen.cmake)
enable_eigen()

add_subdirectory(mesh_interface/headers)
add_subdirectory(mesh_interface/sources)
add_subdirectory(solid/headers)
add_subdirectory(solid/sources)

#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/mesh_interface/headers"
	"${PROJECT_SOURCE_DIR}/mesh_interface/sources"
   "${PROJECT_SOURCE_DIR}/solid/headers"
   "${PROJECT_SOURCE_DIR}/solid/sources"
)

add_library(solid_library ${HEADERS} ${SOURCES} ${MESH_INT} ${MESH_INT_SRC})

set(GMSH_PATH "NONE" CACHE STRING "Gmsh build directory.")
if(${GMSH_PATH} STREQUAL "NONE")
   message(FATAL_ERROR "Please, include Gmsh build directory -- CMake will exit!")
else()
    target_include_directories(solid_library PUBLIC             "${GMSH_PATH}/api/"
                                                    "${GMSH_PATH}/src/common/"
    )

   find_library(GMSH_LIB NAMES libgmsh.so.4.14 PATHS ${GMSH_PATH}/build/)
endif()

target_link_libraries(solid_library PRIVATE Eigen3::Eigen)

link_directories(${PETSC_DIR}/${PETSC_ARCH}/lib)

add_executable(${PROJECT_NAME} main.cpp ${CXX_SOURCE_FILES} ${GMSH_SOURCE_FILES} ${ELEMENT_SOURCE_FILES} ${BOUNDARY_CONDITION_SOURCE_FILES} ${IGA_SOURCE_FILES} ${QUADRATURE_SOURCE_FILES})

target_link_libraries (s metis solid_library Eigen3::Eigen ${GMSH_LIB} ${PETSC_KSP_LIB} ${PETSC_LIBRARIES} ${HDF5_LIBRARIES} ${METIS_LIBRARIES} ${MPI_LIBRARIES})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")