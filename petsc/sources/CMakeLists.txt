# @file FEM/mesh_interface/sources/CMakeLists.txt  -- CMake file for the Sources module

file(GLOB petsc_src *.cpp)

set(PETSC_SCR ${petsc_src} PARENT_SCOPE)