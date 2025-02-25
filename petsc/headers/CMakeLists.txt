# @file FEM/headers/CMakeLists.txt  -- CMake file for the Headers module

file(GLOB petsc_headers *.h)

set(PETSC_HEADERS ${petsc_headers} PARENT_SCOPE)