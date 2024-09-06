static char help[] = "Solves the phase-field problem using Finete Element Method "
                     "Start date: July 22nd, 2024 "
                     "Written by: Deborah Nardi"
                     "Advisor: Prof. Edson Denner Leonel"
                     "Co-advisor: Prof. Ayrton Ferreira"
                     "São Carlos School of Engineering - University of São Paulo";

#include "mesh_interface/headers/Geometry.h"
#include "solid/headers/FEM.h"
#include "solid/headers/Quadrature.h"
#include "solid/headers/ShapeFunction.h"
#include "solid/headers/DenseEigen.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help); // Starts main program invoking PETSc

// #include "examples/Ex01Inclusions.hpp"
#include "examples/square.hpp"
    //   #include "examples/Ex02NumericalIntegration.hpp"
    // #include "examples/Ex03Truss.hpp"

    PetscFinalize(); // Finalize main program
    return 0;
}