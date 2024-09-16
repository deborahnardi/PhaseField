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

#include <petscdmplex.h>  // DM -> Data Management: Deals with mesh and problems related to domain discretization
#include <petscsnes.h>    // SNES -> Scalable Nonlinear Equations Solvers
#include <petscds.h>      // DS -> Data Structures: Deals with data structures and data management
#include <petscbag.h>     // BAG -> Basic Algebraic Graph: Deals with graph theory and algebraic graph theory
#include <petscconvest.h> // CONVEST -> Convergence Estimation: Deals with convergence estimation

int main(int argc, char **argv)
{

    PetscInitialize(&argc, &argv, (char *)0, help); // Starts main program invoking PETSc

    // DISCOMMENT THE ABOVE LINE BEFORE RUNNING THE NON PETSC EXAMPLES
    // #include "examples/pointerAndReference.hpp"
    // #include "examples/Ex01Inclusions.hpp"
    // #include "examples/square.hpp"
    // #include "examples/squareEllipse.hpp"
    //      #include "examples/Ex02NumericalIntegration.hpp"
    //  #include "examples/Ex03Truss.hpp"
#include "examples/Ex04Truss.hpp"

    // =======================================

    // PETSC Examples
    // PetscFunctionBeginUser; // States the beginning of a user-defined function/program
    // PetscCall(PetscInitialize(&argc, &argv, NULL, help));

    // #include "examples/Petsc.cpp"

    PetscFinalize(); // Finalize main program
    return 0;
}