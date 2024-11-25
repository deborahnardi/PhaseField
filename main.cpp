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
#include "petsc/headers/PETScExs.h"
#include "solid/headers/DenseEigen.h"

int main(int argc, char **argv)
{

    PetscInitialize(&argc, &argv, (char *)0, help); // Starts main program invoking PETSc

    // DISCOMMENT THE ABOVE LINE BEFORE RUNNING THE NON PETSC EXAMPLES
    // #include "examples/pointerAndReference.hpp"
    // #include "examples/Ex01Inclusions.hpp"
    // #include "examples/square.hpp"
    //       #include "examples/squareEllipse.hpp"
    //       #include "examples/Ex02NumericalIntegration.hpp"
    //          #include "examples/Ex03Truss.hpp"
    //   #include "examples/Ex04Truss.hpp"

    // =======================================

    // PETSC Examples
    // PetscFunctionBeginUser; // States the beginning of a user-defined function/program
    // PetscErrorCode ierr; // PETSc error code
    // ierr = PetscInitialize(&argc, &argv, (char *)0, help);
    // CHKERRQ(ierr);

    // PETScExs *p = new PETScExs();
    //  p->PETScSequentialTest();
    // p->PETScParallelTest();

    // ================ Phase Field Examples =================
    // #include "examples/phaseField1D.hpp"
#include "examples/phaseField2D.hpp"
    PetscFinalize(); // Finalize main program
    // return ierr;
}