#pragma once

#include <petscsnes.h>
#include <petscksp.h> // KSP -> Krylov Subspace Methods: Deals with Krylov subspace methods

class PETScExs
{
public:
    PETScExs();
    ~PETScExs();

    PetscErrorCode PETScSequentialTest();
};