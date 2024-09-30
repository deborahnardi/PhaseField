#include "../headers/FEM.h"

void FEM::showResults()
{
    PetscPrintf(PETSC_COMM_WORLD, "Exporting data to Paraview...\n");

    std::stringstream text;
    text << "results/" << name << ".vtu";
    std::ofstream file(text.str()); // write the results to a file

}