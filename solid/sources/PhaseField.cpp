#include "../headers/FEM.h"
void FEM::solvePhaseFieldProblem()
{
    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];
    norm = sqrt(norm);

    createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, true);

    for (int iStep = 0; iStep < params->getNSteps(); iStep++)
    {
        staggeredAlgorithm(iStep);
    }
}

void FEM::staggeredAlgorithm(int _iStep)
{
    // Staggered algorithm
    // 1. Solve the displacement problem
    // 2. Solve the phase field problem
    // 3. Repeat until convergence

    int it = 0;
    double resStag = 0.0;
    do
    {
        it++;
        solveDisplacementField(_iStep);
        solvePhaseField();

    } while (resStag > params->getResidStaggered() && it < params->getMaxItStaggered());
}

void FEM::solveDisplacementField(int _iStep)
{
    double lambda = (1. + double(_iStep)) / double(params->getNSteps());
    updateBoundaryValues(lambda);

    if (boundaryFunction) // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime());

    int itNR = 0;
    res = 1.;

    do
    {
        itNR++;
        assembleProblem();
        solveLinearSystem(matrix, rhs, solution);
        updateVariables(solution);
        res = res / norm;

    } while (res > params->getTol() && itNR < params->getMaxNewtonRaphsonIt());

    if (rank == 0)
        showResults(_iStep);
}

void FEM::solvePhaseField()
{
    assemblePhaseFieldProblem();
    solveLinearSystem(matrixPF, rhsPF, solutionPF);
}

PetscErrorCode FEM::assemblePhaseFieldProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling Phase Field problem...\n");

    ierr = MatZeroEntries(matrixPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhsPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solutionPF);
    CHKERRQ(ierr);

    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getPhaseFieldContribution(matrixPF, rhsPF);

    return ierr;
}
