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
        showResults(_iStep); // Paraview
}

void FEM::solvePhaseField()
{
    assemblePhaseFieldProblem();
    solveSystemByPSORPetsc(matrixPF, rhsPF, solutionPF);
    updateFieldVariables(solutionPF);
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

    ierr = VecAssemblyBegin(rhsPF);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhsPF);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    return ierr;
}

PetscErrorCode FEM::solveSystemByPSORPetsc(Mat &A, Vec &b, Vec &x)
{
    /*
        PSOR (Projected Successive Over-Relaxation) algorithm;
        The PSOR method is a modification of the Gauss-Seidel method that allows for a relaxation factor;
        The value of Delta_d is calculated in this algorithm;
    */

    KSP ksp;
    PC pc;
    PetscInt its;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPRICHARDSON); // KSPRICHARDSON is offer used with SOR and PSOR methods
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    CHKERRQ(ierr);

    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = PCSetType(pc, PCSOR);
    CHKERRQ(ierr);
    ierr = PCSetFromOptions(pc);
    CHKERRQ(ierr);

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);

    // Calculated values must be > 0.
    for (PetscInt i = 0; i < numNodes; i++)
    {
        PetscScalar value;
        VecGetValues(x, 1, &i, &value);
        if (value < 0)
            value = 0;
        VecSetValues(x, 1, &i, &value, INSERT_VALUES);
    }

    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    ierr = KSPGetIterationNumber(ksp, &its);
    CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &res);
    CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);

    return ierr;
}

void FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    getPSORVecs(A, b);
}

double FEM::getPSORVecs(Mat &A, Vec &b)
{
    /*
        Here the PA, IR and JC arrays are used to store the matrix A in compressed column format;
        PA stores the non-zero values of the matrix A;
        IR stores the row indices of the non-zero values of the matrix A;
        JC stores the pointers to the entries of the IR array;
    */

    double JC[numNodes + 1] = {};
    



}