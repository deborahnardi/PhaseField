#include "../headers/FEM.h"

FEM::FEM() {}
FEM::FEM(const std::string _name)
    : name(_name)
{
    setResultsPath();
    deleteResults(true);
    createResultsPath();
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
};
FEM::~FEM() {}

/*----------------------------------------------------------------------------------
                            DATA INPUT METHODS
------------------------------------------------------------------------------------*/
void FEM::setLoadingVector(double ubar, int nSteps)
{
    double stepSize = ubar / nSteps;

    // Load ramp
    for (int i = 1; i <= nSteps; ++i)
        load.push_back(stepSize * i);

    // Load unloading
    for (int i = nSteps; i >= 1; --i)
        load.push_back(stepSize * i);

    // Negative load ramp
    for (int i = 1; i <= nSteps; ++i)
        load.push_back(-stepSize * i);

    // Negative load unloading
    for (int i = nSteps; i >= 1; --i)
        load.push_back(-stepSize * i);
}

void FEM::createResultsPath()
{
    if (rank == 0)
    {
        std::string command = "mkdir -p " + resultsPath + "results/hdf5";
        command += " && mkdir -p examples/";
        system(command.c_str());
    }
}

void FEM::deleteResults(bool deleteFiles)
{
    if (rank == 0)
    {
        if (deleteFiles)
        {
            std::cout << "Deleting files\n";
            std::string command = "rm -r " + resultsPath + "*";
            std::cout << command << "\n";
            system(command.c_str());
        }
    }
}

/*----------------------------------------------------------------------------------
                Assembling and solving problem PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode FEM::solveFEMProblem()
{
    int it = 0;
    double norm = 0.;
    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];

    norm = sqrt(norm);

    for (int iStep = 0; iStep < params->getNSteps(); iStep++)
    {
        double lambda = (1. + double(iStep)) / double(params->getNSteps());
        updateBoundaryValues(lambda);

        if (boundaryFunction) // 0 is false, any non zero value is true;
            updateBoundaryFunction(double(iStep) * params->getDeltaTime());

        it = 0;
        res = 1.;

        do
        {
            it++;
            assembleProblem();
            solveLinearSystem(matrix, rhs, solution);
            updateVariables(solution);
            res = res / norm;

        } while (res > params->getTol() && it < params->getMaxNewtonRaphsonIt());

        if (rank == 0)
            showResults(iStep);
    }

    cleanSolution(rhs, solution, matrix);
}

PetscErrorCode FEM::cleanSolution(Vec &x, Vec &b, Mat &A)
{
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);

    return ierr;
}

void FEM::updateBoundaryValues(double _lambda)
{
    for (auto bdElem : bdElements)
        bdElem->updateBoundaryValues(_lambda);
}

void FEM::updateBoundaryFunction(double _time)
{
    for (auto node : nodes)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
                boundaryFunction(node->getInitialCoordinates(), _time, dof, load);
}

PetscErrorCode FEM::assembleProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling FEM problem...\n");

    ierr = MatZeroEntries(matrix);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhs);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solution);
    CHKERRQ(ierr);

    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getContribution(matrix, rhs);

    for (int Ii = IIstart; Ii < IIend; Ii++) // Neumann boundary conditions
        bdElements[Ii]->getContribution(rhs);

    // Assemble the matrix and the right-hand side vector
    ierr = VecAssemblyBegin(rhs);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // ierr = VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);

    ierr = MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs); // Apply Dirichlet boundary conditions
    CHKERRQ(ierr);

    if (showMatrix && rank == 0) // Print the global stiffness matrix on the terminal
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, " --- GLOBAL STIFFNESS MATRIX: ----\n");
        CHKERRQ(ierr);
        ierr = MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRQ(ierr);
        printGlobalMatrix(matrix);
    }

    return ierr;
}

PetscErrorCode FEM::printGlobalMatrix(Mat &A)
{
    PetscInt i, j, rows, cols;
    PetscScalar value;
    const int width = 10; // Columns width

    ierr = MatGetSize(A, &rows, &cols);
    CHKERRQ(ierr);

    std::cout << "Global stiffness matrix:" << std::endl;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            ierr = MatGetValue(A, i, j, &value);
            CHKERRQ(ierr);
            // std::cout << std::setw(width) << std::fixed << std::setprecision(0) << value << " ";
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    return ierr;
}

PetscErrorCode FEM::solveLinearSystem(Mat &A, Vec &b, Vec &x)
{
    KSP ksp;
    PC pc;
    PetscInt its;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    switch (size)
    {
    case 1:
        ierr = KSPGetPC(ksp, &pc);
        CHKERRQ(ierr);
        ierr = PCSetType(pc, PCJACOBI);
        CHKERRQ(ierr);
        break;
    }

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &its); // Gets the number of iterations
    CHKERRQ(ierr);

    // ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // Prints the Krylov subspace method information
    // CHKERRQ(ierr);

    // ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    // ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); // Prints the solution vector
    // CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);

    return ierr;
}

void FEM::updateVariables(Vec &x)
{
    // Set the solution to the final coordinates of the nodes
    finalDisplacements = new double[globalDOFs.size()];
    int rankLocalDOFs = IIIend - IIIstart;
    double *localDisplacements;
    localDisplacements = new double[rankLocalDOFs];

    for (int i = IIIstart; i < IIIend; i++)
    {
        DOF *dof = globalDOFs[i];
        PetscScalar retrieved;
        PetscInt index = dof->getIndex();
        VecGetValues(x, 1, &index, &retrieved);

        localDisplacements[i - IIIstart] = retrieved;
    }

    // Transfering data between processes
    int *numEachProcess = new int[size];
    // Collecting the number of local DOFs from all processes
    MPI_Allgather(&rankLocalDOFs, 1, MPI_INT, numEachProcess, 1, MPI_INT, PETSC_COMM_WORLD);

    int *globalBuffer = new int[size];
    globalBuffer[0] = 0;
    for (int i = 1; i < size; i++)
        globalBuffer[i] = globalBuffer[i - 1] + numEachProcess[i - 1];

    // Collecting the displacements from all processes to the root process
    MPI_Gatherv(localDisplacements, rankLocalDOFs, MPI_DOUBLE, finalDisplacements, numEachProcess, globalBuffer, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    // Broadcasting the final displacements to all processes
    MPI_Bcast(finalDisplacements, nDOFs, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    res = 0.;
    for (auto node : nodes)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
            {
                double value = finalDisplacements[dof->getIndex()];
                dof->incrementValue(value);
                res += value * value;
            }

    // Print Update dof values
    for (auto node : nodes)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
                std::cout << dof->getValue() << std::endl;

    std::cout << "Values have been updated" << std::endl;
}
