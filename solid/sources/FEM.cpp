#include "../headers/FEM.h"

FEM::FEM() {}
FEM::FEM(const std::string _name)
    : name(_name)
{
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    setResultsPath();
    deleteResults(true);
    createResultsPath();
};
FEM::~FEM() {}

double FEM::elapsedTime(std::chrono::_V2::system_clock::time_point t1, std::chrono::_V2::system_clock::time_point t2)
{
    std::chrono::duration<double> elapsed = t2 - t1;
    return elapsed.count();
}

/*----------------------------------------------------------------------------------
                            DATA INPUT METHODS
------------------------------------------------------------------------------------*/
void FEM::setLoadingVector(double ubar, int nSteps)
{
    // From 0 to ubar with x variating from 0 to 19 (20 steps)
    double step1 = ubar / 19; // 20 points -> 19 intervals
    for (int i = 0; i < 20; ++i)
        load.push_back(step1 * i);

    double step2 = 2 * ubar / 40;
    for (int i = 1; i <= 40; i++) // Índices 20 a 59
        load.push_back(ubar - step2 * i);

    double step3 = 2 * ubar / 20;
    for (int i = 1; i <= 20; i++)
        load.push_back(-ubar + step3 * i);

    for (int i = 0; i < load.size(); i++)
        std::cout << i << " " << load[i] << std::endl;

    std::cout << std::endl;
}

// void FEM::setLoadingVector(double ubar, int nSteps)
// {
//     double stepSize = ubar / nSteps;

//     // Load ramp
//     for (int i = 1; i <= nSteps; ++i)
//         load.push_back(stepSize * i);

//     // Load unloading
//     for (int i = nSteps; i >= 1; --i)
//         load.push_back(stepSize * i);

//     // Negative load ramp
//     for (int i = 1; i <= nSteps; ++i)
//         load.push_back(-stepSize * i);

//     // Negative load unloading
//     for (int i = nSteps; i >= 1; --i)
//         load.push_back(-stepSize * i);

//     std::cout << "Loading vector size: " << load.size() << std::endl;
// }

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
double FEM::computeNorm(const double *vec1, const double *vec2, const int &size)
{
    double norm = 0.0;
    for (int i = 0; i < size; i++)
        norm += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);

    return std::sqrt(norm);
}

PetscErrorCode FEM::solveFEMProblem()
{
    auto start_timer = std::chrono::high_resolution_clock::now();
    int it = 0;
    double norm = 0.;
    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];

    norm = sqrt(norm);

    for (int iStep = 0; iStep < params->getNSteps(); iStep++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\n================ Step %d ================\n", iStep);
        double lambda = (1. + double(iStep)) / double(params->getNSteps());
        updateBoundaryValues(lambda);

        if (boundaryFunction) // 0 is false, any non zero value is true;
            updateBoundaryFunction(double(iStep) * params->getDeltaTime());

        it = 0;
        res = 1.;

        do
        {
            it++;
            PetscPrintf(PETSC_COMM_WORLD, "\n------- Iteration %d -------\n", it);
            assembleProblem();
            solveLinearSystem(matrix, rhs, solution);
            updateVariables(solution);
            res = res / norm;
            PetscPrintf(PETSC_COMM_WORLD, "Residual: %e\n", res);
        } while (res > params->getTolNR() && it < params->getMaxNewtonRaphsonIt());

        if (rank == 0)
            showResults(iStep);
    }

    cleanSolution(rhs, solution, matrix);

    auto end_timer = std::chrono::high_resolution_clock::now();
    PetscPrintf(PETSC_COMM_WORLD, "Total elapsed time: %f\n", elapsedTime(start_timer, end_timer));
    return ierr;
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

    ierr = MatZeroEntries(matrix);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhs);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solution);
    CHKERRQ(ierr);

    // ====================== CALCULATING CONTRIBUTIONS ======================
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getContribution(matrix, rhs, negativeLoad);

    for (int Ii = IIstart; Ii < IIend; Ii++) // Neumann boundary conditions
        bdElements[Ii]->getContribution(rhs);
    auto t2 = std::chrono::high_resolution_clock::now();

    // ====================== ASSEMBLING MATRIX AND VECTOR ======================
    ierr = VecAssemblyBegin(rhs);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    auto t3 = std::chrono::high_resolution_clock::now();

    // ====================== APPLYING DIRICHLET BOUNDARY CONDITIONS ======================
    ierr = MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs); // Apply Dirichlet boundary conditions
    CHKERRQ(ierr);
    auto t4 = std::chrono::high_resolution_clock::now();

    if (showMatrix && rank == 0) // Print the global stiffness matrix on the terminal
        printGlobalMatrix(matrix);

    PetscPrintf(PETSC_COMM_WORLD, "Assembling (s) %f %f %f\n", elapsedTime(t1, t2), elapsedTime(t2, t3), elapsedTime(t3, t4));

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

    PetscReal norm1;
    PetscReal residual_norm = 0.;

    ierr = VecNorm(b, NORM_2, &norm1);
    CHKERRQ(ierr);
    // std::cout << "Norm of rhs DISPLACEMENT:" << norm1 << std::endl;
    // std::cout << "--------------------" << std::endl;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);

    switch (params->getSolverType())
    {
    case ESuiteSparse: // Sequential
        ierr = PCSetType(pc, PCLU);
        CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
        CHKERRQ(ierr);
        break;
    case EMumps: // Parallel
        ierr = PCSetType(pc, PCLU);
        CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        CHKERRQ(ierr);
        break;
    case EIterative: // Parallel - faster than MUMPS, however it is harder to converge;
        ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, params->getTolEIterative(), PETSC_DEFAULT, params->getMaxIterEIterative());
        CHKERRQ(ierr);
        ierr = KSPSetType(ksp, KSPFGMRES);
        CHKERRQ(ierr);
        ierr = PCSetType(pc, PCBJACOBI); // PCBJACOBI can be modified to other preconditioners
        CHKERRQ(ierr);
        break;
    }

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &its); // Gets the number of iterations
    CHKERRQ(ierr);

    if (params->getSolverType() == EIterative)
    {
        ierr = KSPGetResidualNorm(ksp, &residual_norm);
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "EIterative(%d) residual norm: %e\n", its, (double)residual_norm);
        CHKERRQ(ierr);
    }

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);

    return ierr;
}

void FEM::updateVariables(Vec &x, bool _hasConverged)
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

    for (auto node : nodes)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
            {
                double value = finalDisplacements[dof->getIndex()];
                dof->incrementValue(value);
            }

    // Erase the memory
    delete[] localDisplacements;
    delete[] numEachProcess;
    delete[] globalBuffer;
    delete[] finalDisplacements;

    // // Print Update dof values
    // std::cout << "--------------------------------------------" << std::endl;
    // for (auto node : nodes)
    //     for (auto dof : node->getDOFs())
    //         if (dof->getDOFType() != D)
    //             std::cout << dof->getValue() << std::endl;

    // // std::cout << "Displacement values have been updated" << std::endl;
    // std::cout << "--------------------------------------------" << std::endl;
}
