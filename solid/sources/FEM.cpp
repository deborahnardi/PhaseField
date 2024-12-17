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

PetscErrorCode FEM::printMemoryUsage(const int &iStep)
{
    ierr = PetscMemoryGetCurrentUsage(&bytes);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes / (1024 * 1024));
    CHKERRQ(ierr);

    if (rank == 0)
    {
        std::ofstream file(resultsPath + "memoryUsage.txt", std::ios::app);
        file << iStep << " " << bytes / (1024 * 1024) << std::endl;
        file.close();
    }

    return ierr;
}

/*----------------------------------------------------------------------------------
                            DATA INPUT METHODS
------------------------------------------------------------------------------------*/
void FEM::setLoadingVector2(double ubar, int nSteps)
{
    double ubar1 = ubar;
    double ubar2 = -ubar; //; / 10;
    double ubar3 = ubar;  // / 10;

    double totalLoad = 0;

    // From 0 to ubar with x variating from 0 to 19 (20 steps)
    for (int i = 0; i < 20; ++i)
        load.push_back(ubar * i);

    for (int i = 20; i < 60; ++i)
        load.push_back(ubar * (40 - i));

    for (int i = 60; i < 80; ++i)
        load.push_back(ubar * (i - 80));

    // double step2 = (ubar1 - ubar2) / 39;
    // for (int i = 1; i <= 40; ++i)
    //     load.push_back(ubar1 - step2 * i);

    // double step3 = (ubar3 - ubar2) / 19;
    // for (int i = 1; i <= 20; ++i)
    //     load.push_back(ubar2 + step3 * i);

    // for (int i = 0; i < load.size(); i++)
    //     std::cout << i << " " << load[i] << " " << std::endl;

    std::cout << std::endl;
}

void FEM::setLoadingVector1(double ubar, int nSteps)
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

void FEM::setLoadingVector3(double ubar, int nSteps)
{
    load.clear();

    double delta1 = 1e-3;
    double delta2 = 3e-4;
    double delta3 = -3e-4;

    double currentVal = 0.;

    for (int i = 0; i < 6; i++)
    {
        currentVal += delta1;
        load.push_back(currentVal);
    }

    for (int i = 6; i < 26; i++)
    {
        currentVal += delta2;
        load.push_back(currentVal);
    }

    for (int i = 26; i < 106; i++)
    {
        currentVal += delta3;
        load.push_back(currentVal);
    }

    for (int i = 106; i < 146; i++)
    {
        currentVal += delta2;
        load.push_back(currentVal);
    }

    // for (size_t i = 0; i < load.size(); ++i)
    //     std::cout << i + 1 << " " << load[i] << std::endl;

    std::cout << std::endl;
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
    for (auto n : discritizedNodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];

    norm = sqrt(norm);

    if (prescribedDamageField)
    {
        params->setCalculateReactionForces(false);
        matrixPreAllocationPF(IstartPF, IendPF);
        createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, true);

        DdkMinus1 = new double[numNodes]{}; // Damage field at the previous iteration
        Ddk = new double[numNodes]{};       // Damage field at the current iteration
        totalMatrixQ = new double *[numNodes] {};
        totalVecq = new double[numNodes]{};

        for (int i = 0; i < numNodes; i++)
            totalMatrixQ[i] = new double[numNodes]{};

        updateFieldDistribution();
        updateFieldVariables(solutionPF);
        params->setCalculateReactionForces(true);

        for (int i = 0; i < numNodes; i++)
            delete[] totalMatrixQ[i];
        delete[] totalMatrixQ;
        delete[] totalVecq;
    }

    for (int iStep = 0; iStep < params->getNSteps(); iStep++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\n================ Step %d ================\n", iStep);
        printMemoryUsage(iStep);

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

        postProc();

        if (rank == 0)
        {
            if (params->getCalculateReactionForces())
                computeReactionForces();
            showResults(iStep);
        }
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

    delete[] Ddk;

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

    if (params->getCalculateReactionForces())
    {
        ierr = VecZeroEntries(nodalForces);
        CHKERRQ(ierr);
    }

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

    if (params->getCalculateReactionForces())
    {
        ierr = VecAssemblyBegin(nodalForces);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(nodalForces);
        CHKERRQ(ierr);
    }

    if (params->getCalculateReactionForces())
    {
        for (auto dof : globalDOFs)
            if (dof->getDOFType() != D)
            {
                PetscScalar value = 0.;
                PetscInt idx = dof->getIndex();
                ierr = VecGetValues(rhs, 1, &idx, &value);
                CHKERRQ(ierr);
                ierr = VecSetValues(nodalForces, 1, &idx, &value, INSERT_VALUES);
                CHKERRQ(ierr);
            }
    }

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
    PetscReal residual_norm = 0.;

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
    double *finalDisplacements = new double[globalDOFs.size()];
    int rankLocalDOFs = IIIend - IIIstart;
    double *localDisplacements = new double[rankLocalDOFs];

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
}

PetscErrorCode FEM::computeReactionForces()
{
    int count = 0;
    double sumDisp = 0.;
    // double sumForces = 0.;

    for (auto dof : globalDOFs)
        if (dof->isControlledDOF())
            if (dof->getDOFType() != D)
            {
                // Print in a txt file force vs displacement
                // PetscInt idx = dof->getIndex();
                // PetscScalar value;
                // ierr = VecGetValues(nodalForces, 1, &idx, &value);
                // CHKERRQ(ierr);
                sumDisp += dof->getValue();
                // sumForces -= value;
                count++;
                // dof->setReactionForce(-value);
            }

    if (count != 0)
        sumDisp /= count;

    std::ofstream file;
    file.open(resultsPath + "results/force_displacement.txt", std::ios::app);
    file << sumDisp << " " << force[0] << std::endl;
    file.close();

    return ierr;
}
