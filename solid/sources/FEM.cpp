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
std::array<Tensor, 3> FEM::computeConstitutiveTensors()
{
    Tensor tensorK = {};
    Tensor tensorI = {};
    Tensor tensorJ = {};

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            tensorK[i][j] = 1.0;

    for (int i = 0; i < 3; i++)
        tensorI[i][i] = 1.0;

    tensorI[2][2] = 0.5;

    // J = I - 1/3 * K
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tensorJ[i][j] = tensorI[i][j] - 1.0 / 3.0 * tensorK[i][j];

    return {tensorK, tensorI, tensorJ};
}

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

    assembleSymmStiffMatrix(matrix);

    auto t2 = std::chrono::high_resolution_clock::now();
    // ====================== ASSEMBLING MATRIX AND VECTOR ======================

    PetscCall(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY));
    PetscCall(MatSetOption(matrix, MAT_SPD, PETSC_TRUE)); // symmetric positive-definite

    // ASSEMBLING THE RHS VECTOR: multiply column i of the stiffness matrix by the prescribed displacement of node i and set it to the right-hand side vector
    updateRHS(matrix, rhs);

    for (int Ii = IstartBD; Ii < IendBD; Ii++) // Neumann boundary conditions
        bdElements[Ii]->getContribution(rhs);

    PetscCall(VecAssemblyBegin(rhs));
    PetscCall(VecAssemblyEnd(rhs));

    PetscCall(VecView(rhs, PETSC_VIEWER_STDOUT_WORLD));
    exit(0);
    PetscCall(MatView(matrix, PETSC_VIEWER_STDOUT_WORLD));

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
    PetscCall(MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs)); // Apply Dirichlet boundary conditions

    auto t4 = std::chrono::high_resolution_clock::now();

    MPI_Barrier(PETSC_COMM_WORLD);
    exit(0);

    if (showMatrix && rank == 0) // Print the global stiffness matrix on the terminal
        printGlobalMatrix(matrix);

    PetscPrintf(PETSC_COMM_WORLD, "Assembling (s) %f %f %f\n", elapsedTime(t1, t2), elapsedTime(t2, t3), elapsedTime(t3, t4));

    return ierr;
}

PetscErrorCode FEM::assembleSymmStiffMatrix(Mat &A)
{
    // -------------------------------------------------------------------------------------------------

    std::array<Tensor, 3> tensors = computeConstitutiveTensors(); // tensor[0] = K, tensor[1] = I, tensor[2] = C;

    // std::vector<int> val(totalNnz, 0), idxRows(totalNnz, 0), idxCols(totalNnz, 0);
    // PetscScalar *val = new PetscScalar[totalNnz]();
    // PetscInt *idxRows = new PetscInt[totalNnz]();
    // PetscInt *idxCols = new PetscInt[totalNnz]();

    /*
        NOTE OF THE PRESENT DEVELOPER FOR THE FUTURE DEVELOPER:

        For some reason, when using BAIJ matrix, the code crashes when setting the values considering pointers for val, idxRows, and idxCols
        The code works fine when using the arrays directly.

        THIS WON'T WORK:
        ierr = MatSetValues(A, totalNnz, idxRows, totalNnz, idxCols, val, INSERT_VALUES);
        CHKERRQ(ierr);

        THIS WILL WORK:
        for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++)
        PetscCall(MatSetValue(A, idxRows[i], idxCols[i], val[i], INSERT_VALUES));

    */

    PetscScalar val[totalNnz] = {0};
    PetscInt idxRows[totalNnz] = {0};
    PetscInt idxCols[totalNnz] = {0};

    int kkn2n = 0, kk = 0;
    for (int iNode1 = 0; iNode1 < n2nCSRUpper.size() - 1; iNode1++)
    {
        const int n1 = nodesForEachRankCSR[rank] + iNode1;
        std::vector<int> friendNodes(n2nUpperMat[iNode1].begin(), n2nUpperMat[iNode1].end());
        const int numFriends = friendNodes.size();
        int iFriendCount = 0;

        for (auto n2 : friendNodes)
        {
            std::vector<int> elemInfo = eSameList[kkn2n];
            const int numLocalElems = elemInfo.size() / 3;

            /*  COMPUTE HERE THE Kglobal COMPONENTS ASSOCIATED TO THE INFLUENCE OF NODE n2 (COLUMN) ON NODE n1 (LINE)
                IF n1 == n2, ONLY 3 DIFFERENT COMPONENTS ARE COMPUTED
                CONSIDERING ONLY A 2D ANALYSIS, THOSE COMPONENTS MUST BE PLACED AT THE FOLLOWING POSITIONS ON VECTOR val:
             */

            if (n1 != n2)
            { // idof = 0,       jdof = 0
                int p1 = kk;
                // idof = 0,       jdof = 1
                int p2 = kk + 1;
                // idof = 1,       jdof = 0
                int p3 = kk + numFriends * nDOF - 1;
                // idof = 1,       jdof = 1
                int p4 = kk + numFriends * nDOF + 0;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];
                    std::vector<double> localStiffValue = elements[elemIndex]->getStiffnessIIOrIJ(tensors, idxLocalNode1, idxLocalNode2);

                    val[p1] += localStiffValue[0];
                    val[p2] += localStiffValue[1];
                    val[p3] += localStiffValue[2];
                    val[p4] += localStiffValue[3];
                }

                idxRows[p1] = nDOF * n1;     // First DOF
                idxRows[p2] = nDOF * n1;     // First DOF
                idxRows[p3] = nDOF * n1 + 1; // + 1 because it is the second DOF
                idxRows[p4] = nDOF * n1 + 1; // + 1 because it is the second DOF

                idxCols[p1] = nDOF * n2;
                idxCols[p2] = nDOF * n2 + 1;
                idxCols[p3] = nDOF * n2;
                idxCols[p4] = nDOF * n2 + 1;
            }
            else
            {
                // idof = 0,       jdof = 0
                int p1 = kk;
                // idof = 0,       jdof = 1
                int p2 = kk + 1;
                // idof = 1,       jdof = 1
                int p3 = kk + numFriends * nDOF;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];
                    std::vector<double> localStiffValue = elements[elemIndex]->getStiffnessIIOrIJ(tensors, idxLocalNode1, idxLocalNode2);

                    val[p1] += localStiffValue[0];
                    val[p2] += localStiffValue[1];
                    val[p3] += localStiffValue[2];
                }

                idxRows[p1] = nDOF * n1;
                idxRows[p2] = nDOF * n1;
                idxRows[p3] = nDOF * n1 + 1; // + 1 because it is the second DOF

                idxCols[p1] = nDOF * n2;
                idxCols[p2] = nDOF * n2 + 1;
                idxCols[p3] = nDOF * n2 + 1;

                // There is no p4 due to symmetry
            }

            iFriendCount++;
            kkn2n++;
            kk += nDOF;
        }

        kk += nDOF * numFriends - (nDOF * nDOF - (nDOF * (nDOF + 1)) / 2.0);
    }

    for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++)                // sizeof(idxRows) = size in bytes of the array; sizeof(idxRows[0]) = size in bytes of the first element of the array
        PetscCall(MatSetValue(A, idxRows[i], idxCols[i], val[i], INSERT_VALUES)); // THIS WORKS

    for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++) // Setting the lower triangular part of the matrix
        PetscCall(MatSetValue(A, idxCols[i], idxRows[i], val[i], INSERT_VALUES));

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
    case EMumps: // Parallel - (Multifrontal Massively Parallel Solver), used for large and sparse linear systems
        ierr = PCSetType(pc, PCLU);
        CHKERRQ(ierr);
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        CHKERRQ(ierr);
        break;
    case EIterative: // Parallel - faster than MUMPS, however it is harder to converge;
        // KSPFGMRES is a generalization of the GMRES algorithm that allows for flexible preconditioning (you can set the preconditioner type and other parameters)
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
    int rankLocalDOFs = Iend - Istart;
    double *localDisplacements = new double[rankLocalDOFs];

    for (int i = Istart; i < Iend; i++)
    {
        DOF *dof = globalDOFs[i];
        PetscScalar retrieved;
        PetscInt index = dof->getIndex();
        VecGetValues(x, 1, &index, &retrieved);

        localDisplacements[i - Istart] = retrieved;
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
    file.open(resultsPath + "force_displacement.txt", std::ios::app);

    if (params->getReactionDir() == "X")
        file << sumDisp << " " << force[0] << std::endl;
    else if (params->getReactionDir() == "Y")
        file << sumDisp << " " << force[1] << std::endl;

    file.close();

    return ierr;
}

PetscErrorCode FEM::updateRHS(Mat &A, Vec &b)
{
    PetscInt n = (Iend - Istart) * nDOF;

    PetscInt Jstart, Jend, Rstart, Rend;
    MatGetOwnershipRange(A, &Rstart, &Rend);
    MatGetOwnershipRangeColumn(A, &Jstart, &Jend);

    // PetscPrintf(PETSC_COMM_SELF, "Processo %d possui as linhas de %d a %d\n", rank, Rstart, Rend - 1);
    // PetscPrintf(PETSC_COMM_SELF, "Processo %d possui as colunas de %d a %d\n", rank, Jstart, Jend - 1);

    for (int line = Rstart; line < Rend; line++)
    {

        for (int j = 0; j < nDOFs; j++)
        {
            DOF *dof = globalDOFs[j];
            PetscInt jdx = dof->getIndex();
            PetscScalar dofValue = dof->getValue();
            PetscScalar value = 0.;

            PetscScalar product = 0.;
            PetscCall(MatGetValue(A, line, jdx, &value)); // MatGetValue access any value in the matrix, since the line idx is in the range of the process. PETSc allows each process to access any column (you can have access to any column of the matrix, but only to the lines that belong to the process)
            product = -value * dofValue;
            PetscCall(VecSetValues(b, 1, &line, &product, ADD_VALUES));
        }
    }

    return ierr;
}