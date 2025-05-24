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
    for (int i = 0; i < 20; i++)
        load.push_back(ubar * (i + 1));

    for (int i = 20; i < 60; i++)
        load.push_back(ubar * (40 - (i + 1)));

    for (int i = 60; i < 80; i++)
        load.push_back(ubar * ((i + 1) - 80));

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
        matrixPreAllocationPF(nDOFPF);
        createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, nDOFPF, true);

        DdkMinus1 = new double[numNodes]{}; // Damage field at the previous iteration
        Ddk = new double[numNodes]{};       // Damage field at the current iteration
        totalVecq = new double[numNodes]{};

        updateFieldDistribution();
        updateFieldVariables(solutionPF);
        params->setCalculateReactionForces(true);

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
        double res = 1.e5;

        stepAUX = iStep;

        do
        {
            it++;
            PetscPrintf(PETSC_COMM_WORLD, "\n------- Iteration %d -------\n", it);
            assembleProblem();
            solveLinearSystem(matrix, rhs, solution); // right hand side is the residual
            updateVariables(matrix, solution, rhs, res);
            // res = res / norm;
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

PetscErrorCode FEM::performLineSearch(Mat &A, Vec &solution, Vec &rhs, Vec &copyRHS, std::vector<PetscScalar> &_backup, double &_res) // solution = delta_u, rhs = Fext - Fint = R(uk)
{
    const int maxBT = 10; // Maximum back-tracking
    double R0 = {}, R1 = {}, alpha = {}, eta = {}, resEta = {};

    // Compute R0 and R1: R0 = delta_u ^ T * R(uk) and R1 = delta_u ^ T * R(uk + delta_u)
    VecDot(solution, copyRHS, &R0); // R0 stands for the iteration at uk (before the line search)
    VecDot(solution, rhs, &R1);     // R1 stands for the iteration at uk + delta_u (during the line search)
    alpha = R0 / R1;

    // Compute the line search parameter eta
    eta = (alpha > 0
               ? 0.5 * alpha
               : 0.5 * alpha + std::sqrt(alpha * alpha / 4.0 - alpha));
    eta = std::clamp(eta, 1e-8, 1.0);

    for (int bt = 0; bt < maxBT; bt++)
    {
        // 1) Gets the current state of the displacement DOFs
        for (auto dof : globalDOFs)
            dof->setValue(_backup[dof->getIndex()]);

        Vec All;
        VecScatter ctx;
        // Gathers the solution vector to the master process
        VecScatterCreateToAll(solution, &ctx, &All);
        VecScatterBegin(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&ctx);

        // x_ls = currentVal + eta * deltaU
        for (auto dof : globalDOFs)
        {
            PetscInt Ii = dof->getIndex();
            PetscScalar currentVal = dof->getValue();
            PetscScalar deltaU = 0.0;
            VecGetValues(All, 1, &Ii, &deltaU);
            dof->setValue(currentVal + eta * deltaU); // Set the displacement DOFs to the new state
        }
        VecDestroy(&All);

        PetscCall(VecZeroEntries(rhs));
        updateRHS(matrix, rhs); // Update the residual with the new displacement values from x_ls

        for (int Ii = IstartBD; Ii < IendBD; Ii++) // Neumann boundary conditions
            bdElements[Ii]->getContribution(rhs);

        PetscCall(VecAssemblyBegin(rhs));
        PetscCall(VecAssemblyEnd(rhs));

        PetscCall(MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs)); // Apply Dirichlet boundary conditions

        // assembleProblem();

        // updateVariables(matrix, solution, rhs, resEta); // Compute the residual norm in case x_ls is the solution
        computeNorm(rhs, resEta); // Compute the residual norm in case x_ls is the solution
        // PetscPrintf(PETSC_COMM_WORLD, "Line search: %d  Residual: %e\n", bt, resEta);
        if (resEta <= 0.5 * _res) // The residual has decreased
        {
            _res = resEta;
            break;
        }

        eta *= 0.5; // Decrease the line search parameter
    }

    _res = resEta;
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

    // delete[] Ddk;

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

    PetscCall(MatZeroEntries(matrix));
    PetscCall(VecZeroEntries(rhs));
    PetscCall(VecZeroEntries(solution));

    // ====================== CALCULATING CONTRIBUTIONS ======================

    assembleSymmStiffMatrix(matrix);
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

    if (params->getCalculateReactionForces())
    {
        PetscCall(VecAssemblyBegin(nodalForces));
        PetscCall(VecAssemblyEnd(nodalForces));
        PetscCall(VecZeroEntries(nodalForces));

        PetscInt localStart = Istart * nDOF;
        PetscInt localEnd = Iend * nDOF;

        for (int i = localStart; i < localEnd; i++)
        {
            DOF *dof = globalDOFs[i];
            PetscScalar retrieved;
            PetscInt index = dof->getIndex();
            PetscCall(VecGetValues(rhs, 1, &index, &retrieved));
            PetscCall(VecSetValues(nodalForces, 1, &index, &retrieved, INSERT_VALUES));
        }
    }

    // ====================== APPLYING DIRICHLET BOUNDARY CONDITIONS ======================
    PetscCall(MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs)); // Apply Dirichlet boundary conditions

    if (showMatrix && rank == 0) // Print the global stiffness matrix on the terminal
        printGlobalMatrix(matrix);

    return ierr;
}

PetscErrorCode FEM::assembleSymmStiffMatrix(Mat &A)
{
    std::array<Tensor, 3> tensors = computeConstitutiveTensors(); // tensor[0] = K, tensor[1] = I, tensor[2] = C;
    int kk = 0;
    SplitModel splitModel = params->getSplitModel();

    // #pragma omp parallel
    //     {
    int thread_id = omp_get_thread_num();
    int num_threads = omp_get_num_threads();

    // #pragma omp for
    for (int iNode1 = 0; iNode1 < n2nCSRUpper.size() - 1; iNode1++)
    {
        const int n1 = nodesForEachRankCSR[rank] + iNode1;
        std::vector<int> friendNodes(n2nUpperMat[iNode1].begin(), n2nUpperMat[iNode1].end());
        const int numFriends = friendNodes.size();
        PetscScalar localVal[nDOF * nDOF * numFriends - 1]{};  // Node 0 has 0, 1, 4 and 7 as friends. 0-0: 3 positions, 0-1: 4 positions, 0-4: 4 positions, 0-7: 4 positions. Total: 15 positions.
        PetscInt idxLocalRows[nDOF * nDOF * numFriends - 1]{}; // Local values belong to the current thread
        PetscInt idxLocalCols[nDOF * nDOF * numFriends - 1]{};

        int kkn2n = 0, kkLocal = 0, kkStart = 0;
        kkStart = n2nCSRUpper[iNode1] * nDOF * nDOF - (iNode1);

        for (auto n2 : friendNodes)
        {
            // std::vector<int> elemInfo = eSameList[kkn2n];
            std::vector<int> elemInfo = eSameListForEachNode[iNode1][kkn2n];
            const int numLocalElems = elemInfo.size() / 3;

            /*  COMPUTE HERE THE Kglobal COMPONENTS ASSOCIATED TO THE INFLUENCE OF NODE n2 (COLUMN) ON NODE n1 (LINE)
                IF n1 == n2, ONLY 3 DIFFERENT COMPONENTS ARE COMPUTED
                CONSIDERING ONLY A 2D ANALYSIS, THOSE COMPONENTS MUST BE PLACED AT THE FOLLOWING POSITIONS ON VECTOR val:
             */

            if (n1 != n2)
            { // idof = 0,       jdof = 0
                int p1 = kkLocal;
                // idof = 0,       jdof = 1
                int p2 = kkLocal + 1;
                // idof = 1,       jdof = 0
                int p3 = kkLocal + numFriends * nDOF - 1;
                // idof = 1,       jdof = 1
                int p4 = kkLocal + numFriends * nDOF + 0;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];

                    std::vector<double> localStiffValue = elements[elemIndex]->getStiffnessIIOrIJ(tensors, idxLocalNode1, idxLocalNode2, splitModel, stepGlobal);

                    localVal[p1] += localStiffValue[0];
                    localVal[p2] += localStiffValue[1];
                    localVal[p3] += localStiffValue[2];
                    localVal[p4] += localStiffValue[3];
                }

                idxLocalRows[p1] = nDOF * n1;     // First DOF
                idxLocalRows[p2] = nDOF * n1;     // First DOF
                idxLocalRows[p3] = nDOF * n1 + 1; // + 1 because it is the second DOF
                idxLocalRows[p4] = nDOF * n1 + 1; // + 1 because it is the second DOF

                idxLocalCols[p1] = nDOF * n2;
                idxLocalCols[p2] = nDOF * n2 + 1;
                idxLocalCols[p3] = nDOF * n2;
                idxLocalCols[p4] = nDOF * n2 + 1;
            }
            else
            {
                // idof = 0,       jdof = 0
                int p1 = kkLocal;
                // idof = 0,       jdof = 1
                int p2 = kkLocal + 1;
                // idof = 1,       jdof = 1
                int p3 = kkLocal + numFriends * nDOF;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];
                    std::vector<double> localStiffValue = elements[elemIndex]->getStiffnessIIOrIJ(tensors, idxLocalNode1, idxLocalNode2, splitModel, stepGlobal);

                    localVal[p1] += localStiffValue[0];
                    localVal[p2] += localStiffValue[1];
                    localVal[p3] += localStiffValue[2];
                }

                idxLocalRows[p1] = nDOF * n1;
                idxLocalRows[p2] = nDOF * n1;
                idxLocalRows[p3] = nDOF * n1 + 1; // + 1 because it is the second DOF

                idxLocalCols[p1] = nDOF * n2;
                idxLocalCols[p2] = nDOF * n2 + 1;
                idxLocalCols[p3] = nDOF * n2 + 1;

                // There is no p4 due to symmetry
            }
            kkn2n++;
            kkLocal += nDOF;
        }

        // #pragma omp critical
        for (int i = 0; i < nDOF * nDOF * numFriends - 1; i++)
        {
            MatSetValue(A, idxLocalRows[i], idxLocalCols[i], localVal[i], INSERT_VALUES);
            MatSetValue(A, idxLocalCols[i], idxLocalRows[i], localVal[i], INSERT_VALUES);
        }
    }
    //} // comment this when OMP is not used

    return 0;
}

PetscErrorCode FEM::updateRHS(Mat &A, Vec &b)
{
    Vec KTimesU, U;
    PetscCall(VecDuplicate(b, &KTimesU));
    PetscCall(VecDuplicate(b, &U));
    PetscCall(VecZeroEntries(KTimesU));
    PetscCall(VecZeroEntries(U));

    for (int j = 0; j < nDOFs; j++)
    {
        DOF *dof = globalDOFs[j];
        PetscInt idx = dof->getIndex();
        PetscScalar dofValue = dof->getValue();
        PetscCall(VecSetValues(U, 1, &idx, &dofValue, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(KTimesU));
    PetscCall(VecAssemblyEnd(KTimesU));
    PetscCall(VecAssemblyBegin(U));
    PetscCall(VecAssemblyEnd(U));

    PetscCall(MatMult(A, U, KTimesU));
    PetscCall(VecAXPY(b, -1.0, KTimesU)); // rhs = rhs - KTimesU

    PetscCall(VecDestroy(&KTimesU));
    PetscCall(VecDestroy(&U));

    return ierr;
}

PetscErrorCode FEM::solveLinearSystem(Mat &A, Vec &b, Vec &x)
{
    // KSP ksp;
    // PC pc;
    PetscInt its = 0;
    PetscReal residual_norm = 0.;

    // PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    // PetscCall(KSPSetReusePreconditioner(ksp, PETSC_FALSE));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
    PetscCall(KSPGetPC(ksp, &pc));

    PetscCall(VecSet(x, 0.0)); // Initialize the solution vector to zero

    switch (params->getSolverType())
    {
    case ESuiteSparse: // Sequential (it is a direct solver)
        PetscCall(PCSetType(pc, PCLU));
        PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK));
        break;
    case EMumps: // Parallel - (Multifrontal Massively Parallel Solver), used for large and sparse linear systems (MUMPS is a direct solver)
        PetscCall(PCSetType(pc, PCLU));
        PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
        break;
    case EIterative: // Parallel - faster than MUMPS, however it is harder to converge;
        // KSPFGMRES is a generalization of the GMRES algorithm that allows for flexible preconditioning (you can set the preconditioner type and other parameters)
        PetscCall(KSPSetTolerances(ksp, PETSC_DEFAULT, params->getTolEIterative(), PETSC_DEFAULT, params->getMaxIterEIterative()));
        PetscCall(KSPSetType(ksp, KSPFGMRES));
        PetscCall(PCSetType(pc, PCBJACOBI)); // PCBJACOBI can be modified to other preconditioners
        break;
    case ECholesky:
        PetscCall(KSPSetType(ksp, KSPPREONLY));
        PetscCall(PCSetType(pc, PCCHOLESKY));
        PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
        break;
    case EPardiso: // Parallel Direct Solver
        if (size == 1)
        {
            PetscCall(KSPSetType(ksp, KSPPREONLY));
            PetscCall(PCSetType(pc, PCLU));
            PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMKL_PARDISO)); // MATSOLVERMKL_CPARDISO for sequential
        }
        else
        {
            /* NOT WORKING YET */
            PetscCall(KSPSetType(ksp, KSPPREONLY));
            PetscCall(PCSetType(pc, PCLU));
            PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMKL_CPARDISO)); // MATSOLVERMKL_CPARDISO for parallel
        }
        break;
    }

    PetscCall(PCFactorSetReuseOrdering(pc, PETSC_TRUE));
    PetscCall(PCFactorSetReuseFill(pc, PETSC_TRUE));

    PetscCall(KSPSolve(ksp, b, x));
    PetscCall(KSPGetIterationNumber(ksp, &its)); // Gets the number of iterations

    if (params->getSolverType() == EIterative)
    {
        PetscCall(KSPGetResidualNorm(ksp, &residual_norm));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "EIterative(%d) residual norm: %e\n", its, (double)residual_norm));
    }

    // PetscCall(KSPDestroy(&ksp));

    return ierr;
}

PetscErrorCode FEM::updateVariables(Mat A, Vec &x, Vec &b, double &_res, bool _hasConverged)
{
    Vec All, AllForces;
    VecScatter ctx;

    _res = 0.;

    // Gathers the solution vector to the master process
    PetscCall(VecScatterCreateToAll(x, &ctx, &All));
    PetscCall(VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&ctx));

    PetscCall(VecScatterCreateToAll(b, &ctx, &AllForces));
    PetscCall(VecScatterBegin(ctx, b, AllForces, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, b, AllForces, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&ctx));

    for (auto dof : globalDOFs)
    {
        PetscInt Ii = dof->getIndex();
        PetscScalar val;
        PetscScalar valForces;
        PetscCall(VecGetValues(AllForces, 1, &Ii, &valForces));
        PetscCall(VecGetValues(All, 1, &Ii, &val));
        dof->incrementValue(val);
        _res += valForces * valForces;
    }
    //_res = sqrt(_res);
    PetscCall(VecDestroy(&All));
    PetscCall(VecDestroy(&AllForces));

    return ierr;
}

PetscErrorCode FEM::computeNorm(Vec &b, double &_res)
{
    Vec All;
    VecScatter ctx;

    PetscCall(VecScatterCreateToAll(b, &ctx, &All));
    PetscCall(VecScatterBegin(ctx, b, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, b, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&ctx));
    PetscCall(VecDot(All, All, &_res));

    PetscCall(VecDestroy(&All));
    return ierr;
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