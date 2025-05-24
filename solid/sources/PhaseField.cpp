#include "../headers/FEM.h"

void FEM::matrixPreAllocationPF(int nDOF)
{
    // TOTAL NUMBER OF NONZERO (UPPER TRIANGULAR PART ONLY) ONLY FOR THE LOCAL PARTITION

    for (int i = 0; i < numNodesForEachRank[rank]; i++)
    {
        std::set<int> neighbours = n2nUpperMat[i];
        int numFriends = neighbours.size();
        totalNnzPF += nDOF * nDOF * numFriends - (nDOF * nDOF - (nDOF * (nDOF + 1)) / 2.0);
    }

    // PREALLOCATION OF THE MPI MATRIX
    int m = numNodesForEachRank[rank] * nDOF; // Number of local lines at the rank. Ex.: 3 nodes at the local partition * 2 DOFs = 6 lines
    int nn = numNodes * nDOF;                 // Number of local columns at the rank. Ex.: 8 nodes * 2 DOFs = 16 columns

    // _nnz: only the nonzero terms of the upper triangular part of the matrix
    d_nnz = new PetscInt[m]();
    o_nnz = new PetscInt[m]();

    // _nz: the nonzero terms of the complete matrix
    d_nz = new PetscInt[m]();
    o_nz = new PetscInt[m]();

    for (int ii = 0; ii < numNodesForEachRank[rank]; ii++)
        for (int iDOF = 0; iDOF < nDOF; iDOF++)
        {
            d_nnz[nDOF * ii + iDOF] = nDOF * n2nUpperMat[ii].size() - n2nDRankUpper[ii] * nDOF - iDOF; // (n2nDRank[ii] * nDOF - iDOF) removes the DOFs that do not belong to the local partition
            o_nnz[nDOF * ii + iDOF] = n2nDRankUpper[ii] * nDOF;

            d_nz[nDOF * ii + iDOF] = (n2nMat[ii].size() - n2nDRankUpper[ii] - n2nDRankLower[ii]) * nDOF;
            o_nz[nDOF * ii + iDOF] = (n2nDRankUpper[ii] + n2nDRankLower[ii]) * nDOF;
        }

    // Total number of non-zero elements in the matrix, stores it in nzQ
    int localNzQ = 0;
    nzQ = 0;

    for (int i = 0; i < numNodesForEachRank[rank]; i++)
        localNzQ += d_nz[i] + o_nz[i];

    MPI_Reduce(&localNzQ, &nzQ, 1, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
    MPI_Allreduce(&localNzQ, &nzQ, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
}
// --------------------------------------------------------------------------------------------------
void FEM::solvePhaseFieldProblem() // Called by the main program
{
    auto start_timer = std::chrono::high_resolution_clock::now();

    for (auto n : discritizedNodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];
    norm = sqrt(norm);

    matrixPreAllocationPF(nDOFPF);

    params->setCalculateReactionForces(false);
    createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, nDOFPF, true);
    params->setCalculateReactionForces(true);

    // Ddk = new double[numNodes]{}; // Damage field at the current iteration

    if (prescribedDamageField)
    {
        Ddk = new double[numNodes]{}; // Damage field at the current iteration
        updateFieldDistribution();
        delete[] Ddk;
    }

    prescribedDamageField = false;

    int nSteps = params->getNSteps();

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "================ BEGINNING PHASE FIELD PROBLEM ================\n");

    for (int iStep = 0; iStep < nSteps; iStep++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\n================ STEP %d ================\n", iStep);
        printMemoryUsage(iStep);
        Ddk = new double[numNodes]{};
        staggeredAlgorithm(iStep); // Returns converged uStag and dStag
        updateFieldVariables(solutionPF);
        delete[] Ddk;
        stepGlobal++;
        postProc();

        if (rank == 0)
        {
            if (params->getCalculateReactionForces())
                computeReactionForces();
            showResults(iStep);
        }
    }

    cleanSolution(rhs, solution, matrix);
    cleanSolution(rhsPF, solutionPF, matrixPF);

    KSPDestroy(&ksp);
    MatDestroy(&seqMatQ);
    VecDestroy(&nodalForces);

    auto end_timer = std::chrono::high_resolution_clock::now();
    PetscPrintf(PETSC_COMM_WORLD, "Total elapsed time: %f\n", elapsedTime(start_timer, end_timer));
    PetscPrintf(PETSC_COMM_WORLD, "================ END OF PHASE FIELD PROBLEM ================\n");
}

void FEM::staggeredAlgorithm(int _iStep)
{
    int it = 0;
    double resStag = 0.0;
    //  double previousU[globalDOFs.size()]{}; // Previous displacement field
    double maxTol = params->getTolStaggered();
    // Ddk = new double[numNodes]{}; // Damage field at the current iteration
    // KSPCreate(PETSC_COMM_WORLD, &ksp);
    do
    {
        it++;
        PetscPrintf(PETSC_COMM_WORLD, "\n------- Stag Iteration %d -------\n", it);
        //================================================================================================
        //                                 STAGE 01 - SOLVE DISPLACEMENT FIELD
        //================================================================================================
        // PetscLogStage stageNum;
        // char stageName[100];
        // snprintf(stageName, sizeof(stageName), "SolveEqProb: %d", globalCounter++, "_", it);

        // PetscLogStageRegister(stageName, &stageNum);
        // PetscLogStagePush(stageNum);

        solveDisplacementField(_iStep); // Obtains ustag

        //  PetscLogStagePop();
        //================================================================================================
        //                                 STAGE 02 - SOLVE PHASE FIELD
        //================================================================================================
        //    snprintf(stageName, sizeof(stageName), "SolvePFProb: %d", globalCounter, "_", it);

        //      PetscLogStageRegister(stageName, &stageNum);
        //        PetscLogStagePush(stageNum);

        solvePhaseField(); // Obtains dstag

        // PetscLogStagePop();
        //================================================================================================

        evalStaggeredRes(resStag); // Computes the residual norm of the phase field

        // double normU = 0.0;
        // for (int i = 0; i < globalDOFs.size(); i++)
        // {
        //     DOF *dof = globalDOFs[i];
        //     if (dof->getDOFType() != D)
        //     {
        //         double value = dof->getValue();
        //         normU += (value - previousU[i]) * (value - previousU[i]);
        //         previousU[i] = value;
        //     }
        // }
        // resStag = sqrt(normU);

        PetscPrintf(PETSC_COMM_WORLD, "Residual stag: %e\n", resStag);

    } while (resStag > params->getTolStaggered() && it < params->getMaxItStaggered());
    // KSPDestroy(&ksp);
}

void FEM::solveDisplacementField(int _iStep)
{
    double lambda = (1. + double(_iStep)) / double(params->getNSteps());
    updateBoundaryValues(lambda);

    double entry = double(_iStep) * params->getDeltaTime();
    if (boundaryFunction)                                                // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime()); //
    // double aux = 0.0;
    //  assembleProblem(0);
    //  solveLinearSystem(matrix, rhs, solution);
    //  updateVariables(matrix, solution, aux);

    int it = 0;
    const int minNewtonLS = 5;       // minimum number of iterations for the line search
    double res0 = 1., resTrial = {}; // res0 = previous residual, res1 = current residual
    Vec copyRHS;
    VecDuplicate(rhs, &copyRHS);
    VecZeroEntries(copyRHS);

    do
    {
        it++;

        if (it <= minNewtonLS) // No need for LS (Line Search)
        {
            assembleProblem();
            solveLinearSystem(matrix, rhs, solution);
            VecCopy(solution, copyRHS);
            // MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
            // VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
            updateVariables(matrix, solution, rhs, res0);

            for (auto dof : globalDOFs)
            {
                std::cout << "id: " << dof->getIndex()
                          << " currentVal: " << dof->getValue() << std::endl;
            }
            throw std::runtime_error("Newton-Raphson iteration diverged!");
            // computeNorm(rhs, res0);
            resTrial = res0;
        }
        else
        {
            // First try the complete Newton Raphson step
            assembleProblem();
            solveLinearSystem(matrix, rhs, solution);

            // Save the current state of the displacement DOFs
            std::vector<double> backup = {};

            for (auto dof : globalDOFs)
            {
                backup.push_back(dof->getValue());
                std::cout << "id: " << dof->getIndex()
                          << " currentVal: " << dof->getValue()
                          << " backUpVal: " << backup[dof->getIndex()] << std::endl;
            }

            updateVariables(matrix, solution, rhs, resTrial);

            computeNorm(rhs, resTrial); // Compute the residual norm in case xTrial is the solution

            if (resTrial > res0) // The residual has increased
            {
                // Perform the line search
                // 1. Go back to the previous state
                for (auto dof : globalDOFs)
                    dof->setValue(backup[dof->getIndex()]);

                // std::cout << "----------------------------------------------" << std::endl;
                // std::cout << "Line search has been triggered" << std::endl;
                performLineSearch(matrix, solution, rhs, copyRHS, backup, res0); // Perform the line search, gets x_ls
                // VecCopy(rhs, copyRHS);
                // std::cout << "Line search has been completed" << std::endl;
                // std::cout << "----------------------------------------------" << std::endl;
                resTrial = res0;
            }
            else
                res0 = resTrial;
        }
        // PetscPrintf(PETSC_COMM_WORLD, "it: %d  Residual: %e\n", it, resTrial);
    } while (res0 > params->getTolNR() && it < params->getMaxNewtonRaphsonIt());
    VecDestroy(&copyRHS);
    std::cout << "Displacement field solved, residual: " << resTrial << " numIts: " << it << std::endl;
}

PetscErrorCode FEM::solvePhaseField()
{
    assemblePhaseFieldProblem();
    solveSystemByPSOR(matrixPF, rhsPF, solutionPF); // Solves Ddk
    updateFieldVariables(solutionPF, false);        // d = dn + delta_d, where dn is the damage field from the PREVIOUS STEP

    return ierr;
}

PetscErrorCode FEM::assemblePhaseFieldProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);

    PetscCall(MatZeroEntries(matrixPF));
    PetscCall(VecZeroEntries(rhsPF));
    PetscCall(VecZeroEntries(solutionPF));

    // #pragma omp parallel
    //     {
    // #pragma omp for
    for (int Ii = DStart; Ii < DEnd; Ii++)
        elements[Ii]->getPhaseFieldContribution(matrixPF, rhsPF);
    //    }

    PetscCall(VecAssemblyBegin(rhsPF));
    PetscCall(VecAssemblyEnd(rhsPF));
    PetscCall(MatAssemblyBegin(matrixPF, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(matrixPF, MAT_FINAL_ASSEMBLY));

    return ierr;
}

PetscErrorCode FEM::assembleQMatrix(Mat &A)
{
    std::array<Tensor, 3> tensors = computeConstitutiveTensors(); // tensor[0] = K, tensor[1] = I, tensor[2] = C;

    PetscScalar val[totalNnzPF] = {0};
    PetscInt idxRows[totalNnzPF] = {0};
    PetscInt idxCols[totalNnzPF] = {0};

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
            {
                int p1 = kk;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];
                    double Qvalue = elements[elemIndex]->getQValue(idxLocalNode1, idxLocalNode2, prescribedDamageField);

                    val[p1] += Qvalue; // localStiffValue[0];
                }

                idxRows[p1] = nDOFPF * n1; // First DOF

                idxCols[p1] = nDOFPF * n2;
            }
            else
            {
                // idof = 0,       jdof = 0
                int p1 = kk;

                for (int iElem = 0; iElem < numLocalElems; iElem++)
                {
                    const int elemIndex = elemInfo[3 * iElem];
                    const int idxLocalNode1 = elemInfo[3 * iElem + 1];
                    const int idxLocalNode2 = elemInfo[3 * iElem + 2];
                    double Qvalue = elements[elemIndex]->getQValue(idxLocalNode1, idxLocalNode2, prescribedDamageField);

                    val[p1] += Qvalue; // localStiffValue[0];
                }

                idxRows[p1] = nDOFPF * n1;

                idxCols[p1] = nDOFPF * n2;
            }

            iFriendCount++;
            kkn2n++;
            kk += nDOFPF;
        }
    }

    for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++)                // sizeof(idxRows) = size in bytes of the array; sizeof(idxRows[0]) = size in bytes of the first element of the array
        PetscCall(MatSetValue(A, idxRows[i], idxCols[i], val[i], INSERT_VALUES)); // THIS WORKS

    for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++) // Setting the lower triangular part of the matrix
        PetscCall(MatSetValue(A, idxCols[i], idxRows[i], val[i], INSERT_VALUES));

    return ierr;
}

PetscErrorCode FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    getCompressedColumnStorage(A, b);
    MPI_Barrier(PETSC_COMM_WORLD);

    if (rank == 0) // Only rank 0 solves the problem
        PSORAlgorithm();

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscCall(MatSeqAIJRestoreArrayRead(seqMatQ, &PA));
    MPI_Bcast(Ddk, numNodes, MPI_DOUBLE, 0, PETSC_COMM_WORLD); // Communicates the solution to all processes

    // PetscCall(PetscFree(seqMatQ));

    return 0;
}

PetscErrorCode FEM::PSORAlgorithm()
{
    double resPSOR = params->getResPSOR();
    int maxPSORIt = params->getMaxItPSOR();
    double tolPSOR = params->getTolPSOR();
    int itPSOR = 0;

    DdkMinus1 = new double[numNodes]{};

    for (int i = 0; i < numNodes; i++)
    {
        DdkMinus1[i] = 0.0;
        Ddk[i] = 0.0;
    }

    do
    {
        itPSOR++;

        for (int i = 0; i < numNodes; i++)
            DdkMinus1[i] = Ddk[i]; // Update the damage field at the previous iteration

        for (int jCol = 0; jCol < numNodes; jCol++)
        {
            double dotQDd = 0.0;
            double dotLDd = 0.0;

            double Dinv = 0.0;

            int ini = JC[jCol];
            int end = JC[jCol + 1];

            for (int h = ini; h < end; h++)
            {
                int iRow = IR[h];
                dotQDd += PA[h] * DdkMinus1[iRow];

                if (iRow == jCol)
                    Dinv = 1.0 / PA[h];

                if (iRow < jCol)
                    dotLDd += PA[h] * (Ddk[iRow] - DdkMinus1[iRow]);
            }

            Ddk[jCol] = DdkMinus1[jCol] - Dinv * (dotQDd + totalVecq[jCol] + dotLDd);

            if (Ddk[jCol] < 0.0)
                Ddk[jCol] = 0.0; // Damage field must be positive
        }
        resPSOR = computeNorm(Ddk, DdkMinus1, numNodes);

    } while (resPSOR > tolPSOR && itPSOR < maxPSORIt);

    delete[] totalVecq;
    delete[] DdkMinus1;

    // if (rank == 0)
    // {
    //     delete[] JC;
    //     delete[] IR;
    // }

    return 0;
}

PetscErrorCode FEM::getCompressedColumnStorage(Mat &A, Vec &b)
{
    totalVecq = new double[numNodes]{};

    // Get the values from vector q and store them in the totalVecq

    Vec q;
    VecScatter ctx;

    // Gathers the solution vector to the master process
    PetscCall(VecScatterCreateToAll(b, &ctx, &q));
    PetscCall(VecScatterBegin(ctx, b, q, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, b, q, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&ctx));

    if (rank == 0)
        for (int i = 0; i < numNodes; i++)
        {
            PetscInt idx = i;
            PetscScalar val;
            PetscCall(VecGetValues(q, 1, &idx, &val));
            totalVecq[i] = val;
        }

    // Get values from the matrix A and store them in the totalQMatrix
    PetscBool done = PETSC_FALSE;
    PetscInt numNodesInt = numNodes;
    PetscInt *nCols = new PetscInt;
    *nCols = numNodesInt;

    if (!reuse)
    {
        PetscCall(MatCreateRedundantMatrix(A, 0, PETSC_COMM_SELF, MAT_INITIAL_MATRIX, &seqMatQ));
        reuse = PETSC_TRUE;
    }
    else
        PetscCall(MatCreateRedundantMatrix(A, 0, PETSC_COMM_SELF, MAT_REUSE_MATRIX, &seqMatQ));

    if (rank == 0) // ATTENTION: Only process 0 has JC, IR, and PA
    {
        PetscCall(MatGetColumnIJ(seqMatQ, 0, PETSC_FALSE, PETSC_FALSE, nCols, &JC, &IR, &done));
        PetscCall(MatSeqAIJGetArrayRead(seqMatQ, &PA));
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscCall(VecDestroy(&q));

    delete[] nCols;

    return 0;
}

PetscErrorCode FEM::updateFieldVariables(Vec &x, bool _hasConverged)
{
    /*
        getValue() returns the value of the DOF at the previous STEP; dn
        getDamageValue() returns the value of the DOF at the current iteration; delta_d^i (i-th iteration)
    */

    double dtol = 0.9999999999;

    if (_hasConverged)
    {
        for (int i = 0; i < numNodes; i++)
        {

            DOF *damageDOF = discritizedNodes[i]->getDOFs()[2];
            double deltaD = Ddk[i];
            double d_stag = damageDOF->getValue() + deltaD; // damageDOF->getValue() = dn
            damageDOF->setDamageValue(d_stag);
            damageDOF->setValue(d_stag); // dn value
        }

        // Setting d = 1 for the nodes with d > 1
        for (auto node : discritizedNodes)
        {
            DOF *damageDOF = node->getDOFs()[2];
            if (damageDOF->getDamageValue() >= dtol)
            {
                damageDOF->setDamageValue(dtol);
                damageDOF->setValue(dtol);
            }
        }
    }
    else
    {
        for (int i = 0; i < numNodes; i++)
        {
            DOF *damageDOF = discritizedNodes[i]->getDOFs()[2];
            double deltaD = Ddk[i];
            double d_stag = damageDOF->getValue() + deltaD;

            damageDOF->setDamageValue(d_stag);
        }
        // Setting d = 1 for the nodes with d > 1
        for (auto node : discritizedNodes)
        {
            DOF *damageDOF = node->getDOFs()[2];
            if (damageDOF->getDamageValue() >= dtol)
                damageDOF->setDamageValue(dtol);
        }
    }

    return ierr;
}

PetscErrorCode FEM::updateFieldDistribution() // Used only in case a prescribed damage field is setted
{

    MPI_Barrier(PETSC_COMM_WORLD);

    PetscCall(MatZeroEntries(matrixPF));
    PetscCall(VecZeroEntries(rhsPF));
    PetscCall(VecZeroEntries(solutionPF));

    // #pragma omp parallel
    //     {
    // #pragma omp for
    for (int Ii = DStart; Ii < DEnd; Ii++)
        elements[Ii]->getPhaseFieldContribution(matrixPF, rhsPF);
    //}

    PetscCall(VecAssemblyBegin(rhsPF));
    PetscCall(VecAssemblyEnd(rhsPF));
    PetscCall(MatAssemblyBegin(matrixPF, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(matrixPF, MAT_FINAL_ASSEMBLY));

    solveSystemByPSOR(matrixPF, rhsPF, solutionPF); // Solves Ddk
    updateFieldVariables(solutionPF, true);         // d = dn + delta_d, where dn is the damage field from the PREVIOUS STEP

    return ierr;
}

PetscErrorCode FEM::updateRHSPF(Mat &A, Vec &b)
{
    // Performing the multiplication of matrixPF by dn and adding to rhsPF
    Vec QTimesDd, Dn;
    PetscCall(VecDuplicate(b, &QTimesDd));
    PetscCall(VecDuplicate(b, &Dn));
    PetscCall(VecZeroEntries(QTimesDd));
    PetscCall(VecZeroEntries(Dn));

    for (int i = 0; i < numNodes; i++)
    {
        DOF *damageDOF = discritizedNodes[i]->getDOFs()[2];
        PetscScalar value = damageDOF->getValue(); // dn, converged value from the last step
        PetscCall(VecSetValues(Dn, 1, &i, &value, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(QTimesDd));
    PetscCall(VecAssemblyEnd(QTimesDd));
    PetscCall(VecAssemblyBegin(Dn));
    PetscCall(VecAssemblyEnd(Dn));

    PetscCall(MatMult(A, Dn, QTimesDd));  // QTimesDd = matrixPF * Dn
    PetscCall(VecAXPY(b, 1.0, QTimesDd)); // rhsPF = rhsPF - QTimesDd

    PetscCall(VecDestroy(&QTimesDd));
    PetscCall(VecDestroy(&Dn));

    return ierr;
}

PetscErrorCode FEM::evalStaggeredRes(double &_res)
{
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

    // ====================== APPLYING DIRICHLET BOUNDARY CONDITIONS ======================
    PetscCall(MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs)); // Apply Dirichlet boundary conditions

    Vec All;
    VecScatter ctx;

    _res = 0.;

    // Gathers the solution vector to the master process
    PetscCall(VecScatterCreateToAll(rhs, &ctx, &All));
    PetscCall(VecScatterBegin(ctx, rhs, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(ctx, rhs, All, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterDestroy(&ctx));

    for (auto dof : globalDOFs)
    {
        PetscInt Ii = dof->getIndex();
        PetscScalar valForces;
        PetscCall(VecGetValues(All, 1, &Ii, &valForces));
        _res += valForces * valForces;
    }
    VecDestroy(&All);
    return ierr;
}
