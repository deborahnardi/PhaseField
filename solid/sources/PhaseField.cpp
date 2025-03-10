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

    // int rankLocalDOFs = end - start; // Number of nodes in the local partition

    // d_nnz = new PetscInt[rankLocalDOFs]();
    // o_nnz = new PetscInt[rankLocalDOFs]();

    // // In the phase-field problem, there is only one DOF

    // for (auto node : discritizedNodes)
    // {
    //     DOF *damageDOF0 = node->getDOFs()[2];
    //     if (damageDOF0->getIndex() >= IstartPF && damageDOF0->getIndex() < IendPF)
    //         for (auto node2 : n2nMatTotal[node->getIndex()])
    //         {
    //             DOF *damageDOF1 = discritizedNodes[node2]->getDOFs()[2];
    //             if (damageDOF1->getIndex() >= IstartPF && damageDOF1->getIndex() < IendPF)
    //                 d_nnz[damageDOF0->getIndex() - IstartPF]++;
    //             else
    //                 o_nnz[damageDOF0->getIndex() - IstartPF]++;
    //         }
    // }

    // // Total number of non-zero elements in the matrix, stores it in nzQ
    // int localNzQ = 0;
    // nzQ = 0;

    // for (int i = 0; i < rankLocalDOFs; i++)
    //     localNzQ += d_nnz[i] + o_nnz[i];

    // MPI_Reduce(&localNzQ, &nzQ, 1, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD);
    // MPI_Allreduce(&localNzQ, &nzQ, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
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

    Ddk = new double[numNodes]{}; // Damage field at the current iteration

    if (prescribedDamageField)
        updateFieldDistribution();

    prescribedDamageField = false;

    int nSteps = params->getNSteps();

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "================ BEGINNING PHASE FIELD PROBLEM ================\n");

    for (int iStep = 0; iStep < nSteps; iStep++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\n================ STEP %d ================\n", iStep);
        printMemoryUsage(iStep);
        staggeredAlgorithm(iStep); // Returns converged uStag and dStag
        updateFieldVariables(solutionPF);

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
    PetscPrintf(PETSC_COMM_WORLD, "================ END OF PHASE FIELD PROBLEM ================\n");
}

void FEM::staggeredAlgorithm(int _iStep)
{
    int it = 0;
    double resStag = 0.0;
    double previousU[globalDOFs.size()]{}; // Previous displacement field
    double maxTol = params->getTolStaggered();
    do
    {
        it++;
        PetscPrintf(PETSC_COMM_WORLD, "\n------- Stag Iteration %d -------\n", it);
        solveDisplacementField(_iStep); // Obtains ustag
        solvePhaseField();              // Obtains dstag

        // Compute the norm of the difference between the previous and current displacement fields
        double normU = 0.0;
        for (int i = 0; i < globalDOFs.size(); i++)
        {
            DOF *dof = globalDOFs[i];
            if (dof->getDOFType() != D)
            {
                double value = dof->getValue();
                normU += (value - previousU[i]) * (value - previousU[i]);
                previousU[i] = value;
            }
        }
        resStag = sqrt(normU);
        PetscPrintf(PETSC_COMM_WORLD, "Residual stag: %e\n", resStag);

    } while (resStag > params->getTolStaggered() && it < params->getMaxItStaggered());
}

void FEM::solveDisplacementField(int _iStep)
{
    double lambda = (1. + double(_iStep)) / double(params->getNSteps());
    updateBoundaryValues(lambda);

    double entry = double(_iStep) * params->getDeltaTime();
    if (boundaryFunction)                                                // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime()); //

    assembleProblem(0);
    solveLinearSystem(matrix, rhs, solution);
    updateVariables(matrix, solution);
}

PetscErrorCode FEM::solvePhaseField()
{
    assemblePhaseFieldProblem();
    solveSystemByPSOR(matrixPF, rhsPF, solutionPF); // Solves Ddk
    updateFieldVariables(solutionPF, false);        // d = dn + delta_d, where dn is the damage field from the PREVIOUS STEP

    return ierr;
}

PetscErrorCode FEM::updateFieldDistribution() // Used only in case a prescribed damage field is setted
{
    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = MatZeroEntries(matrixPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhsPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solutionPF);
    CHKERRQ(ierr);

    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getPhaseFieldContribution(matrixPF, rhsPF, true);

    ierr = VecAssemblyBegin(rhsPF);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhsPF);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    solveSystemByPSOR(matrixPF, rhsPF, solutionPF); // Solves Ddk
    updateFieldVariables(solutionPF, true);         // d = dn + delta_d, where dn is the damage field from the PREVIOUS STEP

    return ierr;
}

PetscErrorCode FEM::assemblePhaseFieldProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);

    ierr = MatZeroEntries(matrixPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhsPF);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solutionPF);
    CHKERRQ(ierr);

    // for (int Ii = Istart; Ii < Iend; Ii++)
    //     elements[Ii]->getPhaseFieldContribution(matrixPF, rhsPF);

    assembleQMatrix(matrixPF, rhsPF);

    ierr = VecAssemblyBegin(rhsPF);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhsPF);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrixPF, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    if (elemDim == 1)
    {
        // Performing the multiplication of matrixPF by dn and adding to rhsPF - solutionPF is equal Ddk
        Vec QDotDd, Dn;
        ierr = VecDuplicate(rhsPF, &QDotDd);
        CHKERRQ(ierr);
        ierr = VecDuplicate(rhsPF, &Dn);
        CHKERRQ(ierr);
        ierr = VecZeroEntries(QDotDd);
        CHKERRQ(ierr);
        ierr = VecZeroEntries(Dn);
        CHKERRQ(ierr);

        for (int i = 0; i < numNodes; i++)
        {
            DOF *damageDOF = discritizedNodes[i]->getDOFs()[2];
            PetscScalar value = damageDOF->getValue(); // dn, converge value from the last step
            ierr = VecSetValues(Dn, 1, &i, &value, INSERT_VALUES);
            CHKERRQ(ierr);
        }

        ierr = VecAssemblyBegin(QDotDd);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(QDotDd);
        CHKERRQ(ierr);

        ierr = VecAssemblyBegin(Dn);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(Dn);
        CHKERRQ(ierr);

        ierr = MatMult(matrixPF, Dn, QDotDd); // QDotDd = matrixPF * Dn
        CHKERRQ(ierr);

        ierr = VecAXPY(rhsPF, 1.0, QDotDd); // rhsPF = rhsPF + QDotDd
        CHKERRQ(ierr);

        ierr = VecDestroy(&QDotDd);
        CHKERRQ(ierr);
        ierr = VecDestroy(&Dn);
        CHKERRQ(ierr);
    }

    return ierr;
}

PetscErrorCode FEM::assembleQMatrix(Mat &A, Vec &b)
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
                    // std::vector<double> localStiffValue = elements[elemIndex]->getStiffnessIIOrIJ(tensors, idxLocalNode1, idxLocalNode2);
                    double Qvalue = elements[elemIndex]->getQValue(idxLocalNode1, idxLocalNode2);

                    val[p1] += 1; // localStiffValue[0];
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
                    double Qvalue = elements[elemIndex]->getQValue(idxLocalNode1, idxLocalNode2);

                    val[p1] += 1; // localStiffValue[0];
                }

                idxRows[p1] = nDOFPF * n1;

                idxCols[p1] = nDOFPF * n2;
            }

            iFriendCount++;
            kkn2n++;
            kk += nDOFPF;
        }

        // kk += nDOFPF * numFriends - (nDOFPF * nDOFPF - (nDOFPF * (nDOFPF + 1)) / 2.0);
    }

    // for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++)                // sizeof(idxRows) = size in bytes of the array; sizeof(idxRows[0]) = size in bytes of the first element of the array
    //     PetscCall(MatSetValue(A, idxRows[i], idxCols[i], val[i], INSERT_VALUES)); // THIS WORKS

    // for (int i = 0; i < sizeof(idxRows) / sizeof(idxRows[0]); i++) // Setting the lower triangular part of the matrix
    //     PetscCall(MatSetValue(A, idxCols[i], idxRows[i], val[i], INSERT_VALUES));

    return ierr;
}

PetscErrorCode FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    assembleBetweenProcesses(A, b);
    getPSORVecs();

    double resPSOR = params->getResPSOR();
    int maxPSORIt = params->getMaxItPSOR();
    double tolPSOR = params->getTolPSOR();
    int itPSOR = 0;

    MPI_Barrier(PETSC_COMM_WORLD);

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

    delete[] JC;
    delete[] IR;
    delete[] PA;
    delete[] totalVecq;
    delete[] DdkMinus1;

    return ierr;
}

PetscErrorCode FEM::assembleBetweenProcesses(Mat &A, Vec &b)
{
    PetscInt startRow, endRow;
    ierr = MatGetOwnershipRange(A, &startRow, &endRow);
    CHKERRQ(ierr);
    int rankLocalRows = endRow - startRow;

    totalMatrixQ = new double *[numNodes] {};
    for (int i = 0; i < numNodes; i++)
        totalMatrixQ[i] = new double[numNodes]{};

    totalVecq = new double[numNodes]{};

    // Get the values from vector q and store them in the totalVecq

    double *localVecq = new double[rankLocalRows]{};
    for (int i = startRow; i < endRow; i++)
    {
        PetscScalar value;
        ierr = VecGetValues(b, 1, &i, &value);
        CHKERRQ(ierr);
        localVecq[i - startRow] = value;
    }

    int *numEachProcess = new int[size];
    MPI_Allgather(&rankLocalRows, 1, MPI_INT, numEachProcess, 1, MPI_INT, PETSC_COMM_WORLD);

    int *globalBuffer0 = new int[size];
    globalBuffer0[0] = 0;
    for (int i = 1; i < size; i++)
        globalBuffer0[i] = globalBuffer0[i - 1] + numEachProcess[i - 1];

    MPI_Gatherv(localVecq, rankLocalRows, MPI_DOUBLE, totalVecq, numEachProcess, globalBuffer0, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(totalVecq, numNodes, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Get values from the matrix A and store them in the totalQMatrix

    // Broadcast the totalMatrixQ to all processes
    double *totalMatrixQVecFormat = new double[numNodes * numNodes]{};
    double *localQMatrix = new double[rankLocalRows * numNodes]{};

    for (PetscInt i = startRow; i < endRow; i++)
    {
        for (int j = 0; j < numNodes; j++)
        {
            PetscScalar value;
            ierr = MatGetValues(A, 1, &i, 1, &j, &value);
            CHKERRQ(ierr);
            localQMatrix[(i - startRow) * numNodes + j] = value;
        }
    }

    int numIJ = rankLocalRows * numNodes;
    MPI_Allgather(&numIJ, 1, MPI_INT, numEachProcess, 1, MPI_INT, PETSC_COMM_WORLD);

    int *globalBuffer = new int[size];
    globalBuffer[0] = 0;
    for (int i = 1; i < size; i++)
        globalBuffer[i] = globalBuffer[i - 1] + numEachProcess[i - 1];

    MPI_Gatherv(localQMatrix, rankLocalRows * numNodes, MPI_DOUBLE, totalMatrixQVecFormat, numEachProcess, globalBuffer, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    MPI_Bcast(totalMatrixQVecFormat, numNodes * numNodes, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    // Convert the totalMatrixQVecFormat to totalMatrixQ
    for (int i = 0; i < numNodes; i++)
        for (int j = 0; j < numNodes; j++)
            totalMatrixQ[i][j] = totalMatrixQVecFormat[i * numNodes + j];

    // Print the totalMatrixQ
    std::cout << "Number of nodes: " << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++)
    {
        for (int j = 0; j < numNodes; j++)
            std::cout << totalMatrixQ[i][j] << " ";
        std::cout << std::endl;
    }

    exit(0);

    // Erase the Memory
    delete[] totalMatrixQVecFormat;
    delete[] localQMatrix;
    delete[] numEachProcess;
    delete[] globalBuffer0;
    delete[] localVecq;
    delete[] globalBuffer;

    return ierr;
}

PetscErrorCode FEM::getPSORVecs()
{
    // Get the matrix A in compressed column format

    JC = new int[numNodes + 1]{};
    IR = new int[nzQ]{};
    PA = new double[nzQ]{};

    int count = -1., var = 0;
    for (int jj = 0; jj < numNodes; jj++)
    {
        var = 0;
        for (int ii = 0; ii < numNodes; ii++)
        {
            if (totalMatrixQ[ii][jj] != 0.0)
            {
                count++;
                PA[count] = totalMatrixQ[ii][jj];
                IR[count] = ii;

                if (var == 0)
                {
                    JC[jj] = count;
                    var = 1;
                }
            }
        }
    }

    JC[numNodes] = nzQ; // nzQ is the total number of non-zero elements in the matrix

    for (int i = 0; i < numNodes; i++)
        delete[] totalMatrixQ[i];

    delete[] totalMatrixQ;

    return ierr;
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