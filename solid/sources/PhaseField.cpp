#include "../headers/FEM.h"

void FEM::matrixPreAllocationPF(PetscInt start, PetscInt end)
{
    int rankLocalDOFs = end - start; // Number of nodes in the local partition

    d_nnz = new PetscInt[rankLocalDOFs]();
    o_nnz = new PetscInt[rankLocalDOFs]();

    // In the phase-field problem, there is only one DOF

    for (auto node : nodes)
    {
        DOF *damageDOF0 = node->getDOFs()[2];
        if (damageDOF0->getIndex() >= IstartPF && damageDOF0->getIndex() < IendPF)
            for (auto node2 : nodeNeighbours[node->getIndex()])
            {
                DOF *damageDOF1 = nodes[node2]->getDOFs()[2];
                if (damageDOF1->getIndex() >= IstartPF && damageDOF1->getIndex() < IendPF)
                    d_nnz[damageDOF0->getIndex() - IstartPF]++;
                else
                    o_nnz[damageDOF0->getIndex() - IstartPF]++;
            }
    }
}

// --------------------------------------------------------------------------------------------------
void FEM::solvePhaseFieldProblem() // Called by the main program
{
    auto start_timer = std::chrono::high_resolution_clock::now();

    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];
    norm = sqrt(norm);

    matrixPreAllocationPF(IstartPF, IendPF);
    createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, true);

    DdkMinus1 = new double[numNodes]{}; // Damage field at the previous iteration
    Ddk = new double[numNodes]{};       // Damage field at the current iteration

    int nSteps = params->getNSteps();

    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "================ BEGINNING PHASE FIELD PROBLEM ================\n");

    for (int iStep = 0; iStep < nSteps; iStep++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\n================ STEP %d ================\n", iStep);

        staggeredAlgorithm(iStep); // Returns converged uStag and dStag
        updateFieldVariables(solutionPF);

        if (rank == 0)
            showResults(iStep); // Paraview
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

    if (boundaryFunction)                                                // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime()); //

    for (auto node : nodes)
    {
        for (auto dof : node->getDOFs())
            if (dof->isControlledDOF() && dof->getValue() < 0.)
            {
                negativeLoad = true;
                break; // Get out of first loop
            }
        if (negativeLoad)
            break; // Get out of second loop
    }

    assembleProblem();
    solveLinearSystem(matrix, rhs, solution);
    updateVariables(solution);
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
            DOF *damageDOF = nodes[i]->getDOFs()[2];
            PetscScalar value = damageDOF->getValue();
            // PetscScalar value = damageDOF->getDamageValue();
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

PetscErrorCode FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    assembleBetweenProcesses(A, b);
    getPSORVecs(A);

    double resPSOR = params->getResPSOR();
    int maxPSORIt = params->getMaxItPSOR();
    double tolPSOR = params->getTolPSOR();
    int itPSOR = 0;

    MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = 0; i < numNodes; i++)
    {
        DdkMinus1[i] = 0.0;
        Ddk[i] = 0.0;
    }

    // Print PA, IR and JC

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

            double q = 0.0;
            ierr = VecGetValues(totalVecq, 1, &jCol, &q);
            CHKERRQ(ierr);
            Ddk[jCol] = DdkMinus1[jCol] - Dinv * (dotQDd + q + dotLDd);

            if (Ddk[jCol] < 0.0)
                Ddk[jCol] = 0.0; // Damage field must be positive
            // PetscPrintf(PETSC_COMM_WORLD, "=========== OVER HERE 1 ===========\n");
        }
        // PetscPrintf(PETSC_COMM_WORLD, "=========== OVER HERE 2 ===========\n");
        resPSOR = computeNorm(Ddk, DdkMinus1, numNodes);

    } while (resPSOR > tolPSOR && itPSOR < maxPSORIt);

    // Update the solution vector (solution vector is Delta_d)

    for (int i = 0; i < numNodes; i++)
    {
        ierr = VecSetValues(x, 1, &i, &Ddk[i], INSERT_VALUES);
        CHKERRQ(ierr);
    }

    return ierr;
}

PetscErrorCode FEM::assembleBetweenProcesses(Mat &A, Vec &b)
{
    VecScatter ctx;

    ierr = VecScatterCreateToAll(b, &ctx, &totalVecq);
    CHKERRQ(ierr);

    ierr = VecScatterBegin(ctx, b, totalVecq, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, b, totalVecq, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    CHKERRQ(ierr);

    ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    // Get values from the matrix A and store them in the totalQMatrix

    // PetscInt startRow, endRow;
    // ierr = MatGetOwnershipRange(A, &startRow, &endRow);
    // CHKERRQ(ierr);

    // for (PetscInt i = startRow; i < endRow; i++)
    // {
    //     PetscInt ncols;
    //     const PetscInt *cols;
    //     const PetscScalar *vals;
    //     ierr = MatGetRow(A, i, &ncols, &cols, &vals);
    //     CHKERRQ(ierr);

    //     for (PetscInt j = 0; j < ncols; j++)
    //     {
    //         ierr = MatSetValue(totalQMatrix, i, cols[j], vals[j], INSERT_VALUES);
    //         CHKERRQ(ierr);
    //     }

    //     ierr = MatRestoreRow(A, i, &ncols, &cols, &vals);
    //     CHKERRQ(ierr);
    // }

    // ierr = MatAssemblyBegin(totalQMatrix, MAT_FINAL_ASSEMBLY);
    // CHKERRQ(ierr);
    // ierr = MatAssemblyEnd(totalQMatrix, MAT_FINAL_ASSEMBLY);
    // CHKERRQ(ierr);

    // // Print the totalQMatrix

    // ierr = MatView(totalQMatrix, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    // PetscPrintf(PETSC_COMM_WORLD, "=========== OVER HERE 3 ===========\n");

    return ierr;
}

PetscErrorCode FEM::getPSORVecs(Mat &A)
{
    /*
        Here the PA, IR and JC arrays are used to store the matrix A in compressed column format;
        PA stores the non-zero values of the matrix A;
        IR stores the row indices of the non-zero values of the matrix A;
        JC stores the pointers to the entries of the IR array;

        CCS and CRS are the same if the matrix is symmetric;
    */
    PetscInt nRows; // Number of rows in the matrix A
    PetscBool done; // Flag to indicate if the matrix A has been completely processed
    const PetscInt *localJC;
    const PetscInt *localIR;
    PetscScalar *localPA;

    // Get the matrix A in compressed column format

    ierr = MatGetRowIJ(A, 0, PETSC_TRUE, PETSC_FALSE, &nRows, &JC, &IR, &done);
    CHKERRQ(ierr);
    localPA = new PetscScalar[JC[nRows]]; // Non-zero values of the matrix A
    PA = new PetscScalar[JC[nRows]];      // Non-zero values of the matrix A
    // Get the non-zero values of the matrix A
    for (int i = 0; i < nRows; i++)
    {
        for (int j = JC[i]; j < JC[i + 1]; j++) // Iterate over the non-zero values
        {
            PetscInt col = IR[j];
            ierr = MatGetValues(A, 1, &i, 1, &col, &PA[j]);
            CHKERRQ(ierr);
        }
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    // print PA
    PetscPrintf(PETSC_COMM_WORLD, "PA array:\n");
    for (int i = 0; i < JC[numNodes]; i++)
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PA[i]);

    // std::cout << std::endl;

    // // print IR
    // std::cout << "IR array:" << std::endl;
    // for (int i = 0; i < JC[numNodes]; i++)
    //     std::cout << IR[i] << " ";

    // std::cout << std::endl;

    // // print JC
    // std::cout << "JC array:" << std::endl;
    // for (int i = 0; i < numNodes + 1; i++)
    //     std::cout << JC[i] << " ";

    return ierr;
}

PetscErrorCode FEM::updateFieldVariables(Vec &x, bool _hasConverged)
{

    /*
        getValue() returns the value of the DOF at the previous STEP; dn
        getDamageValue() returns the value of the DOF at the current iteration; delta_d^i (i-th iteration)
    */

    if (_hasConverged)
    {
        for (int i = 0; i < numNodes; i++)
        {
            DOF *damageDOF = nodes[i]->getDOFs()[2];
            PetscScalar value;
            ierr = VecGetValues(x, 1, &i, &value);
            CHKERRQ(ierr);
            damageDOF->incrementValue(value);
            double convergedValue = damageDOF->getValue();
            damageDOF->setValue(convergedValue);
            damageDOF->setDamageValue(convergedValue);
        }

        // Setting d = 1 for the nodes with d > 1
        for (auto node : nodes)
        {
            DOF *damageDOF = node->getDOFs()[2];
            if (damageDOF->getDamageValue() > 1.0)
            {
                damageDOF->setDamageValue(1.0);
                damageDOF->setValue(1.0);
            }
        }
    }
    else
    {
        for (int i = 0; i < numNodes; i++)
        {
            DOF *damageDOF = nodes[i]->getDOFs()[2];
            PetscScalar value = 0;
            ierr = VecGetValues(x, 1, &i, &value);
            CHKERRQ(ierr);
            PetscScalar d_stag = damageDOF->getValue() + value;
            damageDOF->setDamageValue(d_stag);
        }

        // Setting d = 1 for the nodes with d > 1
        for (auto node : nodes)
        {
            DOF *damageDOF = node->getDOFs()[2];
            if (damageDOF->getDamageValue() > 1.0)
                damageDOF->setDamageValue(1.0);
        }
    }

    // Print the damage field
    // if (_hasConverged)
    // {
    //     for (auto node : nodes)
    //     {
    //         DOF *damageDOF = node->getDOFs()[2];
    //         // std::cout << "Damage field at node " << node->getIndex() << ": " << damageDOF->getDamageValue() << std::endl;
    //         std::cout << node->getX() << " " << damageDOF->getDamageValue() << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    return ierr;
}