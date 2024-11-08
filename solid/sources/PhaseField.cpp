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
    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];
    norm = sqrt(norm);

    matrixPreAllocationPF(IstartPF, IendPF);
    createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, true);

    DdkMinus1 = new double[numNodes]{}; // Damage field at the previous iteration
    Ddk = new double[numNodes]{};       // Damage field at the current iteration

    int nSteps = params->getNSteps();

    std::cout << "BEGINNING OF PHASE FIELD PROBLEM..." << std::endl;
    for (int iStep = 0; iStep < nSteps; iStep++)
    {
        std::cout << "STEP " << iStep << std::endl;
        staggeredAlgorithm(iStep); // Returns converged uStag and dStag
        updateFieldVariables(solutionPF);
        std::cout << "-----------------------------------------------" << std::endl;
    }
}

void FEM::staggeredAlgorithm(int _iStep)
{
    int it = 0;
    double resStag = 0.0;
    double previousU[globalDOFs.size()]; // Previous displacement field

    do
    {
        it++;
        solveDisplacementField(_iStep); // Obtains u^{n+1}
        solvePhaseField();              // Obtains D^{n+1} (d^n)

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

    } while (resStag > params->getTolStaggered() && it < params->getMaxItStaggered());
}

void FEM::solveDisplacementField(int _iStep)
{
    double lambda = (1. + double(_iStep)) / double(params->getNSteps());
    updateBoundaryValues(lambda);

    if (boundaryFunction)                                                // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime()); //

    assembleProblem();
    solveLinearSystem(matrix, rhs, solution);
    updateVariables(solution);

    // if (rank == 0)
    //     showResults(_iStep); // Paraview
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

    // Adding the potential energy evaluated at the previous load step into the rhs
    // Performing the multiplication of matrixPF by dn and adding to rhsPF - solutionPF is equal Ddk
    Vec QDotDd, Dd;
    ierr = VecDuplicate(rhsPF, &QDotDd);
    CHKERRQ(ierr);
    ierr = VecDuplicate(rhsPF, &Dd);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(QDotDd);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(Dd);
    CHKERRQ(ierr);

    for (int i = 0; i < numNodes; i++)
    {
        DOF *damageDOF = nodes[i]->getDOFs()[2];
        PetscScalar value = damageDOF->getValue();
        ierr = VecSetValues(Dd, 1, &i, &value, INSERT_VALUES);
        CHKERRQ(ierr);
    }

    ierr = MatMult(matrixPF, Dd, QDotDd); // QDotDd = matrixPF * Dd
    CHKERRQ(ierr);

    ierr = VecAXPY(rhsPF, 1.0, QDotDd); // rhsPF = rhsPF + QDotDd
    CHKERRQ(ierr);

    std::cout << " q: " << std::endl;
    ierr = VecView(rhsPF, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);
    std::cout << "-----------------------------------------------" << std::endl;

    return ierr;
}

PetscErrorCode FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    getPSORVecs(A, b);

    double resPSOR = params->getResPSOR();
    int maxPSORIt = params->getMaxItPSOR();
    double tolPSOR = params->getTolPSOR();
    int itPSOR = 0;

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

            double q = 0.0;
            ierr = VecGetValues(b, 1, &jCol, &q);
            CHKERRQ(ierr);
            Ddk[jCol] = DdkMinus1[jCol] - Dinv * (dotQDd + q + dotLDd);

            if (Ddk[jCol] < 0.0)
                Ddk[jCol] = 0.0; // Damage field must be positive
        }

        resPSOR = computeNorm(Ddk, DdkMinus1, numNodes);

    } while (resPSOR > tolPSOR && itPSOR < maxPSORIt);

    // print PA
    std::cout << "PA array:" << std::endl;
    for (int i = 0; i < JC[numNodes]; i++)
        std::cout << PA[i] << " ";

    std::cout << std::endl;

    // print Ddk
    std::cout << "Ddk: " << std::endl;
    for (int i = 0; i < numNodes; i++)
        std::cout << Ddk[i] << " ";

    std::cout << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    // Update the solution vector (solution vector is Delta_d)
    for (int i = 0; i < numNodes; i++)
    {
        ierr = VecSetValues(x, 1, &i, &Ddk[i], INSERT_VALUES);
        CHKERRQ(ierr);
    }

    return ierr;
}

PetscErrorCode FEM::getPSORVecs(Mat &A, Vec &b)
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

    // Get the matrix A in compressed column format
    ierr = MatGetRowIJ(A, 0, PETSC_TRUE, PETSC_FALSE, &nRows, &JC, &IR, &done);
    CHKERRQ(ierr);

    PA = new PetscScalar[JC[nRows]]; // Non-zero values of the matrix A

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

    return ierr;
}

PetscErrorCode FEM::updateFieldVariables(Vec &x, bool _hasConverged)
{
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
                damageDOF->setDamageValue(1.0);
        }
    }
    else
    {
        for (int i = 0; i < numNodes; i++)
        {
            DOF *damageDOF = nodes[i]->getDOFs()[2];
            PetscScalar value;
            ierr = VecGetValues(x, 1, &i, &value);
            CHKERRQ(ierr);
            damageDOF->setDamageValue(value);
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
    if (_hasConverged)
    {
        for (auto node : nodes)
        {
            DOF *damageDOF = node->getDOFs()[2];
            std::cout << "Damage field at node " << node->getIndex() << ": " << damageDOF->getDamageValue() << std::endl;
        }
    }

    return ierr;
}