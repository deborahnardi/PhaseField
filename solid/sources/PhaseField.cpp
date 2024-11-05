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

void FEM::solvePhaseFieldProblem()
{
    for (auto n : nodes)
        norm += n->getInitialCoordinates()[0] * n->getInitialCoordinates()[0] + n->getInitialCoordinates()[1] * n->getInitialCoordinates()[1];
    norm = sqrt(norm);

    matrixPreAllocationPF(IstartPF, IendPF);
    createPETScVariables(matrixPF, rhsPF, solutionPF, numNodes, true);

    DdkMinus1 = new double[numNodes]{}; // Damage field at the previous iteration
    Ddk = new double[numNodes]{};       // Damage field at the current iteration

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

    if (boundaryFunction)                                                // 0 is false, any non zero value is true;
        updateBoundaryFunction(double(_iStep) * params->getDeltaTime()); //

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
    solveSystemByPSOR(matrixPF, rhsPF, solutionPF); // Solves Ddk
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

    // Adding the potential energy evaluated at the previous load step into the rhs
    // Performing the multiplication of matrixPF by dn and adding to rhsPF - solutionPF is equal Ddk
    Vec QDotDd;
    ierr = VecDuplicate(rhsPF, &QDotDd);
    CHKERRQ(ierr);

    ierr = MatMult(matrixPF, solutionPF, QDotDd); // QDotDd = matrixPF * solutionPF
    CHKERRQ(ierr);

    double potentialEnergy;
    ierr = VecDot(QDotDd, solutionPF, &potentialEnergy);
    CHKERRQ(ierr);

    

    return ierr;
}

PetscErrorCode FEM::solveSystemByPSOR(Mat &A, Vec &b, Vec &x)
{
    getPSORVecs(A, b);

    // Solve the system using the PSOR method

    double resPSOR = params->getResPSOR();
    int maxPSORIt = params->getMaxItPSOR();
    double tolPSOR = params->getTolPSOR();
    int itPSOR = 0;

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

    // Update the solution vector
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

    // // Print PA
    // std::cout << "PA array:" << std::endl;
    // for (int i = 0; i < JC[nRows]; i++)
    //     std::cout << PA[i] << " ";

    // // Print IR and JC arrays
    // int nnz = JC[nRows];
    // int jcSize = nRows + 1;
    // int irSize = nnz;
    // std::cout << "-------------------------" << std::endl;
    // std::cout << "IR and JC arrays:" << std::endl;
    // for (int i = 0; i < jcSize; i++)
    //     std::cout << JC[i] << " ";
    // std::cout << std::endl;
    // for (int i = 0; i < irSize; i++)
    //     std::cout << IR[i] << " ";
    // std::cout << std::endl;

    return ierr;
}

PetscErrorCode FEM::updateFieldVariables(Vec &x)
{
    // Update the damage field
    for (int i = 0; i < numNodes; i++)
    {
        DOF *damageDOF = nodes[i]->getDOFs()[2];
        double value;
        ierr = VecGetValues(x, 1, &i, &value);
        CHKERRQ(ierr);
        damageDOF->incrementValue(value);
    }

    return ierr;
}