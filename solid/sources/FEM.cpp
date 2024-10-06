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
    std::cout << "Number of nodes: " << numNodes << std::endl;
    assembleProblem();
    solveLinearSystem(matrix, rhs, solution);
    if (rank == 0)
        showResults();
}

PetscErrorCode FEM::assembleProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling problem...\n");

    createPETScVariables(matrix, rhs, solution, nDOFs, true);

    ierr = MatZeroEntries(matrix);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhs);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solution);
    CHKERRQ(ierr);

    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getContribution(matrix);

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

    ierr = MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs); // Apply Dirichlet boundary conditions
    CHKERRQ(ierr);

    if (showMatrix) // Print the global stiffness matrix on the terminal
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, " --- GLOBAL STIFFNESS MATRIX: ----\n");
        CHKERRQ(ierr);
        ierr = MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRQ(ierr);
        printGlobalMatrix(matrix);
    }
}

PetscErrorCode FEM::createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo) // mSize stands for matrix size, mSize = DOFs = rows = cols
{
    PetscLogDouble bytes;

    (size == 1)
        ? ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, mSize, mSize, 1800, NULL, &A)                                        // 1800 is the number of non-zero elements per row, this can be improved by pre allocating the matrix
        : ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, mSize, mSize, 1800, NULL, 3000, NULL, &A); // Again, 1800 is the number of non-zero elements per row, 3000 is the number of non-zero elements per row in the off-diagonal portion of the matrix. Both can be improved by pre allocating the matrix
    CHKERRQ(ierr);

    // (size == 1)
    //     ? ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, mSize, mSize, 1800, NULL, &A) // 1800 is the number of non-zero elements per row, this can be improved by pre allocating the matrix
    //     : ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD, 2, PETSC_DECIDE, PETSC_DECIDE, mSize, mSize, PETSC_DECIDE, PETSC_DECIDE, NULL, NULL, &A);
    // CHKERRQ(ierr);

    // ----------------------------------------------------------------
    // PARTIONING DOMAIN ELEMENTS
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, elements.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    // PARTIONING BOUNDARY ELEMENTS (FOR NEUMANN CONDITIONS)
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, bdElements.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IIstart, &IIend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    // PARTIONING NODES
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, nodes.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IIIstart, &IIIend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------

    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, mSize);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);

    if (showInfo)
    {
        ierr = PetscMemoryGetCurrentUsage(&bytes);
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes / (1024 * 1024));
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Matrix and vectors created...\n");
        CHKERRQ(ierr);
    }
}

PetscErrorCode FEM::setBoundaryConditions()
{
    for (auto bd : bdElements)
        for (auto node : bd->getElemConnectivity())
            for (auto dof : node->getDOFs())
                if (dof->isNeumann())
                {
                    PetscInt pos = dof->getIndex();
                    PetscScalar value = dof->getNeumannValue();
                    ierr = VecSetValues(rhs, 1, &pos, &value, ADD_VALUES);
                    CHKERRQ(ierr);
                    numNeumannDOFs++;
                }
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

    ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // Prints the Krylov subspace method information
    CHKERRQ(ierr);

    ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);

    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); // Prints the solution vector
    CHKERRQ(ierr);

    // Set the solution to the final coordinates of the nodes
    for (int i = IIIstart; i < IIIend; i++)
    {
        Node *node = nodes[i];
        PetscScalar retrieved[2];
        PetscInt indices[2] = {node->getDOFs()[0]->getIndex(), node->getDOFs()[1]->getIndex()};

        VecGetValues(x, 2, indices, retrieved);
        node->setFinalDisplacement({retrieved[0], retrieved[1], 0.});
    }

    /*
         Clean up
    */

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
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
}
/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
----------------------------------------------------------------------------------
*/
void FEM::solveFEMProblemNoPetsc()
{
    K = MatrixXd::Zero(nDOFs, nDOFs);
    F = VectorXd::Zero(nDOFs);
    U = VectorXd::Zero(nDOFs);

    assembleProblemNoPetsc();
    setBoundaryConditionsNoPetsc();
    solveLinearSystemNoPetsc();
}

void FEM::solveLinearSystemNoPetsc()
{
    U = K.fullPivLu().solve(F);
    std::cout << "Displacement vector: " << std::endl;
    std::cout << U << std::endl;
}

void FEM::setBoundaryConditionsNoPetsc()
{
    std::cout << "K:" << std::endl;
    std::cout << K << std::endl;

    std::cout << "----------------------------------" << std::endl;

    for (auto bd : bdElements)
        bd->getContributionNoPetsc(F, K);

    // // Setting NEUMANN boundary conditions

    // for (auto bd : bdElements)
    //     for (auto node : bd->getElemConnectivity())
    //         for (auto dof : node->getDOFs())
    //             if (dof->isNeumann())
    //             {
    //                 F(dof->getIndex()) += dof->getNeumannValue();
    //                 numNeumannDOFs++;
    //             }

    // // Setting DIRICHLET boundary conditions -> Adding 0 to the stiffness matrix and to the force vector: column and row of the DOF

    // for (auto elem : bdElements)
    //     for (auto node : elem->getElemConnectivity())
    //         for (auto dof : node->getDOFs())
    //             if (dof->isDirichlet())
    //             {
    //                 K.row(dof->getIndex()).setZero();
    //                 K.col(dof->getIndex()).setZero();
    //                 K(dof->getIndex(), dof->getIndex()) = 1;       // Setting the diagonal to 1
    //                 F(dof->getIndex()) = dof->getDirichletValue(); // If a prescribed displacement value is given
    //                 numDirichletDOFs++;
    //             }

    std::cout << "K:" << std::endl;
    std::cout << K << std::endl;

    std::cout << "F:" << std::endl;
    std::cout << F << std::endl;
}

void FEM::assembleProblemNoPetsc()
{
    for (auto elem : elements)
    {
        elem->getContribution();
        elem->assembleGlobalStiffnessMatrix(K);
    }
}