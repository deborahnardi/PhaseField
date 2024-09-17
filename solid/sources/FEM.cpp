#include "../headers/FEM.h"

FEM::FEM() {}
FEM::FEM(const std::string _name)
    : name(_name)
{
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
};
FEM::~FEM() {}

/*----------------------------------------------------------------------------------
                Assembling and solving problem PETSc
----------------------------------------------------------------------------------
*/

PetscErrorCode FEM::solveFEMProblem()
{
    assembleProblem();
    solveLinearSystem(matrix, rhs, solution);
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

    for (auto e : partitionedElements)
        e->getContribution(matrix);

    // Neumann boundary conditions
    setBoundaryConditions();

    // Assemble the matrix and the right-hand side vector
    ierr = VecAssemblyBegin(rhs);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // Apply Dirichlet boundary conditions
    ierr = MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs);
    CHKERRQ(ierr);

    // Print the global stiffness matrix on the terminal
    if (showMatrix)
    {
        ierr = PetscPrintf(PETSC_COMM_WORLD, " --- GLOBAL STIFFNESS MATRIX: ----\n");
        CHKERRQ(ierr);
        ierr = MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
        CHKERRQ(ierr);
        printGlobalMatrix(matrix);
    }
}

PetscErrorCode FEM::createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo)
{
    PetscLogDouble bytes;

    (getSolverType() == SEQ)
        ? ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, mSize, mSize, 1800, NULL, &A) // 1800 is the number of non-zero elements per row
        : ierr = MatCreateAIJ(PETSC_COMM_WORLD, size, size, PETSC_DECIDE, PETSC_DECIDE, 0, NULL, 0, NULL, &A);
    CHKERRQ(ierr);

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
        PetscMemoryGetCurrentUsage(&bytes);
        PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes / (1024 * 1024));
        PetscPrintf(PETSC_COMM_WORLD, "Matrix and vectors created...\n");
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
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    switch (getSolverType())
    {
    case SEQ:
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
            std::cout << std::setw(width) << std::fixed << std::setprecision(4) << value << " ";
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
    // Setting NEUMANN boundary conditions

    for (auto bd : bdElements)
        for (auto node : bd->getElemConnectivity())
            for (auto dof : node->getDOFs())
                if (dof->isNeumann())
                {
                    F(dof->getIndex()) += dof->getNeumannValue();
                    numNeumannDOFs++;
                }

    // Setting DIRICHLET boundary conditions -> Adding 0 to the stiffness matrix and to the force vector: column and row of the DOF

    for (auto elem : bdElements)
        for (auto node : elem->getElemConnectivity())
            for (auto dof : node->getDOFs())
                if (dof->isDirichlet())
                {
                    K.row(dof->getIndex()).setZero();
                    K.col(dof->getIndex()).setZero();
                    K(dof->getIndex(), dof->getIndex()) = 1;       // Setting the diagonal to 1
                    F(dof->getIndex()) = dof->getDirichletValue(); // If a prescribed displacement value is given
                    numDirichletDOFs++;
                }

    std::cout << "K:" << std::endl;
    std::cout << K << std::endl;

    std::cout << "F:" << std::endl;
    std::cout << F << std::endl;
}

void FEM::assembleProblemNoPetsc()
{
    int dim = 2;

    for (auto elem : elements)
    {
        elem->getContribution();
        MatrixXd Kelem = elem->getElemStiffnessMatrix();

        // Print the stiffness matrix of each element on the terminal
        std::cout << "Element stiffness matrix: " << std::endl;
        std::cout << Kelem << std::endl;

        int dof1 = elem->getNode(0)->getDOF(0)->getIndex();
        int dof2 = elem->getNode(1)->getDOF(0)->getIndex();

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                K(dof1 + i, dof1 + j) += Kelem(i, j);
                K(dof1 + i, dof2 + j) += Kelem(i, j + 2);
                K(dof2 + i, dof1 + j) += Kelem(i + 2, j);
                K(dof2 + i, dof2 + j) += Kelem(i + 2, j + 2);
            }
    }

    // Print the global stiffness matrix on the terminal
    std::cout << "Global stiffness matrix: " << std::endl;
    std::cout << K << std::endl;
}