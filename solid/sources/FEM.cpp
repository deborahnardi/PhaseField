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

void FEM::solveFEMProblemPETSc()
{
    // Defining matrix and vector using PETSc

    Mat K;
    Vec F, U;

    PetscErrorCode ierr;

    assembleProblemPETSc();
    setBoundaryConditionsPETSc();
    solveLinearSystemPETSc();
}

PetscErrorCode FEM::createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo)
{
    PetscLogDouble bytes;

    (getSolverType() == SolverType::Sequential)
        ? ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, mSize, mSize, 0, NULL, &A)
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

PetscErrorCode FEM::assembleProblemPETSc()
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

    // boundary conditions de neumann aqui
    // monta matriz e vetores

    // MatZeroRowsColumns(tangent, numDirichletBC, dirichletBC, 1., sol, rhs); APLICA CONDICOES DE DIRICHLET

    // Print the global stiffness matrix on the terminal
    PetscPrintf(PETSC_COMM_WORLD, "Global stiffness matrix:\n");
    MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
}
/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
----------------------------------------------------------------------------------
*/

void FEM::solveFEMProblem()
{
    K = MatrixXd::Zero(nDOFs, nDOFs);
    F = VectorXd::Zero(nDOFs);
    U = VectorXd::Zero(nDOFs);

    assembleProblem();
    setBoundaryConditions();
    solveLinearSystem();
}

void FEM::solveLinearSystem()
{
    U = K.fullPivLu().solve(F);
    std::cout << "Displacement vector: " << std::endl;
    std::cout << U << std::endl;
}

void FEM::setBoundaryConditions()
{
    // Setting NEUMANN boundary conditions

    // for (auto bd : bdElements)
    //     for (auto node : bd->getElemConnectivity())
    //         for (auto dof : node->getDOFs())
    //             if (dof->isNeumann())
    //             {
    //                 F(dof->getIndex()) += dof->getNeumannValue();
    //                 numNeumannDOFs++;
    //             }

    for (auto bd : bdElements)
        for (auto node : bd->getElemConnectivity())
            for (auto dof : node->getDOFs())
                if (dof->isNeumann())
                {
                    PetscInt pos = dof->getIndex();
                    PetscScalar value = dof->getNeumannValue();
                    VecSetValues(rhs, 1, &pos, &value, ADD_VALUES);
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

void FEM::assembleProblem()
{
    // int dim = 2;

    // for (auto elem : elements)
    // {
    //     elem->getContribution();
    //     MatrixXd Kelem = elem->getElemStiffnessMatrix();

    //     // Print the stiffness matrix of each element on the terminal
    //     std::cout << "Element stiffness matrix: " << std::endl;
    //     std::cout << Kelem << std::endl;

    //     for (int d = 0; d < dim; d++)
    //         for (int i = 0; i < 2; i++)
    //             for (int j = 0; j < 2; j++)
    //                 K(2 * elem->getNode(i)->getIndex() + d, 2 * elem->getNode(j)->getIndex() + d) += Kelem(2 * i + d, 2 * j + d);
    // }

    // // Print the global stiffness matrix on the terminal
    // std::cout << "Global stiffness matrix: " << std::endl;
    // std::cout << K << std::endl;
}