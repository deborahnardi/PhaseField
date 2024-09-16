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

void FEM::assembleProblemPETSc()
{
    // Defining matrix and vector using PETSc

    PetscInt dim = 2;
    PetscInt nDOFs = nodes.size() * dim;
    PetscErrorCode ierr;
    // Getting each element contribution to the global stiffness matrix with PETSc

    for (auto elem : elements)
    {
        elem->getContribution();
        Mat Kelem = elem->getElemStiffnessMatrix();

        // Print the stiffness matrix of each element on the terminal
        CHKERRQ(ierr);

        for (int d = 0; d < dim; d++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    PetscInt row = 2 * elem->getNode(i)->getIndex() + d;
                    PetscInt col = 2 * elem->getNode(j)->getIndex() + d;
                    ierr = MatSetValue(K, row, col, Kelem(2 * i + d, 2 * j + d), ADD_VALUES);
                    CHKERRQ(ierr);
                }
    }
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

void FEM::assembleProblem()
{
    int dim = 2;

    for (auto elem : elements)
    {
        elem->getContribution();
        MatrixXd Kelem = elem->getElemStiffnessMatrix();

        // Print the stiffness matrix of each element on the terminal
        std::cout << "Element stiffness matrix: " << std::endl;
        std::cout << Kelem << std::endl;

        for (int d = 0; d < dim; d++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    K(2 * elem->getNode(i)->getIndex() + d, 2 * elem->getNode(j)->getIndex() + d) += Kelem(2 * i + d, 2 * j + d);
    }

    // Print the global stiffness matrix on the terminal
    std::cout << "Global stiffness matrix: " << std::endl;
    std::cout << K << std::endl;
}