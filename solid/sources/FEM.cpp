#include "../headers/FEM.h"

FEM::FEM() {}
FEM::FEM(const std::string _name, const int &_problemDimension)
    : name(_name), problemDimension(_problemDimension) {};
FEM::~FEM() {}

void FEM::assembleProblem()
{
    int dim = problemDimension;
    MatrixXd K = MatrixXd::Zero(nDOFs, nDOFs);
    VectorXd F = VectorXd::Zero(nDOFs);
    VectorXd U = VectorXd::Zero(nDOFs);

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
    // Solve the linear system

    U = K.fullPivLu().solve(F);
    std::cout << "Displacement vector: " << std::endl;
    std::cout << U << std::endl;
}