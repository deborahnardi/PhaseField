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

    for (auto elem : trussElements)
    {
        elem->getContribution();
        MatrixXd Kelem = elem->getElemStiffnessMatrix();

        int n1 = elem->getNode1()->getIndex();
        int n2 = elem->getNode2()->getIndex();

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                K(2 * n1 + i, 2 * n1 + j) += Kelem(i, j);
                K(2 * n1 + i, 2 * n2 + j) += Kelem(i, j + 2);
                K(2 * n2 + i, 2 * n1 + j) += Kelem(i + 2, j);
                K(2 * n2 + i, 2 * n2 + j) += Kelem(i + 2, j + 2);
            }
        }
    }
}