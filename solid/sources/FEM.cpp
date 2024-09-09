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

    // for (auto elem : trussElements)
    // {
    //     elem->getContribution();
    //     MatrixXd Kelem = elem->getElemStiffnessMatrix();

    //     for (int d = 0; d < dim; d++)
    //         for (int i = 0; i < 2; i++)
    //             for (int j = 0; j < 2; j++)
    //                 K(2 * elem->getNode(i)->getIndex() + d, 2 * elem->getNode(j)->getIndex() + d) += Kelem(2 * i + d, 2 * j + d);
    // }
}