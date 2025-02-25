#include "../headers/FEM.h"

/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
------------------------------------------------------------------------------------
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