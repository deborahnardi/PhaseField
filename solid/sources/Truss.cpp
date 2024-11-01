#include "../headers/Element.h"

Truss::Truss() {}
Truss::Truss(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, const double &area_)
    : Element(_index, _elemDimension, _elemConnectivity, _material, _physicalEntity), area(area_)
{
    Node *_node1 = elemConnectivity[0];
    Node *_node2 = elemConnectivity[1];
    area = area_;

    length = pow((_node2->getX() - _node1->getX()) * (_node2->getX() - _node1->getX()) +
                     (_node2->getY() - _node1->getY()) * (_node2->getY() - _node1->getY()),
                 0.5);

    theta = atan((_node2->getY() - _node1->getY()) / (_node2->getX() - _node1->getX())); // Theta is in radians

    if (_node2->getX() - _node1->getX() < 0)
        theta += M_PI;
    else if (_node2->getX() - _node1->getX() == 0) // If the element is on the vertical
    {
        if (_node2->getY() - _node1->getY() > 0) // If the element is pointing upwards
            theta = M_PI / 2;
        else if (_node2->getY() - _node1->getY() < 0) // If the element is pointing downwards
            theta = -M_PI / 2;
    }
}
Truss::~Truss() {}

/*----------------------------------------------------------------------------------
                Assembling and solving problem with PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode Truss::getContribution(Mat &matrix, Vec &rhs)
{
    PetscInt numElDOF = 4;
    PetscInt *indx = new PetscInt[4]();
    PetscScalar *localStiff = new PetscScalar[16]();
    PetscInt count = 0;

    PetscScalar d0 = 0., d1 = 0.;

    d0 = elemConnectivity[0]->getDOF(2)->getDamageValue();
    d1 = elemConnectivity[1]->getDOF(2)->getDamageValue();

    PetscScalar damageValue = 1 / 3. * (d0 * d0 + (d1 - 3.) * d0 + d1 * d1 - 3. * d1 + 3.);

    PetscScalar k = damageValue * material->getYoungModulus() * area / length;

    localStiff[0] = k * cos(theta) * cos(theta);
    localStiff[1] = k * cos(theta) * sin(theta);
    localStiff[2] = k * -cos(theta) * cos(theta);
    localStiff[3] = k * -cos(theta) * sin(theta);
    localStiff[4] = k * cos(theta) * sin(theta);
    localStiff[5] = k * sin(theta) * sin(theta);
    localStiff[6] = k * -cos(theta) * sin(theta);
    localStiff[7] = k * -sin(theta) * sin(theta);
    localStiff[8] = k * -cos(theta) * cos(theta);
    localStiff[9] = k * -cos(theta) * sin(theta);
    localStiff[10] = k * cos(theta) * cos(theta);
    localStiff[11] = k * cos(theta) * sin(theta);
    localStiff[12] = k * -cos(theta) * sin(theta);
    localStiff[13] = k * -sin(theta) * sin(theta);
    localStiff[14] = k * cos(theta) * sin(theta);
    localStiff[15] = k * sin(theta) * sin(theta);

    for (auto node : elemConnectivity)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
                indx[count++] = dof->getIndex();

    ierr = MatSetValues(matrix, numElDOF, indx, numElDOF, indx, localStiff, ADD_VALUES);
    CHKERRQ(ierr);

    // Setting prescribed Dirichlet boundary conditions
    for (PetscInt a = 0; a < numElNodes; a++)
        for (PetscInt i = 0; i < 2; i++)
            for (PetscInt b = 0; b < numElNodes; b++)
                for (PetscInt j = 0; j < 2; j++)
                {
                    DOF *dof = elemConnectivity[b]->getDOFs()[j];
                    double value = dof->getDirichletValue();

                    if (value != 0)
                    {
                        int pos = numElDOF * (2 * a + i) + 2 * b + j;
                        int ps2 = 2 * a + i;
                        double fi = -localStiff[numElDOF * (2 * a + i) + 2 * b + j] * value;
                        ierr = VecSetValues(rhs, 1, &indx[2 * a + i], &fi, ADD_VALUES);
                        CHKERRQ(ierr);
                    }
                }

    // for (PetscInt a = 0; a < numElNodes; a++)
    //     for (PetscInt i = 0; i < 2; i++)
    //         for (PetscInt b = 0; b < numElNodes; b++)
    //             for (PetscInt j = 0; j < 2; j++)
    //             {
    //                 double value = elemConnectivity[b]->getDOFs()[j]->getValue();
    //                 double fi = -localStiffnessMatrix[numElDOF * (2 * a + i) + 2 * b + j] * elemConnectivity[b]->getDOF(j)->getValue();
    //                 ierr = VecSetValues(rhs, 1, &idx[2 * a + i], &fi, ADD_VALUES);
    //                 CHKERRQ(ierr);
    //             }

    delete[] indx;
    delete[] localStiff;
}

PetscErrorCode Truss::getPhaseFieldContribution(Mat &A, Vec &rhs)
{
    /*
        Analytical solution for the phase field problem in a truss element
    */
    const PetscInt numNodeDOF = 1;                                  // Number of DOFs per node considering displacements only
    PetscInt numElDOF = numElNodes;                                 // Only one DOF per node when considering only phase field
    PetscReal *localStiff = new PetscScalar[numElDOF * numElDOF](); // Equivalent to matrix Qlocal in the phase field problem
    PetscReal *localRHS = new PetscReal[numElDOF]();
    PetscInt *idx = new PetscInt[numElDOF]();
    PetscScalar l0 = material->getL0();
    PetscScalar elas = material->getYoungModulus();
    PetscScalar Gc = material->getGriffithCriterion();

    PetscInt count = 0;
    for (auto node : elemConnectivity)
        idx[count++] = node->getIndex(); // Phase field DOF has the same index as the node for the local problem

    PetscScalar Qkk = l0 / length + length / (3. * l0);
    PetscScalar Qjk = -l0 / length + length / (6. * l0);

    computeDeformation();

    double def = getDeformation();

    localStiff[0] = Gc * area * Qkk;
    localStiff[1] = Gc * area * Qjk;
    localStiff[2] = Gc * area * Qjk;
    localStiff[3] = Gc * area * Qkk;

    if (def > 0.)
    {
        PetscScalar gammaU = def * def * elas * area * length / 6.;
        PetscScalar Qkk2 = 2 * gammaU;
        PetscScalar Qjk2 = gammaU;

        localStiff[0] += Qkk2;
        localStiff[1] += Qjk2;
        localStiff[2] += Qjk2;
        localStiff[3] += Qkk2;

        localRHS[0] = -3. * gammaU;
        localRHS[1] = -3. * gammaU;
    }

    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localStiff, ADD_VALUES);
    CHKERRQ(ierr);

    ierr = VecSetValues(rhs, numElDOF, idx, localRHS, ADD_VALUES);
    CHKERRQ(ierr);

    // Print idx
    for (int i = 0; i < 2; i++)
        std::cout << idx[i] << " ";

    // Print localStiff
    for (int i = 0; i < 4; i++)
        std::cout << localStiff[i] << " ";

    delete[] idx;
    delete[] localStiff;

    return ierr;
}

void Truss::computeDeformation()
{
    Node *_node0 = elemConnectivity[0];
    Node *_node1 = elemConnectivity[1];

    double currentLenght = pow((_node1->getDOF(0)->getValue() - _node0->getDOF(0)->getValue()) * (_node1->getDOF(0)->getValue() - _node0->getDOF(0)->getValue()) + (_node1->getDOF(1)->getValue() - _node0->getDOF(1)->getValue()) * (_node1->getDOF(1)->getValue() - _node0->getDOF(1)->getValue()),
                               0.5);

    epsilon = (currentLenght - length) / length;
    setDeformation(epsilon);
}

/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
----------------------------------------------------------------------------------
*/

void Truss::getContribution()
{
    localStiffnessMatrix = MatrixXd::Zero(4, 4);
    localStiffnessMatrix << 1, 0, -1, 0,
        0, 0, 0, 0,
        -1, 0, 1, 0,
        0, 0, 0, 0;
    localStiffnessMatrix *= material->getYoungModulus() * area / length;

    rotationMatrix = MatrixXd::Zero(4, 4);
    rotationMatrix << cos(theta), sin(theta), 0, 0,
        -sin(theta), cos(theta), 0, 0,
        0, 0, cos(theta), sin(theta),
        0, 0, -sin(theta), cos(theta);
    localStiff = rotationMatrix.transpose() * localStiffnessMatrix * rotationMatrix;
}

void Truss::assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff)
{
    /*
        dof1: first degree of freedom of the first node
        dof2: second degree of freedom of the first node

        Note that in the loop below, dof1 + i gets the degrees of freedom of the first node (x and y);
        dof2 + i gets the degrees of freedom of the second node (x and y).
    */
    int dof1 = getNode1()->getDOF(0)->getIndex();
    int dof2 = getNode2()->getDOF(0)->getIndex();

    std::cout << "Kelem:" << std::endl;
    std::cout << localStiff << std::endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            GlobalStiff(dof1 + i, dof1 + j) += localStiff(i, j);
            GlobalStiff(dof1 + i, dof2 + j) += localStiff(i, j + 2);
            GlobalStiff(dof2 + i, dof1 + j) += localStiff(i + 2, j);
            GlobalStiff(dof2 + i, dof2 + j) += localStiff(i + 2, j + 2);
        }
}