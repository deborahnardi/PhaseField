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

PetscErrorCode Truss::getContribution(Mat &matrix)
{
    PetscInt numElDOF = 4;
    PetscInt *indx = new PetscInt[4]();
    PetscScalar *localStiff = new PetscScalar[16]();
    PetscInt count = 0;
    PetscScalar k = material->getYoungModulus() * area / length;

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
        {
            indx[count] = dof->getIndex();
            count++;
        }

    ierr = MatSetValues(matrix, numElDOF, indx, numElDOF, indx, localStiff, ADD_VALUES);
    CHKERRQ(ierr);

    delete[] indx;
    delete[] localStiff;
}

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