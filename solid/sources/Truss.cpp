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

    for (auto node : elemConnectivity)
        for (auto dof : node->getDOFs())
        {
            indx[count] = dof->getIndex();
            count++;
        }

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

    K = rotationMatrix.transpose() * localStiffnessMatrix * rotationMatrix;
}