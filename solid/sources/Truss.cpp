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

void Truss::getContribution(Mat &matrix, Vec &rhs)
{
    // localStiffnessMatrix = MatrixXd::Zero(4, 4);
    // localStiffnessMatrix << 1, 0, -1, 0,
    //     0, 0, 0, 0,
    //     -1, 0, 1, 0,
    //     0, 0, 0, 0;
    // localStiffnessMatrix *= material->getYoungModulus() * area / length;

    // rotationMatrix = MatrixXd::Zero(4, 4);
    // rotationMatrix << cos(theta), sin(theta), 0, 0,
    //     -sin(theta), cos(theta), 0, 0,
    //     0, 0, cos(theta), sin(theta),
    //     0, 0, -sin(theta), cos(theta);

    // K = rotationMatrix.transpose() * localStiffnessMatrix * rotationMatrix;

    // Translating the above code to PETSc
    PetscScalar localStiffnessMatrix[16] = {1, 0, -1, 0,
                                            0, 0, 0, 0,
                                            -1, 0, 1, 0,
                                            0, 0, 0, 0};
    for (int i = 0; i < 16; i++)
        localStiffnessMatrix[i] *= material->getYoungModulus() * area / length;

    PetscScalar rotationMatrix[16] = {
        cos(theta), sin(theta), 0, 0,
        -sin(theta), cos(theta), 0, 0,
        0, 0, cos(theta), sin(theta),
        0, 0, -sin(theta), cos(theta)};

    PetscScalar rotationMatrix[16] = {0};

    // Multiplying the rotation matrix by the local stiffness matrix
}