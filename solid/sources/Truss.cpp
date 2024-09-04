#include "../headers/Element.h"

Truss::Truss() {}
Truss::Truss(const int &index, Node *_node1, Node *_node2)
    : Element(index, elemConnectivity, 1)
{
    node1 = _node1;
    node2 = _node2;
    elemConnectivity.push_back(_node1);
    elemConnectivity.push_back(_node2);

    length = pow((_node2->getX() - _node1->getX()) * (_node2->getX() - _node1->getX()) +
                     (_node2->getY() - _node1->getY()) * (_node2->getY() - _node1->getY()),
                 0.5);

    theta = atan((_node2->getY() - _node1->getY()) / (_node2->getX() - _node1->getX()));

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