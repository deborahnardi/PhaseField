#include "../headers/Node.h"

Node::Node() {}
Node::Node(const int &_index, const std::vector<double> &_initialCoordinates)
    : index(_index), initialCoordinates(_initialCoordinates) {}
Node::~Node() {}

void Node::addDOF(DOF *_dof)
{
    for (const auto &existingDOF : dofs)
        if (existingDOF->getDOFType() == _dof->getDOFType()) // For avoiding the addition of the same DOF twice
            return;

    dofs.push_back(_dof);
}

void Node::addInverseIncidence(int el)
{
    if (std::any_of(inverseIncidence.begin(), inverseIncidence.end(), [el](int existingEl)
                    { return existingEl == el; })) // Check if element already exists
        return;                                    // element already exists, no need to add it again

    inverseIncidence.push_back(el);
}