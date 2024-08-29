#include "../headers/NodeSet.h"

NodeSet::NodeSet() {}
NodeSet::NodeSet(const std::string _name, const std::vector<Node *> _nodes)
    : name(_name), nodes(_nodes) {}
NodeSet::~NodeSet() {}

void NodeSet::addCondition(BoundaryType _bType, DOFType _dofType, double _value) // Add the boundary condition to all the nodes of the bundary element
{
    std::vector<DOF *> dofVec;
    for (auto node : nodes)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() == _dofType) // If the DOF is the same as the one in the boundary condition
            {
                if (_bType == DIRICHLET)
                {
                    dof->setDirichlet();
                    dof->setDirichletValue(_value);
                }
                else if (_bType == NEUMANN)
                {
                    dof->setNeumann();
                    dof->setNeumannValue(_value);
                }
                dofVec.push_back(dof);
            }

    // appliedBCs.push_back({_bType, dofVec});
}