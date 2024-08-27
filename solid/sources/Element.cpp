#include "../headers/Element.h"

Element::Element() {}
Element::Element(const int &_index, const std::vector<Node *> &_elemConnectivity, const int &_elemDimension)
    : index(_index), elemConnectivity(_elemConnectivity), elemDimension(_elemDimension) {}
Element::~Element() {}

Solid2D::Solid2D() {}
Solid2D::Solid2D(const int &index, const std::vector<Node *> &elemConnectivity, const int &elemDimension)
    : Element(index, elemConnectivity, 2) {}
Solid2D::~Solid2D() {}

BoundaryElement::BoundaryElement() {}
BoundaryElement::BoundaryElement(const int &index, const std::vector<Node *> &elemConnectivity, const int &elemDimension)
    : Element(index, elemConnectivity, 1) {}
BoundaryElement::~BoundaryElement() {}

void BoundaryElement::addCondition(BoundaryType _bType, DOFType _dofType, double _value) // Add the boundary condition to all the nodes of the bundary element
{
    std::vector<DOF *> dofVec;
    for (auto node : elemConnectivity)
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

    appliedBCs.push_back({_bType, dofVec});
}