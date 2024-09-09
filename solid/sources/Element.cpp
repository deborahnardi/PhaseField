#include "../headers/Element.h"

Element::Element() {}
Element::Element(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : index(_index), elemDimension(_elemDimension), elemConnectivity(_elemConnectivity), material(_material), physicalEntity(_physicalEntity) {}
Element::~Element() {}

// **************************

BoundaryElement::BoundaryElement() {}
BoundaryElement::BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : Element(_index, _elemDimension, _elemConnectivity, _material, _physicalEntity)
{
    for (auto n : _elemConnectivity)
    {
        n->setIsDiscritized();
        n->addDOF(new DOF(X, 0.));
        n->addDOF(new DOF(Y, 0.));
    }
}
BoundaryElement::~BoundaryElement() {}

void BoundaryElement::addCondition(BoundaryType _bdType, DOFType _type, double _value)
{
    bdType = _bdType;
    type = _type;
    value = _value;

    for (auto n : elemConnectivity)
        for (auto dof : n->getDOFs())
            if (dof->getDOFType() == type)
            {
                if (bdType = NEUMANN)
                {
                    dof->setNeumann();
                    dof->setNeumannValue(value);
                }
                else if (bdType == DIRICHLET)
                {
                    dof->setDirichlet();
                    dof->setDirichletValue(value);
                }
            }
}

// **************************

Solid2D::Solid2D() {}
Solid2D::Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : Element(_index, _elemDimension, _elemConnectivity, _material, _physicalEntity)
{
    for (auto n : _elemConnectivity)
    {
        n->setIsDiscritized();
        n->addDOF(new DOF(X, 0.));
        n->addDOF(new DOF(Y, 0.));
    }
}
Solid2D::~Solid2D() {}
