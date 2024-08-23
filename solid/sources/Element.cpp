#include "../headers/Element.h"

Element::Element() {}
Element::Element(const int &_index, const std::vector<Node *> &_elemConnectivity, const int &_elemDimension)
    : index(_index), elemConnectivity(_elemConnectivity), elemDimension(_elemDimension) {}
Element::~Element() {}

Solid2D::Solid2D() {}
Solid2D::Solid2D(const int &index, const std::vector<Node *> &elemConnectivity)
    : Element(index, elemConnectivity, 2) {}
Solid2D::~Solid2D() {}

BoundaryElement::BoundaryElement() {}
BoundaryElement::BoundaryElement(const int &index, const std::vector<Node *> &elemConnectivity)
    : Element(index, elemConnectivity, 1) {}
BoundaryElement::~BoundaryElement() {}