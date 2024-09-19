#include "../headers/Element.h"

Element::Element() {}
Element::Element(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : index(_index), elemDimension(_elemDimension), elemConnectivity(_elemConnectivity), material(_material), physicalEntity(_physicalEntity) {}
Element::~Element() {}
