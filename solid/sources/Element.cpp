#include "../headers/Element.h"

Element::Element() {}
Element::Element(const int &_index, const std::vector<Node *> &_elemConnectivity, const int &_elemDimension)
    : index(_index), elemConnectivity(_elemConnectivity), elemDimension(_elemDimension) {}
Element::~Element() {}