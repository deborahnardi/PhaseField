#include "../headers/ElementSet.h"

ElementSet::ElementSet() {};
ElementSet::ElementSet(const std::string &_name, const std::vector<Element *> &_elements)
    : name(_name), elements(_elements) {};
ElementSet::~ElementSet() {};