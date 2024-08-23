#pragma once

#include "Element.h"

class ElementSet
{
private:
    std::string name;
    std::vector<Element *> elements;

public:
    ElementSet();
    ElementSet(const std::string &_name, const std::vector<Element *> &_elements);
    ~ElementSet();
};