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

    std::string getName() const { return name; }
    std::vector<Element *> getElements() const { return elements; }
    Element *getElement(const int &_index) const { return elements[_index]; }

    void setName(const std::string &_name) { name = _name; }
    void setElements(const std::vector<Element *> &_elements) { elements = _elements; }
    void setElement(const int &_index, Element *_element) { elements[_index] = _element; }
};