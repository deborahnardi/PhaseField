#pragma once

#include "Node.h"

class Element
{
protected:
    int index, elemDimension;
    std::vector<Node *> elemConnectivity;

public:
    Element();
    Element(const int &_index, const std::vector<Node *> &_elemConnectivity, const int &_elemDimension);
    ~Element();

    int getIndex() const { return index; }
    int getPhysicalEntity() const { return elemDimension; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &_index) const { return elemConnectivity[_index]; }

    void setIndex(const int &_index) { index = _index; }
    void setPhysicalEntity(const int &_elemDimension) { elemDimension = _elemDimension; }
    void setElemConnectivity(const std::vector<Node *> &_elemConnectivity) { elemConnectivity = _elemConnectivity; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }
};