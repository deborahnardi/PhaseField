#include "../headers/Node.h"

Node::Node() {}
Node::Node(const int &_index, const std::vector<double> &_initialCoordinates)
    : index(_index), initialCoordinates(_initialCoordinates) {}
Node::~Node() {}