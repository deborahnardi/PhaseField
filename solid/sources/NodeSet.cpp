#include "../headers/NodeSet.h"

NodeSet::NodeSet() {}
NodeSet::NodeSet(const std::string _name, const std::vector<Node *> _nodes)
    : name(_name), nodes(_nodes) {}
NodeSet::~NodeSet() {}