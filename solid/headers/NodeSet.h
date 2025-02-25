#pragma once

#include "Node.h"

class NodeSet
{
private:
    int index;
    std::string name;
    std::vector<Node *> nodes;

public:
    NodeSet();
    NodeSet(const std::string _name, const std::vector<Node *> _nodes);
    ~NodeSet();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }
    Node *getNode(const int &_index) const { return nodes[_index]; }

    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }
    void setNode(const int &_index, Node *_node) { nodes[_index] = _node; }

    void addCondition(BoundaryType _bType, DOFType _dofType, double _value);
};