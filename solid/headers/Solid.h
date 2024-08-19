#pragma once

#include "hdf5.h"
#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <metis.h>

#include "Node.h"

class Solid
{
private:
    std::string name;
    std::string filename;
    int numNodes, num2DElements;
    std::vector<Node *> nodes;

public:
    Solid();
    Solid(const std::string _name);
    ~Solid();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }

    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }

    void readGeometry(const std::string &_filename);
};