#pragma once

#include "hdf5.h"
#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <metis.h>

#include "Node.h"
#include "Element.h"
#include "NodeSet.h"
#include "ElementSet.h"

#include "../../enumclass.hpp"

class Solid
{
private:
    std::string name, nSetName, elSetName;
    std::string filename;
    int numNodes, num2DElements, numBoundaryElements = 0;
    std::vector<Node *> nodes;
    std::vector<Solid2D *> elements;
    std::vector<BoundaryElement *> boundaryElements;
    std::vector<NodeSet *> nodeSets;
    std::vector<ElementSet *> elementSets;

public:
    Solid();
    Solid(const std::string _name);
    ~Solid();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }
    std::vector<Solid2D *> getElements() const { return elements; }

    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }
    void setElements(const std::vector<Solid2D *> &_elements) { elements = _elements; }

    void readGeometry(const std::string &_filename);
};