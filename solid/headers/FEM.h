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
#include "DOF.h"

#include "../../enumclass.hpp"

class FEM
{
private:
    std::string name, nSetName, elSetName, abaqusElementType;
    std::string filename;
    int numNodes, num2DElements, numBoundaryElements = 0, nDOFs, numDirichletDOFs = 0, numNeumannDOFs = 0;
    int problemDimension;
    std::vector<Node *> nodes;
    std::vector<Truss *> trussElements;
    std::vector<Solid2D *> elements;
    std::vector<BoundaryElement *> boundaryElements;
    std::vector<NodeSet *> nodeSets;
    std::vector<ElementSet *> elementSets;
    std::vector<DOF *> globalDOFs;

public:
    FEM();
    FEM(const std::string _name, const int &_problemDimension = 2); // Only 2D problems are supported
    ~FEM();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }
    std::vector<Solid2D *> getElements() const { return elements; }

    /*
                        DATA INPUT METHODS
    */
    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }
    void setElements(const std::vector<Solid2D *> &_elements) { elements = _elements; }

    void readGeometry(const std::string &_filename);
    void readGeometryCPS3(std::ifstream &file);
    void readGeometryT3D2(std::ifstream &file);
    void removeNonDiscritizedNodes(std::vector<Node *> &_nodes);
    void renumberNodesIndexes(std::vector<Node *> &_nodes);

    /*
                        SOLVE FEM PROBLEM METHODS
    */

    void assembleProblem();
};