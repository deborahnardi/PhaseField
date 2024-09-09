#pragma once

#include "hdf5.h"
#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <metis.h>

#include "Node.h"
#include "Element.h"
#include "DOF.h"
#include "Material.h"

#include "../../enumclass.hpp"

class FEM
{
private:
    std::string name, filename;
    int numNodes, num2DElements, numBoundaryElements = 0, nDOFs, numDirichletDOFs = 0, numNeumannDOFs = 0;
    int problemDimension;
    std::vector<Material *> materials;
    std::vector<Node *> nodes;
    std::vector<Element *> elements, bdElements;
    std::vector<DOF *> globalDOFs;

public:
    FEM();
    FEM(const std::string _name, const int &_problemDimension = 2); // Only 2D problems are supported
    ~FEM();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }

    /*
                        DATA INPUT METHODS
    */
    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }

    void readGeometry(const std::string &_filename);
    void removeNonDiscritizedNodes(std::vector<Node *> &_nodes);
    void renumberNodesIndexes(std::vector<Node *> &_nodes);

    /*
                        SOLVE FEM PROBLEM METHODS
    */

    void assembleProblem();
};