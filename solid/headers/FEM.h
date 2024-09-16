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
    int numNodes, nDOFs, numDirichletDOFs = 0, numNeumannDOFs = 0;
    int rank, size;
    std::string name, filename;
    std::vector<Material *> materials;
    std::vector<Node *> nodes, partitionedNodes;
    std::vector<Element *> elements, bdElements, partitionedElements, partitionedBoundaryElements;
    std::vector<DOF *> globalDOFs, partitionedDOFs;

    MatrixXd K;
    VectorXd F;
    VectorXd U;

public:
    FEM();
    FEM(const std::string _name); // Only 2D problems are supported
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
    void decomposeElements();
    void matrixPreAllocation();

    /*
                        SOLVE FEM PROBLEM METHODS
    */
    void solveFEMProblem();
    void assembleProblem();
    void setBoundaryConditions();
    void solveLinearSystem();

    void solveFEMProblemPETSc();
    void assembleProblemPETSc();
    void setBoundaryConditionsPETSc() {};
    void solveLinearSystemPETSc() {};
};