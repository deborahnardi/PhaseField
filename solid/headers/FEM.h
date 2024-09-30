#pragma once

#include <iomanip> // Para std::setw e std::fixed
#include "hdf5.h"
#include <petscksp.h>
#include <metis.h>

#include "Node.h"
#include "Element.h"
#include "DOF.h"
#include "Material.h"
#include "BoundaryElement.h"

#include "../../enumclass.hpp"

class FEM
{
private:
    int numNodes, nDOFs, numDirichletDOFs = 0, numNeumannDOFs = 0;
    int rank, size;
    std::string name, filename;
    std::vector<Material *> materials;
    std::vector<Node *> nodes, partitionedNodes;
    std::vector<Element *> elements, partitionedElements, partitionedBoundaryElements;
    std::vector<BoundaryElement *> bdElements;
    std::vector<DOF *> globalDOFs;

    MatrixXd K;
    VectorXd F;
    VectorXd U;

    Mat matrix;
    Vec rhs, solution;
    PetscInt Istart, Iend, IIstart, IIend;
    PetscInt *dirichletBC;
    PetscErrorCode ierr;
    bool showMatrix = false;

public:
    FEM();
    FEM(const std::string _name); // Only 2D problems are supported
    ~FEM();

    std::string getName() const { return name; }
    std::vector<Node *> getNodes() const { return nodes; }
    std::vector<Element *> getElements() const { return elements; }
    std::vector<DOF *> getGlobalDOFs() const { return globalDOFs; }

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
    void showResults();

    /*
                        SOLVE FEM PROBLEM METHODS
    */
    void solveFEMProblemNoPetsc();
    void assembleProblemNoPetsc();
    void setBoundaryConditionsNoPetsc();
    void solveLinearSystemNoPetsc();
    /*----------------------------------------------------------------------------------
                                    PETSc Methods
    ------------------------------------------------------------------------------------
    */
    PetscErrorCode solveFEMProblem();
    PetscErrorCode assembleProblem();
    PetscErrorCode createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo);
    PetscErrorCode setBoundaryConditions();
    PetscErrorCode solveLinearSystem(Mat &A, Vec &b, Vec &x);
    PetscErrorCode printGlobalMatrix(Mat &A);
    void setPrintMatrix(const bool &_showMatrix) { showMatrix = _showMatrix; }
};