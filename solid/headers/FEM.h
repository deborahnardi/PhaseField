#pragma once

#include <set>
#include <iomanip> // Para std::setw e std::fixed
#include "hdf5.h"
#include <petscksp.h>
#include <petscmat.h>
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
    int numNodes = 0, numElements = 0, nDOFs = 0, numDirichletDOFs = 0, numNeumannDOFs = 0, numElNodes = 0, elemDim = 0;
    int rank, size;
    std::string name, filename, resultsPath;
    std::vector<double> finalDisplacements;
    std::vector<std::set<int>> nodeNeighbours;
    std::vector<Material *> materials;
    std::vector<Node *> nodes, partitionedNodes, discritizedNodes;
    std::vector<Element *> elements, partitionedElements, partitionedBoundaryElements;
    std::vector<BoundaryElement *> bdElements;
    std::vector<DOF *> globalDOFs;

    MatrixXd K;
    VectorXd F;
    VectorXd U;

    Mat matrix;
    Vec rhs, solution;
    PetscInt Istart, Iend, IIstart, IIend, IIIstart, IIIend;
    PetscInt *d_nnz, *o_nnz;
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
    std::string getResultsPath() const { return resultsPath; }

    /*----------------------------------------------------------------------------------
                                    DATA INPUT METHODS
    ------------------------------------------------------------------------------------
    */
    void setName(const std::string _name) { name = _name; }
    void setNodes(const std::vector<Node *> &_nodes) { nodes = _nodes; }
    void setResultsPath() { resultsPath = "./../output/" + name + "/"; }

    void readGeometry(const std::string &_filename);
    void removeNonDiscritizedNodes(std::vector<Node *> &_nodes);
    void renumberNodesIndexes(std::vector<Node *> &_nodes);
    void showResults();
    void createResultsPath();
    void deleteResults(bool _deleteResults);

    /*----------------------------------------------------------------------------------
                            SOLVE FEM PROBLEM METHODS
    ------------------------------------------------------------------------------------
    */
    void findNeighbours();
    void solveFEMProblemNoPetsc();
    void assembleProblemNoPetsc();
    void setBoundaryConditionsNoPetsc();
    void solveLinearSystemNoPetsc();
    /*----------------------------------------------------------------------------------
                                    PETSc Methods
    ------------------------------------------------------------------------------------
    */
    PetscErrorCode decomposeElements(Vec &b, Vec &x);
    PetscErrorCode matrixPreAllocation();
    PetscErrorCode solveFEMProblem();
    PetscErrorCode assembleProblem();
    PetscErrorCode createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo);
    PetscErrorCode solveLinearSystem(Mat &A, Vec &b, Vec &x);
    PetscErrorCode printGlobalMatrix(Mat &A);
    /*----------------------------------------------------------------------------------
                                    OUTPUT METHODS
    ------------------------------------------------------------------------------------
    */
    void setPrintMatrix(const bool &_showMatrix) { showMatrix = _showMatrix; }
    void deleteFromString(std::string &fullStr, std::string removeStr);
    void writeInHDF5(hid_t &file, std::fstream &output_v, herr_t &status, hid_t &dataset, hid_t &dataspace, std::string AttributeName, std::string AttributeType, double *valueVector, hsize_t valueVectorDims[], std::string s1);
};