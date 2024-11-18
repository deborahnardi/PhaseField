#pragma once

#include <set>
#include <iomanip> // Para std::setw e std::fixed
#include "hdf5.h"
#include <petscksp.h>
#include <petscmat.h>
#include <metis.h>
#include <functional>

#include "Node.h"
#include "Element.h"
#include "DOF.h"
#include "Material.h"
#include "BoundaryElement.h"
#include "AnalysisParameters.h"

#include "../../enumclass.hpp"

class FEM
{
private:
    int numNodes = 0, numElements = 0, nDOFs = 0, numDirichletDOFs = 0, numNeumannDOFs = 0, numElNodes = 0, elemDim = 0, numOfPrescribedDisp = 0;
    int rank, size;
    std::string name, filename, resultsPath;
    double *finalDisplacements, norm = 0., res = 0.;
    std::vector<std::set<int>> nodeNeighbours;
    std::vector<Material *> materials;
    std::vector<Node *> nodes, partitionedNodes, discritizedNodes;
    std::vector<Element *> elements, partitionedElements, partitionedBoundaryElements;
    std::vector<BoundaryElement *> bdElements;
    std::vector<DOF *> globalDOFs;
    AnalysisParameters *params;
    bool negativeLoad = false;

    MatrixXd K;
    VectorXd F;
    VectorXd U;

    Mat matrix, matrixPF;
    const PetscInt *JC, *IR; // ia contais the row pointer (equivalent to JC), ja contais the column index (equivalent to IR)
    PetscScalar *PA;
    double *DdkMinus1, *Ddk;
    Vec rhs, solution, rhsPF, solutionPF;
    PetscInt Istart, Iend, IIstart, IIend, IIIstart, IIIend, IstartPF, IendPF;
    PetscInt *d_nnz, *o_nnz;
    PetscInt *dirichletBC;
    PetscErrorCode ierr;
    bool showMatrix = false;

    std::vector<double> load;

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
    void setAnalysisParameters(AnalysisParameters *_params) { params = _params; }

    void readGeometry(const std::string &_filename);
    void removeNonDiscritizedNodes(std::vector<Node *> &_nodes);
    void renumberNodesIndexes(std::vector<Node *> &_nodes);
    void createResultsPath();
    void deleteResults(bool _deleteResults);
    void setLoadingVector(double ubar, int nSteps);

    std::vector<double> getLoadingVector() { return load; }
    /*----------------------------------------------------------------------------------
                                 FEM PROBLEM METHODS
    ------------------------------------------------------------------------------------
    */
    void findNeighbours();
    void updateBoundaryValues(double _lambda);
    void solveFEMProblemNoPetsc();
    void assembleProblemNoPetsc();
    void setBoundaryConditionsNoPetsc();
    void solveLinearSystemNoPetsc();
    double computeNorm(const double *vec1, const double *vec2, const int &size);
    /*----------------------------------------------------------------------------------
                                    Phase Field Methods
    ------------------------------------------------------------------------------------
    */
    void matrixPreAllocationPF(PetscInt start, PetscInt end);
    void solveDisplacementField(int _iStep);
    PetscErrorCode solvePhaseField();
    PetscErrorCode assemblePhaseFieldProblem();
    void solvePhaseFieldProblem();
    PetscErrorCode solveSystemByPSOR(Mat &A, Vec &b, Vec &x);
    PetscErrorCode getPSORVecs(Mat &A, Vec &b);
    void staggeredAlgorithm(int _iStep);
    PetscErrorCode updateFieldVariables(Vec &x, bool _hasConverged = true);
    /*----------------------------------------------------------------------------------
                                    PETSc Methods
    ------------------------------------------------------------------------------------
    */
    PetscErrorCode decomposeElements(Vec &b, Vec &x);
    PetscErrorCode matrixPreAllocation(PetscInt start, PetscInt end);
    PetscErrorCode solveFEMProblem();
    PetscErrorCode assembleProblem();
    PetscErrorCode createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo);
    PetscErrorCode solveLinearSystem(Mat &A, Vec &b, Vec &x);
    PetscErrorCode printGlobalMatrix(Mat &A);
    PetscErrorCode cleanSolution(Vec &x, Vec &b, Mat &A);
    /*----------------------------------------------------------------------------------
                                    OUTPUT METHODS
    ------------------------------------------------------------------------------------
    */
    void showResults(int _nStep);
    void setPrintMatrix(const bool &_showMatrix) { showMatrix = _showMatrix; }
    void updateVariables(Vec &x, bool _hasConverged = true);
    void deleteFromString(std::string &fullStr, std::string removeStr);
    void writeInHDF5(hid_t &file, std::fstream &output_v, herr_t &status, hid_t &dataset, hid_t &dataspace, std::string AttributeName, std::string AttributeType, double *valueVector, hsize_t valueVectorDims[], std::string s1);

    /*----------------------------------------------------------------------------------
                                        FUNCTIONS
    ------------------------------------------------------------------------------------
    */

    /*
         In setBoundaryFunction, the function is passed as a parameter.
         In updateBoundaryFunction, the function is called, and the arguments are passed to the function.
    */
public:
    void setBoundaryFunction(std::function<void(const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)> bFunc) { boundaryFunction = bFunc; }
    std::function<void(const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)> &getBoundaryFunction() { return boundaryFunction; }

    void updateBoundaryFunction(double _time);

private:
    std::function<void(const std::vector<double> &coord, const double &pseudoTime, DOF *dof, const std::vector<double> &load)> boundaryFunction = 0;
};