#pragma once

#include <chrono>

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

using Tensor = std::array<std::array<double, 3>, 3>;

class FEM
{
private:
    int numNodes = 0, numElements = 0, nDOFs = 0, numDirichletDOFs = 0, numNeumannDOFs = 0, numElNodes = 0, elemDim = 0, numOfPrescribedDisp = 0;
    int rank, size, nDOF = 2;
    double norm = 0., res = 0.;
    std::string name, filename, resultsPath;
    std::vector<Material *> materials;
    std::vector<Node *> nodes, partitionedNodes, discritizedNodes;
    std::vector<Element *> elements, partitionedElements, partitionedBoundaryElements;
    std::vector<BoundaryElement *> bdElements, tractionBd;
    std::vector<DOF *> globalDOFs;
    std::vector<double> load;
    AnalysisParameters *params;
    bool negativeLoad = false, prescribedDamageField = false, showMatrix = false;
    PetscLogDouble bytes = 0.0;
    PetscInt totalNnz = 0;

    Mat matrix, matrixPF, matrixCopy;
    Vec rhs, solution, rhsPF, solutionPF, disp, nodalForces, reactionForces;
    PetscInt Istart, Iend, IstartBD, IendBD, IstartPF, IendPF;
    PetscInt *d_nnz, *o_nnz, *d_nz, *o_nz, *dirichletBC;
    PetscErrorCode ierr;
    int *JC, *IR, nzQ = 0;
    double *PA;
    double *DdkMinus1, *Ddk, *totalVecq, **totalMatrixQ;

    std::vector<int> n2nUpperTotal, n2nCSRUpperTotal, nodesForEachRankCSR, numNodesForEachRank, n2nDRankUpper, n2nDRankLower, n2nCSRUpper;
    std::vector<std::set<int>> n2nUpperMatTotal, n2e, n2nMat, n2nUpperMat, n2nLowerMat;
    std::vector<std::vector<int>> eSameList;

    std::vector<std::set<int>> n2nMatTotal, n2nLowerMatTotal;
    std::vector<int> n2nTotal, n2nLowerTotal;

    MatrixXd K;
    VectorXd F;
    VectorXd U;

    double force[2] = {};

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
    void setLoadingVector1(double ubar, int nSteps);
    void setLoadingVector2(double ubar, int nSteps);
    void setLoadingVector3(double ubar, int nSteps);

    std::vector<double> getLoadingVector() { return load; }
    void setPrescribedDamageField(bool _prescribedDamageField) { prescribedDamageField = _prescribedDamageField; }
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
    PetscErrorCode printMemoryUsage(const int &iStep);
    PetscErrorCode computeReactionForces();
    double computeNorm(const double *vec1, const double *vec2, const int &size);
    std::array<Tensor, 3> computeConstitutiveTensors();
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
    PetscErrorCode getPSORVecs();
    void staggeredAlgorithm(int _iStep);
    PetscErrorCode updateFieldVariables(Vec &x, bool _hasConverged = true);
    PetscErrorCode updateFieldDistribution();
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
    PetscErrorCode assembleBetweenProcesses(Mat &A, Vec &b);
    PetscErrorCode assembleSymmStiffMatrix(Mat &A);
    PetscErrorCode updateRHS(Mat &A, Vec &b);
    /*----------------------------------------------------------------------------------
                                    OUTPUT METHODS
    ------------------------------------------------------------------------------------
    */
    void showResults(int _nStep);
    void setPrintMatrix(const bool &_showMatrix) { showMatrix = _showMatrix; }
    void updateVariables(Vec &x, bool _hasConverged = true);
    void deleteFromString(std::string &fullStr, std::string removeStr);
    void writeInHDF5(hid_t &file, std::fstream &output_v, herr_t &status, hid_t &dataset, hid_t &dataspace, std::string AttributeName, std::string AttributeType, double *valueVector, hsize_t valueVectorDims[], std::string s1);

    double elapsedTime(std::chrono::_V2::system_clock::time_point t1, std::chrono::_V2::system_clock::time_point t2);

    /*----------------------------------------------------------------------------------
                                POST PROCESSING METHODS
  ------------------------------------------------------------------------------------
  */
    void postProc();
    void computeNodalStress();
    void computeTractions();
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