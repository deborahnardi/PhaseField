#pragma once

#include "Node.h"
#include "DenseEigen.h"
#include "Material.h"
#include "ShapeFunction.h"
#include "Quadrature.h"
#include "../../enumclass.hpp"
#include "AnalysisParameters.h"

using Tensor = std::array<std::array<double, 3>, 3>;

class Element
{
protected:
    int index, elemDimension, physicalEntity;
    std::vector<Node *> elemConnectivity;
    Material *material;
    BoundaryType bdType;
    DOFType type;
    double value, epsilon = 0.;
    double **coords, *weights;
    MatrixXd K, localStiff;
    PetscErrorCode ierr;
    AnalysisParameters *params;

public:
    Element();
    Element(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, AnalysisParameters *_params);
    ~Element();

    int getIndex() const { return index; }
    int getElemDimension() const { return elemDimension; }
    int getPhysicalEntity() const { return physicalEntity; }
    double getDeformation() { return epsilon; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &_index) const { return elemConnectivity[_index]; }

    void setIndex(const int &_index) { index = _index; }
    void setPhysicalEntity(const int &_elemDimension) { elemDimension = _elemDimension; }
    void setElemConnectivity(const std::vector<Node *> &_elemConnectivity) { elemConnectivity = _elemConnectivity; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }
    void setDeformation(const double &_epsilon) { epsilon = _epsilon; }

    virtual MatrixXd getElemStiffnessMatrix() const = 0; // const 0 means that the function has no implementation in the base class
    virtual void assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff) {};
    virtual void addCondition(BoundaryType _bdType, DOFType _type, double _value) {};

    virtual PetscErrorCode getContribution(Mat &matrix, Vec &rhs, bool negativeLoad = false, bool _PrescribedDamageField = false) {};
    virtual std::vector<double> getStiffnessIIOrIJ(std::array<Tensor, 3> tensors, const int idxLocalNode1, const int idxLocalNode2, bool _PrescribedDamageField = false) {};
    virtual double getQValue(const int idxLocalNode1, const int idxLocalNode2) {};
    virtual double stiffnessValue(const int localPos1, const int localPos2, Tensor &tensorC, const PetscReal B[3][6]) {};
    virtual PetscErrorCode getPhaseFieldContribution(Mat &A, Vec &rhs, bool _PrescribedDamageField = false) {};
    virtual void getContribution() {};
    virtual void Test(PetscScalar &integral) {};
    virtual void computeDeformation() {};

    // =======================================
    // =========== POST PROCESSING ===========
    // =======================================
    virtual PetscErrorCode calculateStress() {};
};

class Truss : public Element
{
private:
    int numHammerPoints = 2, numElNodes = 2;
    double length, area, theta;
    MatrixXd localStiffnessMatrix, rotationMatrix;
    ShapeFunction *sF;
    Quadrature *q;

public:
    Truss();
    Truss(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, const double &area_, AnalysisParameters *_params);
    ~Truss();

    int getLength() const { return length; }
    int getArea() const { return area; }
    int getTheta() const { return theta; }
    Node *getNode1() const { return elemConnectivity[0]; }
    Node *getNode2() const { return elemConnectivity[1]; }
    MatrixXd getLocalStiffnessMatrix() const { return localStiffnessMatrix; } // Without rotation
    MatrixXd getRotationMatrix() const { return rotationMatrix; }
    MatrixXd getElemStiffnessMatrix() const override { return localStiff; }
    void assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff) override;

    void setIndex(const int &_index) { index = _index; }
    void setLength(const double &_length) { length = _length; }
    void setArea(const double &_area) { area = _area; }
    void setTheta(const double &_theta) { theta = _theta; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }

    PetscErrorCode getContribution(Mat &matrix, Vec &rhs, bool negativeLoad = false, bool _PrescribedDamageField = false) override;
    PetscErrorCode getPhaseFieldContribution(Mat &matrix, Vec &rhs, bool _PrescribedDamageField = false) override;
    void getContribution() override;
    void computeDeformation() override;
};

class Solid2D : public Element
{
private:
    int numHammerPoints = 3, numElNodes = 3;
    double area;
    ShapeFunction *sF;
    Quadrature *q;

public:
    Solid2D();
    Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, AnalysisParameters *_params);
    ~Solid2D();

    double getArea() const { return area; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }
    void setArea(const double &_area) { area = _area; }
    void assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff) override;

    MatrixXd getElemStiffnessMatrix() const override { return localStiff; }
    PetscErrorCode getContribution(Mat &matrix, Vec &rhs, bool negativeLoad = false, bool _PrescribedDamageField = false) override;
    std::vector<double> getStiffnessIIOrIJ(std::array<Tensor, 3> tensors, const int idxLocalNode1, const int idxLocalNode2, bool _PrescribedDamageField = false) override;
    double getQValue(const int idxLocalNode1, const int idxLocalNode2) override;
    double stiffnessValue(const int localPos1, const int localPos2, Tensor &tensorC, const PetscReal B[3][6]) override;
    PetscErrorCode getPhaseFieldContribution(Mat &matrix, Vec &rhs, bool _PrescribedDamageField = false) override;
    void getContribution() override;
    void Test(PetscScalar &integral) override;

    // =======================================
    // =========== POST PROCESSING ===========
    // =======================================
    PetscErrorCode calculateStress() override;
};
