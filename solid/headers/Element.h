#pragma once

#include "Node.h"
#include "DenseEigen.h"
#include "Material.h"
#include "ShapeFunction.h"
#include "Quadrature.h"
#include "../../enumclass.hpp"

class Element
{
protected:
    int index, elemDimension, physicalEntity;
    std::vector<Node *> elemConnectivity;
    Material *material;
    BoundaryType bdType;
    DOFType type;
    double value;
    MatrixXd K, localStiff;
    PetscErrorCode ierr;

public:
    Element();
    Element(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~Element();

    int getIndex() const { return index; }
    int getElemDimension() const { return elemDimension; }
    int getPhysicalEntity() const { return physicalEntity; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &_index) const { return elemConnectivity[_index]; }

    void setIndex(const int &_index) { index = _index; }
    void setPhysicalEntity(const int &_elemDimension) { elemDimension = _elemDimension; }
    void setElemConnectivity(const std::vector<Node *> &_elemConnectivity) { elemConnectivity = _elemConnectivity; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }

    virtual MatrixXd getElemStiffnessMatrix() const = 0; // const 0 means that the function has no implementation in the base class
    virtual void assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff) {};
    virtual void addCondition(BoundaryType _bdType, DOFType _type, double _value) {};

    virtual PetscErrorCode getContribution(Mat &matrix, Vec &rhs) {};
    virtual void getContribution() {};
    virtual void Test(PetscScalar &integral) {};
};

class Truss : public Element
{
private:
    double length, area, theta;
    MatrixXd localStiffnessMatrix, rotationMatrix;

public:
    Truss();
    Truss(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, const double &area_);
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

    PetscErrorCode getContribution(Mat &matrix, Vec &rhs) override;
    void getContribution() override;
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
    Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~Solid2D();

    double getArea() const { return area; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }
    void setArea(const double &_area) { area = _area; }
    void assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff) override;

    MatrixXd getElemStiffnessMatrix() const override { return localStiff; }
    PetscErrorCode getContribution(Mat &matrix, Vec &rhs) override;
    void getContribution() override;
    void Test(PetscScalar &integral) override;
};
