#pragma once

#include "Node.h"
#include "DenseEigen.h"

class Element
{
protected:
    int index, elemDimension;
    std::vector<Node *> elemConnectivity;

public:
    Element();
    Element(const int &_index, const std::vector<Node *> &_elemConnectivity, const int &_elemDimension);
    ~Element();

    int getIndex() const { return index; }
    int getElemDimension() const { return elemDimension; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &_index) const { return elemConnectivity[_index]; }

    void setIndex(const int &_index) { index = _index; }
    void setPhysicalEntity(const int &_elemDimension) { elemDimension = _elemDimension; }
    void setElemConnectivity(const std::vector<Node *> &_elemConnectivity) { elemConnectivity = _elemConnectivity; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }

    virtual void getContribution() {};
};

class BoundaryElement : public Element
{
private:
    std::string entityName;

    struct AppliedBoundaryCondition
    {
        BoundaryType bType;
        std::vector<DOF *> dofs;
    };

    std::vector<AppliedBoundaryCondition> appliedBCs;

public:
    BoundaryElement();
    BoundaryElement(const int &index, const std::vector<Node *> &elemConnectivity, const int &elemDimension = 1);
    ~BoundaryElement();

    int getIndex() const { return index; }
    int getElemDimension() const { return 1; }
    std::string getEntityName() const { return "Boundary_" + std::to_string(index); }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }

    void getContribution() override {};
};

class Truss : public Element
{
private:
    double length, area, theta;
    MatrixXd localStiffnessMatrix, rotationMatrix, K;
    Node *node1, *node2;

public:
    Truss();
    Truss(const int &index, Node *_node1, Node *_node2);
    ~Truss();

    int getIndex() const { return index; }
    int getLength() const { return length; }
    int getArea() const { return area; }
    int getTheta() const { return theta; }
    Node *getNode1() const { return elemConnectivity[0]; }
    Node *getNode2() const { return elemConnectivity[1]; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }
    MatrixXd getLocalStiffnessMatrix() const { return localStiffnessMatrix; }
    MatrixXd getRotationMatrix() const { return rotationMatrix; }
    MatrixXd getElemStiffnessMatrix() const { return K; }

    void setIndex(const int &_index) { index = _index; }
    void setLength(const double &_length) { length = _length; }
    void setArea(const double &_area) { area = _area; }
    void setTheta(const double &_theta) { theta = _theta; }
    void setNode(const int &_index, Node *_node) { elemConnectivity[_index] = _node; }
    void setLocalStiffnessMatrix(const MatrixXd &_localStiffnessMatrix) { localStiffnessMatrix = _localStiffnessMatrix; }
    void setRotationMatrix(const MatrixXd &_rotationMatrix) { rotationMatrix = _rotationMatrix; }
    void setElemStiffnessMatrix(const MatrixXd &_K) { K = _K; }

    void getContribution() override;
};

class Solid2D : public Element
{
private:
    double area;
    MatrixXd K;

public:
    Solid2D();
    Solid2D(const int &index, const std::vector<Node *> &elemConnectivity, const int &elemDimension = 2);
    ~Solid2D();

    int getIndex() const { return index; }
    int getElemDimension() const { return 2; }
    double getArea() const { return area; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }
    std::string getEntityName() const { return "Solid2D_" + std::to_string(index); }
    void setArea(const double &_area) { area = _area; }

    MatrixXd getElemStiffnessMatrix() const { return K; }
    void getContribution() override {};
};
