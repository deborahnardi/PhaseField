#pragma once

#include "Node.h"
#include "DenseEigen.h"
#include "Material.h"

class Element
{
protected:
    int index, elemDimension, physicalEntity;
    std::vector<Node *> elemConnectivity;
    Material *material;

public:
    Element();
    Element(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
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
public:
    BoundaryElement();
    BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~BoundaryElement();

    void getContribution() override {};
};

class Truss : public Element
{
private:
    double length, area, theta;
    MatrixXd localStiffnessMatrix, rotationMatrix, K;

public:
    Truss();
    Truss(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, const double &area_);
    ~Truss();

    int getLength() const { return length; }
    int getArea() const { return area; }
    int getTheta() const { return theta; }
    Node *getNode1() const { return elemConnectivity[0]; }
    Node *getNode2() const { return elemConnectivity[1]; }
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
    Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~Solid2D();

    double getArea() const { return area; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    Node *getNode(const int &index) const { return elemConnectivity[index]; }
    void setArea(const double &_area) { area = _area; }

    MatrixXd getElemStiffnessMatrix() const { return K; }
    void getContribution() override {};
};
