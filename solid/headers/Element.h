#pragma once

#include "Node.h"

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

    void addCondition(BoundaryType _bType, DOFType _dofType, double _value);
};

class Solid2D : public Element
{
private:
    double area;

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
};
