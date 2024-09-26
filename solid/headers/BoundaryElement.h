#pragma once

#include "Node.h"
#include "DenseEigen.h"
#include "Material.h"
#include "../../enumclass.hpp"

class BoundaryElement
{
private:
    int index, elemDimension, physicalEntity, numBdNodes, numNeumannDOFs = 0;
    std::vector<Node *> elemConnectivity;
    Material *material;

    struct BoundaryCondition
    {
        BoundaryType bdType;
        std::vector<DOF *> dofs;
        double value;
    };

    std::vector<BoundaryCondition> conditions;

public:
    BoundaryElement();
    BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity);
    ~BoundaryElement();

    int getIndex() const { return index; }
    int getElemDimension() const { return elemDimension; }
    int getNumOfConditions() const { return conditions.size(); }
    int getPhysicalEntity() const { return physicalEntity; }
    std::vector<Node *> getElemConnectivity() const { return elemConnectivity; }
    BoundaryType getConditionType(const int &_index) const { return conditions[_index].bdType; }
    std::vector<DOF *> getDOFs(const int &_index) const { return conditions[_index].dofs; }

    void addCondition(BoundaryType _bdType, DOFType _type, double _value);
    void getContributionNoPetsc(VectorXd &F, MatrixXd &K);
    /*
        PETSc Methods
    */
    void getContribution(Vec &rhs);
};