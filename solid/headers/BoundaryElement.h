#pragma once

#include "Node.h"
#include "DenseEigen.h"
#include "Material.h"
#include "ShapeFunction.h"
#include "Quadrature.h"
#include "../../enumclass.hpp"
#include "AnalysisParameters.h"

class BoundaryElement
{
private:
    int index, elemDimension, physicalEntity, numQuadraturePoints, numBdNodes;
    std::vector<Node *> elemConnectivity;
    Material *material;
    ShapeFunction *sF;
    Quadrature *q;
    PetscErrorCode ierr;
    AnalysisParameters *params;

    struct BoundaryCondition
    {
        BoundaryType bdType;
        std::vector<DOF *> dofs;
        double value;
    };

    std::vector<BoundaryCondition> conditions;

public:
    BoundaryElement();
    BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, AnalysisParameters *_params);
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
    void updateBoundaryValues(double _lambda);

    /*
        PETSc Methods
    */
    PetscErrorCode getContribution(Vec &rhs);
};