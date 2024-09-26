#include "../headers/BoundaryElement.h"

BoundaryElement::BoundaryElement() {}
BoundaryElement::BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : index(_index), elemDimension(_elemDimension), elemConnectivity(_elemConnectivity), material(_material), physicalEntity(_physicalEntity)
{
    numBdNodes = _elemConnectivity.size();

    for (auto n : _elemConnectivity)
    {
        n->setIsDiscritized();
        n->addDOF(new DOF(X, 0.));
        n->addDOF(new DOF(Y, 0.));
    }
}
BoundaryElement::~BoundaryElement() {}

void BoundaryElement::addCondition(BoundaryType _bdType, DOFType _type, double _value)
{
    std::vector<DOF *> dofVec;
    for (auto n : elemConnectivity)
        for (auto dof : n->getDOFs())
            if (dof->getDOFType() == _type)
            {
                if (_bdType == NEUMANN)
                {
                    dof->setNeumann();
                    dof->setNeumannValue(_value);
                }
                else if (_bdType == DIRICHLET)
                {
                    dof->setDirichlet();
                    dof->setDirichletValue(_value);
                }
                dofVec.push_back(dof);
            }
    conditions.push_back({_bdType, dofVec, _value}); // Each boundary element has its own conditions
}

void BoundaryElement::getContributionNoPetsc(VectorXd &F, MatrixXd &K)
{
    for (auto c : conditions)
        if (c.bdType == NEUMANN)
            if (elemDimension == 0)
            {
                F(c.dofs[0]->getIndex()) += c.value;
                numNeumannDOFs++;
            }
            else
            {
                for (auto dof : c.dofs)
                    if (dof->getIndex() == -1)
                        return;

                for (auto dof : c.dofs)
                    F(dof->getIndex()) += c.value;
            }
        else if (c.bdType == DIRICHLET)
            for (auto dof : c.dofs)
            {
                K.row(dof->getIndex()).setZero();
                K.col(dof->getIndex()).setZero();
                K(dof->getIndex(), dof->getIndex()) = 1; // Setting the diagonal to 1
                F(dof->getIndex()) = c.value;            // If a prescribed displacement value is given
            }
}

void BoundaryElement::getContribution(Vec &rhs)
{
    for (auto c : conditions)
        if (c.bdType == NEUMANN)
            if (elemDimension == 0)
            {
                PetscInt dofA = c.dofs[0]->getIndex(); // c.dofs[0] because if elemDimension == 0, it is only a node -> only one DOF
                PetscScalar value = c.value;
                VecSetValues(rhs, 1, &dofA, &value, ADD_VALUES);
            }
            else
            {
                for (auto dof : c.dofs)
                    if (dof->getIndex() == -1)
                        return;

                PetscInt *idx = new PetscInt[numBdNodes]();
                PetscScalar *values = new PetscScalar[numBdNodes]();

                int count = 0;
                for (auto dof : c.dofs)
                {
                    idx[count] = dof->getIndex();
                    values[count] = c.value;
                    count++;
                }

                VecSetValues(rhs, numBdNodes, idx, values, ADD_VALUES);
                delete[] idx;
                delete[] values;
            }
}