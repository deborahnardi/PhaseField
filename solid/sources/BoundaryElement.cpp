#include "../headers/BoundaryElement.h"

BoundaryElement::BoundaryElement() {}
BoundaryElement::BoundaryElement(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, AnalysisParameters *_params)
    : index(_index), elemDimension(_elemDimension), elemConnectivity(_elemConnectivity), material(_material), physicalEntity(_physicalEntity), params(_params)
{
    numBdNodes = _elemConnectivity.size();

    // for (auto n : _elemConnectivity)
    // {
    //     // n->setIsDiscritized();
    //     n->addDOF(new DOF(X, 0.));
    //     n->addDOF(new DOF(Y, 0.));
    //     n->addDOF(new DOF(D, 0.));
    // }

    if (elemDimension == 1)
    {
        numQuadraturePoints = 1;
        sF = new S2ShapeFunction();
        q = new LineQuadrature(numQuadraturePoints);
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
                else if (_bdType == DAMAGE)
                {
                    dof->setPrescribedDamage();
                    dof->setDamageValue(_value);
                    dof->setValue(_value);
                }
                dofVec.push_back(dof);
            }
    conditions.push_back({_bdType, dofVec, _value}); // Each boundary element has its own conditions"
}

void BoundaryElement::updateBoundaryValues(double _lambda)
{
    for (auto c : conditions)
        if (c.bdType == DIRICHLET)
            for (auto dof : c.dofs)
                dof->setValue(_lambda * dof->getDirichletValue());
}
/*----------------------------------------------------------------------------------
                Assembling and solving problem with PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode BoundaryElement::getContribution(Vec &rhs) // entra com lambda aq
{
    // Neumann is applied in weak form
    for (auto c : conditions)
        if (c.bdType == NEUMANN)
            if (elemDimension == 0)
            {
                PetscInt dofA = c.dofs[0]->getIndex(); // c.dofs[0] because if elemDimension == 0, it is only a node -> only one DOF
                PetscScalar value = c.value;           // entra com lambda aqui
                VecSetValues(rhs, 1, &dofA, &value, ADD_VALUES);
            }
            else
            {
                double **coords = q->getQuadratureCoordinates();
                double *weights = q->getQuadratureWeights();

                PetscInt *idx = new PetscInt[numBdNodes]();
                PetscScalar *localRhs = new PetscScalar[numBdNodes]();

                for (int a = 0; a < numBdNodes; a++)
                    idx[a] = c.dofs[a]->getIndex();

                for (int ig = 0; ig < numQuadraturePoints; ig++) // ig stands for gauss points
                {
                    double *xi = coords[ig];
                    double weight = weights[ig];

                    double *N = sF->evaluateShapeFunction(xi);
                    double **dN = sF->getShapeFunctionDerivative(xi);

                    // jac for 1D element in 2D space is the same as the length of the element

                    double tangent[2] = {};
                    for (int a = 0; a < numBdNodes; a++)
                        for (int i = 0; i < 2; i++)
                            tangent[i] += dN[a][0] * elemConnectivity[a]->getInitialCoordinates()[i];

                    double jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
                    double wJac = weight * jac;

                    for (int a = 0; a < numBdNodes; a++)
                    {
                        PetscScalar ti = N[a] * c.value * wJac; // entra com lambda aqui
                        localRhs[a] += ti;
                    }

                    delete[] N;
                    for (int a = 0; a < numBdNodes; a++)
                        delete[] dN[a];
                    delete[] dN;
                }
                ierr = VecSetValues(rhs, numBdNodes, idx, localRhs, ADD_VALUES);
                CHKERRQ(ierr);

                delete[] idx;
                delete[] localRhs;
            }
    return ierr;
}

void BoundaryElement::calculateNormal(double **dN, double normal[], double &djac)
{
    double dX_dxi[2] = {};
    for (int i = 0; i < 2; i++)
        for (int a = 0; a < numBdNodes; a++)
            dX_dxi[i] += elemConnectivity[a]->getInitialCoordinates()[i] * dN[a][0];

    normal[0] = dX_dxi[1];
    normal[1] = -dX_dxi[0];

    djac = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    normal[0] /= djac;
    normal[1] /= djac;
}

void BoundaryElement::calculateTraction(double force[])
{
    double **coords = q->getQuadratureCoordinates();
    double *weights = q->getQuadratureWeights();

    for (int ip = 0; ip < numQuadraturePoints; ip++)
    {
        double *xi = coords[ip];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        double _stress[2][2] = {};
        for (int a = 0; a < numBdNodes; a++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    _stress[i][j] += elemConnectivity[a]->getStress(i, j) * N[a];

        double normal[2] = {}, djac = 0.;
        calculateNormal(dN, normal, djac);
        const double wjac = weights[ip] * djac;

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                force[i] += _stress[i][j] * normal[j] * wjac; // It is actually a force, not a traction

        delete[] N;
        for (int i = 0; i < numBdNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }
}

/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
----------------------------------------------------------------------------------
*/

void BoundaryElement::getContributionNoPetsc(VectorXd &F, MatrixXd &K)
{
    for (auto c : conditions)
        if (c.bdType == NEUMANN)
            if (elemDimension == 0)
                F(c.dofs[0]->getIndex()) += c.value;
            else
            {
                double **coords = q->getQuadratureCoordinates();
                double *weights = q->getQuadratureWeights();

                for (int ig = 0; ig < numQuadraturePoints; ig++) // ig stands for gauss points
                {
                    double *xi = coords[ig];
                    double weight = weights[ig];

                    double *N = sF->evaluateShapeFunction(xi);
                    double **dN = sF->getShapeFunctionDerivative(xi);

                    // jac for 1D element in 2D space is the same as the length of the element

                    double tangent[2] = {};
                    for (int a = 0; a < numBdNodes; a++)
                        for (int i = 0; i < 2; i++)
                            tangent[i] += dN[a][0] * elemConnectivity[a]->getInitialCoordinates()[i];

                    double jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
                    double wJac = weight * jac;

                    for (int a = 0; a < numBdNodes; a++)
                    {
                        int pos = c.dofs[a]->getIndex();
                        double ti = N[a] * c.value * wJac;
                        F(pos) += ti;
                    }
                }
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