#include "../headers/Element.h"

Solid2D::Solid2D() {}
Solid2D::Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity, AnalysisParameters *_params)
    : Element(_index, _elemDimension, _elemConnectivity, _material, _physicalEntity, _params)
{
    for (auto n : _elemConnectivity)
    {
        n->setIsDiscritized();
        n->addDOF(new DOF(X, 0.));
        n->addDOF(new DOF(Y, 0.));
        n->addDOF(new DOF(D, 0.));
    }

    sF = new T3ShapeFunction();
    q = new TriangularQuadrature(numHammerPoints);
    coords = q->getQuadratureCoordinates();
    weights = q->getQuadratureWeights();
}
Solid2D::~Solid2D() {}

/*----------------------------------------------------------------------------------
                Assembling and solving problem with PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode Solid2D::getContribution(Mat &A, Vec &rhs, bool negativeLoad)
{
    PetscInt numElDOF = numElNodes * 2;
    PetscReal *localStiffnessMatrix = new PetscScalar[numElDOF * numElDOF]();
    PetscReal *localRHS = new PetscScalar[numElDOF]();
    PetscInt *idx = new PetscInt[numElDOF]();

    PetscReal kroen[2][2] = {{1., 0.}, {0., 1.}};

    const double G = material->getShearModulus();
    const double lame = material->getLameConstant();

    PetscInt count = 0;
    for (auto node : elemConnectivity)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
                idx[count++] = dof->getIndex();

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        PetscReal dX_dXsi[2][2] = {};

        for (PetscInt a = 0; a < 3; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dX_dXsi[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];

        PetscReal jac = dX_dXsi[0][0] * dX_dXsi[1][1] - dX_dXsi[0][1] * dX_dXsi[1][0];
        PetscReal wJac = weight * jac;

        PetscReal dX_dXsiInv[2][2] = {};
        dX_dXsiInv[0][0] = dX_dXsi[1][1] / jac;
        dX_dXsiInv[0][1] = -dX_dXsi[0][1] / jac;
        dX_dXsiInv[1][0] = -dX_dXsi[1][0] / jac;
        dX_dXsiInv[1][1] = dX_dXsi[0][0] / jac;

        PetscReal dN_dX[numElNodes][2] = {}; // Derivative of shape functions with respect to global coordinates; number of nodes x number of dimensions
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * dX_dXsiInv[j][i];

        double damageValue = 0.;

        // if (negativeLoad)
        //     damageValue = 0.;
        // else
        for (PetscInt a = 0; a < numElNodes; a++)
            damageValue += N[a] * elemConnectivity[a]->getDOFs()[2]->getDamageValue(); // DamageValue -> dstag = dn + delta_d^i

        PetscReal dCoeff = pow(1 - damageValue, 2);

        if (dCoeff > 1)
        {
            std::cout << " ================= ATENTION =================" << std::endl;
            std::cout << "Damage value is greater than 1:" << dCoeff << std::endl;
            std::cout << " ==========================================" << std::endl;
            // std::exit(EXIT_FAILURE);
        }

        if (damageValue < 0)
        {
            std::cout << " ================= ATENTION =================" << std::endl;
            std::cout << "Damage value is negative:" << dCoeff << std::endl;
            std::cout << " ==========================================" << std::endl;
            // std::exit(EXIT_FAILURE);
        }

        for (PetscInt a = 0; a < numElNodes; a++)

            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt b = 0; b < numElNodes; b++)
                {
                    PetscReal contraction = 0.;
                    for (PetscInt k = 0; k < 2; k++)
                        contraction += dN_dX[a][k] * dN_dX[b][k];

                    for (PetscInt j = 0; j < 2; j++)
                    {
                        PetscInt pos = numElDOF * (2 * a + i) + 2 * b + j;
                        localStiffnessMatrix[pos] += dCoeff * (G * contraction * wJac * kroen[i][j] + G * dN_dX[a][j] * dN_dX[b][i] * wJac + lame * dN_dX[a][i] * dN_dX[b][j] * wJac);
                    }
                }
        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    /*
        Applying prescribed Dirichlet boundary conditions
    */

    for (PetscInt a = 0; a < numElNodes; a++)
        for (PetscInt i = 0; i < 2; i++)
            for (PetscInt b = 0; b < numElNodes; b++)
                for (PetscInt j = 0; j < 2; j++)
                {
                    double value = elemConnectivity[b]->getDOFs()[j]->getValue();
                    double fi = -localStiffnessMatrix[numElDOF * (2 * a + i) + 2 * b + j] * elemConnectivity[b]->getDOF(j)->getValue();
                    ierr = VecSetValues(rhs, 1, &idx[2 * a + i], &fi, ADD_VALUES);
                    CHKERRQ(ierr);
                }

    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localStiffnessMatrix, ADD_VALUES);
    CHKERRQ(ierr);

    delete[] idx;
    delete[] localRHS;
    delete[] localStiffnessMatrix;

    return ierr;
}

PetscErrorCode Solid2D::getPhaseFieldContribution(Mat &A, Vec &rhs)
{
    const PetscInt numNodeDOF = 2;                              // Number of DOFs per node considering displacements only
    PetscInt numElDOF = numElNodes;                             // Only one DOF per node when considering only phase field
    PetscReal *localQ = new PetscScalar[numElDOF * numElDOF](); // Equivalent to matrix Qlocal in the phase field problem
    PetscReal *localRHS = new PetscScalar[numElDOF]();          // Equivalent to vector RHSlocal in the phase field problem
    PetscInt *idx = new PetscInt[numElDOF]();
    PetscScalar l0 = material->getL0();
    PetscScalar Gc = material->getGriffithCriterion();
    PetscScalar lame = material->getLameConstant();
    PetscScalar G = material->getShearModulus();
    std::string PFmodel = params->getPFModel();

    PetscInt count = 0;
    for (auto node : elemConnectivity)
        idx[count++] = node->getIndex(); // Phase field DOF has the same index as the node for the local problem

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        /*
            COMPUTING THE JACOBIAN AND ITS INVERSE
        */
        PetscReal dX_dXsi[2][2] = {};
        for (PetscInt a = 0; a < 3; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dX_dXsi[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];

        PetscReal jac = dX_dXsi[0][0] * dX_dXsi[1][1] - dX_dXsi[0][1] * dX_dXsi[1][0];
        PetscReal wJac = weight * jac;

        PetscReal dX_dXsiInv[2][2] = {};
        dX_dXsiInv[0][0] = dX_dXsi[1][1] / jac;
        dX_dXsiInv[0][1] = -dX_dXsi[0][1] / jac;
        dX_dXsiInv[1][0] = -dX_dXsi[1][0] / jac;
        dX_dXsiInv[1][1] = dX_dXsi[0][0] / jac;

        PetscReal dN_dX[numElNodes][2] = {}; // Derivative of shape functions with respect to global coordinates; number of nodes x number of dimensions
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * dX_dXsiInv[j][i];

        // Computing uk,l and ul,k
        PetscScalar gradU[2][2] = {};
        for (PetscInt c = 0; c < numElNodes; c++)
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    gradU[k][l] += elemConnectivity[c]->getDOFs()[k]->getValue() * dN_dX[c][l];

        PetscScalar divU = gradU[0][0] + gradU[1][1]; // uk,k and ul,l are the same

        // ======================= FIRST DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE ========================
        PetscScalar damageValue = 0.;
        for (PetscInt c = 0; c < numElNodes; c++)
            damageValue += N[c] * elemConnectivity[c]->getDOFs()[2]->getValue(); // dn

        PetscScalar firstInt[numElNodes] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    firstInt[a] += (1 - damageValue) * N[a] * (G * 0.5 * (gradU[k][l] + gradU[l][k]) * (gradU[k][l] + gradU[l][k])) * wJac;

            firstInt[a] += (1 - damageValue) * N[a] * lame * divU * divU * wJac;
        }

        /*
            AT1 OR AT2 PHASE FIELD MODEL:
            The second integral of the first derivative with respect to the field variable is different for AT1 and AT2 models.
        */
        PetscScalar secondInt[numElNodes] = {};
        if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
            {
                for (PetscInt c = 0; c < numElNodes; c++)
                    for (PetscInt k = 0; k < 2; k++)
                        secondInt[a] += Gc * (l0 * elemConnectivity[c]->getDOFs()[2]->getValue() * dN_dX[a][k] * dN_dX[c][k]) * wJac;

                secondInt[a] += Gc * (1 / l0 * damageValue * N[a]) * wJac;
            }
        }
        else if (PFmodel == "AT1")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
            {
                for (PetscInt c = 0; c < numElNodes; c++)
                    for (PetscInt k = 0; k < 2; k++)
                        secondInt[a] += Gc * (0.75 * l0 * elemConnectivity[c]->getDOFs()[2]->getValue() * dN_dX[a][k] * dN_dX[c][k]) * wJac;

                secondInt[a] += Gc * ((0.375 / l0) * N[a]) * wJac;
            }
        }

        for (PetscInt a = 0; a < numElNodes; a++)
            localRHS[a] = -firstInt[a] + secondInt[a];

        ierr = VecSetValues(rhs, numElDOF, idx, localRHS, ADD_VALUES);
        CHKERRQ(ierr);
        // ======================= SECOND DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE =======================
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt b = 0; b < numElNodes; b++)
            {
                for (PetscInt k = 0; k < 2; k++)
                    for (PetscInt l = 0; l < 2; l++)
                    {
                        PetscInt pos = numElDOF * a + b;
                        // PetscScalar value = N[a] * N[b] * (G * 0.5 * (gradU[k][l] + gradU[l][k]) * (gradU[k][l] + gradU[l][k]) + lame * divU * divU) * wJac;
                        PetscScalar value = N[a] * N[b] * (G * 0.5 * (gradU[k][l] + gradU[l][k]) * (gradU[k][l] + gradU[l][k])) * wJac;
                        localQ[pos] += value; // Integral 1
                    }

                PetscInt pos = numElDOF * a + b;
                localQ[pos] += N[a] * N[b] * lame * divU * divU * wJac;
            }

        /*
            AT1 OR AT2 PHASE FIELD MODEL:
            The second integral of the first derivative with respect to the field variable is different for AT1 and AT2 models.
        */
        if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt b = 0; b < numElNodes; b++)
                {
                    PetscScalar contraction = 0.;
                    for (PetscInt k = 0; k < 2; k++)
                        contraction += dN_dX[a][k] * dN_dX[b][k];

                    PetscInt pos = numElDOF * a + b;
                    double value = Gc * (1 / l0 * N[a] * N[b] + l0 * contraction) * wJac;
                    localQ[pos] += value; // Integral 2 for AT2
                }
        }
        else if (PFmodel == "AT1")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt b = 0; b < numElNodes; b++)
                {
                    PetscScalar contraction = 0.;
                    for (PetscInt k = 0; k < 2; k++)
                        contraction += dN_dX[a][k] * dN_dX[b][k];

                    PetscInt pos = numElDOF * a + b;
                    double value = Gc * (0.75 * l0 * contraction) * wJac;
                    localQ[pos] += value; // Integral 2 for AT1
                }
        }

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localQ, ADD_VALUES);
    CHKERRQ(ierr);

    delete[] idx;
    delete[] localQ;
    delete[] localRHS;

    return ierr;
}

// For an initial prescribed damage field

PetscErrorCode Solid2D::getPhaseFieldContribution(Mat &A, Vec &rhs, bool _PrescribedDamageField)
{
    const PetscInt numNodeDOF = 2;                              // Number of DOFs per node considering displacements only
    PetscInt numElDOF = numElNodes;                             // Only one DOF per node when considering only phase field
    PetscReal *localQ = new PetscScalar[numElDOF * numElDOF](); // Equivalent to matrix Qlocal in the phase field problem
    PetscReal *localRHS = new PetscScalar[numElDOF]();          // Equivalent to vector RHSlocal in the phase field problem
    PetscInt *idx = new PetscInt[numElDOF]();
    PetscScalar l0 = material->getL0();
    PetscScalar Gc = material->getGriffithCriterion();
    PetscScalar lame = material->getLameConstant();
    PetscScalar G = material->getShearModulus();

    PetscInt count = 0;
    for (auto node : elemConnectivity)
        idx[count++] = node->getIndex(); // Phase field DOF has the same index as the node for the local problem

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        /*
            COMPUTING THE JACOBIAN AND ITS INVERSE
        */
        PetscReal dX_dXsi[2][2] = {};
        for (PetscInt a = 0; a < 3; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dX_dXsi[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];

        PetscReal jac = dX_dXsi[0][0] * dX_dXsi[1][1] - dX_dXsi[0][1] * dX_dXsi[1][0];
        PetscReal wJac = weight * jac;

        PetscReal dX_dXsiInv[2][2] = {};
        dX_dXsiInv[0][0] = dX_dXsi[1][1] / jac;
        dX_dXsiInv[0][1] = -dX_dXsi[0][1] / jac;
        dX_dXsiInv[1][0] = -dX_dXsi[1][0] / jac;
        dX_dXsiInv[1][1] = dX_dXsi[0][0] / jac;

        PetscReal dN_dX[numElNodes][2] = {}; // Derivative of shape functions with respect to global coordinates; number of nodes x number of dimensions
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * dX_dXsiInv[j][i];

        // Computing uk,l and ul,k
        PetscScalar gradU[2][2] = {};
        for (PetscInt c = 0; c < numElNodes; c++)
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    gradU[k][l] += 0.0 * dN_dX[c][l];

        PetscScalar divU = gradU[0][0] + gradU[1][1]; // uk,k and ul,l are the same

        // ======================= FIRST DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE ========================
        PetscScalar damageValue = 0.;
        for (PetscInt c = 0; c < numElNodes; c++)
            damageValue += N[c] * elemConnectivity[c]->getDOFs()[2]->getValue();

        PetscScalar firstInt[numElNodes] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    firstInt[a] += (1 - damageValue) * N[a] * (G * 0.5 * (gradU[k][l] + gradU[l][k]) * (gradU[k][l] + gradU[l][k])) * wJac;

            firstInt[a] += (1 - damageValue) * N[a] * lame * divU * divU * wJac;
        }

        PetscScalar secondInt[numElNodes] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            for (PetscInt c = 0; c < numElNodes; c++)
                for (PetscInt k = 0; k < 2; k++)
                    secondInt[a] += Gc * (l0 * elemConnectivity[c]->getDOFs()[2]->getValue() * dN_dX[a][k] * dN_dX[c][k]) * wJac;

            secondInt[a] += Gc * (1 / l0 * damageValue * N[a]) * wJac;
        }

        for (PetscInt a = 0; a < numElNodes; a++)
            localRHS[a] = -firstInt[a] + secondInt[a];

        ierr = VecSetValues(rhs, numElDOF, idx, localRHS, ADD_VALUES);
        CHKERRQ(ierr);
        // ======================= SECOND DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE =======================
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt b = 0; b < numElNodes; b++)
                for (PetscInt k = 0; k < 2; k++)
                    for (PetscInt l = 0; l < 2; l++)
                    {
                        PetscInt pos = numElDOF * a + b;
                        PetscScalar value = N[a] * N[b] * (G * 0.5 * (gradU[k][l] + gradU[l][k]) * (gradU[k][l] + gradU[l][k]) + lame * divU * divU) * wJac;
                        localQ[pos] += value; // Integral 1
                    }

        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt b = 0; b < numElNodes; b++)
            {
                PetscScalar contraction = 0.;
                for (PetscInt k = 0; k < 2; k++)
                    contraction += dN_dX[a][k] * dN_dX[b][k];

                PetscInt pos = numElDOF * a + b;
                double value = Gc * (1 / l0 * N[a] * N[b] + l0 * contraction) * wJac;
                localQ[pos] += Gc * (1 / l0 * N[a] * N[b] + l0 * contraction) * wJac; // Integral 2
            }

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localQ, ADD_VALUES);
    CHKERRQ(ierr);

    delete[] idx;
    delete[] localQ;
    delete[] localRHS;

    return ierr;
}

void Solid2D::Test(PetscScalar &integral)
{
    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        double x[2] = {};
        double dX[2][2] = {};

        for (int a = 0; a < 3; a++)
            for (int i = 0; i < 2; i++)
            {
                x[i] += N[a] * elemConnectivity[a]->getInitialCoordinates()[i];

                for (int j = 0; j < 2; j++)
                    dX[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];
            }
        double jac = (dX[0][0] * dX[1][1] - dX[0][1] * dX[1][0]);
        integral += (3. * x[0] + 7. * x[1] - 2.) * weight * jac;
        // try with a quadratic function -> refine
    }
}

/*----------------------------------------------------------------------------------
                Assembling and solving problem without PETSc
----------------------------------------------------------------------------------
*/
void Solid2D::getContribution()
{
    const int numElDOF = numElNodes * 2;
    localStiff = MatrixXd::Zero(6, 6);

    const double G = material->getShearModulus();
    const double lame = material->getLameConstant();

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        double dX_dXsi[2][2] = {};
        for (int a = 0; a < numElNodes; a++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    dX_dXsi[i][j] += elemConnectivity[a]->getInitialCoordinates()[i] * dN[a][j];

        double jac = dX_dXsi[0][0] * dX_dXsi[1][1] - dX_dXsi[0][1] * dX_dXsi[1][0];
        double wJac = jac * weight;

        double dX_dXsiInv[2][2] = {};
        dX_dXsiInv[0][0] = dX_dXsi[1][1] / jac;
        dX_dXsiInv[0][1] = -dX_dXsi[0][1] / jac;
        dX_dXsiInv[1][0] = -dX_dXsi[1][0] / jac;
        dX_dXsiInv[1][1] = dX_dXsi[0][0] / jac;

        double dN_dX[numElNodes][2] = {}; // Derivative of shape functions with respect to global coordinates; number of nodes x number of dimensions
        for (int a = 0; a < numElNodes; a++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * dX_dXsiInv[j][i];

        /*
        For a solid element, four loops are needed to compute the element stiffness matrix.
        The first loop is over the number of nodes in the element.
        The second loop is over the number of degrees of freedom per node.
        The third loop is over the number of nodes in the element.
        The fourth loop is over the number of degrees of freedom per node.

        K_ij ^ {ab} = G * \delta_{ij} * N,k ^ a * N,k ^ b + G * N,i ^ b * N,j ^ a + \lambda * N,i ^ a * N,j ^ b
        */

        for (int a = 0; a < numElNodes; a++)
        {
            for (int b = 0; b < numElNodes; b++)
            {
                double contraction = 0.;
                for (int k = 0; k < 2; k++)
                    contraction += dN_dX[a][k] * dN_dX[b][k];

                for (int i = 0; i < 2; i++)
                {
                    localStiff(2 * a + i, 2 * b + i) += G * contraction * wJac; // Due to Kroenecker delta

                    for (int j = 0; j < 2; j++)
                        localStiff(2 * a + i, 2 * b + j) += (G * dN_dX[a][j] * dN_dX[b][i] + lame * dN_dX[a][i] * dN_dX[b][j]) * wJac;
                }
            }
        }
    }
    std::cout << "Local stiffness matrix: " << " for element " << index << std::endl;
    std::cout << localStiff << std::endl;
}

void Solid2D::assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff)
{
    for (int a = 0; a < numElNodes; a++)
        for (int i = 0; i < 2; i++)
        {
            int dof1 = getNode(a)->getDOFs()[i]->getIndex();
            for (int b = 0; b < numElNodes; b++)
                for (int j = 0; j < 2; j++)
                {
                    int dof2 = getNode(b)->getDOFs()[j]->getIndex();
                    double value = localStiff(2 * a + i, 2 * b + j);
                    GlobalStiff(dof1, dof2) += localStiff(2 * a + i, 2 * b + j);
                }
        }

    // for (int i = 0; i < 2; i++)
    //     for (int j = 0; j < 2; j++)
    //     {
    //         GlobalStiff(dof1 + i, dof1 + j) += localStiff(i, j);
    //         GlobalStiff(dof1 + i, dof2 + j) += localStiff(i, j + 2);
    //         GlobalStiff(dof1 + i, dof3 + j) += localStiff(i, j + 4);

    //         GlobalStiff(dof2 + i, dof1 + j) += localStiff(i + 2, j);
    //         GlobalStiff(dof2 + i, dof2 + j) += localStiff(i + 2, j + 2);
    //         GlobalStiff(dof2 + i, dof3 + j) += localStiff(i + 2, j + 4);

    //         GlobalStiff(dof3 + i, dof1 + j) += localStiff(i + 4, j);
    //         GlobalStiff(dof3 + i, dof2 + j) += localStiff(i + 4, j + 2);
    //         GlobalStiff(dof3 + i, dof3 + j) += localStiff(i + 4, j + 4);
    //     }
}

// =========================================================
// ==================== POST PROCESSING ====================
// =========================================================
PetscErrorCode Solid2D::calculateStress()
{
    const PetscInt &numElDOF = numElNodes * 2;
    const double &G = material->getShearModulus();
    const double &lame = material->getLameConstant();

    double **coord;
    sF->getNodalXi(coord);

    for (int a = 0; a < numElNodes; a++)
    {
        double *xi = coords[a];
        double **dN = sF->getShapeFunctionDerivative(xi);

        PetscReal dX_dXsi[2][2] = {}, dX_dXsiInv[2][2] = {}, gradU[2][2] = {}, strain[2][2] = {}, stress[2][2] = {};
        for (PetscInt a = 0; a < 3; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dX_dXsi[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];
        PetscReal jac = dX_dXsi[0][0] * dX_dXsi[1][1] - dX_dXsi[0][1] * dX_dXsi[1][0];

        dX_dXsiInv[0][0] = dX_dXsi[1][1] / jac;
        dX_dXsiInv[0][1] = -dX_dXsi[0][1] / jac;
        dX_dXsiInv[1][0] = -dX_dXsi[1][0] / jac;
        dX_dXsiInv[1][1] = dX_dXsi[0][0] / jac;

        PetscReal dN_dX[numElNodes][2] = {}; // Derivative of shape functions with respect to global coordinates; number of nodes x number of dimensions
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * dX_dXsiInv[j][i];

        for (PetscInt a = 0; a < 3; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    gradU[i][j] += elemConnectivity[a]->getDOFs()[i]->getValue() * dN_dX[a][j];

        for (PetscInt i = 0; i < 2; i++)
            for (PetscInt j = 0; j < 2; j++)
                strain[i][j] = 0.5 * (gradU[i][j] + gradU[j][i]);

        material->Lame(strain, stress);

        const int &recurrence = elemConnectivity[a]->getInverseIncidence().size();
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                elemConnectivity[a]->incrementStress(i, j, stress[i][j] / double(recurrence));

        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }
    for (int i = 0; i < numElNodes; i++)
        delete[] coord[i];
    delete[] coord;

    return ierr;
}