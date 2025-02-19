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
PetscErrorCode Solid2D::getContribution(Mat &A, Vec &rhs, bool negativeLoad, bool _PrescribedDamageField)
{
    /*
         Ke = int{B^T C B}dOmega_e
         B^T C B IS HEREIN DEFINED THROUGH AN ANALYTICAL EXPRESSION
         ONLY THE UPPER TRIANGULAR PART OF THE MATRIX IS COMPUTED (21 ELEMENTS DUE TO SYMMETRY)
         ONLY ISOTROPIC MATERIALS ARE CONSIDERED (SEE CONSTITUTIVE TENSOR)

         B = [dN1/dx, 0, dN2/dx, 0, dN3/dx, 0]
                [0, dN1/dy, 0, dN2/dy, 0, dN3/dy]
                [dN1/dy, dN1/dx, dN2/dy, dN2/dx, dN3/dy, dN3/dx]

            dN/dx = dN/dxi * dxi/dx + dN/deta * deta/dx
            dN/dy = dN/dxi * dxi/dy + dN/deta * deta/dy

        Jacobian matrix:
            J = [dx/dxi, dy/dxi]
                [dx/deta, dy/deta]

        Inverse of the Jacobian matrix:
            J^-1 = 1/det(J) * [dy/deta, -dy/dxi]
                              [-dx/deta, dx/dxi]

        Derivative of the shape functions with respect to global coordinates:
            [dN1/dx, dN1/dy]             [dN1/dxi, dN1/deta]
            [dN2/dx, dN2/dy](3x2)   =    [dN2/dxi, dN2/deta](3x2)  * [J^-1](2x2)
            [dN3/dx, dN3/dy]             [dN3/dxi, dN3/deta]
    */
    PetscInt numElDOF = numElNodes * 2;                      // = 6
    PetscReal *localStiffnessMatrix = new PetscScalar[36](); // 6x6 matrix, 21 is the number of elements in the upper triangular part of the matrix
    PetscReal *localRHS = new PetscScalar[numElDOF]();
    PetscInt *idx = new PetscInt[numElDOF]();
    PetscScalar mu = material->getShearModulus();
    PetscScalar kappa = 2.0 / 3.0 * mu + material->getLameConstant();

    // SETTING idx VALUES CONSIDERING ONLY THE SYMMETRIC PART OF THE MATRIX
    PetscInt count = 0;
    for (auto node : elemConnectivity)
        for (auto dof : node->getDOFs())
            if (dof->getDOFType() != D)
                idx[count++] = dof->getIndex();

    double tensorK[3][3] = {};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            tensorK[i][j] = 1.0;

    double tensorI[3][3] = {};
    for (int i = 0; i < 3; i++)
        tensorI[i][i] = 1.0;

    tensorI[2][2] = 0.5;

    double tensorJ[3][3] = {}; // J = I - 1/3 * K
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            tensorJ[i][j] = tensorI[i][j] - 1.0 / 3.0 * tensorK[i][j];

    // const double c00 = material->getLameConstant() + 2.0 * material->getShearModulus();
    // const double c01 = material->getLameConstant();
    // const double c10 = c01;
    // const double c11 = c00;
    // const double c22 = material->getShearModulus();

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        // COMPUTING THE JACOBIAN MATRIX
        PetscReal jacobianMatrix[2][2] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    jacobianMatrix[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];

        // COMPUTING THE JACOBIAN FOR 2D SPACE
        PetscReal jac = jacobianMatrix[0][0] * jacobianMatrix[1][1] - jacobianMatrix[0][1] * jacobianMatrix[1][0];
        PetscReal wJac = weight * jac;

        // COMPUTING THE INVERSE OF THE JACOBIAN MATRIX
        PetscReal jacobianMatrixInv[2][2] = {};
        jacobianMatrixInv[0][0] = jacobianMatrix[1][1] / jac;
        jacobianMatrixInv[0][1] = -jacobianMatrix[0][1] / jac;
        jacobianMatrixInv[1][0] = -jacobianMatrix[1][0] / jac;
        jacobianMatrixInv[1][1] = jacobianMatrix[0][0] / jac;

        // COMPUTING THE DERIVATIVE OF THE SHAPE FUNCTIONS WITH RESPECT TO GLOBAL COORDINATES (3x2 matrix - nNodes x nDirections)
        PetscReal dN_dX[numElNodes][2] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * jacobianMatrixInv[j][i];

        // COMPUTING uk,l AND ul,k
        PetscScalar gradU[2][2] = {};
        for (PetscInt c = 0; c < numElNodes; c++)
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    gradU[k][l] += elemConnectivity[c]->getDOFs()[k]->getValue() * dN_dX[c][l];

        // if (_PrescribedDamageField)
        //     for (PetscInt c = 0; c < numElNodes; c++)
        //         for (PetscInt k = 0; k < 2; k++)
        //             for (PetscInt l = 0; l < 2; l++)
        //                 gradU[k][l] += 0.0 * dN_dX[c][l];

        PetscScalar divU = gradU[0][0] + gradU[1][1]; // uk,k and ul,l are the same

        // COMPUTING THE DEGRADATION FUNCTION (1-d)^2
        double damageValue = 0.;
        for (PetscInt a = 0; a < numElNodes; a++)
            damageValue += N[a] * elemConnectivity[a]->getDOFs()[2]->getDamageValue(); // DamageValue -> dstag = dn + delta_d^i

        PetscReal dCoeff = pow(1 - damageValue, 2);

        // ASSEMBLING CONSTITUTIVE MATRIX CONSIDERING THE ENERGY SPLIT

        // COMPUTING THE CONSTITUTIVE TENSOR
        double tensorCPlus[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                tensorCPlus[i][j] = 2.0 * mu * tensorJ[i][j];

        double tensorCMinus[3][3] = {};

        if (divU >= 0)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tensorCPlus[i][j] += kappa * tensorK[i][j];
                    tensorCPlus[i][j] *= dCoeff;
                }
        else if (divU < 0)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tensorCMinus[i][j] += kappa * tensorK[i][j];

        double tensorC[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];

        // RELATING dN_dX TO THE B MATRIX
        PetscReal B[3][6] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            B[0][2 * a] = dN_dX[a][0];
            B[1][2 * a + 1] = dN_dX[a][1];
            B[2][2 * a] = dN_dX[a][1];
            B[2][2 * a + 1] = dN_dX[a][0];
        }

        // COMPUTING THE LOCAL STIFFNESS MATRIX

        localStiffnessMatrix[0] += dCoeff * (B[0][0] * B[0][0] * tensorC[0][0] + B[1][1] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[1] += dCoeff * (B[0][0] * B[1][1] * tensorC[0][1] + B[0][0] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[2] += dCoeff * (B[0][0] * B[0][2] * tensorC[0][0] + B[1][1] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[3] += dCoeff * (B[0][0] * B[1][3] * tensorC[0][1] + B[0][2] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[4] += dCoeff * (B[0][0] * B[0][4] * tensorC[0][0] + B[1][1] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[5] += dCoeff * (B[0][0] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[6] = localStiffnessMatrix[1];
        localStiffnessMatrix[7] += dCoeff * (B[0][0] * B[0][0] * tensorC[2][2] + B[1][1] * B[1][1] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[8] += dCoeff * (B[0][0] * B[1][3] * tensorC[2][2] + B[0][2] * B[1][1] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[9] += dCoeff * (B[0][0] * B[0][2] * tensorC[2][2] + B[1][1] * B[1][3] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[10] += dCoeff * (B[0][0] * B[1][5] * tensorC[2][2] + B[0][4] * B[1][1] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[11] += dCoeff * (B[0][0] * B[0][4] * tensorC[2][2] + B[1][1] * B[1][5] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[12] = localStiffnessMatrix[2];
        localStiffnessMatrix[13] = localStiffnessMatrix[8];
        localStiffnessMatrix[14] += dCoeff * (B[0][2] * B[0][2] * tensorC[0][0] + B[1][3] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[15] += dCoeff * (B[0][2] * B[1][3] * tensorC[0][1] + B[0][2] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[16] += dCoeff * (B[0][2] * B[0][4] * tensorC[0][0] + B[1][3] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[17] += dCoeff * (B[0][2] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[18] = localStiffnessMatrix[3];
        localStiffnessMatrix[19] = localStiffnessMatrix[9];
        localStiffnessMatrix[20] = localStiffnessMatrix[15];
        localStiffnessMatrix[21] += dCoeff * (B[0][2] * B[0][2] * tensorC[2][2] + B[1][3] * B[1][3] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[22] += dCoeff * (B[0][2] * B[1][5] * tensorC[2][2] + B[0][4] * B[1][3] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[23] += dCoeff * (B[0][2] * B[0][4] * tensorC[2][2] + B[1][3] * B[1][5] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[24] = localStiffnessMatrix[4];
        localStiffnessMatrix[25] = localStiffnessMatrix[10];
        localStiffnessMatrix[26] = localStiffnessMatrix[16];
        localStiffnessMatrix[27] = localStiffnessMatrix[22];
        localStiffnessMatrix[28] += dCoeff * (B[0][4] * B[0][4] * tensorC[0][0] + B[1][5] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[29] += dCoeff * (B[0][4] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[30] = localStiffnessMatrix[5];
        localStiffnessMatrix[31] = localStiffnessMatrix[11];
        localStiffnessMatrix[32] = localStiffnessMatrix[17];
        localStiffnessMatrix[33] = localStiffnessMatrix[23];
        localStiffnessMatrix[34] = localStiffnessMatrix[29];
        localStiffnessMatrix[35] += dCoeff * (B[0][4] * B[0][4] * tensorC[2][2] + B[1][5] * B[1][5] * tensorC[1][1]) * wJac;

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    // ASSEMBLING THE RHS VECTOR

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
    PetscScalar kappa = 2.0 / 3.0 * G + lame;
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

        // COMPUTING THE JACOBIAN MATRIX
        PetscReal jacobianMatrix[2][2] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    jacobianMatrix[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];

        // COMPUTING THE JACOBIAN FOR 2D SPACE
        PetscReal jac = jacobianMatrix[0][0] * jacobianMatrix[1][1] - jacobianMatrix[0][1] * jacobianMatrix[1][0];
        PetscReal wJac = weight * jac;

        // COMPUTING THE INVERSE OF THE JACOBIAN MATRIX
        PetscReal jacobianMatrixInv[2][2] = {};
        jacobianMatrixInv[0][0] = jacobianMatrix[1][1] / jac;
        jacobianMatrixInv[0][1] = -jacobianMatrix[0][1] / jac;
        jacobianMatrixInv[1][0] = -jacobianMatrix[1][0] / jac;
        jacobianMatrixInv[1][1] = jacobianMatrix[0][0] / jac;

        // COMPUTING THE DERIVATIVE OF THE SHAPE FUNCTIONS WITH RESPECT TO GLOBAL COORDINATES (3x2 matrix - nNodes x nDirections)
        PetscReal dN_dX[numElNodes][2] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt i = 0; i < 2; i++)
                for (PetscInt j = 0; j < 2; j++)
                    dN_dX[a][i] += dN[a][j] * jacobianMatrixInv[j][i];

        // ======================= FIRST DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE ========================

        // COMPUTING uk,l AND ul,k
        PetscScalar gradU[2][2] = {};
        for (PetscInt c = 0; c < numElNodes; c++)
            for (PetscInt k = 0; k < 2; k++)
                for (PetscInt l = 0; l < 2; l++)
                    gradU[k][l] += elemConnectivity[c]->getDOFs()[k]->getValue() * dN_dX[c][l];

        if (_PrescribedDamageField)
            for (PetscInt c = 0; c < numElNodes; c++)
                for (PetscInt k = 0; k < 2; k++)
                    for (PetscInt l = 0; l < 2; l++)
                        gradU[k][l] += 0.0 * dN_dX[c][l];

        PetscScalar divU = gradU[0][0] + gradU[1][1]; // uk,k and ul,l are the same

        // COMPUTING \psi^+(\epsilon): ELASTIC ENERGY DENSITY (POSITIVE PART) - DRIVING FORCE
        double psiPlus = 0.;
        for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
                psiPlus += G * ((gradU[k][l] * gradU[k][l] + gradU[k][l] * gradU[l][k]) * 0.5 - 1.0 / 3.0 * divU * divU) * wJac;

        if (divU >= 0)
            psiPlus += kappa * 0.5 * divU * divU * wJac;

        PetscScalar damageValue = 0.;
        for (PetscInt c = 0; c < numElNodes; c++)
            damageValue += N[c] * elemConnectivity[c]->getDOFs()[2]->getValue(); // dn

        const PetscScalar dCoeff = -2 * (1 - damageValue);

        PetscScalar firstInt[numElNodes] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
            firstInt[a] += dCoeff * N[a] * psiPlus;

        /*
            AT1 OR AT2 PHASE FIELD MODEL:
            The second integral of the first derivative with respect to the field variable is different for AT1 and AT2 models.

            RELATING dN_dX TO THE B_d MATRIX (B_d is the derivative of the B matrix with respect to the field variable = B_d(2x3))
            B_d = [dN1/dx, dN2/dx, dN3/dx]
                  [dN1/dy, dN2/dy, dN3/dy]
        */

        // ASSEMBLING B_d MATRIX
        PetscReal B_d[2][3] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            B_d[0][a] = dN_dX[a][0];
            B_d[1][a] = dN_dX[a][1];
        }

        // ASSEMBLING B_d^T
        PetscReal B_dT[3][2] = {};
        for (PetscInt i = 0; i < 3; i++)
            for (PetscInt j = 0; j < 2; j++)
                B_dT[i][j] = B_d[j][i];

        PetscScalar secondInt[numElNodes] = {};
        if (PFmodel == "AT1")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                        secondInt[a] += Gc * (0.375 / l0 * N[a] + 0.75 * l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;
        }
        else if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
            {
                double valueD = elemConnectivity[a]->getDOFs()[2]->getValue();
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                    {
                        double valueSecondInt = Gc * (l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;
                        secondInt[a] += Gc * (l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;
                    }

                secondInt[a] += Gc * (1 / l0 * damageValue * N[a]) * wJac;
            }
        }

        for (PetscInt a = 0; a < numElNodes; a++)
            localRHS[a] = firstInt[a] + secondInt[a];

        ierr = VecSetValues(rhs, numElDOF, idx, localRHS, ADD_VALUES);
        CHKERRQ(ierr);
        // ======================= SECOND DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE =======================
        for (PetscInt a = 0; a < numElNodes; a++)
            for (PetscInt b = 0; b < numElNodes; b++)
            {
                PetscInt pos = numElDOF * a + b;
                localQ[pos] += 2 * N[a] * N[b] * psiPlus; // Integral 1 is the same for AT1 and AT2
            }

        /*
            AT1 OR AT2 PHASE FIELD MODEL:
            The second integral of the second derivative with respect to the field variable is different for AT1 and AT2 models.
        */
        if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt b = 0; b < numElNodes; b++)
                {
                    PetscInt pos = numElDOF * a + b;
                    localQ[pos] += 1.0 / l0 * Gc * N[a] * N[b] * wJac;
                }

            // PERFORMING Gc * c_d * B_d^T * B_d
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt b = 0; b < numElNodes; b++)
                    for (PetscInt i = 0; i < 2; i++)
                    {
                        PetscInt pos = numElDOF * a + b;
                        localQ[pos] += Gc * l0 * B_dT[a][i] * B_d[i][b] * wJac;
                    }
        }
        else if (PFmodel == "AT1")
        {
            // PERFORMING Gc * c_d * B_d^T * B_d
            for (PetscInt a = 0; a < numElNodes; a++)
                for (PetscInt b = 0; b < numElNodes; b++)
                    for (PetscInt i = 0; i < 2; i++)
                    {
                        PetscInt pos = numElDOF * a + b;
                        localQ[pos] += Gc * 0.75 * l0 * B_dT[a][i] * B_d[i][b] * wJac;
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