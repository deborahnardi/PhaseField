#include "../headers/Element.h"
#include <iomanip>

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
PetscErrorCode Solid2D::getContribution(std::array<Tensor, 3> tensors, Mat &A, Vec &rhs, bool negativeLoad, bool _PrescribedDamageField)
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
        // double tensorCPlus[3][3] = {};
        // double tensorCMinus[3][3] = {};
        // double tensorC[3][3] = {};

        // if (divU > 0)
        // {
        //     for (int i = 0; i < 3; i++)
        //         for (int j = 0; j < 3; j++)
        //         {
        //             tensorCPlus[i][j] = dCoeff * (kappa * tensorK[i][j] + 2.0 * mu * tensorJ[i][j]);
        //             tensorCMinus[i][j] = 0.0;
        //             tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];
        //         }
        // }
        // else if (divU <= 0)
        // {
        //     for (int i = 0; i < 3; i++)
        //         for (int j = 0; j < 3; j++)
        //         {
        //             tensorCPlus[i][j] = dCoeff * 2.0 * mu * tensorJ[i][j];
        //             tensorCMinus[i][j] = kappa * tensorK[i][j];
        //             tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];
        //         }
        // }
        SplitModel splitModel = params->getSplitModel();
        std::array<std::array<double, 3>, 3> tensorC = computeConstitutiveTensor(splitModel, tensors, gradU, divU, kappa, mu, dCoeff);

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

        localStiffnessMatrix[0] += (B[0][0] * B[0][0] * tensorC[0][0] + B[1][1] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[1] += (B[0][0] * B[1][1] * tensorC[0][1] + B[0][0] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[2] += (B[0][0] * B[0][2] * tensorC[0][0] + B[1][1] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[3] += (B[0][0] * B[1][3] * tensorC[0][1] + B[0][2] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[4] += (B[0][0] * B[0][4] * tensorC[0][0] + B[1][1] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[5] += (B[0][0] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][1] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[6] = localStiffnessMatrix[1];
        localStiffnessMatrix[7] += (B[0][0] * B[0][0] * tensorC[2][2] + B[1][1] * B[1][1] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[8] += (B[0][0] * B[1][3] * tensorC[2][2] + B[0][2] * B[1][1] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[9] += (B[0][0] * B[0][2] * tensorC[2][2] + B[1][1] * B[1][3] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[10] += (B[0][0] * B[1][5] * tensorC[2][2] + B[0][4] * B[1][1] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[11] += (B[0][0] * B[0][4] * tensorC[2][2] + B[1][1] * B[1][5] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[12] = localStiffnessMatrix[2];
        localStiffnessMatrix[13] = localStiffnessMatrix[8];
        localStiffnessMatrix[14] += (B[0][2] * B[0][2] * tensorC[0][0] + B[1][3] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[15] += (B[0][2] * B[1][3] * tensorC[0][1] + B[0][2] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[16] += (B[0][2] * B[0][4] * tensorC[0][0] + B[1][3] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[17] += (B[0][2] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][3] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[18] = localStiffnessMatrix[3];
        localStiffnessMatrix[19] = localStiffnessMatrix[9];
        localStiffnessMatrix[20] = localStiffnessMatrix[15];
        localStiffnessMatrix[21] += (B[0][2] * B[0][2] * tensorC[2][2] + B[1][3] * B[1][3] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[22] += (B[0][2] * B[1][5] * tensorC[2][2] + B[0][4] * B[1][3] * tensorC[1][0]) * wJac;
        localStiffnessMatrix[23] += (B[0][2] * B[0][4] * tensorC[2][2] + B[1][3] * B[1][5] * tensorC[1][1]) * wJac;
        localStiffnessMatrix[24] = localStiffnessMatrix[4];
        localStiffnessMatrix[25] = localStiffnessMatrix[10];
        localStiffnessMatrix[26] = localStiffnessMatrix[16];
        localStiffnessMatrix[27] = localStiffnessMatrix[22];
        localStiffnessMatrix[28] += (B[0][4] * B[0][4] * tensorC[0][0] + B[1][5] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[29] += (B[0][4] * B[1][5] * tensorC[0][1] + B[0][4] * B[1][5] * tensorC[2][2]) * wJac;
        localStiffnessMatrix[30] = localStiffnessMatrix[5];
        localStiffnessMatrix[31] = localStiffnessMatrix[11];
        localStiffnessMatrix[32] = localStiffnessMatrix[17];
        localStiffnessMatrix[33] = localStiffnessMatrix[23];
        localStiffnessMatrix[34] = localStiffnessMatrix[29];
        localStiffnessMatrix[35] += (B[0][4] * B[0][4] * tensorC[2][2] + B[1][5] * B[1][5] * tensorC[1][1]) * wJac;

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

std::vector<double> Solid2D::getStiffnessIIOrIJ(std::array<Tensor, 3> tensors, const int idxLocalNode1, const int idxLocalNode2, SplitModel splitModel, int _stepGlobal, bool _PrescribedDamageField)
{
    /*
        THIS METHOD IS USED TO COMPUTE THE CONTRIBUTION OF A SPECIFIC ELEMENT TO THE GLOBAL MATRIX AND RHS VECTOR

        elemInfo: Vector containing the information of the element to be computed. All the elements that one node is connected to are stored in elemInfo.
                  elemInfo.size() = 3 * numSharedElements = 3 * [elemIndex, localNode1Index, localNode2Index];
                  For example, if a node is connected to 4 elements, elemInfo.size() = 12.
                  Consider the example below: [0, 1, 1, 1, 0, 0, 2, 0, 0, 3, 0, 0] -> 0, 1, 1: element 0, local node 1, local node 2
                                                                                   1, 0, 0: element 1, local node 0, local node 0
                                                                                   2, 0, 0: element 2, local node 0, local node 0
                                                                                   3, 0, 0: element 3, local node 0, local node 0

        THE METHOD RETURNS THE STIFFNESS CONTRIBUTION VALUE OF THE GLOBAL NODE, ALREADY SUMMING ALL THE ELEMENTS THAT THE NODE IS CONNECTED TO, I.E, ALL THE LOCAL CONTRIBUTIONS

        -----------------------------------------------------------------------------------------------------------------------------------------------------------

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

    PetscScalar mu = material->getShearModulus();
    PetscScalar kappa = 2.0 / 3.0 * mu + material->getLameConstant();

    int nDOF = 2;

    std::vector<double> values;
    if (idxLocalNode1 == idxLocalNode2)
        values.resize(3, 0.0);
    else
        values.resize(4, 0.0);

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

        PetscScalar divU = gradU[0][0] + gradU[1][1]; // uk,k and ul,l are the same

        // COMPUTING THE DEGRADATION FUNCTION (1-d)^2
        double damageValue = 0.;
        for (PetscInt a = 0; a < numElNodes; a++)
            damageValue += N[a] * elemConnectivity[a]->getDOFs()[2]->getDamageValue(); // DamageValue -> dstag = dn + delta_d^i

        PetscReal dCoeff = pow(1 - damageValue, 2);

        // ASSEMBLING CONSTITUTIVE MATRIX CONSIDERING THE ENERGY SPLIT
        // COMPUTING THE CONSTITUTIVE TENSOR

        std::array<std::array<double, 3>, 3> tensorC = computeConstitutiveTensor(splitModel, tensors, gradU, divU, kappa, mu, dCoeff);

        // RELATING dN_dX TO THE B MATRIX
        PetscReal B[3][6] = {};
        for (PetscInt a = 0; a < numElNodes; a++)
        {
            B[0][2 * a] = dN_dX[a][0];
            B[1][2 * a + 1] = dN_dX[a][1];
            B[2][2 * a] = dN_dX[a][1];
            B[2][2 * a + 1] = dN_dX[a][0];
        }

        /*
            COMPUTE ONLY THE NEEDED VALUES - THIS INFORMATION COMES FROM idxLocalNode1 AND idxLocalNode2
            For II only 3 values are needed;
        */

        int count = 0;

        /*
                Note that if idxLocalNode1 = 0 and idxLocalNode2 = 0, for example, we need to compute Ke_00, Ke_01, Ke_11
                if idxLocalNode1 = 0 and idxLocalNode2 = 1, for example, we need to compute Ke_02, Ke_03, Ke_12 and Ke_13
        */

        // COMPUTE B^T C B -> Ke_kl = B_ki C_ij B_lj
        if (idxLocalNode1 == idxLocalNode2) // Only symmetric part is assembled
            for (int iDir = 0; iDir < nDOF; iDir++)
                for (int jDir = iDir; jDir < nDOF; jDir++)
                {
                    int localRow = nDOF * idxLocalNode1 + iDir; // k
                    int localCol = nDOF * idxLocalNode2 + jDir; // l

                    for (int i = 0; i < 3; i++)
                        for (int j = i; j < 3; j++)
                            values[count] += B[i][localRow] * tensorC[i][j] * B[j][localCol] * wJac;

                    count++;
                }
        else if (idxLocalNode1 != idxLocalNode2)
            for (int iDir = 0; iDir < nDOF; iDir++)
                for (int jDir = 0; jDir < nDOF; jDir++)
                {
                    int localRow = nDOF * idxLocalNode1 + iDir;
                    int localCol = nDOF * idxLocalNode2 + jDir;

                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            values[count] += B[i][localRow] * tensorC[i][j] * B[j][localCol] * wJac;

                    count++;
                }

        // for (std::size_t i = 0; i < values.size(); ++i)
        // {
        //     if (std::isnan(values[i]))
        //         std::cout << "⚠️  NaN detectado em values[" << i << "] dentro de getStiffnessIIOrIJ!" << std::endl;
        // }

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    return values;
}

std::array<std::array<double, 3>, 3> Solid2D::computeConstitutiveTensor(SplitModel splitModel, std::array<Tensor, 3> tensors, PetscScalar gradU[2][2], double divU, double kappa, double mu, double dCoeff)
{
    /*
        Compute Constitutive Tensor based on the energy split model. Implemented split energy models:tensorCPlus
        volDev: stands for volumetric deviatoric split model;
        spectral: stands for spectral split model;
    */

    Tensor tensorK = tensors[0];
    Tensor tensorI = tensors[1];
    Tensor tensorJ = tensors[2];

    std::array<std::array<double, 3>, 3> tensorCPlus = {};
    std::array<std::array<double, 3>, 3> tensorCMinus = {};
    std::array<std::array<double, 3>, 3> tensorC = {};

    switch (splitModel)
    {
    case volDev:
    {
        // std::ofstream teste;
        // teste.open("testevoldev.txt", std::ios::out | std::ios::app);
        // teste << 1 << "\n";
        // teste.close();

        if (divU > 0)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tensorCPlus[i][j] = dCoeff * (kappa * tensorK[i][j] + 2.0 * mu * tensorJ[i][j]);
                    tensorCMinus[i][j] = 0.0;
                    tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];
                }
        }
        else if (divU <= 0)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tensorCPlus[i][j] = dCoeff * 2.0 * mu * tensorJ[i][j];
                    tensorCMinus[i][j] = kappa * tensorK[i][j];
                    tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];
                }
        }

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                tensorC[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];

        break;
    }

    case spectral:
    {
        double strain[3][3] = {};
        for (PetscInt i = 0; i < 2; i++)
            for (PetscInt j = 0; j < 2; j++)
                strain[i][j] = 0.5 * (gradU[i][j] + gradU[j][i]);

        switch (material->getPlaneAnalysis())
        {
        case PLANE_STRESS:
        {
            double _poisson = material->getPoisson();
            strain[2][2] = -(_poisson / (1.0 - _poisson)) * (strain[0][0] + strain[1][1]);
            break;
        }
        default:
            strain[2][2] = 0.0; // PLANE_STRAIN
            break;
        }

        strain[0][0] = 4.201433961241174e-004;
        strain[1][1] = 1.803278624210739e-004;
        strain[0][1] = 1.479491857090474e-003;
        strain[1][0] = strain[0][1];

        /*
            delta < 0 -> there are 3 distinct real eigenvalues; delta = 0 -> there are 2 equal eigenvalues; delta > 0 -> there is one real eigenvalue and two complex conjugate eigenvalues; we only work with cases where delta < 0 and delta = 0
        */

        double lame = material->getLameConstant(); // 1.0
        double mu = material->getShearModulus();   // 0.25

        double volStrain = (1.0 / 3.0) * (strain[0][0] + strain[1][1] + strain[2][2]);

        double d11 = strain[0][0] - volStrain;
        double d22 = strain[1][1] - volStrain;
        double d33 = strain[2][2] - volStrain;
        double d12 = strain[0][1];
        double d21 = strain[1][0];
        double d23 = strain[1][2];
        double d32 = strain[2][1];
        double d31 = strain[2][0];
        double d13 = strain[0][2];

        double eeq = sqrt(2.0 / 3.0 * (d11 * d11 + d22 * d22 + d33 * d33 + 2.0 * (d12 * d12 + d13 * d13 + d23 * d23))); // Equivalent strain

        double detdd =
            d11 * (d22 * d33 - d23 * d32) -
            d12 * (d21 * d33 - d23 * d31) +
            d13 * (d21 * d32 - d22 * d31);

        double const tol = 1e-12;
        double etaReal = 0.0, cos3eta = 0.0;

        if (std::abs(eeq) < tol)
        {
            etaReal = 0.0;
        }
        else
        {
            cos3eta = 4.0 * detdd / (eeq * eeq * eeq);
            cos3eta = std::clamp(cos3eta, -1.0, 1.0);
            // cos3eta = std::round(cos3eta * 1e10) / 1e10;
            etaReal = (1.0 / 3.0) * std::acos(cos3eta);
        }

        double etaStar[3] = {};
        etaStar[0] = etaReal;
        etaStar[1] = etaReal - 2.0 * M_PI / 3.0;
        etaStar[2] = etaReal + 2.0 * M_PI / 3.0;

        double d_ep_de[3][3] = {};      // First derivative of the energy functional with respect to the strain tensor
        double d2_ep_de2[3][3][3] = {}; // Second derivative of the energy functional with respect to the strain tensor

        double ident[3] = {};
        ident[0] = 1.0;
        ident[1] = 1.0;
        ident[2] = 0.0;

        // Cofactor matrix of the deviatoric strain tensor
        double dcof11 = d22 * d33 - d23 * d23;
        double dcof12 = d23 * d31 - d12 * d33;
        double dcof13 = d12 * d23 - d22 * d31;
        double dcof21 = dcof12;
        double dcof22 = d11 * d33 - d31 * d31;
        double dcof23 = d12 * d31 - d11 * d23;
        double dcof31 = dcof13;
        double dcof32 = dcof23;
        double dcof33 = d11 * d22 - d12 * d12;

        double dcofTrace = dcof11 + dcof22 + dcof33;

        // COMPUTING THE PRINCIPAL VALUES VIA LODE ANGLE

        double principalValues[3] = {};
        principalValues[0] = volStrain + eeq * std::cos(etaStar[0]); // ep1
        principalValues[1] = volStrain + eeq * std::cos(etaStar[1]); // ep2
        principalValues[2] = volStrain + eeq * std::cos(etaStar[2]); // ep3

        // for (int ip = 0; ip < 3; ip++)
        // {
        //     double cosEta = std::cos(etaStar[ip]);
        //     double val = volStrain + eeq * cosEta;
        //     principalValues[ip] = std::round(val * 1e10) / 1e10;
        // }

        bool strainIsNull = true;
        const double pTol = 1e-12;

        if (std::abs(principalValues[0]) > pTol ||
            std::abs(principalValues[1]) > pTol ||
            std::abs(principalValues[2]) > pTol)
        {
            strainIsNull = false;
        }

        if (strainIsNull)
        {

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                {
                    double outerProductIdent = ident[i] * ident[j];
                    tensorC[i][j] = lame * outerProductIdent; // C = lambda * I
                }
            break;
        }

        // If the code goes over this point, the strain is not null, so we can compute an equivalent measure of the strain tensor

        double a11 = strain[0][0] / eeq;
        double a22 = strain[1][1] / eeq;
        double a12 = strain[0][1] / eeq;
        double a33 = strain[2][2] / eeq;
        double a21 = strain[1][0] / eeq;
        double a23 = strain[1][2] / eeq;
        double a32 = strain[2][1] / eeq;
        double a31 = strain[2][0] / eeq;
        double a13 = strain[0][2] / eeq;

        double volStrainA = (1.0 / 3.0) * (a11 + a22 + a33);

        double ad11 = a11 - volStrainA;
        double ad22 = a22 - volStrainA;
        double ad33 = a33 - volStrainA;
        double ad12 = a12;
        double ad21 = a21;
        double ad23 = a23;
        double ad32 = a32;
        double ad31 = a31;
        double ad13 = a13;

        double adCof11 = ad22 * ad33 - ad23 * ad23;
        double adCof12 = ad23 * ad31 - ad12 * ad33;
        double adCof13 = ad12 * ad23 - ad22 * ad31;
        double adCof21 = adCof12;
        double adCof22 = ad11 * ad33 - ad31 * ad31;
        double adCof23 = ad12 * ad31 - ad11 * ad23;
        double adCof31 = adCof13;
        double adCof32 = adCof23;
        double adCof33 = ad11 * ad22 - ad12 * ad12;

        double adCofTrace = adCof11 + adCof22 + adCof33;
        double detddA =
            ad11 * (ad22 * ad33 - ad23 * ad32) -
            ad12 * (ad21 * ad33 - ad23 * ad31) +
            ad13 * (ad21 * ad32 - ad22 * ad31);

        // d_eeq_deij--------------------------------------------
        double d_eeq_de[3] = {};

        // d_eeq_de[0] = 2.0 * d11 / (3.0 * eeq);
        // d_eeq_de[1] = 2.0 * d22 / (3.0 * eeq);
        // d_eeq_de[2] = 2.0 * d12 / (3.0 * eeq);

        d_eeq_de[0] = 2.0 / 3.0 * ad11;
        d_eeq_de[1] = 2.0 / 3.0 * ad22;
        d_eeq_de[2] = 2.0 / 3.0 * ad12;

        // d2_eeq_deij_dekl -------------------------------------------------

        double d2_eeq_de2[3][3] = {};

        // d2_eeq_de2[0][0] = 4.0 / (9.0 * eeq) * (1.0 - d11 * d11 / (eeq * eeq));
        // d2_eeq_de2[0][1] = -2.0 / (9.0 * eeq) * (1.0 + 2.0 * d11 * d22 / (eeq * eeq));
        // d2_eeq_de2[0][2] = -4.0 / 9.0 * (d11 * d12 / (eeq * eeq * eeq));

        // d2_eeq_de2[1][0] = d2_eeq_de2[0][1];
        // d2_eeq_de2[1][1] = 4.0 / (9.0 * eeq) * (1.0 - d22 * d22 / (eeq * eeq));
        // d2_eeq_de2[1][2] = -4.0 / 9.0 * (d22 * d12 / (eeq * eeq * eeq));

        // d2_eeq_de2[2][0] = d2_eeq_de2[0][2];
        // d2_eeq_de2[2][1] = d2_eeq_de2[1][2];
        // d2_eeq_de2[2][2] = 1.0 / (3.0 * eeq) - 4.0 / 9.0 * d12 * d12 / (eeq * eeq * eeq);

        d2_eeq_de2[0][0] = 4.0 / 9.0 / eeq * (1.0 - ad11 * ad11);
        d2_eeq_de2[0][1] = -2.0 / 9.0 / eeq * (1.0 + 2.0 * ad11 * ad22);
        d2_eeq_de2[0][2] = -4.0 / 9.0 / eeq * (ad11 * ad12);

        d2_eeq_de2[1][0] = d2_eeq_de2[0][1];
        d2_eeq_de2[1][1] = 4.0 / 9.0 / eeq * (1.0 - ad22 * ad22);
        d2_eeq_de2[1][2] = -4.0 / 9.0 / eeq * (ad22 * ad12);

        d2_eeq_de2[2][0] = d2_eeq_de2[0][2];
        d2_eeq_de2[2][1] = d2_eeq_de2[1][2];
        d2_eeq_de2[2][2] = 1.0 / 3.0 / eeq - 4.0 / 9.0 / eeq * (ad12 * ad12);

        //--
        double d_cos3eta_de[3] = {};
        // d_cos3eta_de[0] = 4.0 * ((eeq * eeq * dcof11 - 2.0 * detdd * d11) / pow(eeq, 5.0) - 1.0 / 3.0 * dcofTrace / pow(eeq, 3.0));
        // d_cos3eta_de[1] = 4.0 * ((eeq * eeq * dcof22 - 2.0 * detdd * d22) / pow(eeq, 5.0) - 1.0 / 3.0 * dcofTrace / pow(eeq, 3.0));
        // d_cos3eta_de[2] = 4.0 * ((eeq * eeq * dcof12 - 2.0 * detdd * d12) / (pow(eeq, 5.0)));

        d_cos3eta_de[0] = 4.0 / eeq * (adCof11 - 2.0 * ad11 * detddA - 1.0 / 3.0 * adCofTrace);
        d_cos3eta_de[1] = 4.0 / eeq * (adCof22 - 2.0 * ad22 * detddA - 1.0 / 3.0 * adCofTrace);
        d_cos3eta_de[2] = 4.0 / eeq * (adCof12 - 2.0 * ad12 * detddA);

        double d2_cos3eta_de2[3][3] = {};
        // d2_cos3eta_de2[0][0] = 4.0 * (-2.0 / 3.0 * dcof11 * d11 / pow(eeq, 5.0) - 4.0 * detdd / (3.0 * pow(eeq, 5.0)) + 2.0 / 3.0 * d11 * dcofTrace / pow(eeq, 5.0) - 10.0 / 3.0 * d11 / pow(eeq, 7.0) * (pow(eeq, 2.0) * dcof11 - 2.0 * d11 * detdd) + 2.0 / 3.0 * d11 / pow(eeq, 3.0) * (1.0 + dcofTrace / pow(eeq, 2.0)));
        // d2_cos3eta_de2[0][1] = 4.0 * (1.0 / pow(eeq, 3.0) * (d33 + d11 / 3.0) + 4.0 / 3.0 * dcof11 * d22 / pow(eeq, 5.0) + 2.0 / 3.0 * detdd / pow(eeq, 5.0) - 2.0 * d11 / pow(eeq, 5.0) * (dcof22 - dcofTrace / 3.0) - 10.0 / 3.0 * d22 / pow(eeq, 7.0) * (eeq * eeq * dcof11 - 2.0 * d11 * detdd) + d22 / pow(eeq, 3.0) / 3.0 * (1.0 + 2.0 * dcofTrace / pow(eeq, 2.0)));
        // d2_cos3eta_de2[0][2] = 4.0 * (4.0 / 3.0 * dcof11 * d12 / pow(eeq, 5.0) - 2.0 * d11 * dcof12 / pow(eeq, 5.0) - 10.0 * d12 / (3.0 * pow(eeq, 7)) * (eeq * eeq * dcof11 - 2.0 * d11 * detdd) + d12 / (3.0 * pow(eeq, 3.0)) * (1.0 + 2.0 * dcofTrace / pow(eeq, 2.0)));

        // d2_cos3eta_de2[1][0] = d2_cos3eta_de2[0][1];
        // d2_cos3eta_de2[1][1] = 4.0 * (-2.0 / 3.0 * dcof22 * d22 / pow(eeq, 5.0) - 4.0 / 3.0 * detdd / pow(eeq, 5.0) + 2.0 / 3.0 * d22 * dcofTrace / pow(eeq, 5.0) - 10.0 / 3.0 * d22 / pow(eeq, 7.0) * (eeq * eeq * dcof22 - 2.0 * d22 * detdd) + 2.0 / 3.0 * d22 / pow(eeq, 3.0) * (1.0 + dcofTrace / pow(eeq, 2.0)));
        // d2_cos3eta_de2[1][2] = 4.0 * (4.0 / 3.0 * dcof22 * d12 / pow(eeq, 5.0) - 2.0 * d22 * dcof12 / pow(eeq, 5.0) - 10.0 / 3.0 * d12 / pow(eeq, 7.0) * (eeq * eeq * dcof22 - 2.0 * d22 * detdd) + d12 / (3.0 * pow(eeq, 3.0)) * (1.0 + 2.0 * dcofTrace / pow(eeq, 2.0)));

        // d2_cos3eta_de2[2][1] = d2_cos3eta_de2[1][2];
        // d2_cos3eta_de2[2][0] = d2_cos3eta_de2[0][2];
        // d2_cos3eta_de2[2][2] = 4.0 * (-d33 / 2.0 / pow(eeq, 3.0) - 2.0 / 3.0 * dcof12 * d12 / pow(eeq, 5.0) - detdd / pow(eeq, 5.0) - 10.0 / 3.0 * d12 / pow(eeq, 7.0) * (eeq * eeq * dcof12 - 2.0 * d12 * detdd));

        d2_cos3eta_de2[0][0] = 4.0 / (eeq * eeq) * (-2.0 / 3.0 * adCof11 * ad11 - 4.0 * detddA / 3.0 + 2.0 / 3.0 * ad11 * adCofTrace - 10.0 / 3.0 * ad11 * (adCof11 - 2.0 * ad11 * detddA) + 2.0 / 3.0 * ad11 * (1.0 + adCofTrace));
        d2_cos3eta_de2[0][1] = 4.0 / (eeq * eeq) * (ad33 + ad11 / 3.0 + 4.0 / 3.0 * adCof11 * ad22 + 2.0 / 3.0 * detddA - 2.0 * ad11 * (adCof22 - adCofTrace / 3.0) - 10.0 / 3.0 * ad22 * (adCof11 - 2.0 * ad11 * detddA) + ad22 / 3.0 * (1.0 + 2.0 * adCofTrace));
        d2_cos3eta_de2[0][2] = 4.0 / (eeq * eeq) * (4.0 / 3.0 * adCof11 * ad12 - 2 * ad11 * adCof12 - 10.0 * ad12 / 3.0 * (adCof11 - 2.0 * ad11 * detddA) + ad12 / 3.0 * (1.0 + 2.0 * adCofTrace));

        d2_cos3eta_de2[1][0] = d2_cos3eta_de2[0][1];
        d2_cos3eta_de2[1][1] = 4.0 / (eeq * eeq) * (-2.0 / 3.0 * adCof22 * ad22 - 4.0 / 3.0 * detddA + 2.0 / 3.0 * ad22 * adCofTrace - 10.0 / 3.0 * ad22 * (adCof22 - 2.0 * ad22 * detddA) + 2.0 / 3.0 * ad22 * (1.0 + adCofTrace));
        d2_cos3eta_de2[1][2] = 4.0 / (eeq * eeq) * (4.0 / 3.0 * adCof22 * ad12 - 2.0 * ad22 * adCof12 - 10.0 / 3.0 * ad12 * (adCof22 - 2.0 * ad22 * detddA) + ad12 / 3.0 * (1.0 + 2.0 * adCofTrace));

        d2_cos3eta_de2[2][0] = d2_cos3eta_de2[0][2];
        d2_cos3eta_de2[2][1] = d2_cos3eta_de2[1][2];
        d2_cos3eta_de2[2][2] = 4.0 / (eeq * eeq) * (-ad33 / 2.0 - 2.0 / 3.0 * adCof12 * ad12 - detddA - 10.0 / 3.0 * ad12 * (adCof12 - 2.0 * ad12 * detddA));

        //--

        // -------------------------------------------------

        if ((abs(principalValues[1] - principalValues[2]) > tol) &&
            (abs(principalValues[0] - principalValues[1]) > tol) &&
            (abs(principalValues[0] - principalValues[2]) > tol)) // THERE ARE 3 DISTINCT REAL EIGENVALUES
            for (int ip = 0; ip < 3; ip++)                        // ip stands for the principal value
            {
                // std::cout << "Entered the case 3 distinct real eigenvalues! " << std::endl;
                double eta = etaStar[ip];
                double cosEta = std::cos(eta);
                cosEta = std::round(cosEta * 1e10) / 1e10;

                // -------------------------------------------------

                double d_coseta_de[3] = {};
                for (int i = 0; i < 3; i++)
                    d_coseta_de[i] = d_cos3eta_de[i] / (12.0 * cosEta * cosEta - 3.0);

                // -------------------------------------------------

                double d2_coseta_de2[3][3] = {};

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double outerProduct1 = d_cos3eta_de[i] * d_coseta_de[j];
                        d2_coseta_de2[i][j] = (d2_cos3eta_de2[i][j] * (12.0 * cosEta * cosEta - 3.0) - 24.0 * cosEta * outerProduct1) / pow((12.0 * cosEta * cosEta - 3.0), 2.0);
                    }

                // -------------------------------------------------

                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double outerProduct2 = d_eeq_de[i] * d_coseta_de[j];
                        double outerProduct3 = d_coseta_de[i] * d_eeq_de[j];
                        d2_ep_de2[ip][i][j] = cosEta * d2_eeq_de2[i][j] + outerProduct2 + outerProduct3 + eeq * d2_coseta_de2[i][j];
                    }

                // -------------------------------------------------

                for (int i = 0; i < 3; i++)
                    d_ep_de[ip][i] = 1.0 / 3.0 * ident[i] + cosEta * d_eeq_de[i] + eeq * d_coseta_de[i];
            }
        else // THERE ARE 2 EQUAL EIGENVALUES
        {
            // std::cout << "Enterted the case of 2 equal eigenvalues. Delta: " << delta << std::endl;
            if (abs(principalValues[1] - principalValues[2]) < tol) // cos(3eta)=1; eta =0
            {
                d_ep_de[0][0] = 1.0 / 3.0 + d_eeq_de[0];
                d_ep_de[0][1] = 1.0 / 3.0 + d_eeq_de[1];
                d_ep_de[0][2] = d_eeq_de[2];

                d_ep_de[2][0] = 2.0 / 3.0 - d_eeq_de[0];
                d_ep_de[2][1] = 2.0 / 3.0 - d_eeq_de[1];
                d_ep_de[2][2] = -d_eeq_de[2];

                d_ep_de[1][0] = 0.0;
                d_ep_de[1][1] = 0.0;
                d_ep_de[1][2] = 0.0;

                // ------------------------------------------------
                double d_eta_de[3] = {};
                d_eta_de[0] = -2.0 / (3.0 * sqrt(3.0) * eeq) + 1.0 / (sqrt(3.0) * eeq) * d_eeq_de[0];
                d_eta_de[1] = -2.0 / (3.0 * sqrt(3.0) * eeq) + 1.0 / (sqrt(3.0) * eeq) * d_eeq_de[1];
                d_eta_de[2] = 1.0 / (sqrt(3.0) * eeq) * d_eeq_de[2];
                // ------------------------------------------------

                d2_ep_de2[0][0][0] = d2_eeq_de2[0][0] - eeq * d_eta_de[0] * d_eta_de[0];
                d2_ep_de2[0][0][1] = d2_eeq_de2[0][1] - eeq * d_eta_de[1] * d_eta_de[1];
                d2_ep_de2[0][0][2] = d2_eeq_de2[0][2] - eeq * d_eta_de[2] * d_eta_de[0];

                d2_ep_de2[0][1][0] = d2_eeq_de2[1][0] - eeq * d_eta_de[0] * d_eta_de[0];
                d2_ep_de2[0][1][1] = d2_eeq_de2[1][1] - eeq * d_eta_de[1] * d_eta_de[1];
                d2_ep_de2[0][1][2] = d2_eeq_de2[1][2] - eeq * d_eta_de[2] * d_eta_de[1];

                d2_ep_de2[0][2][0] = d2_eeq_de2[2][0] - eeq * d_eta_de[2] * d_eta_de[0];
                d2_ep_de2[0][2][1] = d2_eeq_de2[2][1] - eeq * d_eta_de[2] * d_eta_de[1];
                d2_ep_de2[0][2][2] = d2_eeq_de2[2][2] - eeq * d_eta_de[2] * d_eta_de[2];
                //
                d2_ep_de2[2][0][0] = -d2_ep_de2[0][0][0];
                d2_ep_de2[2][0][1] = -d2_ep_de2[0][0][1];
                d2_ep_de2[2][0][2] = -d2_ep_de2[0][0][2];
                d2_ep_de2[2][1][0] = -d2_ep_de2[0][1][0];
                d2_ep_de2[2][1][1] = -d2_ep_de2[0][1][1];
                d2_ep_de2[2][1][2] = -d2_ep_de2[0][1][2];
                d2_ep_de2[2][2][0] = -d2_ep_de2[0][2][0];
                d2_ep_de2[2][2][1] = -d2_ep_de2[0][2][1];
                d2_ep_de2[2][2][2] = -d2_ep_de2[0][2][2];
                //
                d2_ep_de2[1][0][0] = 0.0;
                d2_ep_de2[1][0][1] = 0.0;
                d2_ep_de2[1][0][2] = 0.0;
                d2_ep_de2[1][1][0] = 0.0;
                d2_ep_de2[1][1][1] = 0.0;
                d2_ep_de2[1][1][2] = 0.0;
                d2_ep_de2[1][2][0] = 0.0;
                d2_ep_de2[1][2][1] = 0.0;
                d2_ep_de2[1][2][2] = 0.0;
            }
            else // eta = n Pi / 3, where n is an integer (multiple of 60 degrees)
            {
                // del ep1 del eij
                d_ep_de[0][0] = 2.0 / 3.0 + d_eeq_de[0];
                d_ep_de[0][1] = 2.0 / 3.0 + d_eeq_de[1];
                d_ep_de[0][2] = d_eeq_de[2];
                // del ep2 del eij
                d_ep_de[2][0] = 1.0 / 3.0 - d_eeq_de[0];
                d_ep_de[2][1] = 1.0 / 3.0 - d_eeq_de[1];
                d_ep_de[2][2] = -d_eeq_de[2];
                // del ep3 del eij
                d_ep_de[1][0] = 0.0;
                d_ep_de[1][1] = 0.0;
                d_ep_de[1][2] = 0.0;
                // ------------------------------------------------

                double d_eta_de[3] = {};
                d_eta_de[0] = -2.0 / (3.0 * sqrt(3.0) * eeq) - 1.0 / (sqrt(3) * eeq) * d_eeq_de[0];
                d_eta_de[1] = -2.0 / (3.0 * sqrt(3.0) * eeq) - 1.0 / (sqrt(3) * eeq) * d_eeq_de[1];
                d_eta_de[2] = -1.0 / (sqrt(3.0) * eeq) * d_eeq_de[2];
                // ------------------------------------------------

                // # d2_ep_de2
                d2_ep_de2[0][0][0] = d2_eeq_de2[0][0] - eeq * d_eta_de[0] * d_eta_de[0];
                d2_ep_de2[0][0][1] = d2_eeq_de2[0][1] - eeq * d_eta_de[1] * d_eta_de[0];
                d2_ep_de2[0][0][2] = d2_eeq_de2[0][2] - eeq * d_eta_de[2] * d_eta_de[0];

                d2_ep_de2[0][1][0] = d2_eeq_de2[1][0] - eeq * d_eta_de[0] * d_eta_de[1];
                d2_ep_de2[0][1][1] = d2_eeq_de2[1][1] - eeq * d_eta_de[1] * d_eta_de[1];
                d2_ep_de2[0][1][2] = d2_eeq_de2[1][2] - eeq * d_eta_de[2] * d_eta_de[1];

                d2_ep_de2[0][2][0] = d2_eeq_de2[2][0] - eeq * d_eta_de[2] * d_eta_de[0];
                d2_ep_de2[0][2][1] = d2_eeq_de2[2][1] - eeq * d_eta_de[2] * d_eta_de[1];
                d2_ep_de2[0][2][2] = d2_eeq_de2[2][2] - eeq * d_eta_de[2] * d_eta_de[2];
                //
                d2_ep_de2[2][0][0] = -d2_ep_de2[0][0][0];
                d2_ep_de2[2][0][1] = -d2_ep_de2[0][0][1];
                d2_ep_de2[2][0][2] = -d2_ep_de2[0][0][2];
                d2_ep_de2[2][1][0] = -d2_ep_de2[0][1][0];
                d2_ep_de2[2][1][1] = -d2_ep_de2[0][1][1];
                d2_ep_de2[2][1][2] = -d2_ep_de2[0][1][2];
                d2_ep_de2[2][2][0] = -d2_ep_de2[0][2][0];
                d2_ep_de2[2][2][1] = -d2_ep_de2[0][2][1];
                d2_ep_de2[2][2][2] = -d2_ep_de2[0][2][2];
                //
                d2_ep_de2[1][0][0] = 0.0;
                d2_ep_de2[1][0][1] = 0.0;
                d2_ep_de2[1][0][2] = 0.0;
                d2_ep_de2[1][1][0] = 0.0;
                d2_ep_de2[1][1][1] = 0.0;
                d2_ep_de2[1][1][2] = 0.0;
                d2_ep_de2[1][2][0] = 0.0;
                d2_ep_de2[1][2][1] = 0.0;
                d2_ep_de2[1][2][2] = 0.0;
            }
        }
        // -------------------------------------------------

        // COMPUTING THE CONSTITUTIVE TENSOR D

        double d_ep_de_pos[3][3] = {};
        double d_ep_de_neg[3][3] = {};

        for (int ip = 0; ip < 3; ip++)
        {
            if (principalValues[ip] > pTol)
            {
                for (int j = 0; j < 3; j++)
                {
                    d_ep_de_pos[ip][j] = d_ep_de[ip][j];
                    d_ep_de_neg[ip][j] = 0.0;
                }
            }
            else
            {
                for (int j = 0; j < 3; j++)
                {
                    d_ep_de_pos[ip][j] = 0.0;
                    d_ep_de_neg[ip][j] = d_ep_de[ip][j];
                }
            }
        }

        double d2_ep_de2_pos[3][3][3] = {};
        double d2_ep_de2_neg[3][3][3] = {};

        for (int ip = 0; ip < 3; ip++)
        {
            if (principalValues[ip] > pTol)
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        d2_ep_de2_pos[ip][i][j] = d2_ep_de2[ip][i][j];
                        d2_ep_de2_neg[ip][i][j] = 0.0;
                    }
            }
            else
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        d2_ep_de2_pos[ip][i][j] = 0.0;
                        d2_ep_de2_neg[ip][i][j] = d2_ep_de2[ip][i][j];
                    }
            }
        }

        double outerProductEpEpPos[3][3][3] = {};
        double outerProductEpEpNeg[3][3][3] = {};
        for (int ip = 0; ip < 3; ip++)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    outerProductEpEpPos[ip][i][j] = d_ep_de_pos[ip][i] * d_ep_de_pos[ip][j];
                    outerProductEpEpNeg[ip][i][j] = d_ep_de_neg[ip][i] * d_ep_de_neg[ip][j];
                }

        // COMPUTING THE POSITIVE AND NEGATIVE PART OF THE TENSOR
        for (int ip = 0; ip < 3; ip++)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tensorCPlus[i][j] += 2 * mu * (outerProductEpEpPos[ip][i][j] + principalValues[ip] * d2_ep_de2_pos[ip][i][j]);
                    tensorCMinus[i][j] += 2 * mu * (outerProductEpEpNeg[ip][i][j] + principalValues[ip] * d2_ep_de2_neg[ip][i][j]);
                }
        }

        if (divU > pTol)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    double outerProductIdent = ident[i] * ident[j];
                    tensorCPlus[i][j] += lame * outerProductIdent;
                }
        else if (divU < -pTol)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    double outerProductIdent = ident[i] * ident[j];
                    tensorCMinus[i][j] += lame * outerProductIdent;
                }
        // double tensorCaux[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                tensorC[i][j] = dCoeff * tensorCPlus[i][j] + tensorCMinus[i][j];
                // tensorCaux[i][j] = tensorCPlus[i][j] + tensorCMinus[i][j];
            }

        // COMPUTE DETERMINANT OF THE TENSORCPLUS + TENSORCMINUS
        // double detCPlusMinus = 0.0;

        // -------------------------------------------------
        // DIVIDING BY THE COEFFICIENTS OF THE TENSOR

        // tensorC[0][2] *= 0.5;
        // tensorC[1][2] *= 0.5;
        // tensorC[2][0] *= 0.5;
        // tensorC[2][1] *= 0.5;
        // tensorC[2][2] *= 0.25;

        // for (int i = 0; i < 3; i++)
        // {
        //     for (int j = 0; j < 3; j++)
        //         std::cout << tensorC[i][j] << " ";
        //     std::cout << std::endl;
        // }
        break;
    }
    }

    return tensorC;
}

double Solid2D::getQValue(const int idxLocalNode1, const int idxLocalNode2, bool _PrescribedDamageField)
{
    const PetscInt numNodeDOF = 2;  // Number of DOFs per node considering displacements only
    const PetscInt nDOFPF = 1;      // Number of DOFs per node considering phase field only
    PetscInt numElDOF = numElNodes; // Only one DOF per node when considering only phase field
    PetscScalar l0 = material->getL0();
    PetscScalar Gc = material->getGriffithCriterion();
    PetscScalar lame = material->getLameConstant();
    PetscScalar G = material->getShearModulus();
    PetscScalar kappa = 2.0 / 3.0 * G + lame;
    std::string PFmodel = params->getPFModel();
    double QValue = 0.0;

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
                psiPlus += G * ((gradU[k][l] * gradU[k][l] + gradU[k][l] * gradU[l][k]) * 0.5) * wJac;

        psiPlus += -1.0 / 3.0 * G * divU * divU * wJac;

        if (divU > 0)
            psiPlus += kappa * 0.5 * divU * divU * wJac;

        PetscScalar damageValue = 0.;
        for (PetscInt c = 0; c < numElNodes; c++)
            damageValue += N[c] * elemConnectivity[c]->getDOFs()[2]->getValue(); // dn

        const PetscScalar dCoeff = -2 * (1 - damageValue);

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

        // ======================= SECOND DERIVATIVE WITH RESPECT TO THE FIELD VARIABLE =======================
        QValue += 2 * N[idxLocalNode1] * N[idxLocalNode2] * psiPlus; // Integral 1 is the same for AT1 and AT2

        /*
            AT1 OR AT2 PHASE FIELD MODEL:
            The second integral of the second derivative with respect to the field variable is different for AT1 and AT2 models.
        */

        if (PFmodel == "AT2")
        {
            QValue += Gc * (1.0 / l0) * N[idxLocalNode1] * N[idxLocalNode2] * wJac;

            for (int i = 0; i < 2; i++)
                QValue += Gc * l0 * B_d[i][idxLocalNode1] * B_d[i][idxLocalNode2] * wJac;
        }
        else if (PFmodel == "AT1")
        {
            for (int i = 0; i < 2; i++)
                QValue += Gc * 0.75 * l0 * B_d[i][idxLocalNode1] * B_d[i][idxLocalNode2] * wJac;
        }

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    return QValue;
}

PetscErrorCode Solid2D::getqContribution(Vec &rhs, bool _PrescribedDamageField)
{
    const PetscScalar lame = material->getLameConstant();
    const PetscScalar G = material->getShearModulus();
    const PetscScalar kappa = 2.0 / 3.0 * G + lame;
    PetscScalar l0 = material->getL0();
    PetscScalar Gc = material->getGriffithCriterion();
    std::string PFmodel = params->getPFModel();

    PetscScalar *localRHS = new PetscScalar[numElNodes]();
    PetscInt *idx = new PetscInt[numElNodes]();

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
                psiPlus += G * ((gradU[k][l] * gradU[k][l] + gradU[k][l] * gradU[l][k]) * 0.5) * wJac;

        psiPlus += -1.0 / 3.0 * G * divU * divU * wJac;

        if (divU > 0)
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
            {
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                        secondInt[a] += Gc * (0.75 * l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;

                secondInt[a] += Gc * ((0.375 / l0) * N[a]) * wJac;
            }
        }
        else if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
            {
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                        secondInt[a] += Gc * (l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;

                secondInt[a] += Gc * (1 / l0 * damageValue * N[a]) * wJac;
            }
        }

        for (PetscInt a = 0; a < numElNodes; a++)
            localRHS[a] = firstInt[a] + secondInt[a];

        PetscCall(VecSetValues(rhs, numElNodes, idx, localRHS, ADD_VALUES));

        delete[] N;
        for (int i = 0; i < numElNodes; i++)
            delete[] dN[i];
        delete[] dN;
    }

    delete[] localRHS;
    delete[] idx;

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
        PetscScalar firstInt[numElNodes] = {};

        switch (params->getSplitModel())
        {
        case volDev:
        {
            for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++)
                    psiPlus += G * ((gradU[k][l] * gradU[k][l] + gradU[k][l] * gradU[l][k]) * 0.5) * wJac;

            psiPlus += -1.0 / 3.0 * G * divU * divU * wJac;

            if (divU > 0)
                psiPlus += kappa * 0.5 * divU * divU * wJac;

            break;
        }

        case spectral:
        {
            double strain[3][3] = {};
            for (PetscInt i = 0; i < 2; ++i)
                for (PetscInt j = 0; j < 2; ++j)
                    strain[i][j] = 0.5 * (gradU[i][j] + gradU[j][i]);

            switch (material->getPlaneAnalysis())
            {
            case PLANE_STRESS:
                double _poisson = material->getPoisson();
                strain[2][2] = -(_poisson / (1 - _poisson)) * (strain[0][0] + strain[1][1]);
                break;
            }

            double volStrain = (1.0 / 3.0) * (strain[0][0] + strain[1][1] + strain[2][2]);

            double d11 = strain[0][0] - volStrain;
            double d22 = strain[1][1] - volStrain;
            double d33 = strain[2][2] - volStrain;
            double d12 = strain[0][1];
            double d21 = strain[1][0];
            double d23 = strain[1][2];
            double d32 = strain[2][1];
            double d31 = strain[2][0];
            double d13 = strain[0][2];

            double eeq = sqrt(2.0 / 3.0 * (d11 * d11 + d22 * d22 + d33 * d33 + 2.0 * (d12 * d12 + d13 * d13 + d23 * d23))); // Equivalent strain

            double detdd =
                d11 * (d22 * d33 - d23 * d32) -
                d12 * (d21 * d33 - d23 * d31) +
                d13 * (d21 * d32 - d22 * d31);

            double const tol = 1e-12;
            double etaReal = 0.0, cos3eta = 0.0;
            if (std::abs(eeq) < tol)
            {
                etaReal = 0.0;
            }
            else
            {
                cos3eta = 4.0 * detdd / (eeq * eeq * eeq);
                cos3eta = std::clamp(cos3eta, -1.0, 1.0);
                cos3eta = std::round(cos3eta * 1e10) / 1e10;
                etaReal = (1.0 / 3.0) * std::acos(cos3eta);
            }

            double etaStar[3] = {};
            etaStar[0] = etaReal;
            etaStar[1] = etaReal - 2.0 * M_PI / 3.0;
            etaStar[2] = etaReal + 2.0 * M_PI / 3.0;

            // COMPUTING THE PRINCIPAL VALUES VIA LODE ANGLE
            double principalValues[3] = {};
            principalValues[0] = volStrain + eeq * std::cos(etaStar[0]); // ep1
            principalValues[1] = volStrain + eeq * std::cos(etaStar[1]); // ep2
            principalValues[2] = volStrain + eeq * std::cos(etaStar[2]); // ep3
            // for (int ip = 0; ip < 3; ip++)
            // {
            //     double cosEta = std::cos(etaStar[ip]);
            //     double val = volStrain + eeq * cosEta;
            //     principalValues[ip] = std::round(val * 1e10) / 1e10;
            // }

            const double pTol = 1e-12;
            for (int ip = 0; ip < 3; ip++)
                if (principalValues[ip] > pTol)
                    psiPlus += G * (principalValues[ip] * principalValues[ip]) * wJac;

            if (divU > pTol)
                psiPlus += lame * 0.5 * divU * divU * wJac;
        }
        }

        PetscScalar damageValue = 0.;
        for (PetscInt c = 0; c < numElNodes; c++)
            damageValue += N[c] * elemConnectivity[c]->getDOFs()[2]->getValue(); // dn

        const PetscScalar dCoeff = -2 * (1 - damageValue);

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
            {
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                        secondInt[a] += Gc * (0.75 * l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;

                secondInt[a] += Gc * ((0.375 / l0) * N[a]) * wJac;
            }
        }
        else if (PFmodel == "AT2")
        {
            for (PetscInt a = 0; a < numElNodes; a++)
            {
                for (PetscInt i = 0; i < 3; i++)
                    for (PetscInt j = 0; j < 2; j++)
                        secondInt[a] += Gc * (l0 * B_dT[a][j] * B_d[j][i] * elemConnectivity[i]->getDOFs()[2]->getValue()) * wJac;

                secondInt[a] += Gc * (1 / l0 * damageValue * N[a]) * wJac;
            }
        }

        for (PetscInt a = 0; a < numElNodes; a++)
            localRHS[a] = firstInt[a] + secondInt[a];
        // #pragma omp critical
        //  {
        ierr = VecSetValues(rhs, numElDOF, idx, localRHS, ADD_VALUES);
        //  }
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

    // #pragma omp critical
    //{
    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localQ, ADD_VALUES);
    // }
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