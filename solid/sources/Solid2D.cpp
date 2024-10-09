#include "../headers/Element.h"

Solid2D::Solid2D() {}
Solid2D::Solid2D(const int &_index, const int &_elemDimension, const std::vector<Node *> &_elemConnectivity, Material *_material, const int &_physicalEntity)
    : Element(_index, _elemDimension, _elemConnectivity, _material, _physicalEntity)
{
    for (auto n : _elemConnectivity)
    {
        n->setIsDiscritized();
        n->addDOF(new DOF(X, 0.));
        n->addDOF(new DOF(Y, 0.));
    }

    sF = new T3ShapeFunction();
    q = new TriangularQuadrature(numHammerPoints);
}
Solid2D::~Solid2D() {}

/*----------------------------------------------------------------------------------
                Assembling and solving problem with PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode Solid2D::getContribution(Mat &A, Vec &rhs)
{
    PetscInt numElDOF = numElNodes * 2;
    PetscReal *localStiffnessMatrix = new PetscScalar[numElDOF * numElDOF]();
    PetscReal *localRHS = new PetscScalar[numElDOF]();
    PetscInt *idx = new PetscInt[numElDOF]();

    double **coords = q->getQuadratureCoordinates();
    double *weights = q->getQuadratureWeights();
    PetscReal kroen[2][2] = {{1., 0.}, {0., 1.}};

    const double G = material->getShearModulus();
    const double lame = material->getLameConstant();

    PetscInt count = 0;
    for (auto node : elemConnectivity)
        for (auto dof : node->getDOFs())
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
                        localStiffnessMatrix[pos] += (G * contraction * wJac * kroen[i][j] + G * dN_dX[a][j] * dN_dX[b][i] * wJac + lame * dN_dX[a][i] * dN_dX[b][j] * wJac);
                    }
                }
        delete[] N;
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
                    DOF *dof = elemConnectivity[b]->getDOFs()[j];
                    double value = dof->getDirichletValue();
                    if (value != 0)
                    {
                        double fi = -localStiffnessMatrix[numElDOF * (2 * a + i) + 2 * b + j] * value;
                        ierr = VecSetValues(rhs, 1, &idx[2 * a + i], &fi, ADD_VALUES);
                        CHKERRQ(ierr);
                    }
                }

    ierr = MatSetValues(A, numElDOF, idx, numElDOF, idx, localStiffnessMatrix, ADD_VALUES);
    CHKERRQ(ierr);

    delete[] coords;
    delete[] weights;
    delete[] idx;
    delete[] localStiffnessMatrix;

    return ierr;
}

void Solid2D::Test(PetscScalar &integral)
{
    double **coords = q->getQuadratureCoordinates();
    double *weights = q->getQuadratureWeights();

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
    double **coords = q->getQuadratureCoordinates();
    double *weights = q->getQuadratureWeights();
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