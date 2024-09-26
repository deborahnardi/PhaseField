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

PetscErrorCode Solid2D::getContribution(Mat &A)
{
    const int numElDOF = numElNodes * 2;
    double **coords = q->getQuadratureCoordinates();
    double *weights = q->getQuadratureWeights();
    double localStiffnessMatrix[numElDOF][numElDOF] = {}; // {} initializes all elements to 0

    const double G = material->getShearModulus();
    const double lame = material->getLameConstant();

    for (int ih = 0; ih < numHammerPoints; ih++)
    {
        double *xi = coords[ih];
        double weight = weights[ih];

        double *N = sF->evaluateShapeFunction(xi);
        double **dN = sF->getShapeFunctionDerivative(xi);

        double dX_dXi[2][2] = {};

        for (int a = 0; a < 3; a++)
        {
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    dX_dXi[i][j] += dN[a][j] * elemConnectivity[a]->getInitialCoordinates()[i];
                }
            }
        }

        double jac = dX_dXi[0][0] * dX_dXi[1][1] - dX_dXi[0][1] * dX_dXi[1][0];
        double wJac = weight * jac;

        /*
        For a solid element, four loops are needed to compute the element stiffness matrix.
        The first loop is over the number of nodes in the element.
        The second loop is over the number of degrees of freedom per node.
        The third loop is over the number of nodes in the element.
        The fourth loop is over the number of degrees of freedom per node.
        */

        for (int a = 0; a < 3; a++)
        {
            for (int b = 0; b < 3; b++)
            {
                double aux = 0.;

                for (int k = 0; k < 2; k++)
                {
                    aux += dN[a][k] * dN[b][k];
                }

                for (int i = 0; i < 2; i++)
                {
                    localStiffnessMatrix[2 * a + i][2 * b + i] += G * aux * wJac; // Due to Kronnecker delta

                    for (int j = 0; j < 2; j++)
                    {
                        localStiffnessMatrix[2 * a + i][2 * b + j] += (G * dN[a][j] * dN[b][i] + lame * dN[a][i] * dN[b][j]) * wJac;
                    }
                }
            }
        }
    }
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

        /*
        For a solid element, four loops are needed to compute the element stiffness matrix.
        The first loop is over the number of nodes in the element.
        The second loop is over the number of degrees of freedom per node.
        The third loop is over the number of nodes in the element.
        The fourth loop is over the number of degrees of freedom per node.
        */

        for (int a = 0; a < numElNodes; a++)
        {
            for (int b = 0; b < numElNodes; b++)
            {
                double contraction = 0.;
                for (int k = 0; k < 2; k++)
                    contraction += dN[a][k] * dN[b][k];

                for (int i = 0; i < 2; i++)
                {
                    localStiff(2 * a + i, 2 * a + i) += G * contraction * wJac; // Due to Kroenecker delta

                    for (int j = 0; j < 2; j++)
                        localStiff(2 * a + i, 2 * b + j) += G * dN[a][j] * dN[b][i] + lame * dN[a][i] * dN[b][j] * wJac;
                }
            }
        }
    }
}

void Solid2D::assembleGlobalStiffnessMatrix(MatrixXd &GlobalStiff)
{
    int dof1 = getNode(0)->getDOF(0)->getIndex();
    int dof2 = getNode(1)->getDOF(0)->getIndex();
    int dof3 = getNode(2)->getDOF(0)->getIndex();

    std::cout << "Kelem:" << std::endl;
    std::cout << localStiff << std::endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
        {
            GlobalStiff(dof1 + i, dof1 + j) += localStiff(i, j);
            GlobalStiff(dof1 + i, dof2 + j) += localStiff(i, j + 2);
            GlobalStiff(dof1 + i, dof3 + j) += localStiff(i, j + 4);

            GlobalStiff(dof2 + i, dof1 + j) += localStiff(i + 2, j);
            GlobalStiff(dof2 + i, dof2 + j) += localStiff(i + 2, j + 2);
            GlobalStiff(dof2 + i, dof3 + j) += localStiff(i + 2, j + 4);

            GlobalStiff(dof3 + i, dof1 + j) += localStiff(i + 4, j);
            GlobalStiff(dof3 + i, dof2 + j) += localStiff(i + 4, j + 2);
            GlobalStiff(dof3 + i, dof3 + j) += localStiff(i + 4, j + 4);
        }
}