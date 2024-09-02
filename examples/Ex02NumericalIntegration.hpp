int intPoints = 3;
ShapeFunction *sF = new T3ShapeFunction();
Quadrature *q = new TriangularQuadrature(intPoints);

double **coords = q->getQuadratureCoordinates();
double *weights = q->getQuadratureWeights();

std::vector<std::vector<double>> X(3, std::vector<double>(2, 0.));
X[0] = {5., -1.};
X[1] = {5., 5.};
X[2] = {10., 0.};

std::vector<double> x(2, 0.);
MatDouble dX = MatrixXd::Zero(2, 2);

double integral = 0.;
double jac = 0.;

for (int ip = 0; ip < intPoints; ip++)
{
    x = {0., 0.};
    dX.setZero();

    double *N = sF->evaluateShapeFunction(coords[ip]);
    double **dN = sF->getShapeFunctionDerivative(coords[ip]);

    for (int a = 0; a < 3; a++)
    {
        for (int i = 0; i < 2; i++)
        {
            x[i] += N[a] * X[a][i];

            for (int j = 0; j < 2; j++)
            {
                dX(i, j) += dN[a][j] * X[a][i];
            }
        }
    }
    jac = std::abs(dX.determinant());
    integral += (3. * x[0] + 7. * x[1] - 2.) * weights[ip] * jac;
}

std::cout << "Integral: " << integral << std::endl;