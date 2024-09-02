#include "../headers/ShapeFunction.h"

ShapeFunction::ShapeFunction() {};
ShapeFunction::~ShapeFunction() {};

T3ShapeFunction::T3ShapeFunction() { numNodes = 3; };
T3ShapeFunction::~T3ShapeFunction() {};

double *T3ShapeFunction::evaluateShapeFunction(double *&xi)
{
    /*
    For T3 element, there are 3 nodes and 3 shape functions;
    xi is a 2D array of size 2x1, containing the natural (parametric) coordinates of the point at which the shape functions are to be evaluated -> they vary from 0 to 1;
    */
    double *N = new double[numNodes];
    double xi1 = xi[0];
    double xi2 = xi[1];
    double xi3 = 1 - xi1 - xi2;

    N[0] = xi1;
    N[1] = xi2;
    N[2] = xi3;

    return N;
}

double **T3ShapeFunction::getShapeFunctionDerivative(double *&xi)
{
    /*
        dN -> [array address 1, array address 2, array address 3]
        dN[0] -> [value 1, value 2] -> [dN1/dxi1, dN1/dxi2]
        dN[1] -> [value 3, value 4] -> [dN2/dxi1, dN2/dxi2]
        dN[2] -> [value 5, value 6] -> [dN3/dxi1, dN3/dxi2]
    */
    double **dN = new double *[numNodes];
    for (int i = 0; i < numNodes; i++)
        dN[i] = new double[2];

    dN[0][0] = -1.;
    dN[1][0] = 1.;
    dN[2][0] = 0.;

    dN[0][1] = -1.;
    dN[1][1] = 0.;
    dN[2][1] = 1.;

    return dN;
}

void T3ShapeFunction::getNodalXi(double **&coord)
{
    coord = new double *[numNodes];
    for (int i = 0; i < numNodes; i++)
        coord[i] = new double[2];

    coord[0][0] = 0.;
    coord[0][1] = 0.;

    coord[1][0] = 1.;
    coord[1][1] = 0.;

    coord[2][0] = 0.;
    coord[2][1] = 1.;
}