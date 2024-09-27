#pragma once

#include <math.h>

class Quadrature
{
protected:
    int numIntPoints;
    double **coordinates;
    double *weights;

public:
    Quadrature();
    Quadrature(int _numIntPoints);
    ~Quadrature();

    virtual void getQuadrature() {};
    double **getQuadratureCoordinates() { return coordinates; }
    double *getQuadratureWeights() { return weights; }
};

class LineQuadrature : public Quadrature
{
public:
    LineQuadrature();
    LineQuadrature(int numIntPoints);
    ~LineQuadrature();

    void getQuadrature() override;
};

class TriangularQuadrature : public Quadrature
{
public:
    TriangularQuadrature();
    TriangularQuadrature(int numIntPoints);
    ~TriangularQuadrature();

    void getQuadrature() override;
};
