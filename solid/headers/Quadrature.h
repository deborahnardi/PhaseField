#pragma once

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

class TriangularQuadrature : public Quadrature
{
public:
    TriangularQuadrature();
    TriangularQuadrature(int numIntPoints);
    ~TriangularQuadrature();

    void getQuadrature() override;
};
