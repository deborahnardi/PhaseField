#pragma once

#include <iostream>

class ShapeFunction
{
protected:
    int numNodes;

public:
    ShapeFunction();
    ~ShapeFunction();

    virtual double *evaluateShapeFunction(double *&xi) {};
    virtual double **getShapeFunctionDerivative(double *&xi) {};
    virtual void getNodalXi(double **&coord) {};
    int getNumNodes() { return numNodes; };
};

class T3ShapeFunction : public ShapeFunction
{
public:
    T3ShapeFunction();
    ~T3ShapeFunction();

    double *evaluateShapeFunction(double *&xi) override;
    double **getShapeFunctionDerivative(double *&xi) override;
    void getNodalXi(double **&coord) override;
};