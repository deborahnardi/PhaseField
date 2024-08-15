#pragma once

/*
This class creates the Input Data header file for providing the necessary input data for the mesh generation in Gmsh.
*/

#include <vector>
#include <iostream>

class Inclusion
{
protected:
    int index;
    double a, b, alpha, xc, yc, lc;

public:
    Inclusion();
    Inclusion(int _index, const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc, const double &_lc);
    ~Inclusion();

    int getIndex() const { return index; }
    double getA() const { return a; }
    double getB() const { return b; }
    double getAlpha() const { return alpha; }
    double getXc() const { return xc; }
    double getYc() const { return yc; }
    double getLc() const { return lc; }

    void setIndex(const int &_index) { index = _index; }
    void setA(const double &_a) { a = _a; }
    void setB(const double &_b) { b = _b; }
    void setAlpha(const double &_alpha) { alpha = _alpha; }
    void setXc(const double &_xc) { xc = _xc; }
    void setYc(const double &_yc) { yc = _yc; }
    void setLc(const double &_lc) { lc = _lc; }
};

/*
    meshMinFac: Minimum factor for mesh size
    meshMaxFac: Maximum factor for mesh size
    meshDistFac: Factor for mesh size distribution
    meshMinSize: Minimum mesh size
    meshMaxSize: Maximum mesh size
*/

class MeshFactor : public Inclusion
{
private:
    int index;
    double meshMinFac, meshMaxFac, meshDistFac, meshMinSize, meshMaxSize;

public:
    MeshFactor(const int &_index, const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize);
    MeshFactor();
    ~MeshFactor();

    double getMeshMinFac() const { return meshMinFac; }
    double getMeshMaxFac() const { return meshMaxFac; }
    double getMeshDistFac() const { return meshDistFac; }
    double getMeshMinSize() const { return meshMinSize; }
    double getMeshMaxSize() const { return meshMaxSize; }

    void setMeshMinFac(const double &_meshMinFac) { meshMinFac = _meshMinFac; }
    void setMeshMaxFac(const double &_meshMaxFac) { meshMaxFac = _meshMaxFac; }
    void setMeshDistFac(const double &_meshDistFac) { meshDistFac = _meshDistFac; }
    void setMeshMinSize(const double &_meshMinSize) { meshMinSize = _meshMinSize; }
    void setMeshMaxSize(const double &_meshMaxSize) { meshMaxSize = _meshMaxSize; }
};