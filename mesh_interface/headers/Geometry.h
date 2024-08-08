#pragma once

#include "hdf5.h"
#include <string>
#include <iostream>
#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>

#include "Inclusion.h"
#include "../../enumclass.hpp"

class Geometry
{
private:
    double edgeLength;
    double meshMinSizeIncl, meshMaxSizeIncl, meshDistMin, meshDistMax, meshMinSizeGlobal, meshMaxSizeGlobal;
    std::string name;
    std::vector<Inclusion *> inclusions;
    std::vector<MeshFactor *> meshFactors;
    MeshAlgorithm algorithm;

public:
    Geometry();
    Geometry(const std::string _name);
    ~Geometry();

    Inclusion *addInclusion(const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc);
    MeshFactor *addMeshFactor(const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize);

    int getEdgeLength() const { return edgeLength; }
    MeshAlgorithm getAlgorithm() const { return algorithm; }
    void setEdgeLength(const int &_edgeLength) { edgeLength = _edgeLength; }
    void setAlgorithm(const MeshAlgorithm &_algorithm) { algorithm = _algorithm; }
};