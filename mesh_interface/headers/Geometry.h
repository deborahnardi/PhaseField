#pragma once

#include "hdf5.h"
#include <string>
#include <iostream>
#include <fstream>

#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <set>

#include <gmsh.h>

#include "Inclusion.h"
#include "Point.h"
#include "Line.h"
#include "LineLoop.h"
#include "Surface.h"
#include "PlaneSurface.h"
#include "../../enumclass.hpp"

class Geometry
{
private:
    int dim;
    int elpPoints = 5;
    double edgeLength;
    double xc, yc, xo, yo, xm, ym, xa, ya, xb, yb;
    double meshMinSizeIncl, meshMaxSizeIncl, meshDistMin, meshDistMax, meshMinSizeGlobal, meshMaxSizeGlobal;
    double **ellipseCoordinates;
    std::string name;
    std::vector<int> linesIndexes;
    std::vector<int> tags, ellipseArcs, ellipseCurves, ellipseSurfaces;
    std::vector<std::pair<int, int>> entities;
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    std::vector<Point *> points;
    std::vector<Line *> lines;
    std::vector<LineLoop *> lineLoops;
    std::vector<PlaneSurface *> planeSurfaces;
    std::vector<Inclusion *> inclusions;
    std::vector<MeshFactor *> meshFactors;
    MeshAlgorithm algorithm;

public:
    Geometry();
    Geometry(const std::string _name);
    ~Geometry();

    double getEdgeLength() const { return edgeLength; }
    MeshAlgorithm getAlgorithm() const { return algorithm; }
    int getDimention() const { return dim; }

    void setDimention(const int &_dim) { dim = _dim; }
    void setEdgeLength(const double &_edgeLength) { edgeLength = _edgeLength; }
    void setAlgorithm(const MeshAlgorithm &_algorithm) { algorithm = _algorithm; }

    Point *addPoint(const std::vector<double> &_coordinates, const double &_lc);
    Line *addLine(const std::vector<Point *> &_points);
    LineLoop *addLineLoop(const std::vector<Line *> &_lines);
    PlaneSurface *addPlaneSurface(LineLoop *_lineLoop);
    Inclusion *addInclusion(const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc, const double &_lc);
    MeshFactor *addMeshFactor(const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize);

    void InitializeGmshAPI(const bool &showInterface = false);
    void generateInclusions();
    void getMeshInfo();
    void writeMeshInfo();
};