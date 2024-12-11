#pragma once

#include "hdf5.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <petscsnes.h>
#include <petscksp.h>
#include <petscdraw.h>
#include <petscmat.h>
#include <set>

#include <gmsh.h>

#include "Point.h"
#include "Line.h"
#include "Circle.h"
#include "Wire.h"
#include "Surface.h"
#include "PlaneSurface.h"
#include "BoundaryCondition.h"
#include "Material.h"

#include "../../enumclass.hpp"

class Geometry
{
private:
    int dim, elemDim = 2;
    std::string name;
    std::vector<Point *> points;
    std::vector<Line *> lines;
    std::vector<Circle *> circles;
    std::vector<Wire *> wires;
    std::vector<Ellipse *> ellipses;
    std::vector<LineLoop *> lineLoops;
    std::vector<PlaneSurface *> planeSurfaces;
    std::vector<Material *> materials;
    std::vector<BoundaryCondition *> boundaryConditions;
    MeshAlgorithm algorithm;

public:
    Geometry();
    Geometry(const std::string _name);
    ~Geometry();

    MeshAlgorithm getAlgorithm() const { return algorithm; }
    int getDimension() const { return dim; }
    int getElemDimension() const { return elemDim; }

    void setDimension(const int &_dim) { dim = _dim; }
    void setElemDimension(const int &_elemDim) { elemDim = _elemDim; }
    void setAlgorithm(const MeshAlgorithm &_algorithm) { algorithm = _algorithm; }
    void setMeshSizeFactors(const double &_minSize, const double &_maxSize, const double &_factor) {};

    Point *addPoint(const std::vector<double> &_coordinates, const double &_lc = 0.);
    Line *addLine(const std::vector<Point *> &_points);
    Line *addSpline(std::vector<Point *> points);
    Circle *addCircle(const std::vector<Point *> &_points);
    LineLoop *addLineLoop(const std::vector<Line *> &_lines);
    PlaneSurface *addPlaneSurface(std::vector<Wire *> _wire);
    Ellipse *addEllipse(double _a, double _b, const double _alpha, std::vector<double> center);
    BoundaryCondition *addBoundaryCondition(Point *point, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues);
    BoundaryCondition *addBoundaryCondition(Line *line, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues);
    Material *addMaterial(const double &_E, const double &_nu, const PlaneAnalysis &_plane = PLANE_STRESS);
    void addTransfiniteLine(const std::vector<Line *> &_lines, const int &_divisions, const double &_progression = 1.0);

    void GenerateMeshAPI(const bool &showInterface = false);
    void setSurfaceRefinement(std::vector<Wire *> _elipseWires, double _meshMinSizeIncl, double _meshMaxSizeIncl, double _meshDistMin, double _meshDistMax);
    void setGlobalMeshSize(double meshMinSizeGlobal, double meshMaxSizeGlobal, double meshSizeFactorGlobal);
    // void generateInclusions();
    void writeMeshInfo();

    void setRefiningFieldCurves(std::vector<Line *> lines, const int &tag);
    void setThresholdRefinement(const double &_meshMinSize, const double &_meshMaxSize, const double &_meshDistMin, const double &_meshDistMax, const int &tagInField, const int &tag);
    void setBoxRefinement(const double &_meshMinSize, const double &_meshMaxSize, const double &xmin, const double &xmax, const double &ymin, const double &ymax, const double &thickness_, const int &tag);
    void setBackgroundMesh(std::vector<double> FieldsList, const int &tag);
};