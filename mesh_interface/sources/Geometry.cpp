#include "../headers/Geometry.h"

Geometry::Geometry() {}
Geometry::Geometry(const std::string _name, const bool &_hasInclusions)
    : name(_name), hasInclusions(_hasInclusions) {}
Geometry::~Geometry() {}

Point *Geometry::addPoint(const std::vector<double> &_coordinates, const double &_lc)
{
    Point *point = new Point(_coordinates, _lc, points.size());
    point->setEntityName("p" + std::to_string(point->getIndex() + 1));
    points.push_back(point);

    return point;
}

Line *Geometry::addLine(const std::vector<Point *> &_points)
{
    Line *line = new Line(_points, lines.size());
    line->setEntityName("l" + std::to_string(line->getIndex() + 1));
    lines.push_back(line);
    return line;
}

LineLoop *Geometry::addLineLoop(const std::vector<Line *> &_lines)
{
    LineLoop *lineLoop = new LineLoop(_lines, lineLoops.size());
    lineLoop->setEntityName("ll" + std::to_string(lineLoop->getIndex() + 1));
    lineLoops.push_back(lineLoop);
    return lineLoop;
}

PlaneSurface *Geometry::addPlaneSurface(LineLoop *_lineLoop)
{
    PlaneSurface *planeSurface = new PlaneSurface(_lineLoop, planeSurfaces.size());
    planeSurface->setEntityName("s" + std::to_string(planeSurface->getIndex() + 1));
    planeSurfaces.push_back(planeSurface);
    return planeSurface;
}

Inclusion *Geometry::addInclusion(const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc, const double &_lc)
{
    Inclusion *incl = new Inclusion(inclusions.size(), _a, _b, _alpha, _xc, _yc, _lc);
    inclusions.push_back(incl);
    return incl;
}

Material *Geometry::addMaterial(const double &_E, const double &_nu, const ApplyMaterial &_whereToApply, const PlaneAnalysis &_plane, const std::string &_matType)
{
    Material *material = new Material(materials.size(), _nu, _E, _plane, _matType, _whereToApply);
    materials.push_back(material);
    return material;
} // remove where to apply

void Geometry::addTransfiniteLine(const std::vector<Line *> &_lines, const int &_divisions, const double &_progression)
{
    for (int i = 0; i < _lines.size(); i++)
        transfiniteLines.push_back(std::make_pair(_lines[i], _divisions));
}

MeshFactor *Geometry::addMeshFactor(const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize)

{
    MeshFactor *meshFac = new MeshFactor(meshFactors.size(), _meshMinFac, _meshMaxFac, _meshDistFac, _meshMinSize, _meshMaxSize);
    meshFactors.push_back(meshFac);

    meshMinSizeIncl = _meshMinFac * inclusions[0]->getA() * edgeLength;
    meshMaxSizeIncl = _meshMaxFac * inclusions[0]->getA() * edgeLength;
    meshDistMin = inclusions[0]->getA() * edgeLength;
    meshDistMax = _meshDistFac * inclusions[0]->getA() * edgeLength;
    meshMinSizeGlobal = _meshMinSize * edgeLength;
    meshMaxSizeGlobal = _meshMaxSize * edgeLength;

    return meshFac;
}

BoundaryCondition *Geometry::addBoundaryCondition(Point *point, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues)
{
    BoundaryCondition *bCondition = new BoundaryCondition(boundaryConditions.size(), point->getEntityName(), _dofValues, _bType);
    boundaryConditions.push_back(bCondition);
    return bCondition;
}

BoundaryCondition *Geometry::addBoundaryCondition(Line *line, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues)
{
    BoundaryCondition *bCondition = new BoundaryCondition(boundaryConditions.size(), line->getEntityName(), _dofValues, _bType);
    boundaryConditions.push_back(bCondition);
    return bCondition;
}

void Geometry::generateInclusions()
{
    double xc, yc, xo, yo, xm, ym, xa, ya, xb, yb;
    double **ellipseCoordinates;
    std::vector<int> tags;
    int elpPoints = 5; // Number of points in the ellipse
    tags.resize(elpPoints);

    ellipseCoordinates = new double *[elpPoints];

    for (int i = 0; i < elpPoints; i++)
        ellipseCoordinates[i] = new double[2];

    for (auto inclusion : inclusions)
    {
        // Generating coordinates for the ellipse
        double aaux = inclusion->getA() * edgeLength;
        double baux = inclusion->getB() * aaux;
        double rad = inclusion->getAlpha() * M_PI / 180.; // Converting to radians

        xc = inclusion->getXc() * edgeLength;
        yc = inclusion->getYc() * edgeLength;

        xo = xc - aaux * sin(rad);
        yo = yc + aaux * cos(rad);

        xm = xc + aaux * sin(rad);
        ym = yc - aaux * cos(rad);

        xa = xc + baux * cos(rad);
        ya = yc + baux * sin(rad);

        xb = xc - baux * cos(rad);
        yb = yc - baux * sin(rad);

        ellipseCoordinates[0][0] = xo;
        ellipseCoordinates[0][1] = yo;

        ellipseCoordinates[1][0] = xa;
        ellipseCoordinates[1][1] = ya;

        ellipseCoordinates[2][0] = xm;
        ellipseCoordinates[2][1] = ym;

        ellipseCoordinates[3][0] = xb;
        ellipseCoordinates[3][1] = yb;

        ellipseCoordinates[4][0] = xc;
        ellipseCoordinates[4][1] = yc;

        for (int i = 0; i < elpPoints; i++) // Each ellipse contains 5 points
        {
            Point *point = new Point({ellipseCoordinates[i][0], ellipseCoordinates[i][1], 0.0}, inclusion->getLc(), points.size());
            points.push_back(point);

            gmsh::model::occ::addPoint(point->getX(), point->getY(), point->getZ());

            tags[i] = point->getIndex() + 1;
        }

        for (int i = 0; i < 4; i++) // Each ellipse contains 4 arc lines
        {
            int centerTag = tags[4];
            int majorTag = tags[2];
            int startTag = tags[i];
            int endTag = tags[(i + 1) % 4]; // Cyclic form of the array

            int elp = gmsh::model::occ::addEllipseArc(startTag, centerTag, majorTag, endTag);
            ellipseArcs.push_back(elp);
            lines.push_back(new Line({points[startTag - 1], points[endTag - 1]}, lines.size()));
            lines[lines.size() - 1]->setEntityName("Inclusion_" + std::to_string(lines[lines.size() - 1]->getIndex() + 1));
        }

        int cInc = gmsh::model::occ::addCurveLoop({ellipseArcs[4 * inclusion->getIndex()], ellipseArcs[4 * inclusion->getIndex() + 1], ellipseArcs[4 * inclusion->getIndex() + 2], ellipseArcs[4 * inclusion->getIndex() + 3]}, lines.size() + ellipseCurves.size() + 1);
        lineLoops.push_back(new LineLoop({lines[lines.size() - 4], lines[lines.size() - 3], lines[lines.size() - 2], lines[lines.size() - 1]}, lineLoops.size()));
        ellipseCurves.push_back(lineLoops.size());

        int pInc = gmsh::model::occ::addPlaneSurface({cInc});
        ellipseSurfaces.push_back(pInc);

        gmsh::model::occ::synchronize();
    }

    delete[] ellipseCoordinates;
};

void Geometry::setMeshInclusionProperties()
{
    generateInclusions();
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize();
    gmsh::model::removePhysicalGroups();

    for (auto line : lines)
    {
        if (line->getEntityName().find("Boundary") != std::string::npos)
        {
            gmsh::model::addPhysicalGroup(1, {line->getIndex() + 1}, line->getIndex() + 1, "Boundary_" + std::to_string(line->getIndex() + 1));
            gmsh::model::occ::synchronize();
        }
    }

    for (auto eC : ellipseCurves)
    {
        // get 4 lines of the ellipse
        gmsh::model::addPhysicalGroup(1, {lines[4 * eC - 4]->getIndex() + 1, lines[4 * eC - 3]->getIndex() + 1, lines[4 * eC - 2]->getIndex() + 1, lines[4 * eC - 1]->getIndex() + 1}, -1, "Inclusion_" + std::to_string(eC - 1));
        gmsh::model::occ::synchronize();
    }

    // for (int i = 0; i < ellipseSurfaces.size() + 1; i++)
    // {
    //     gmsh::model::addPhysicalGroup(2, {ellipseSurfaces[i]}, -1, "DomainInc_" + std::to_string(i + 1));
    //     gmsh::model::occ::synchronize();

    //     // if (i == ellipseSurfaces.size() - 1)
    //     // {
    //     //     gmsh::model::addPhysicalGroup(2, {ellipseSurfaces[i] + 1}, -1, "Host"); // Contain all the boundaries
    //     //     gmsh::model::occ::synchronize();
    //     // }
    // }

    // Global definition for mesh size generation
    gmsh::option::setNumber("Mesh.MeshSizeMin", meshMinSizeGlobal);      // Defines the minimum mesh size
    gmsh::option::setNumber("Mesh.MeshSizeMax", meshMaxSizeGlobal);      // Defines the maximum mesh size
    gmsh::option::setNumber("Mesh.MeshSizeFactor", getMeshSizeFactor()); // Defines the mesh size factor

    // 0 -> Deactivated; 1 -> Activated
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    // Refining the region around and inside the inclusions
    gmsh::model::mesh::field::add("Distance", 1);
    std::vector<double> doubleEllipseSurfaces(ellipseSurfaces.begin(), ellipseSurfaces.end());
    gmsh::model::mesh::field::setNumbers(1, "SurfacesList", doubleEllipseSurfaces); // List of surfaces to be refined
    gmsh::model::mesh::field::setNumber(1, "Sampling", 1000);                       // Number of points to be sampled

    // We then define a `Threshold' field, which uses the return value of the
    // `Distance' field 1 in order to define a simple change in element size
    // depending on the computed distances
    //
    // SizeMax -                     /------------------
    //                              /
    //                             /
    //                            /
    // SizeMin -o----------------/
    //          |                |    |
    //        Point         DistMin  DistMax

    gmsh::model::mesh::field::add("Threshold", 2); // Threshold field allows to refine the mesh in a specific region
    gmsh::model::mesh::field::setNumber(2, "IField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", meshMinSizeIncl);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", meshMaxSizeIncl);
    gmsh::model::mesh::field::setNumber(2, "DistMin", meshDistMin);
    gmsh::model::mesh::field::setNumber(2, "DistMax", meshDistMax);

    gmsh::model::mesh::field::add("Min", 3);
    gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});

    gmsh::model::mesh::field::setAsBackgroundMesh(3);
}

void Geometry::InitializeGmshAPI(const bool &showInterface)
{
    gmsh::initialize();
    gmsh::model::add(name);

    for (auto point : points)
    {
        gmsh::model::occ::addPoint(point->getX(), point->getY(), point->getZ(), point->getLC(), point->getIndex() + 1);
        gmsh::model::occ::synchronize();
        gmsh::model::addPhysicalGroup(0, {point->getIndex() + 1}, -1, "p" + std::to_string(point->getIndex() + 1));
    }

    for (auto line : lines)
    {
        gmsh::model::occ::addLine(line->getPoint(0)->getIndex() + 1, line->getPoint(1)->getIndex() + 1, line->getIndex() + 1);
        gmsh::model::occ::synchronize();
        gmsh::model::addPhysicalGroup(1, {line->getIndex() + 1}, -1, "l" + std::to_string(line->getIndex() + 1));
    }

    for (auto lineLoop : lineLoops)
    {
        linesIndexes.clear();
        for (int i = 0; i < lineLoop->getNumLines(); i++)
            linesIndexes.push_back(lineLoop->getLine(i)->getIndex() + 1);

        gmsh::model::occ::addCurveLoop(linesIndexes, lineLoop->getIndex() + 1);
        gmsh::model::occ::synchronize();
        lineLoop->setEntityName("ll" + std::to_string(lineLoop->getIndex() + 1));
        gmsh::model::occ::synchronize();
    }

    for (auto planeSurface : planeSurfaces)
    {
        int surfaceTag = gmsh::model::occ::addPlaneSurface({planeSurface->getLineLoop()->getIndex() + 1}, -1);
        gmsh::model::occ::synchronize();
        planeSurface->setEntityName("s" + std::to_string(planeSurface->getIndex() + 1));
        gmsh::model::addPhysicalGroup(2, {surfaceTag}, -1, "s" + std::to_string(planeSurface->getIndex() + 1));
    }

    for (auto l : transfiniteLines)
    {
        gmsh::model::mesh::setTransfiniteCurve(l.first->getIndex() + 1, l.second + 1, "Progression");
        gmsh::model::occ::synchronize();
    }

    gmsh::option::setNumber("Mesh.Algorithm", algorithm);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(dim);
    gmsh::write(name + ".msh");

    writeMeshInfo();

    // if (elemDim == 2)
    //     writeMeshInfo2D();
    // else
    //     writeMeshInfo1D();
    gmsh::write("modelo_oc.geo_unrolled");
    if (showInterface)
        gmsh::fltk::run();

    gmsh::clear();
    gmsh::finalize();
}