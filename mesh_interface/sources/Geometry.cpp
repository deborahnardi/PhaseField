#include "../headers/Geometry.h"

Geometry::Geometry() {}
Geometry::Geometry(const std::string _name)
    : name(_name) {}
Geometry::~Geometry() {}

Point *Geometry::addPoint(const std::vector<double> &_coordinates, const double &_lc)
{
    Point *point = new Point(_coordinates, _lc, points.size());
    points.push_back(point);
    return point;
}

Line *Geometry::addLine(const std::vector<Point *> &_points)
{
    Line *line = new Line(_points, lines.size());
    lines.push_back(line);
    return line;
}

LineLoop *Geometry::addLineLoop(const std::vector<Line *> &_lines)
{
    LineLoop *lineLoop = new LineLoop(_lines, lineLoops.size());
    lineLoops.push_back(lineLoop);
    return lineLoop;
}

PlaneSurface *Geometry::addPlaneSurface(LineLoop *_lineLoop)
{
    PlaneSurface *planeSurface = new PlaneSurface(_lineLoop, planeSurfaces.size());
    planeSurfaces.push_back(planeSurface);
    return planeSurface;
}

Inclusion *Geometry::addInclusion(const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc, const double &_lc)
{
    Inclusion *incl = new Inclusion(inclusions.size(), _a, _b, _alpha, _xc, _yc, _lc);
    inclusions.push_back(incl);
    return incl;
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

void Geometry::generateInclusions()
{
    std::cout << "Generating inclusions coordinates..." << std::endl;

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

            gmsh::model::occ::addPoint(point->getX(), point->getY(), point->getZ(), point->getLC(), point->getIndex() + 1);

            tags[i] = point->getIndex() + 1;
        }

        for (int i = 0; i < 4; i++) // Each ellipse contains 4 arc lines
        {
            int centerTag = tags[4];
            int majorTag = tags[2];
            int startTag = tags[i];
            int endTag = tags[(i + 1) % 4]; // Cyclic form of the array

            std::cout << "gmsh::model::occ::addEllipseArc(" << startTag << ", " << centerTag << ", " << majorTag << ", " << endTag << ", " << lines.size() + ellipseArcs.size() + 1 << ")" << std::endl;
            int elp = gmsh::model::occ::addEllipseArc(startTag, centerTag, majorTag, endTag, lines.size() + ellipseArcs.size() + 1);
            ellipseArcs.push_back(elp);
        }

        std::cout << "gmsh::model::occ::addCurveLoop({" << ellipseArcs[4 * inclusion->getIndex()] << ", " << ellipseArcs[4 * inclusion->getIndex() + 1] << ", " << ellipseArcs[4 * inclusion->getIndex() + 2] << ", " << ellipseArcs[4 * inclusion->getIndex() + 3] << "}, " << lineLoops.size() + inclusion->getIndex() + 1 << ")" << std::endl;

        int cInc = gmsh::model::occ::addCurveLoop({ellipseArcs[4 * inclusion->getIndex()], ellipseArcs[4 * inclusion->getIndex() + 1], ellipseArcs[4 * inclusion->getIndex() + 2], ellipseArcs[4 * inclusion->getIndex() + 3]}, lineLoops.size() + inclusion->getIndex() + 1);
        ellipseCurves.push_back(cInc);

        int pInc = gmsh::model::occ::addPlaneSurface({cInc}, planeSurfaces.size() + inclusion->getIndex() + 1);
        ellipseSurfaces.push_back(pInc);

        std::cout << "gmsh::model::occ::addPlaneSurface({" << cInc << "}, " << planeSurfaces.size() + inclusion->getIndex() + 1 << ")" << std::endl;
        std::cout << "--------------------------------" << std::endl;

        gmsh::model::occ::synchronize();

        gmsh::model::addPhysicalGroup(2, {pInc}, planeSurfaces.size() + inclusion->getIndex() + 1, "Inclusions");

        gmsh::model::occ::synchronize();
    }
};

void Geometry::InitializeGmshAPI()
{
    gmsh::initialize();
    gmsh::model::add(name);

    // Generating the geometry
    // Adding the HOST geometry -> HOST is the main domain

    for (auto point : points)
    {
        gmsh::model::occ::addPoint(point->getX(), point->getY(), point->getZ(), point->getLC(), point->getIndex() + 1);
        int aux = point->getIndex() + 1;
        std::cout << "gmsh::model::occ::addPoint(" << point->getX() << ", " << point->getY() << ", " << point->getZ() << ", " << point->getLC() << ", " << aux << ")" << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;

    for (auto line : lines)
    {
        gmsh::model::occ::addLine(line->getPoint(0)->getIndex() + 1, line->getPoint(1)->getIndex() + 1, line->getIndex() + 1);
        std::cout << "gmsh::model::occ::addLine(" << line->getPoint(0)->getIndex() + 1 << ", " << line->getPoint(1)->getIndex() + 1 << ", " << line->getIndex() + 1 << ")" << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    // Valid only for lines with 2 points

    for (auto lineLoop : lineLoops)
    {
        for (int i = 0; i < lineLoop->getNumLines(); i++)
            linesIndexes.push_back(lineLoop->getLine(i)->getIndex() + 1);

        gmsh::model::occ::addCurveLoop(linesIndexes, lineLoop->getIndex() + 1);
        std::cout << "gmsh::model::occ::addCurveLoop({";
        for (int i = 0; i < lineLoop->getNumLines(); i++)
        {
            std::cout << lineLoop->getLine(i)->getIndex() + 1;
            if (i < lineLoop->getNumLines() - 1)
                std::cout << ", ";
        }
        std::cout << "}, " << lineLoop->getIndex() + 1 << ")" << std::endl;

        gmsh::model::addPhysicalGroup(1, linesIndexes, lineLoop->getIndex() + 1, "External Boundaries");
        gmsh::model::occ::synchronize();
    }

    for (auto planeSurface : planeSurfaces)
    {
        gmsh::model::occ::addPlaneSurface({planeSurface->getLineLoop()->getIndex() + 1}, planeSurface->getIndex() + 1);
        std::cout << "gmsh::model::occ::addPlaneSurface({" << planeSurface->getLineLoop()->getIndex() + 1 << "}, " << planeSurface->getIndex() + 1 << ")" << std::endl;
        std::cout << "--------------------------------" << std::endl;

        gmsh::model::addPhysicalGroup(2, {planeSurface->getIndex() + 1}, planeSurface->getIndex() + 1, "Host");
        gmsh::model::occ::synchronize();
    }

    generateInclusions();

    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", meshMinSizeGlobal);
    gmsh::option::setNumber("Mesh.MeshSizeMax", meshMaxSizeGlobal);
    gmsh::option::setNumber("Mesh.MeshSizeFactor", 1.0);

    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
    gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

    // Refining the region around and inside the inclusions
    gmsh::model::mesh::field::add("Distance", 1);
    std::vector<double> doubleEllipseSurfaces(ellipseSurfaces.begin(), ellipseSurfaces.end());
    gmsh::model::mesh::field::setNumbers(1, "SurfacesList", doubleEllipseSurfaces);
    gmsh::model::mesh::field::setNumber(1, "Sampling", 1000);

    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "IField", 1);
    gmsh::model::mesh::field::setNumber(2, "SizeMin", meshMinSizeIncl);
    gmsh::model::mesh::field::setNumber(2, "SizeMax", meshMaxSizeIncl);
    gmsh::model::mesh::field::setNumber(2, "DistMin", meshDistMin);
    gmsh::model::mesh::field::setNumber(2, "DistMax", meshDistMax);

    gmsh::model::mesh::field::add("Min", 3);
    gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});

    gmsh::model::mesh::field::setAsBackgroundMesh(3);

    gmsh::option::setNumber("Mesh.Algorithm", algorithm);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(dim);
    gmsh::write(name + ".msh");
    gmsh::fltk::run();
    gmsh::finalize();
}