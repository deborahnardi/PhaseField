#include "../headers/Geometry.h"

Geometry::Geometry() {}
Geometry::Geometry(const std::string _name) : name(_name)
{
    gmsh::initialize();
    gmsh::model::add(name);
}
Geometry::~Geometry() {}

Point *Geometry::addPoint(const std::vector<double> &_coordinates, const double &_lc)
{
    Point *point = new Point(_coordinates, _lc, points.size());
    points.push_back(point);

    gmsh::model::occ::addPoint(_coordinates[0], _coordinates[1], _coordinates[2], _lc, point->getIndex() + 1);
    gmsh::model::occ::synchronize();
    gmsh::model::addPhysicalGroup(0, {point->getIndex() + 1}, -1, point->getEntityName());

    return point;
}

Line *Geometry::addLine(const std::vector<Point *> &_points)
{
    Line *line = new Line(_points, lines.size());
    lines.push_back(line);

    gmsh::model::occ::addLine(line->getPoint(0)->getIndex() + 1, line->getPoint(1)->getIndex() + 1, line->getIndex() + 1);
    gmsh::model::occ::synchronize();
    gmsh::model::addPhysicalGroup(1, {line->getIndex() + 1}, -1, line->getEntityName());

    return line;
}

Circle *Geometry::addCircle(const std::vector<Point *> &_points)
{
    Circle *c = new Circle(_points, lines.size());
    c->setEntityName("c" + std::to_string(c->getIndex() + 1));
    circles.push_back(c);
    return c;
}

LineLoop *Geometry::addLineLoop(const std::vector<Line *> &_lines)
{
    LineLoop *lineLoop = new LineLoop(_lines, wires.size());
    wires.push_back(lineLoop);
    lineLoops.push_back(lineLoop);

    std::vector<int> linesIndexes;
    for (auto l : _lines)
        linesIndexes.push_back(l->getIndex() + 1);

    gmsh::model::occ::addWire(linesIndexes, lineLoop->getIndex() + 1);
    gmsh::model::occ::synchronize();

    return lineLoop;
}

PlaneSurface *Geometry::addPlaneSurface(std::vector<Wire *> _wire)
{
    std::vector<int> wireTags, wireTagsAux;
    for (auto w : _wire)
    {
        wireTags.push_back(w->getIndex());
        wireTagsAux.push_back(w->getIndex() + 1);
    }
    PlaneSurface *planeSurface = new PlaneSurface(wireTags, planeSurfaces.size());
    planeSurfaces.push_back(planeSurface);

    int surfaceTag = gmsh::model::occ::addPlaneSurface(wireTagsAux, -1);
    gmsh::model::occ::synchronize();
    gmsh::model::addPhysicalGroup(2, {surfaceTag}, -1, planeSurface->getEntityName());

    return planeSurface;
}

Ellipse *Geometry::addEllipse(double _a, double _b, const double _alpha, std::vector<double> center, const double _lc)
{
    // Generating coordinates for the ellipse
    double rad = _alpha * M_PI / 180.;                     // Converting to radians
    std::vector<double> xAxis = {cos(rad), sin(rad), 0.0}; // Positive anticlockwise

    Ellipse *el = new Ellipse(center, _a, _b, wires.size(), 0., 2 * M_PI, xAxis);
    wires.push_back(el);
    ellipses.push_back(el);

    gmsh::model::occ::addEllipse(center[0], center[1], center[2], _a, _b, -1, 0., 2. * M_PI, {0., 0., 1.}, xAxis);
    gmsh::model::occ::synchronize();
    int curveTag = ellipses.size() + lines.size();
    gmsh::model::occ::addWire({curveTag}, el->getIndex() + 1);
    gmsh::model::addPhysicalGroup(1, {el->getIndex() + 1}, -1, el->getEntityName());

    return el;
}

Material *Geometry::addMaterial(const double &_E, const double &_nu, const PlaneAnalysis &_plane)
{
    Material *material = new Material(materials.size(), _nu, _E, _plane);
    materials.push_back(material);
    return material;
}

void Geometry::addTransfiniteLine(const std::vector<Line *> &_lines, const int &_divisions, const double &_progression)
{
    for (int i = 0; i < _lines.size(); i++)
    {
        transfiniteLines.push_back(std::make_pair(_lines[i], _divisions));

        gmsh::model::mesh::setTransfiniteCurve(_lines[i]->getIndex() + 1, _divisions + 1, "Progression");
        gmsh::model::occ::synchronize();
    }
}

// MeshFactor *Geometry::addMeshFactor(const double &_meshMinFac, const double &_meshMaxFac, const double &_meshDistFac, const double &_meshMinSize, const double &_meshMaxSize)

// {
//     MeshFactor *meshFac = new MeshFactor(meshFactors.size(), _meshMinFac, _meshMaxFac, _meshDistFac, _meshMinSize, _meshMaxSize);
//     meshFactors.push_back(meshFac);

//     meshMinSizeIncl = _meshMinFac * inclusions[0]->getA() * edgeLength;
//     meshMaxSizeIncl = _meshMaxFac * inclusions[0]->getA() * edgeLength;
//     meshDistMin = inclusions[0]->getA() * edgeLength;
//     meshDistMax = _meshDistFac * inclusions[0]->getA() * edgeLength;
//     meshMinSizeGlobal = _meshMinSize * edgeLength;
//     meshMaxSizeGlobal = _meshMaxSize * edgeLength;

//     return meshFac;
// }

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

// void Geometry::setMeshInclusionProperties()
// {

//     // Global definition for mesh size generation
//     gmsh::option::setNumber("Mesh.MeshSizeMin", meshMinSizeGlobal);      // Defines the minimum mesh size
//     gmsh::option::setNumber("Mesh.MeshSizeMax", meshMaxSizeGlobal);      // Defines the maximum mesh size
//     gmsh::option::setNumber("Mesh.MeshSizeFactor", getMeshSizeFactor()); // Defines the mesh size factor

//     // 0 -> Deactivated; 1 -> Activated
//     gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
//     gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 0);
//     gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", 0);

//     // Refining the region around and inside the inclusions
//     gmsh::model::mesh::field::add("Distance", 1);
//     std::vector<double> doubleEllipseSurfaces(ellipseSurfaces.begin(), ellipseSurfaces.end());
//     gmsh::model::mesh::field::setNumbers(1, "SurfacesList", doubleEllipseSurfaces); // List of surfaces to be refined
//     gmsh::model::mesh::field::setNumber(1, "Sampling", 1000);                       // Number of points to be sampled

//     // We then define a `Threshold' field, which uses the return value of the
//     // `Distance' field 1 in order to define a simple change in element size
//     // depending on the computed distances
//     //
//     // SizeMax -                     /------------------
//     //                              /
//     //                             /
//     //                            /
//     // SizeMin -o----------------/
//     //          |                |    |
//     //        Point         DistMin  DistMax

//     gmsh::model::mesh::field::add("Threshold", 2); // Threshold field allows to refine the mesh in a specific region
//     gmsh::model::mesh::field::setNumber(2, "IField", 1);
//     gmsh::model::mesh::field::setNumber(2, "SizeMin", meshMinSizeIncl);
//     gmsh::model::mesh::field::setNumber(2, "SizeMax", meshMaxSizeIncl);
//     gmsh::model::mesh::field::setNumber(2, "DistMin", meshDistMin);
//     gmsh::model::mesh::field::setNumber(2, "DistMax", meshDistMax);

//     gmsh::model::mesh::field::add("Min", 3);
//     gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});

//     gmsh::model::mesh::field::setAsBackgroundMesh(3);
// }

void Geometry::InitializeGmshAPI(const bool &showInterface)
{
    gmsh::option::setNumber("Mesh.Algorithm", algorithm);
    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(dim);
    // gmsh::write(name + ".msh");

    writeMeshInfo();
    // (elemDim == 2) ? writeMeshInfo2D() : writeMeshInfo1D();

    gmsh::write("modelo_oc.geo_unrolled");
    if (showInterface)
        gmsh::fltk::run();

    gmsh::clear();
    gmsh::finalize();
}