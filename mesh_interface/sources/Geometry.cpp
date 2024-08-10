#include "../headers/Geometry.h"

Geometry::Geometry() {}
Geometry::Geometry(const std::string _name)
    : name(_name) {}
Geometry::~Geometry() {}

Point *Geometry::addPoint(const std::vector<double> &_coordinates, const double &_lc)
{
    Point *point = new Point(_coordinates, _lc, points.size());
    std::cout << "Lc: " << point->getLC() << std::endl;
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

PlaneSurface *Geometry::addPlaneSurface(LineLoop * _lineLoop)
{
    PlaneSurface *planeSurface = new PlaneSurface(_lineLoop, planeSurfaces.size());
    planeSurfaces.push_back(planeSurface);
    return planeSurface;
}

Inclusion *Geometry::addInclusion(const double &_a, const double &_b, const double &_alpha, const double &_xc, const double &_yc)
{
    Inclusion *incl = new Inclusion(inclusions.size(), _a, _b, _alpha, _xc, _yc);
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


void Geometry::InitializeGmshAPI()
{
    gmsh::initialize();
    gmsh::model::add(name);

    //Generating the geometry
    //Adding the HOST geometry -> HOST is the main domain

    for (auto point : points){
        gmsh::model::geo::addPoint(point->getX(), point->getY(), point->getZ(), point->getLC(), point->getIndex()+1);
        int aux = point->getIndex()+1;
        std::cout << "gmsh::model::geo::addPoint(" << point->getX() << ", " << point->getY() << ", " << point->getZ() << ", " << point->getLC() << ", " << aux << ")" << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
    for (auto line : lines){
        gmsh::model::geo::addLine(line->getPoint(0)->getIndex()+1, line->getPoint(1)->getIndex()+1, line->getIndex()+1);
        std::cout << "gmsh::model::geo::addLine(" << line->getPoint(0)->getIndex()+1 << ", " << line->getPoint(1)->getIndex()+1 << ", " << line->getIndex()+1 << ")" << std::endl;}
    std::cout << "--------------------------------" << std::endl;
        //Valid only for lines with 2 points

    for (auto lineLoop : lineLoops)
    {
        for (int i = 0; i < lineLoop->getNumLines(); i++)
            linesIndexes.push_back(lineLoop->getLine(i)->getIndex()+1);

        gmsh::model::geo::addCurveLoop(linesIndexes, lineLoop->getIndex()+1);
        std::cout << "gmsh::model::geo::addCurveLoop({";
        for (int i = 0; i < lineLoop->getNumLines(); i++)
        {
            std::cout << lineLoop->getLine(i)->getIndex()+1;
            if (i < lineLoop->getNumLines() - 1)
                std::cout << ", ";
        }
        std::cout << "}, " << lineLoop->getIndex()+1 << ")" << std::endl;
    }

    for (auto planeSurface : planeSurfaces){
        gmsh::model::geo::addPlaneSurface({planeSurface->getLineLoop()->getIndex()+1}, planeSurface->getIndex()+1);
        std::cout << "gmsh::model::geo::addPlaneSurface({" << planeSurface->getLineLoop()->getIndex()+1 << "}, " << planeSurface->getIndex()+1 << ")" << std::endl;
    }

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(dim);
    gmsh::write(name + ".msh");
    gmsh::fltk::run();
    gmsh::finalize();
}