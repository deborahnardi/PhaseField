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
    line->setEntityName("Boundary_" + std::to_string(line->getIndex() + 1));
    lines.push_back(line);
    return line;
}

LineLoop *Geometry::addLineLoop(const std::vector<Line *> &_lines)
{
    LineLoop *lineLoop = new LineLoop(_lines, lineLoops.size());
    lineLoop->setEntityName("LineLoop_" + std::to_string(lineLoop->getIndex() + 1));
    lineLoops.push_back(lineLoop);
    return lineLoop;
}

PlaneSurface *Geometry::addPlaneSurface(LineLoop *_lineLoop)
{
    PlaneSurface *planeSurface = new PlaneSurface(_lineLoop, planeSurfaces.size());
    planeSurface->setEntityName("PlaneSurface_" + std::to_string(planeSurface->getIndex() + 1));
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

BoundaryCondition *Geometry::addBoundaryCondition(Point *point, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues)
{
    BoundaryCondition *bCondition = new BoundaryCondition(boundaryConditions.size(), point->getEntityName(), _dofValues, _bType);
    boundaryConditions.push_back(bCondition);
    return bCondition;
}

BoundaryCondition *Geometry::addBoundaryCondition(Line *line, const BoundaryType &_bType, const std::vector<std::pair<DOFType, double>> &_dofValues)
{
    BoundaryCondition *bCondition = new BoundaryCondition(boundaryConditions.size(), line->getEntityName(), _dofValues, _bType);
    std::string lineName = line->getEntityName();
    std::cout << "Boundary Condition has been added to: " << lineName << std::endl;
    boundaryConditions.push_back(bCondition);
    return bCondition;
}

void Geometry::generateInclusions()
{
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
        }

        int cInc = gmsh::model::occ::addCurveLoop({ellipseArcs[4 * inclusion->getIndex()], ellipseArcs[4 * inclusion->getIndex() + 1], ellipseArcs[4 * inclusion->getIndex() + 2], ellipseArcs[4 * inclusion->getIndex() + 3]}, lines.size() + ellipseCurves.size() + 1);
        ellipseCurves.push_back(cInc);

        int pInc = gmsh::model::occ::addPlaneSurface({cInc});
        ellipseSurfaces.push_back(pInc);

        gmsh::model::occ::synchronize();
    }

    delete[] ellipseCoordinates;
};

void Geometry::writeMeshInfo()
{
    // dim = -1 -> all dimensions; tag = -1 -> all tags (get all nodes), includeBoundary = true, returnParametricCoordinates = false

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true, false);

    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

    std::vector<int> lineElems;
    std::vector<std::vector<std::size_t>> lineElemsTags, lineElemsNodeTags;

    gmsh::model::mesh::getElements(lineElems, lineElemsTags, lineElemsNodeTags, 1, -1);
    PetscPrintf(PETSC_COMM_WORLD, "Generating Output file...\n");

    std::string fileName = name + ".mir";
    std::ofstream file(fileName.c_str());

    file << "*Heading" << std::endl;
    file << fileName.c_str() << std::endl;

    std::set<int> writtenNodes;

    // ********************************************************************************************************************

    file << "*NODES" << std::endl;
    std::vector<std::tuple<int, double, double, double>> auxPrintNodes;

    for (int i = 0; i < nodeTags.size(); i++)
    {
        // If the element is found, the iterator points to the location where the element is in the container.
        // If the element is not found, the iterator returned is equal to the iterator returned by end().

        if (writtenNodes.find(nodeTags[i]) == writtenNodes.end()) // If the node is not written yet
        {
            auxPrintNodes.push_back(std::make_tuple(nodeTags[i], nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]));
            // file << nodeTags[i] << " " << nodeCoords[3 * i] << " " << nodeCoords[3 * i + 1] << " " << nodeCoords[3 * i + 2] << std::endl;
            writtenNodes.insert(nodeTags[i]);
        }
    }

    std::sort(auxPrintNodes.begin(), auxPrintNodes.end(), [](const std::tuple<int, double, double, double> &a, const std::tuple<int, double, double, double> &b)
              { return std::get<0>(a) < std::get<0>(b); }); // Reordering the nodes by their tags, from the smallest to the largest
                                                            // get<0>(a) gets the first element of the tuple a; get<0>(b) gets the first element of the tuple b

    for (int i = 0; i < auxPrintNodes.size(); i++)
        file << std::get<0>(auxPrintNodes[i]) << " " << std::get<1>(auxPrintNodes[i]) << " " << std::get<2>(auxPrintNodes[i]) << " " << std::get<3>(auxPrintNodes[i]) << std::endl;

    auxPrintNodes.clear();
    // ********************************************************************************************************************

    file << "*Element, type=CPS3, 2DElements_Only" << std::endl; // CPS3 -> 3-node triangular element, for gmsh -> elemType = 2

    for (int i = 0; i < elemTypes.size(); i++)
    {
        if (elemTypes[i] == 2)
        {
            for (int j = 0; j < elemTags[i].size(); j++)
                file << j + 1 << " " << elemNodeTags[i][3 * j] << " " << elemNodeTags[i][3 * j + 1] << " " << elemNodeTags[i][3 * j + 2] << std::endl;
        }
    }

    // ********************************************************************************************************************

    writtenNodes.clear();
    std::set<std::string> addedGroups;

    int nodesPerLine = 8, elemsPerLine = 8; // Number of nodes per line
    int count = 0, countE = 0, bdElems = 0;
    std::set<int> writtenBoundaryTags;

    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::cout << "*PHYSICAL GROUPS" << std::endl;

    std::vector<std::pair<double, double>> coordsPoints;
    std::vector<double> firstNodeCoords, lastNodeCoords, parametricCoords;
    int firstNode, lastNode, dimNode, tagNode;

    for (auto line : lines)
    {
        double x1 = line->getPoint(0)->getX();
        double y1 = line->getPoint(0)->getY();
        double x2 = line->getPoint(1)->getX();
        double y2 = line->getPoint(1)->getY();

        coordsPoints.push_back(std::make_pair(x1, y1));
        coordsPoints.push_back(std::make_pair(x2, y2));
    }

    for (const auto &group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

        for (const auto &entity : entities)
        {
            std::vector<int> elementTypes;
            std::vector<std::vector<std::size_t>> elementTags, nodeTags;
            gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);
            firstNode = nodeTags[0][0];
            lastNode = nodeTags[0][nodeTags[0].size() - 1];
            // std::cout << "First Node: " << firstNode << ", Last Node: " << lastNode << std::endl;
            gmsh::model::mesh::getNode(firstNode, firstNodeCoords, parametricCoords, dimNode, tagNode);
            gmsh::model::mesh::getNode(lastNode, lastNodeCoords, parametricCoords, dimNode, tagNode);

            // Identifying the lines with the boundary created by Gmsh

            for (int i = 0; i < coordsPoints.size(); i++)
            {
                if (coordsPoints[2 * i].first == firstNodeCoords[0] && coordsPoints[2 * i].second == firstNodeCoords[1] && coordsPoints[2 * i + 1].first == lastNodeCoords[0] && coordsPoints[2 * i + 1].second == lastNodeCoords[1])
                {
                    // Line is on this boundary
                    std::cout << "Boundary: " << groupName << ", Line: " << (2 * i / 2) + 1 << std::endl;
                    std::cout << "----------------------------------" << std::endl;
                    count = 0;
                    writtenNodes.clear();
                    file << "*Nset, nset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;
                    addedGroups.insert(groupName);

                    for (std::size_t i = 0; i < elementTypes.size(); i++)
                    {
                        for (const auto &node : nodeTags[i])
                        {
                            if (writtenNodes.find(node) == writtenNodes.end())
                            {
                                file << node;
                                count++;
                                writtenNodes.insert(node);

                                if (count == nodesPerLine)
                                {
                                    file << std::endl;
                                    count = 0;
                                }
                                else
                                    file << " ";
                            }
                        }
                        file << std::endl;
                    }

                    file << "*Elset, elset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;

                    for (std::size_t i = 0; i < elementTypes.size(); i++)
                    {
                        for (const auto &elem : elementTags[i])
                        {
                            bdElems++;
                            file << bdElems << " ";
                            std::size_t elementTag = elem;
                            int elementType, dim, tag;
                            std::vector<std::size_t> elemNodeTags;
                            gmsh::model::mesh::getElement(elementTag, elementType, elemNodeTags, dim, tag);

                            for (int i = 0; i < elemNodeTags.size(); i++)
                            {
                                file << elemNodeTags[i] << " ";
                            }

                            file << std::endl;
                        }
                    }
                }
            }
        }
    }

    // ********************************************************************************************************************
    // /*
    //     WRITING INCLUSIONS INFORMATIONS
    //     **ATTENTION**: Please note that the element numbering for the boundary elements is the one created here.
    //                    Similarly, the 2D elements present in the .mir file are numbered from 1 to the number of elements,
    //                     i.e, they are renumbered.
    // */

    count = 0;

    for (const auto &group : physicalGroups)
    {

        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);

        if (addedGroups.find(groupName) == addedGroups.end())
        {
            std::vector<int> entities;
            gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

            for (const auto &entity : entities)
            {
                std::vector<std::pair<int, int>> boundaries;
                std::vector<std::pair<int, int>> dimTags = {{2, entity}};
                gmsh::model::getBoundary(dimTags, boundaries);

                if (groupName != "Host")
                {
                    writtenNodes.clear();
                    count = 0;
                    file << "*Nset, nset=" << groupName << std::endl;

                    for (const auto boundary : boundaries)
                    {
                        std::vector<int> elementTypes;
                        std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

                        for (std::size_t i = 0; i < elementTypes.size(); i++)
                        {
                            for (const auto &node : nodeTags[i])
                            {
                                if (writtenNodes.find(node) == writtenNodes.end())
                                {
                                    file << node;
                                    count++;
                                    writtenNodes.insert(node);

                                    if (count == nodesPerLine)
                                    {
                                        file << std::endl;
                                        count = 0;
                                    }
                                    else
                                        file << " ";
                                }
                            }
                        }
                    }

                    if (count > 0 && count < nodesPerLine)
                        file << std::endl;
                    file << "*Elset, elset=" << groupName << std::endl;

                    for (const auto boundary : boundaries)
                    {
                        std::vector<int> elementTypes;
                        std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

                        bdElems++;
                        file << bdElems << " ";
                        writtenNodes.clear();

                        for (std::size_t i = 0; i < elementTypes.size(); i++)
                        {
                            for (const auto &node : nodeTags[i])
                            {
                                if (writtenNodes.find(node) == writtenNodes.end())
                                {
                                    file << node << " ";
                                    writtenNodes.insert(node);
                                }
                            }
                        }
                        file << std::endl;
                    }
                }
            }
        }
    }
    // ********************************************************************************************************************

    file << "*BOUNDARY" << std::endl;

    for (auto *bc : boundaryConditions)
    {
        int dirichletOrNeumann = bc->getBType();

        if (dirichletOrNeumann == 0) // Dirichlet
        {
            for (auto dofValues : bc->getDOFValues())
            {
                file << bc->getEntityname() << " ";
                file << dofValues.first << " " << dofValues.second;
                file << std::endl;
            }
        }
    }

    file << "*CLOAD" << std::endl;

    for (auto *bc : boundaryConditions)
    {
        int dirichletOrNeumann = bc->getBType();

        if (dirichletOrNeumann == 1) // Neumann
        {
            for (auto dofValues : bc->getDOFValues())
            {
                file << bc->getEntityname() << " ";
                file << dofValues.first << " " << dofValues.second;
                file << std::endl;
            }
        }
    }

    file << "*END" << std::endl;
    file.close();
};

void Geometry::InitializeGmshAPI(const bool &showInterface)
{
    gmsh::initialize();
    gmsh::model::add(name);

    // Generating the geometry
    // Adding the HOST geometry -> HOST is the main domain

    for (auto point : points)
        gmsh::model::occ::addPoint(point->getX(), point->getY(), point->getZ());

    for (auto line : lines) // Valid only for lines with 2 points
        gmsh::model::occ::addLine(line->getPoint(0)->getIndex() + 1, line->getPoint(1)->getIndex() + 1, line->getIndex() + 1);

    for (auto lineLoop : lineLoops)
    {
        linesIndexes.clear();
        for (int i = 0; i < lineLoop->getNumLines(); i++)
            linesIndexes.push_back(lineLoop->getLine(i)->getIndex() + 1);

        gmsh::model::occ::addCurveLoop(linesIndexes, lineLoop->getIndex() + 1);
        lineLoop->setEntityName("LineLoop_" + std::to_string(lineLoop->getIndex() + 1));
        gmsh::model::occ::synchronize();
    }

    for (auto planeSurface : planeSurfaces)
    {
        int surfaceTag = gmsh::model::occ::addPlaneSurface({planeSurface->getLineLoop()->getIndex() + 1}, -1);
        planeSurface->setEntityName("PlaneSurface_" + std::to_string(planeSurface->getIndex() + 1));
        gmsh::model::occ::synchronize();
    }

    generateInclusions();

    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize();

    for (auto line : lines)
    {
        gmsh::model::addPhysicalGroup(1, {line->getIndex() + 1}, -1, "Boundary_" + std::to_string(line->getIndex() + 1));
        // std::cout << "Boundary_" + std::to_string(line->getIndex() + 1) << ", of dimension: 1, and tag: " << line->getIndex() + 1 << std::endl;
    }

    for (int i = 0; i < ellipseSurfaces.size() + 1; i++)
    {
        gmsh::model::addPhysicalGroup(2, {ellipseSurfaces[i]}, -1, "Inclusion_" + std::to_string(i + 1));
        gmsh::model::occ::synchronize();

        if (i == ellipseSurfaces.size() - 1)
        {
            gmsh::model::addPhysicalGroup(2, {ellipseSurfaces[i] + 1}, -1, "Host"); // Contain all the boundaries
            gmsh::model::occ::synchronize();
        }
    }

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

    gmsh::option::setNumber("Mesh.Algorithm", algorithm);

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(dim);
    gmsh::write(name + ".inp");

    writeMeshInfo();

    if (showInterface)
        gmsh::fltk::run();

    gmsh::clear();
    gmsh::finalize();
}