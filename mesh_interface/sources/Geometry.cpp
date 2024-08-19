#include "../headers/Geometry.h"

Geometry::Geometry() {}
Geometry::Geometry(const std::string _name)
    : name(_name) {}
Geometry::~Geometry() {}

Point *Geometry::addPoint(const std::vector<double> &_coordinates, const double &_lc)
{
    Point *point = new Point(_coordinates, _lc, points.size());
    points.push_back(point);

    int x = point->getX();
    int y = point->getY();

    if (x > maxX)
        maxX = x;
    if (x < minX)
        minX = x;
    if (y > maxY)
        maxY = y;
    if (y < minY)
        minY = y;

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

        int cInc = gmsh::model::occ::addCurveLoop({ellipseArcs[4 * inclusion->getIndex()], ellipseArcs[4 * inclusion->getIndex() + 1], ellipseArcs[4 * inclusion->getIndex() + 2], ellipseArcs[4 * inclusion->getIndex() + 3]});
        ellipseCurves.push_back(cInc);

        int pInc = gmsh::model::occ::addPlaneSurface({cInc});
        ellipseSurfaces.push_back(pInc);

        gmsh::model::occ::synchronize();
    }

    delete[] ellipseCoordinates;
};

void Geometry::getMeshInfo()
{ // dim = -1 -> all dimensions; tag = -1 -> all tags (get all nodes), includeBoundary = true, returnParametricCoordinates = false

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true, false);

    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

    std::vector<int> lineElems;
    std::vector<std::vector<std::size_t>> lineElemsTags, lineElemsNodeTags;

    gmsh::model::mesh::getElements(lineElems, lineElemsTags, lineElemsNodeTags, 1, -1);
};

void Geometry::writeMeshInfo()
{
    PetscPrintf(PETSC_COMM_WORLD, "Generating Output file...\n");

    std::string fileName = name + ".mir";
    std::ofstream file(fileName.c_str());

    file << "*Heading" << std::endl;
    file << fileName.c_str() << std::endl;

    std::set<int> writtenNodes;

    // ********************************************************************************************************************

    file << "*NODES" << std::endl;
    std::vector<int> nodeTagsMir;
    std::vector<double> nodeCoordsMir;

    for (int i = 0; i < nodeTags.size(); i++)
    {
        // If the element is found, the iterator points to the location where the element is in the container.
        // If the element is not found, the iterator returned is equal to the iterator returned by end().

        if (writtenNodes.find(nodeTags[i]) == writtenNodes.end()) // If the node is not written yet
        {
            file << nodeTags[i] << " " << nodeCoords[3 * i] << " " << nodeCoords[3 * i + 1] << " " << nodeCoords[3 * i + 2] << std::endl;
            nodeTagsMir.push_back(nodeTags[i]);
            nodeCoordsMir.push_back(nodeCoords[3 * i]);
            nodeCoordsMir.push_back(nodeCoords[3 * i + 1]);
            nodeCoordsMir.push_back(nodeCoords[3 * i + 2]);
            writtenNodes.insert(nodeTags[i]);
        }
    }
    // ********************************************************************************************************************

    file << "*ELEMENTS" << std::endl;
    file << "*Element, type=CPS3, 2DElements_Only" << std::endl; // CPS3 -> 3-node triangular element, for gmsh -> elemType = 2

    std::vector<int> elemTagsMir;     // Saving only the triangular elements
    std::vector<int> elemNodeTagsMir; // 3 nodes per element

    for (int i = 0; i < elemTypes.size(); i++)
    {
        if (elemTypes[i] == 2)
        {
            for (int j = 0; j < elemTags[i].size(); j++)
            {
                file << j + 1 << " " << elemNodeTags[i][3 * j] << " " << elemNodeTags[i][3 * j + 1] << " " << elemNodeTags[i][3 * j + 2] << std::endl;
                elemTagsMir.push_back(j + 1);
                elemNodeTagsMir.push_back(elemNodeTags[i][3 * j]);
                elemNodeTagsMir.push_back(elemNodeTags[i][3 * j + 1]);
                elemNodeTagsMir.push_back(elemNodeTags[i][3 * j + 2]);
            }
        }
    }

    // ********************************************************************************************************************

    writtenNodes.clear();

    int nodesPerLine = 8, elemsPerLine = 8; // Number of nodes per line
    int count = 0, countE = 0;
    std::set<int> writtenBoundaryTags;

    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::cout << "PHYSICAL GROUPS" << std::endl;
    for (const auto &group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);
        for (const auto &entity : entities)
        {
            std::cout << "Entidade: " << entity << std::endl;
            std::vector<std::pair<int, int>> boundaries;
            std::vector<std::pair<int, int>> dimTags = {{2, entity}};
            gmsh::model::getBoundary(dimTags, boundaries);
            std::cout << "BOUNDARY INFO: " << std::endl;
            for (const auto &boundary : boundaries)
            {
                if (writtenBoundaryTags.find(boundary.second) == writtenBoundaryTags.end())
                {
                    writtenBoundaryTags.insert(boundary.second);
                    std::cout << "Dimension: " << boundary.first << ", Tag: " << boundary.second << std::endl;
                    if (groupName == "Host")
                    {
                        std::vector<int> elementTypes;
                        std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

                        for (std::size_t i = 0; i < elementTypes.size(); i++)
                        {
                            writtenNodes.clear();
                            count = 0;
                            std::vector<double> coords;
                            std::vector<double> paramCoords;
                            int dim, tag;
                            int node = nodeTags[i][0]; // Always checking on the first node of the element
                            gmsh::model::mesh::getNode(node, coords, paramCoords, dim, tag);

                            if (coords[0] == minX && coords[1] == minY) // If (x = 0.0 and y = 0.0)
                            {
                                file << "*Nset, nset=inferiorBoundaryNodes" << std::endl;
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
                                file << "*Elset, elset=inferiorBoundaryElements" << std::endl;

                                for (const auto &elem : elementTags[i])
                                {
                                    file << elem;
                                    countE++;
                                    if (count == elemsPerLine)
                                    {
                                        file << std::endl;
                                        countE = 0;
                                    }
                                    else
                                        file << " ";
                                }
                                file << std::endl;
                            }
                            else if (coords[0] == minX && coords[1] == maxY) // If (x = 0.0 and y = edgeLength)
                            {
                                file << "*Nset, nset=leftBoundaryNodes" << std::endl;
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
                                file << "*Elset, elset=leftBoundaryElements" << std::endl;
                                countE = 0;
                                for (const auto &elem : elementTags[i])
                                {
                                    file << elem;
                                    countE++;
                                    if (count == elemsPerLine)
                                    {
                                        file << std::endl;
                                        countE = 0;
                                    }
                                    else
                                        file << " ";
                                }
                                file << std::endl;
                            }
                            else if (coords[0] = maxX && coords[1] == minY) // If (x = edgeLength and y = 0.0)
                            {
                                file << "*Nset, nset=rightBoundaryNodes" << std::endl;
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
                                file << "*Elset, elset=rightBoundaryElements" << std::endl;
                                countE = 0;
                                for (const auto &elem : elementTags[i])
                                {
                                    file << elem;
                                    countE++;
                                    if (count == elemsPerLine)
                                    {
                                        file << std::endl;
                                        countE = 0;
                                    }
                                    else
                                        file << " ";
                                }
                                file << std::endl;
                            }
                            else // If (x = edgeLength and y = edgeLength)
                            {
                                file << "*Nset, nset=superiorBoundaryNodes" << std::endl;
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
                                file << "*Elset, elset=superiorBoundaryElements" << std::endl;
                                countE = 0;
                                for (const auto &elem : elementTags[i])
                                {
                                    file << elem;
                                    countE++;
                                    if (count == elemsPerLine)
                                    {
                                        file << std::endl;
                                        countE = 0;
                                    }
                                    else
                                        file << " ";
                                }
                                file << std::endl;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "---------------------------" << std::endl;
    }
    // ********************************************************************************************************************
    /*
        WRITING INCLUSIONS INFORMATIONS
        **ATTENTION**: Please note that the element numbering for the boundary elements is the one that comes from Gmsh.
                       Differently, the 2D elements present in the .mir file are numbered from 1 to the number of elements,
                        i.e, they are renumbered.
    */

    count = 0;

    for (const auto &group : physicalGroups)
    {

        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);

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
                file << std::endl;
                countE = 0;
                file << "*Elset, elset=BoundaryElements_" << groupName << std::endl;

                for (const auto boundary : boundaries)
                {
                    std::vector<int> elementTypes;
                    std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

                    for (std::size_t i = 0; i < elementTypes.size(); i++)
                    {
                        for (const auto &elem : elementTags[i])
                        {
                            file << elem;
                            countE++;
                            if (countE == elemsPerLine)
                            {
                                file << std::endl;
                                countE = 0;
                            }
                            else
                                file << " ";
                        }
                    }
                }
                file << std::endl;
            }
        }
    }

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
        gmsh::model::occ::synchronize();
    }

    for (auto planeSurface : planeSurfaces)
    {
        gmsh::model::occ::addPlaneSurface({planeSurface->getLineLoop()->getIndex() + 1});
        gmsh::model::occ::synchronize();
    }

    generateInclusions();

    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize();

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

    getMeshInfo();
    writeMeshInfo();

    if (showInterface)
        gmsh::fltk::run();

    gmsh::clear();
    gmsh::finalize();
}