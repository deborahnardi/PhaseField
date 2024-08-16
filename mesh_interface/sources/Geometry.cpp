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

            int elp = gmsh::model::occ::addEllipseArc(startTag, centerTag, majorTag, endTag, lines.size() + ellipseArcs.size() + 1);
            ellipseArcs.push_back(elp);
        }

        int cInc = gmsh::model::occ::addCurveLoop({ellipseArcs[4 * inclusion->getIndex()], ellipseArcs[4 * inclusion->getIndex() + 1], ellipseArcs[4 * inclusion->getIndex() + 2], ellipseArcs[4 * inclusion->getIndex() + 3]}, lineLoops.size() + inclusion->getIndex() + 1);
        ellipseCurves.push_back(cInc);

        int pInc = gmsh::model::occ::addPlaneSurface({cInc}, planeSurfaces.size() + inclusion->getIndex() + 1);
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

    file << "*NODE" << std::endl;
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
    file << "*Element, type=CPS3" << std::endl; // CPS3 -> 3-node triangular element, for gmsh -> elemType = 2

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

    file << "*Nset, nset=inferiorBoundaryNodes" << std::endl; // Nodes on the inferior boundary

    writtenNodes.clear();
    int nodesPerLine = 8; // Number of nodes per line
    int count = 0;

    std::vector<std::pair<int, int>> boundaries;
    std::vector<std::pair<int, int>> dimTags = {{2, 1}};
    gmsh::model::occ::synchronize();
    gmsh::model::getBoundary(dimTags, boundaries);

    try
    {
        // Obter o contorno das superfícies
        gmsh::model::getBoundary(dimTags, boundaries, true, true, false);

        // Imprimir o contorno
        std::cout << "Contorno das superfícies:" << std::endl;
        for (const auto &boundary : boundaries)
        {
            std::cout << "Dimensão: " << boundary.first << ", Tag: " << boundary.second << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Erro: " << e.what() << std::endl;
    }

    //    for (int i = 0; i < nodeTags.size(); i++)
    //    {
    //        if (nodeCoords[3 * i + 1] == 0.0) // If y = 0
    //        {
    //            if (writtenNodes.find(nodeTags[i]) == writtenNodes.end())
    //            {
    //                file << nodeTags[i];
    //                writtenNodes.insert(nodeTags[i]);
    //                count++;
    //
    //                if (count == nodesPerLine)
    //                {
    //                    file << std::endl;
    //                    count = 0;
    //                }
    //                else
    //                    file << " ";
    //            }
    //        }
    //    }

    file << std::endl;

    // ********************************************************************************************************************

    file << "*Elset, elset=inferiorBoundaryElements" << std::endl; // Elements on the inferior boundary

    count = 0;

    std::set<int> writtenElems;

    for (int i = 0; i < elemTagsMir.size(); i++)
    {
        int numOfNodes = 0;
        for (int j = 0; j < 3; j++)
        {
            if (nodeCoordsMir[3 * (elemNodeTagsMir[3 * i + j]) - 2] == 0.0)
                numOfNodes++;
        }

        if (numOfNodes == 2) // 2 nodes of the element must be on the inferior boundary
        {
            if (writtenElems.find(elemTagsMir[i]) == writtenElems.end())
            {
                file << elemTagsMir[i];
                count++;
                writtenElems.insert(elemTagsMir[i]);

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

    file << std::endl;
    // ********************************************************************************************************************

    file << "*Nset, nset=superiorBoundaryNodes" << std::endl; // Nodes on the superior boundary

    writtenNodes.clear();

    // Identifying maximum y coordinate
    double maxY = 0.0;
    count = 0;
    for (int i = 0; i < nodeTagsMir.size(); i++)
    {
        if (nodeCoords[3 * i + 1] >= maxY)
            maxY = nodeCoords[3 * i + 1];
    }

    for (int i = 0; i < nodeTags.size(); i++)
    {
        if (nodeCoords[3 * i + 1] == maxY) // If y = maxY
        {
            if (writtenNodes.find(nodeTags[i]) == writtenNodes.end())
            {
                file << nodeTags[i];
                writtenNodes.insert(nodeTags[i]);
                count++;

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

    file << std::endl;

    // ********************************************************************************************************************
    file << "*Elset, superiorBoundaryElements" << std::endl; // Elements on the superior boundary

    count = 0;
    writtenElems.clear();

    for (int i = 0; i < elemTagsMir.size(); i++)
    {
        int numOfNodes = 0;
        for (int j = 0; j < 3; j++)
        {
            if (nodeCoordsMir[3 * (elemNodeTagsMir[3 * i + j]) - 2] == maxY)
                numOfNodes++;
        }

        if (numOfNodes == 2)
        {
            if (writtenElems.find(elemTagsMir[i]) == writtenElems.end())
            {
                file << elemTagsMir[i];
                count++;
                writtenElems.insert(elemTagsMir[i]);

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
    file << std::endl;

    // ********************************************************************************************************************

    file << "*Nset, nset=rightBoundaryNodes" << std::endl; // Nodes on the right side boundary

    writtenNodes.clear();
    double maxX = 0.0;
    count = 0;
    for (int i = 0; i < nodeTagsMir.size(); i++)
    {
        if (nodeCoords[3 * i] >= maxX)
            maxX = nodeCoords[3 * i];
    }

    for (int i = 0; i < nodeTags.size(); i++)
    {
        if (nodeCoords[3 * i] == maxX) // If x = maxX
        {
            if (writtenNodes.find(nodeTags[i]) == writtenNodes.end())
            {
                file << nodeTags[i];
                writtenNodes.insert(nodeTags[i]);
                count++;

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

    file << std::endl;

    // ********************************************************************************************************************

    file << "*Elset, rightBoundaryElements" << std::endl; // Elements on the right side boundary

    count = 0;
    writtenElems.clear();

    for (int i = 0; i < elemTagsMir.size(); i++)
    {
        int numOfNodes = 0;
        for (int j = 0; j < 3; j++)
        {
            if (nodeCoordsMir[3 * (elemNodeTagsMir[3 * i + j]) - 3] == maxX)
                numOfNodes++;
        }

        if (numOfNodes == 2)
        {
            if (writtenElems.find(elemTagsMir[i]) == writtenElems.end())
            {
                file << elemTagsMir[i];
                count++;
                writtenElems.insert(elemTagsMir[i]);

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
    file << std::endl;

    // ********************************************************************************************************************
    file << "*Nset, nset=leftBoundaryNodes" << std::endl; // Nodes on the left side boundary

    writtenNodes.clear();
    double minX = 0.0;
    count = 0;

    for (int i = 0; i < nodeTagsMir.size(); i++)
    {
        if (nodeCoords[3 * i] <= minX)
            minX = nodeCoords[3 * i];
    }

    for (int i = 0; i < nodeTags.size(); i++)
    {
        if (nodeCoords[3 * i] == minX) // If x = minX
        {
            if (writtenNodes.find(nodeTags[i]) == writtenNodes.end())
            {
                file << nodeTags[i];
                writtenNodes.insert(nodeTags[i]);
                count++;

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

    file << std::endl;

    // ********************************************************************************************************************

    file << "*Elset, leftBoundaryElements" << std::endl; // Elements on the left side boundary

    count = 0;
    writtenElems.clear();

    for (int i = 0; i < elemTagsMir.size(); i++)
    {
        int numOfNodes = 0;
        for (int j = 0; j < 3; j++)
        {
            if (nodeCoordsMir[3 * (elemNodeTagsMir[3 * i + j]) - 3] == minX)
                numOfNodes++;
        }

        if (numOfNodes == 2)
        {
            if (writtenElems.find(elemTagsMir[i]) == writtenElems.end())
            {
                file << elemTagsMir[i];
                count++;
                writtenElems.insert(elemTagsMir[i]);

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
    file << std::endl;

    // ********************************************************************************************************************
    std::string entityType;
    std::string expectedName = "Ellipse";
    std::vector<std::size_t> entityNodeTags;
    std::vector<double> entityNodeCoords;
    int nIncls = 0;

    gmsh::model::getEntities(entities);
    nodeTags.clear();
    nodeCoords.clear();
    writtenNodes.clear();
    writtenElems.clear();

    //    for (auto e : entities)
    //    {
    //        int entityDim = e.first, entityTag = e.second;
    //        gmsh::model::getType(entityDim, entityTag, entityType);
    //
    //        if (entityType == expectedName)
    //        {
    //            nIncls++;
    //
    //            file << "*Nset inclusion_" << nIncls << std::endl;
    //
    //            gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, entityDim, entityTag, true, false);
    //
    //            for (int i = 0; i < nodeTags.size(); i++)
    //            {
    //                entityNodeTags.push_back(nodeTags[i]);
    //                entityNodeCoords.push_back(nodeCoords[3 * i]);
    //                entityNodeCoords.push_back(nodeCoords[3 * i + 1]);
    //                entityNodeCoords.push_back(nodeCoords[3 * i + 2]);
    //            }
    //
    //            for (int i = 0; i < entityNodeTags.size(); i++)
    //            {
    //                if (writtenNodes.find(entityNodeTags[i]) == writtenNodes.end())
    //                {
    //                    file << entityNodeTags[i] << std::endl;
    //                    writtenNodes.insert(entityNodeTags[i]);
    //                }
    //            }
    //        }
    //    }
    //
    //    for (int i = 0; i < entityNodeTags.size(); i++)
    //        std::cout << "Node " << entityNodeTags[i] << " -> " << entityNodeCoords[3 * i] << " " << entityNodeCoords[3 * i + 1] << " " << entityNodeCoords[3 * i + 2] << std::endl;
    //    std::cout << "---------------------------" << std::endl;
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

    gmsh::model::addPhysicalGroup(2, {planeSurfaces[0]->getIndex() + 2}, -1, "Host"); // For some reason the index for host is 2
    gmsh::model::occ::synchronize();

    for (int i = 0; i < ellipseSurfaces.size(); i++)
    {
        gmsh::model::addPhysicalGroup(2, {ellipseSurfaces[i] + 1}, -1, "Inclusion" + std::to_string(i + 1));
        gmsh::model::occ::synchronize();
    }

    std::vector<std::pair<int, int>> physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::cout << "Grupos Físicos Criados:" << std::endl;
    for (const auto &group : physicalGroups)
    {
        std::cout << "Dimensão: " << group.first << std::endl;
        std::cout << "Tags: " << group.second << std::endl;
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