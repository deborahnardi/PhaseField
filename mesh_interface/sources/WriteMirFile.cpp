#include "../headers/Geometry.h"

void Geometry::writeMeshInfo()
{
    // CREATING A GENERAL OUT PUT FILE, REGARDLESS OF THE DIMENSION OF THE GEOMETRY

    PetscPrintf(PETSC_COMM_WORLD, "Generating Output file...\n");
    std::string fileName = name + ".mir";
    std::ofstream file(fileName.c_str()); // Opening the output file

    file << "*HEADING" << std::endl;
    file << fileName.c_str() << std::endl;
    file << "*MATERIALS" << std::endl;
    file << materials.size() << std::endl;

    for (auto m : materials)
    {
        file << m->getIndex() << " " << m->getPoisson() << " " << m->getYoungModulus() << " " << m->getMatType() << " " << m->getPlaneAnalysis() << std::endl;
    }

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, nodeParams;
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true, false);

    std::set<int> writtenNodes;
    file << "*NODES" << std::endl;
    std::vector<std::tuple<int, double, double, double>> auxPrintNodes;
    int bdElems = 0;

    for (int i = 0; i < nodeTags.size(); i++)
    {
        // If the element is found, the iterator points to the location where the element is in the container.
        // If the element is not found, the iterator returned is equal to the iterator returned by end().

        if (writtenNodes.find(nodeTags[i]) == writtenNodes.end()) // If the node is not written yet
        {
            auxPrintNodes.push_back(std::make_tuple(nodeTags[i], nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]));
            writtenNodes.insert(nodeTags[i]);
        }
    }

    std::sort(auxPrintNodes.begin(), auxPrintNodes.end(), [](const std::tuple<int, double, double, double> &a, const std::tuple<int, double, double, double> &b)
              { return std::get<0>(a) < std::get<0>(b); }); // Reordering the nodes by their tags, from the smallest to the largest
                                                            // get<0>(a) gets the first element of the tuple a; get<0>(b) gets the first element of the tuple b
    file << auxPrintNodes.size() << std::endl;
    for (int i = 0; i < auxPrintNodes.size(); i++)
        file << std::get<0>(auxPrintNodes[i]) << " " << std::get<1>(auxPrintNodes[i]) << " " << std::get<2>(auxPrintNodes[i]) << " " << std::get<3>(auxPrintNodes[i]) << std::endl;

    auxPrintNodes.clear();
    int count = 0;

    for (int i = 0; i < elemTypes.size(); i++)
    {
        file << count << std::endl;
        for (int j = 0; j < elemTags[i].size(); j++)
            file << j + 1 << " " << elemNodeTags[i][2 * j] << " " << elemNodeTags[i][2 * j + 1] << std::endl;
    }

    file << "*Element, type=CPS3" << std::endl; // CPS3 -> 3-node triangular element, for gmsh -> elemType = 2
    gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2, -1);
    count = 0;
    for (int i = 0; i < elemTypes.size(); i++)
    {
        for (int j = 0; j < elemTags[i].size(); j++)
            count++;
    }

    for (int i = 0; i < elemTypes.size(); i++)
    {
        file << count << std::endl;
        for (int j = 0; j < elemTags[i].size(); j++)
            file << j + 1 << " " << elemNodeTags[i][3 * j] << " " << elemNodeTags[i][3 * j + 1] << " " << elemNodeTags[i][3 * j + 2] << std::endl;
    }

    // ********************************************************************************************************************

    int nodesPerLine = 8, elemsPerLine = 8; // Number of nodes per line
    count = 0;
    std::vector<std::pair<double, double>> coordsPoints;

    for (auto line : bdLines)
    {
        double x1 = line->getPoint(0)->getX();
        double y1 = line->getPoint(0)->getY();
        double x2 = line->getPoint(1)->getX();
        double y2 = line->getPoint(1)->getY();

        coordsPoints.push_back(std::make_pair(x1, y1));
        coordsPoints.push_back(std::make_pair(x2, y2));
    }

    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::cout << "*PHYSICAL GROUPS" << std::endl;

    std::vector<double> firstNodeCoords, lastNodeCoords, parametricCoords;
    int firstNode, lastNode, dimNode, tagNode;

    file << "*NSETS" << std::endl;
    file << physicalGroups.size() << std::endl;

    for (const auto &group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        // std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);
        for (const auto &entity : entities)
        {
            std::vector<int> elementTypes;
            std::vector<std::vector<std::size_t>> elementTags, nodeTags;
            gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);
            firstNode = nodeTags[0][0];
            lastNode = nodeTags[0][nodeTags[0].size() - 1];

            std::vector<double> firstNodeCoords, lastNodeCoords, parametricCoords;
            gmsh::model::mesh::getNode(firstNode, firstNodeCoords, parametricCoords, dimNode, tagNode);
            gmsh::model::mesh::getNode(lastNode, lastNodeCoords, parametricCoords, dimNode, tagNode);

            for (int i = 0; i < coordsPoints.size(); i++)
            {
                if (coordsPoints[2 * i].first == firstNodeCoords[0] && coordsPoints[2 * i].second == firstNodeCoords[1] && coordsPoints[2 * i + 1].first == lastNodeCoords[0] && coordsPoints[2 * i + 1].second == lastNodeCoords[1])
                {
                    // Line is on this boundary
                    std::cout << "Boundary: " << groupName << ", Line: " << (2 * i / 2) + 1 << std::endl;
                    std::cout << "----------------------------------" << std::endl;
                    count = 0;
                    writtenNodes.clear();
                    file << "*NSET, nset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;

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
                        file << "*ELSET, elset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;

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
    }
    // ********************************************************************************************************************
    for (const auto &group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

        if (groupName.find("Inclusion") != std::string::npos)
        {
            file << "*NSET, n=set" << groupName << std::endl;

            writtenNodes.clear();
            count = 0;

            for (const auto &entity : entities)
            {
                std::vector<int> elementTypes;
                std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);

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
            file << "*ELSET, elset=" << groupName << std::endl;
            for (const auto &entity : entities)
            {
                std::vector<int> elementTypes;
                std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);

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
    //**********************************************************************************************************************
    bdElems = 0;
    for (const auto &group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

        if (groupName.find("DomainInc") != std::string::npos)
        {
            file << "*ELSET, elset=" << groupName << std::endl;
            // Printing now the 2D elements of each inclusion
            for (const auto &entity : entities)
            {
                std::vector<int> elementTypes;
                std::vector<std::vector<std::size_t>> elementTags, nodeTags;
                gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);

                writtenNodes.clear();

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

    // ********************************************************************************************************************

    file << "*BOUNDARY" << std::endl;
    count = 0;
    for (auto *bc : boundaryConditions)
    {
        int dirichletOrNeumann = bc->getBType();

        if (dirichletOrNeumann == 0) // Dirichlet
        {
            for (auto dofValues : bc->getDOFValues())
            {
                count++;
            }
        }
    }

    file << count++ << std::endl;

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

    count = 0;
    file << "*CLOAD" << std::endl;

    for (auto *bc : boundaryConditions)
    {
        int dirichletOrNeumann = bc->getBType();

        if (dirichletOrNeumann == 1) // Neumann
        {
            for (auto dofValues : bc->getDOFValues())
            {
                count++;
            }
        }
    }

    file << count << std::endl;

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
}

// void Geometry::writeMeshInfo2D()
// {
//     std::vector<std::pair<int, int>> entities;
//     std::vector<std::size_t> nodeTags;
//     std::vector<double> nodeCoords, nodeParams;
//     std::vector<int> elemTypes;
//     std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
//     gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, -1, -1, true, false);

//     gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, -1, -1);

//     std::vector<int> lineElems;
//     std::vector<std::vector<std::size_t>> lineElemsTags, lineElemsNodeTags;

//     gmsh::model::mesh::getElements(lineElems, lineElemsTags, lineElemsNodeTags, 1, -1);
//     PetscPrintf(PETSC_COMM_WORLD, "Generating Output file...\n");

//     std::string fileName = name + ".mir";
//     std::ofstream file(fileName.c_str());

//     file << "*Heading" << std::endl;
//     file << fileName.c_str() << std::endl;

//     std::set<int> writtenNodes;

//     file << "*NODES" << std::endl;
//     std::vector<std::tuple<int, double, double, double>> auxPrintNodes;

//     for (int i = 0; i < nodeTags.size(); i++)
//     {
//         if (writtenNodes.find(nodeTags[i]) == writtenNodes.end())
//         {
//             auxPrintNodes.push_back(std::make_tuple(nodeTags[i], nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]));
//             writtenNodes.insert(nodeTags[i]);
//         }
//     }

//     std::sort(auxPrintNodes.begin(), auxPrintNodes.end(), [](const std::tuple<int, double, double, double> &a, const std::tuple<int, double, double, double> &b)
//               { return std::get<0>(a) < std::get<0>(b); });

//     for (int i = 0; i < auxPrintNodes.size(); i++)
//         file << std::get<0>(auxPrintNodes[i]) << " " << std::get<1>(auxPrintNodes[i]) << " " << std::get<2>(auxPrintNodes[i]) << " " << std::get<3>(auxPrintNodes[i]) << std::endl;

//     auxPrintNodes.clear();
//     // ********************************************************************************************************************

//     file << "*Element, type=CPS3" << std::endl; // CPS3 -> 3-node triangular element, for gmsh -> elemType = 2

//     for (int i = 0; i < elemTypes.size(); i++)
//     {
//         if (elemTypes[i] == 2)
//         {
//             for (int j = 0; j < elemTags[i].size(); j++)
//                 file << j + 1 << " " << elemNodeTags[i][3 * j] << " " << elemNodeTags[i][3 * j + 1] << " " << elemNodeTags[i][3 * j + 2] << std::endl;
//         }
//     }

//     // ********************************************************************************************************************

//     writtenNodes.clear();
//     std::set<std::string> addedGroups;

//     int nodesPerLine = 8, elemsPerLine = 8; // Number of nodes per line
//     int count = 0, countE = 0, bdElems = 0;
//     std::set<int> writtenBoundaryTags;

//     gmsh::vectorpair physicalGroups;
//     gmsh::model::getPhysicalGroups(physicalGroups);

//     std::cout << "*PHYSICAL GROUPS" << std::endl;

//     std::vector<std::pair<double, double>> coordsPoints;
//     std::vector<double> firstNodeCoords, lastNodeCoords, parametricCoords;
//     int firstNode, lastNode, dimNode, tagNode;

//     for (auto line : lines)
//     {
//         double x1 = line->getPoint(0)->getX();
//         double y1 = line->getPoint(0)->getY();
//         double x2 = line->getPoint(1)->getX();
//         double y2 = line->getPoint(1)->getY();

//         coordsPoints.push_back(std::make_pair(x1, y1));
//         coordsPoints.push_back(std::make_pair(x2, y2));
//     }

//     for (const auto &group : physicalGroups)
//     {
//         std::string groupName;
//         gmsh::model::getPhysicalName(group.first, group.second, groupName);
//         std::cout << "Group name: " << groupName << ", of dimension: " << group.first << ", and tag: " << group.second << std::endl;
//         std::vector<int> entities;
//         gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

//         for (const auto &entity : entities)
//         {
//             std::vector<int> elementTypes;
//             std::vector<std::vector<std::size_t>> elementTags, nodeTags;
//             gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, group.first, entity);
//             firstNode = nodeTags[0][0];
//             lastNode = nodeTags[0][nodeTags[0].size() - 1];
//             // std::cout << "First Node: " << firstNode << ", Last Node: " << lastNode << std::endl;
//             gmsh::model::mesh::getNode(firstNode, firstNodeCoords, parametricCoords, dimNode, tagNode);
//             gmsh::model::mesh::getNode(lastNode, lastNodeCoords, parametricCoords, dimNode, tagNode);

//             // Identifying the lines with the boundary created by Gmsh

//             for (int i = 0; i < coordsPoints.size(); i++)
//             {
//                 if (coordsPoints[2 * i].first == firstNodeCoords[0] && coordsPoints[2 * i].second == firstNodeCoords[1] && coordsPoints[2 * i + 1].first == lastNodeCoords[0] && coordsPoints[2 * i + 1].second == lastNodeCoords[1])
//                 {
//                     // Line is on this boundary
//                     std::cout << "Boundary: " << groupName << ", Line: " << (2 * i / 2) + 1 << std::endl;
//                     std::cout << "----------------------------------" << std::endl;
//                     count = 0;
//                     writtenNodes.clear();
//                     file << "*Nset, nset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;
//                     addedGroups.insert(groupName);

//                     for (std::size_t i = 0; i < elementTypes.size(); i++)
//                     {
//                         for (const auto &node : nodeTags[i])
//                         {
//                             if (writtenNodes.find(node) == writtenNodes.end())
//                             {
//                                 file << node;
//                                 count++;
//                                 writtenNodes.insert(node);

//                                 if (count == nodesPerLine)
//                                 {
//                                     file << std::endl;
//                                     count = 0;
//                                 }
//                                 else
//                                     file << " ";
//                             }
//                         }
//                         file << std::endl;
//                     }

//                     file << "*Elset, elset=" << "Boundary_" + std::to_string((2 * i / 2) + 1) << std::endl;

//                     for (std::size_t i = 0; i < elementTypes.size(); i++)
//                     {
//                         for (const auto &elem : elementTags[i])
//                         {
//                             bdElems++;
//                             file << bdElems << " ";
//                             std::size_t elementTag = elem;
//                             int elementType, dim, tag;
//                             std::vector<std::size_t> elemNodeTags;
//                             gmsh::model::mesh::getElement(elementTag, elementType, elemNodeTags, dim, tag);

//                             for (int i = 0; i < elemNodeTags.size(); i++)
//                             {
//                                 file << elemNodeTags[i] << " ";
//                             }

//                             file << std::endl;
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     // ********************************************************************************************************************
//     // /*
//     //     WRITING INCLUSIONS INFORMATIONS
//     //     **ATTENTION**: Please note that the element numbering for the boundary elements is the one created here.
//     //                    Similarly, the 2D elements present in the .mir file are numbered from 1 to the number of elements,
//     //                     i.e, they are renumbered.
//     // */

//     count = 0;

//     for (const auto &group : physicalGroups)
//     {

//         std::string groupName;
//         gmsh::model::getPhysicalName(group.first, group.second, groupName);

//         if (addedGroups.find(groupName) == addedGroups.end())
//         {
//             std::vector<int> entities;
//             gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

//             for (const auto &entity : entities)
//             {
//                 std::vector<std::pair<int, int>> boundaries;
//                 std::vector<std::pair<int, int>> dimTags = {{2, entity}};
//                 gmsh::model::getBoundary(dimTags, boundaries);

//                 if (groupName != "Host")
//                 {
//                     writtenNodes.clear();
//                     count = 0;
//                     file << "*Nset, nset=" << groupName << std::endl;

//                     for (const auto boundary : boundaries)
//                     {
//                         std::vector<int> elementTypes;
//                         std::vector<std::vector<std::size_t>> elementTags, nodeTags;
//                         gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

//                         for (std::size_t i = 0; i < elementTypes.size(); i++)
//                         {
//                             for (const auto &node : nodeTags[i])
//                             {
//                                 if (writtenNodes.find(node) == writtenNodes.end())
//                                 {
//                                     file << node;
//                                     count++;
//                                     writtenNodes.insert(node);

//                                     if (count == nodesPerLine)
//                                     {
//                                         file << std::endl;
//                                         count = 0;
//                                     }
//                                     else
//                                         file << " ";
//                                 }
//                             }
//                         }
//                     }

//                     if (count > 0 && count < nodesPerLine)
//                         file << std::endl;
//                     file << "*Elset, elset=" << groupName << std::endl;

//                     for (const auto boundary : boundaries)
//                     {
//                         std::vector<int> elementTypes;
//                         std::vector<std::vector<std::size_t>> elementTags, nodeTags;
//                         gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, boundary.first, boundary.second);

//                         bdElems++;
//                         file << bdElems << " ";
//                         writtenNodes.clear();

//                         for (std::size_t i = 0; i < elementTypes.size(); i++)
//                         {
//                             for (const auto &node : nodeTags[i])
//                             {
//                                 if (writtenNodes.find(node) == writtenNodes.end())
//                                 {
//                                     file << node << " ";
//                                     writtenNodes.insert(node);
//                                 }
//                             }
//                         }
//                         file << std::endl;
//                     }
//                 }
//             }
//         }
//     }
//     // ********************************************************************************************************************

//     file << "*BOUNDARY" << std::endl;

//     for (auto *bc : boundaryConditions)
//     {
//         int dirichletOrNeumann = bc->getBType();

//         if (dirichletOrNeumann == 0) // Dirichlet
//         {
//             for (auto dofValues : bc->getDOFValues())
//             {
//                 file << bc->getEntityname() << " ";
//                 file << dofValues.first << " " << dofValues.second;
//                 file << std::endl;
//             }
//         }
//     }

//     file << "*CLOAD" << std::endl;

//     for (auto *bc : boundaryConditions)
//     {
//         int dirichletOrNeumann = bc->getBType();

//         if (dirichletOrNeumann == 1) // Neumann
//         {
//             for (auto dofValues : bc->getDOFValues())
//             {
//                 file << bc->getEntityname() << " ";
//                 file << dofValues.first << " " << dofValues.second;
//                 file << std::endl;
//             }
//         }
//     }

//     file << "*END" << std::endl;
//     file.close();
// };