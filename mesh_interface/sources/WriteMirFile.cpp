#include "../headers/Geometry.h"

void Geometry::writeMeshInfo()
{
    int numElNodes[] = {2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56, 64, 125};
    PetscPrintf(PETSC_COMM_WORLD, "Generating Output file...\n");
    std::string fileName = name + ".mir";
    std::ofstream file(fileName.c_str()); // Opening the output file
    file << std::fixed;                   // Setting the output to fixed point notation

    file << "*HEADING" << std::endl;
    file << fileName.c_str() << std::endl;

    file << "*MATERIALS" << std::endl;
    file << materials.size() << std::endl;

    for (auto m : materials)
        file << m->getIndex() + 1 << " " << m->getPoisson() << " " << m->getYoungModulus() << " " << m->getPlaneAnalysis() << std::endl;

    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);
    std::unordered_map<std::string, int> physicalEntities;

    file << "*PHYSICALNAMES" << std::endl;
    file << physicalGroups.size() << std::endl;

    for (const auto group : physicalGroups)
    {
        std::string groupName;
        gmsh::model::getPhysicalName(group.first, group.second, groupName);
        file << group.first << " " << group.second << " " << groupName << " ";
        physicalEntities[groupName] = group.second;

        if (groupName[0] == 'l')
        {
            for (auto l : lines)
                if ((l->getEntityName() == groupName) && (l->getMaterial() != nullptr))
                    file << l->getMaterial()->getIndex() + 1 << " " << l->getArea() << " " << l->getElementType() << " ";
        }
        else if (groupName[0] == 's')
        {

            for (auto s : planeSurfaces)
                if ((s->getEntityName() == groupName) && (s->getMaterial() != nullptr))
                    file << s->getMaterial()->getIndex() + 1 << " " << s->getThickness() << " " << s->getElementType() << " ";
        }
        file << std::endl;
    }

    file << "*NODES" << std::endl;
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

    std::size_t numNodes = nodeTags.size();
    file << numNodes << std::endl;
    for (int i = 0; i < numNodes; i++)
        file << i + 1 << " " << coord[3 * i] << " " << coord[3 * i + 1] << " " << coord[3 * i + 2] << std::endl;

    file << "*ELEMENTS" << std::endl;

    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    std::size_t numElems = 0, count = 1;
    std::vector<std::string> elem;

    for (const auto group : physicalGroups)
    {
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(group.first, group.second, entities);

        for (const auto entity : entities)
        {
            gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, group.first, entity);

            for (auto elemTags : elemTags)
                numElems += elemTags.size();

            for (int i = 0; i < elemTypes.size(); i++)
                for (int j = 0; j < elemTags[i].size(); j++)
                {
                    std::string auxString = std::to_string(count) + " " + std::to_string(elemTypes[i]) + " " + std::to_string(group.second) + " ";
                    for (int k = 0; k < numElNodes[elemTypes[i] - 1]; k++)
                        auxString += std::to_string(elemNodeTags[i][numElNodes[elemTypes[i] - 1] * j + k]) + " ";
                    elem.push_back(auxString);
                    count++;
                }
        }
    }

    file << numElems << std::endl;

    for (const auto e : elem)
        file << e << std::endl;

    file << "*BOUNDARY" << std::endl;
    file << boundaryConditions.size() << std::endl;

    for (auto bc : boundaryConditions) // In physicalEntities[bc->getEntityname()] we give the physical name and the function returns the physical tag
    {
        file << bc->getIndex() + 1 << " " << bc->getBType() << " " << physicalEntities[bc->getEntityName()] << " ";
        for (auto dofValues : bc->getDOFValues())
            file << dofValues.first << " " << dofValues.second << " ";
        file << std::endl;
    }
}