#include "../headers/Solid.h"

Solid::Solid() {}
Solid::Solid(const std::string _name) : name(_name) {}
Solid::~Solid() {}

void Solid::readGeometry(const std::string &_filename)
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Reading geometry from file: %s\n", _filename.c_str());

    std::ifstream file(_filename); // Open file
    std::string line;

    while (line != "*NODES")
        std::getline(file, line);

    while (line != "*ELEMENTS")
    {
        std::getline(file, line);
        if (line != "*ELEMENTS")
        {
            int index;
            std::vector<double> coords(3);

            std::istringstream iss(line); // Uses istringstream to parse/analyze the line
            iss >> index >> coords[0] >> coords[1] >> coords[2];
            nodes.push_back(new Node(index - 1, coords));
        }
    }

    numNodes = nodes.size();
    std::cout << "Total number of nodes: " << numNodes << std::endl;

    while (line != "*Element, type=CPS3, 2DElements_Only")
    {
        std::getline(file, line);
        if (line != "*Nset, nset=inferiorBoundaryNodes")
        {
            int index;
            std::vector<int> connec(3);

            std::istringstream iss(line);
            iss >> index >> connec[0] >> connec[1] >> connec[2];
            std::cout << index << " " << connec[0] << " " << connec[1] << " " << connec[2] << std::endl;
        }
    }

    file.close();
}