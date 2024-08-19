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

    while (std::getline(file, line))
    {
        // != std::string::nos -> It is used to indicate that a search or operation has failed. For example, if you use the find function to search for a substring and the substring is not found, find returns std::string::npos to indicate this failure.
        if (line.find("*Element") != std::string::npos) // If the line contains the string "*Element"
        {
            break;
        }
        std::cout << line << std::endl;
    }

    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.find("*") != std::string::npos)
        {
            break;
        }

        int index;
        std::vector<int> elemConnectivity(3);
        std::istringstream iss(line);
        iss >> index >> elemConnectivity[0] >> elemConnectivity[1] >> elemConnectivity[2];
    }

    file.close();
}