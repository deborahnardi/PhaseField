#include "../headers/Solid.h"

Solid::Solid() {}
Solid::Solid(const std::string _name) : name(_name) {}
Solid::~Solid() {}

std::vector<std::string> split(std::string str, char delim)
{
    std::istringstream input(str);
    std::vector<std::string> results;
    std::string token; // Token is a substring of the input string

    while (getline(input, token, delim))
        results.push_back(token);

    return results;
}

void Solid::readGeometry(const std::string &_filename)
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Reading geometry from file: %s\n", _filename.c_str());

    std::ifstream file(_filename); // Open file
    std::string line;

    while (line != "*NODES")
        std::getline(file, line);

    while (std::getline(file, line))
    {
        // != std::string::nos -> It is used to indicate that a search or operation has failed. For example, if you use the find function to search for a substring and the substring is not found, find returns std::string::npos to indicate this failure.
        if (line.find("*Element") != std::string::npos) // If the line contains the string "*Element"
        {
            std::cout << "The code has found the string: " << line << std::endl;
            break;
        }

        int index;
        std::vector<double> coords(3);

        std::istringstream iss(line); // Uses istringstream to parse/analyze the line
        iss >> index >> coords[0] >> coords[1] >> coords[2];
        nodes.push_back(new Node(index - 1, coords));
    }

    numNodes = nodes.size();

    std::cout << "There are: " << numNodes << " nodes" << std::endl;

    while (std::getline(file, line))
    {
        if (line.find("*") != std::string::npos) // If "*" is not found, line.find() returns std::string::npos
        {
            std::cout << "The code has found the string: " << line << std::endl;
            std::vector<std::string> result = split(line, ',');

            const auto &lastString = result.back();
            std::vector<std::string> result2 = split(lastString, '=');
            std::string lastString2 = result2.back();

            std::getline(file, line);
            std::vector<Node *> nodesFromSet;
            std::istringstream iss(line);

            for (int i = 0; i < iss.str().size(); i++)
            {
                int index;
                iss >> index;
                std::cout << "Index: " << index << std::endl;
                nodesFromSet.push_back(nodes[index - 1]);
            }

            nodeSets.push_back(new NodeSet(lastString2, nodesFromSet));
            break;
        }

        int index;
        std::vector<int> elemConnectivity(3);
        std::istringstream iss(line);
        iss >> index >> elemConnectivity[0] >> elemConnectivity[1] >> elemConnectivity[2];

        std::vector<Node *> elemNodes;
        for (int i = 0; i < 3; i++)
            elemNodes.push_back(nodes[elemConnectivity[i] - 1]);

        elements.push_back(new Solid2D(index - 1, elemNodes));
    }

    std::cout << "There are: " << elements.size() << " elements" << std::endl;

    num2DElements = elements.size();

    file.close();
}