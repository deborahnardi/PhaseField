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

    for (; std::getline(file, line);)
    {
        // != std::string::nos -> It is used to indicate that a search or operation has failed. For example, if you use the find function to search for a substring and the substring is not found, find returns std::string::npos to indicate this failure.
        if (line.find("*Element") != std::string::npos) // If the line contains the string "*Element"
            break;

        int index;
        std::vector<double> coords(3);

        std::istringstream iss(line); // Uses istringstream to parse/analyze the line
        iss >> index >> coords[0] >> coords[1] >> coords[2];
        nodes.push_back(new Node(index - 1, coords));
    }

    numNodes = nodes.size();

    std::cout << "There are: " << numNodes << " nodes" << std::endl;
    std::string currentName;

    for (; std::getline(file, line);)
    {
        if (line.find("*") != std::string::npos) // If "*" is not found, line.find() returns std::string::npos
        {
            std::vector<std::string> result = split(line, ',');
            const auto &lastString = result.back();
            std::vector<std::string> result2 = split(lastString, '=');
            const auto lastString2 = result2.back();
            nSetName = lastString2;
            currentName = "nSet";
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

    std::cout << "There are: " << elements.size() << " 2D elements" << std::endl;

    num2DElements = elements.size();

    std::vector<Node *> elemConnectivity;
    std::vector<Element *> elemSets;
    bool flagNode = false, flagElem = false;

    for (; std::getline(file, line);)
    {
        // std::cout << line << std::endl;
        if (line.find("*END") != std::string::npos)
        {
            elementSets.push_back(new ElementSet(elSetName, elemSets));
            elemSets.clear();

            break;
        }
        else if (line.find("*Elset") != std::string::npos)
        {
            std::vector<std::string> result = split(line, ',');
            const auto &lastString = result.back();
            std::vector<std::string> result2 = split(lastString, '=');
            const auto lastString2 = result2.back();
            elSetName = lastString2;
            currentName = "elSet";
            flagElem = false;

            if (flagNode == true)
            {
                nodeSets.push_back(new NodeSet(nSetName, elemConnectivity));
                elemConnectivity.clear();
            }
        }
        else if (line.find("*Nset") != std::string::npos)
        {
            std::vector<std::string> result = split(line, ',');
            const auto &lastString = result.back();
            std::vector<std::string> result2 = split(lastString, '=');
            const auto lastString2 = result2.back();
            nSetName = lastString2;
            currentName = "nSet";
            flagNode = false;

            if (flagElem == true)
            {

                elementSets.push_back(new ElementSet(elSetName, elemSets));
                elemSets.clear();
            }
        }
        else
        {
            if (currentName == "nSet")
            {
                // std::cout << "The code has found the string: " << line << std::endl;
                std::vector<std::string> result = split(line, ' ');

                for (const auto &token : result)
                {
                    int index = std::stoi(token);
                    elemConnectivity.push_back(nodes[index - 1]);
                }
                flagNode = true;
            }
            else
            {
                // std::cout << "The code has found the string: " << line << std::endl;
                std::vector<std::string> result = split(line, ' ');

                for (size_t i = 1; i < result.size(); i++)
                {
                    int index = std::stoi(result[i]);
                    elemConnectivity.push_back(nodes[index - 1]);
                }
                flagElem = true;
                elements.push_back(new BoundaryElement(elements.size(), elemConnectivity));
                elemSets.push_back(elements.back());
                elemConnectivity.clear();
                numBoundaryElements++;
            }
        }
    }

    numBoundaryElements = elements.size() - num2DElements;

    std::cout << "There are: " << numBoundaryElements << " boundary elements" << std::endl;

    // for (int i = 0; i < nodeSets.size(); i++)
    //     std::cout << "Node set: " << nodeSets[i]->getName() << " has " << nodeSets[i]->getNodes().size() << " nodes" << std::endl;
    // for (int i = 0; i < elementSets.size(); i++)
    //     std::cout << "Element set: " << elementSets[i]->getName() << " has " << elementSets[i]->getElements().size() << " elements" << std::endl;
    file.close();
}