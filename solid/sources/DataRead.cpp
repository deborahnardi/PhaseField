#include "../headers/FEM.h"

std::vector<std::string> split(std::string str, char delim)
{
    std::istringstream input(str);
    std::vector<std::string> results;
    std::string token; // Token is a substring of the input string

    while (getline(input, token, delim))
        results.push_back(token);

    return results;
}

void FEM::removeNonDiscritizedNodes(std::vector<Node *> &_nodes)
{
    /*
    remove_if :: (where it begins, where it ends, condition)
    lambda function :: [](const Node *node) { return node->getIsDiscritized() == false; } :: [](input parameters) { return output; }
    lambda function iterates over the vector of nodes and returns the nodes that are not discretized
    ! :: inverts the condition, if getIsDiscritized() == true, ! makes it false for the remove_if function, so it is not removed;
    on the other hand, if getIsDiscritized() == false, ! makes it true for the remove_if function, so it is removed
    */
    auto newEnd = std::remove_if(_nodes.begin(), _nodes.end(),
                                 [](const Node *node)
                                 { return !node->getIsDiscritized(); });
    _nodes.erase(newEnd, _nodes.end());
}

void FEM::renumberNodesIndexes(std::vector<Node *> &_nodes)
{
    /*
    [&] :: captures all variables by reference, this is necessary for lambda function to access the variables and modify them
    */
    int newIndex = 0;
    std::for_each(_nodes.begin(), _nodes.end(), [&](Node *node)
                  { node->setIndex(newIndex++); });
}

void FEM::readGeometry(const std::string &_filename)
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
        {
            std::vector<std::string> result = split(line, ',');        // Split the line by the character ','
            const auto &lastString = result.back();                    // Get the last string from the vector
            std::vector<std::string> result2 = split(lastString, '='); // Split the last string by the character '='
            const auto lastString2 = result2.back();                   // Get the last string from the vector
            abaqusElementType = lastString2;                           // Set the elementType to the last string
            break;
        }

        int index;
        std::vector<double> coords(3);

        std::istringstream iss(line); // Uses istringstream to parse/analyze the line
        iss >> index >> coords[0] >> coords[1] >> coords[2];
        nodes.push_back(new Node(index - 1, coords));
    }

    if (abaqusElementType == "CPS3")
        readGeometryCPS3(file);
    else if (abaqusElementType == "T3D2")
        readGeometryT3D2(file);
}

void FEM::readGeometryCPS3(std::ifstream &file)
{
    std::string currentName;
    std::string line;

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

        for (auto node : elemNodes) // Adding DOFs to the nodes
        {
            node->setIsDiscritized();
            node->addDOF(new DOF(X));
            node->addDOF(new DOF(Y));
        }
    }

    // std::cout << "There are: " << elements.size() << " 2D elements" << std::endl;

    num2DElements = elements.size();

    std::vector<Node *> elemConnectivity;
    std::vector<Element *> elemSets;
    bool flagNode = false, flagElem = false;

    for (; std::getline(file, line);)
    {
        if (line.find("*BOUNDARY") != std::string::npos)
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
                boundaryElements.push_back(new BoundaryElement(boundaryElements.size(), elemConnectivity));
                elemSets.push_back(elements.back());
                for (auto node : elemConnectivity)
                {
                    node->setIsDiscritized();
                    node->addDOF(new DOF(X));
                    node->addDOF(new DOF(Y));
                }
                elemConnectivity.clear();
            }
        }
    }

    removeNonDiscritizedNodes(nodes);
    renumberNodesIndexes(nodes);

    numNodes = nodes.size();
    numBoundaryElements = boundaryElements.size();
    // std::cout << "There are: " << numBoundaryElements << " boundary elements" << std::endl;

    //***** SETTING DOFs

    nDOFs = 0;
    for (auto node : nodes)
        for (auto dof : node->getDOFs())
        {
            dof->setIndex(nDOFs++); // Correspondent index in the global DOF vector
            // std::cout << "Node index: " << node->getIndex() << " DOF index: " << dof->getIndex() << std::endl;
            globalDOFs.push_back(dof); // Adding the DOF to the global DOF vector
        }

    // std::cout << "There are: " << nDOFs << " global DOFs" << std::endl;

    for (; std::getline(file, line);)
    {
        if (line.find("*BOUNDARY") != std::string::npos)
            continue;
        else if (line.find("*CLOAD") != std::string::npos)
            break;
        else
        {
            std::vector<std::string> result = split(line, ' ');
            std::string entityName = result[0];
            int dof = std::stoi(result[1]);
            int value = std::stoi(result[2]);
            BoundaryType bdType = DIRICHLET;

            for (auto nSet : nodeSets)
                if (entityName == nSet->getName())
                    nSet->addCondition(bdType, static_cast<DOFType>(dof), value);
        }
    }

    for (; std::getline(file, line);)
    {
        if (line.find("*END") != std::string::npos)
            break;
        else
        {
            std::vector<std::string> result = split(line, ' ');
            std::string entityName = result[0];
            int dof = std::stoi(result[1]);
            int value = std::stoi(result[2]);
            BoundaryType bdType = NEUMANN;

            for (auto nSet : nodeSets)
                if (entityName == nSet->getName())
                    nSet->addCondition(bdType, static_cast<DOFType>(dof), value);
        }
    }

    for (auto dof : globalDOFs)
        if (dof->isDirichlet())
            numDirichletDOFs++;
        else if (dof->isNeumann())
            numNeumannDOFs++;

    // std::cout << "There are: " << numDirichletDOFs << " Dirichlet DOFs" << std::endl;
    // std::cout << "There are: " << numNeumannDOFs << " Neumann DOFs" << std::endl;
    file.close();

    PetscPrintf(PETSC_COMM_WORLD, "Geometry read successfully!\n");
    PetscPrintf(PETSC_COMM_WORLD, "Number of nodes: %d\n", numNodes);
    PetscPrintf(PETSC_COMM_WORLD, "Number of global DOFs: %d\n", nDOFs);
    PetscPrintf(PETSC_COMM_WORLD, "Number of 2D elements: %d\n", num2DElements);
    PetscPrintf(PETSC_COMM_WORLD, "Number of boundary elements: %d\n", numBoundaryElements);
}

void FEM::readGeometryT3D2(std::ifstream &file)
{
    std::string line;
    std::vector<Node *> connec;

    for (; std::getline(file, line);)
    {
        if (line.find("*") != std::string::npos)
        {
            std::cout << "The code has found the string: " << line << std::endl;
            break;
        }

        connec.clear();
        std::vector<std::string> result = split(line, ' ');

        int index = std::stoi(result[0]);
        Node *node1 = nodes[std::stoi(result[1]) - 1];
        Node *node2 = nodes[std::stoi(result[2]) - 1];

        trussElements.push_back(new Truss(index - 1, node1, node2));

        // for (int i = 1; i < result.size(); i++)
        // {
        //     int index = std::stoi(result[i]);
        //     connec.push_back(nodes[index - 1]);
        // }

        // trussElements.push_back(new Truss(trussElements.size(), connec));
    }
    file.close();
}