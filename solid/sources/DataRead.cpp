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

    // ********** MATERIALS **********

    while (line != "*MATERIALS")
        std::getline(file, line);

    int numMaterials;
    std::getline(file, line);
    numMaterials = std::stoi(line);
    for (int i = 0; i < numMaterials; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        int index = std::stoi(result[0]);
        double poisson = std::stod(result[1]);
        double youngModulus = std::stod(result[2]);
        PlaneAnalysis planeAnalysis = static_cast<PlaneAnalysis>(std::stoi(result[3]));
        materials.push_back(new Material(index, poisson, youngModulus, planeAnalysis));
    }

    // ********** PHYSICALNAMES **********

    while (line != "*PHYSICALNAMES")
        std::getline(file, line);

    int numPhysicalNames;
    std::getline(file, line);
    numPhysicalNames = std::stoi(line);

    struct PhysicalEntity // A strucuct is a user-defined data type that groups related data under one name
    {
        int indexType, material = -1, dimension;
        std::string name;
        double value = 0.0;
        ElementType elementType = NONE;
    };

    PhysicalEntity physicalEntities[numPhysicalNames]; // Array of PhysicalEntity

    for (int i = 0; i < numPhysicalNames; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        physicalEntities[i].dimension = std::stoi(result[0]);
        physicalEntities[i].indexType = std::stoi(result[1]);
        physicalEntities[i].name = result[2];

        if (result.size() > 3)
        {
            physicalEntities[i].material = std::stoi(result[3]) - 1;
            physicalEntities[i].value = std::stod(result[4]);
            physicalEntities[i].elementType = static_cast<ElementType>(std::stoi(result[5]));
        }
    }

    // ********** NODES **********

    while (line != "*NODES")
        std::getline(file, line);

    int numNodes;
    std::getline(file, line);
    numNodes = std::stoi(line);

    for (int i = 0; i < numNodes; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        int index = std::stoi(result[0]);
        double x = std::stod(result[1]);
        double y = std::stod(result[2]);
        double z = std::stod(result[3]);
        nodes.push_back(new Node(index - 1, {x, y, z}));
    }

    // ********** ELEMENTS **********

    while (line != "*ELEMENTS")
        std::getline(file, line);

    int numElements;
    std::getline(file, line);
    numElements = std::stoi(line);

    for (int i = 0; i < numElements; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        int index = std::stoi(result[0]);
        int gmshElemType = std::stoi(result[1]);                  // 2: 3-node triangle; 1: 2-node line; 15: 1-node point;
        int physicalEntity = std::stoi(result[2]) - 1;            // Physical entity is the number of the physical entity, i.e, 1 for p1, 2 for p2, etc.
        std::string name = physicalEntities[physicalEntity].name; // -1 because the physical entity vector starts at 0
        std::vector<Node *> connectivity;

        if (gmshElemType == 15)
            nodes[stoi(result[3]) - 1]->setPhysicalEntity(physicalEntity);

        for (int j = 3; j < result.size(); j++)
            connectivity.push_back(nodes[std::stoi(result[j]) - 1]);

        Material *material = nullptr;

        double value = physicalEntities[physicalEntity].value;
        int elemDim = physicalEntities[physicalEntity].dimension;

        switch (physicalEntities[physicalEntity].elementType)
        {
        case TRUSS_ELEMENT:
            material = materials[physicalEntities[physicalEntity].material];
            elements.push_back(new Truss(elements.size(), elemDim, connectivity, material, physicalEntity, value));
            break;
        case SOLID_ELEMENT:
            material = materials[physicalEntities[physicalEntity].material];
            elements.push_back(new Solid2D(elements.size(), elemDim, connectivity, material, physicalEntity));
            break;
        default:
            material = materials[physicalEntities[physicalEntity].material];
            bdElements.push_back(new BoundaryElement(bdElements.size(), elemDim, connectivity, material, physicalEntity));
            break;
        }
    }

    // ********** SETTING DOFS **********

    for (auto n : nodes)
        for (auto dof : n->getDOFs())
        {
            dof->setIndex(globalDOFs.size());
            globalDOFs.push_back(dof);
        }

    numElements = elements.size();
    nDOFs = globalDOFs.size();
    // ********** BOUNDARY CONDITIONS **********

    while (line != "*BOUNDARY")
        std::getline(file, line);

    int numBCs;
    std::getline(file, line);
    numBCs = std::stoi(line);

    for (int i = 0; i < numBCs; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        int index = std::stoi(result[0]);
        BoundaryType bdType = static_cast<BoundaryType>(std::stoi(result[1])); // 0: Dirichlet; 1: Neumann
        int physicalEntity = std::stoi(result[2]) - 1;
        int numAppliedBCs = (result.size() - 3) / 2;

        for (auto b : bdElements)
            if (b->getPhysicalEntity() == physicalEntity)
                for (int j = 0; j < numAppliedBCs; j++)
                    b->addCondition(bdType, static_cast<DOFType>(std::stoi(result[2 * j + 3])), std::stod(result[2 * j + 4]));
    }

    for (auto dof : globalDOFs)
        if (dof->isDirichlet())
            numDirichletDOFs++;

    dirichletBC = new PetscInt[numDirichletDOFs](); // Array that holds the dofs that have Dirichlet boundary conditions
    for (int i = 0, j = 0; i < nDOFs; i++)
        if (globalDOFs[i]->isDirichlet())
            dirichletBC[j++] = i;

    file.close();

    PetscPrintf(PETSC_COMM_WORLD, "Geometry read successfully!\n");
    PetscPrintf(PETSC_COMM_WORLD, "Number of nodes: %d\n", numNodes);
    PetscPrintf(PETSC_COMM_WORLD, "Number of elements: %d\n", numElements);
    PetscPrintf(PETSC_COMM_WORLD, "Number of boundary elements: %d\n", bdElements.size());
    PetscPrintf(PETSC_COMM_WORLD, "Number of Dirichlet DOFs: %d\n", numDirichletDOFs);
    PetscPrintf(PETSC_COMM_WORLD, "Number of DOFs: %d\n", globalDOFs.size());

    // decomposeElements();
    // matrixPreAllocation();
}

void FEM::decomposeElements()
{
    MPI_Barrier(PETSC_COMM_WORLD); // Synchronizes all processes with PETSc communicator
    PetscPrintf(PETSC_COMM_WORLD, "Decomposing elements...\n");

    // PARTIONING NODES
    Vec x;               // Declaring a vector x in PETSc
    PetscInt start, end; // Array of integers: start and end of the global DOFs for each process

    VecCreate(PETSC_COMM_WORLD, &x);            // Creates a vector x in PETSc
    VecSetSizes(x, PETSC_DECIDE, nodes.size()); // Sets the size of the vector x; PETSC_DECIDE means that PETSc will decide the size of the vector for each process
    VecSetFromOptions(x);                       // Good practice
    VecGetOwnershipRange(x, &start, &end);      // Gets the range of the global DOFs for each process

    VecDestroy(&x); // Destroys the vector x

    // for (int i = start; i < end; i++)
    //     partitionedNodes.push_back(nodes[i]);

    // std::cout << "Rank: " << rank << " start: " << start << " end: " << end << " Number of nodes: " << partitionedNodes.size() << " " << std::endl;

    // // PARTIONING DOMAIN ELEMENTS
    // VecCreate(PETSC_COMM_WORLD, &x);
    // VecSetSizes(x, PETSC_DECIDE, elements.size());
    // VecSetFromOptions(x);
    // VecGetOwnershipRange(x, &start, &end);

    // VecDestroy(&x);

    for (int i = start; i < end; i++)
        partitionedElements.push_back(elements[i]);

    // std::cout << "Rank: " << rank << " start: " << start << " end: " << end << " Number of elements: " << partitionedElements.size() << std::endl;

    // PARTIONING GLOBAL DOF
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, nDOFs);
    VecSetFromOptions(x);
    VecGetOwnershipRange(x, &start, &end);

    VecDestroy(&x);

    // for (int i = start; i < end; i++)
    //     partitionedDOFs.push_back(globalDOFs[i]);

    // std::cout << "Rank: " << rank << "start: " << start << " end: " << end << " Number of DOFs: " << partitionedDOFs.size() << std::endl;

    PetscPrintf(PETSC_COMM_WORLD, "Elements decomposed successfully!\n");

    // MISSING BOUNDARY ELEMENTS YET
}

void FEM::matrixPreAllocation()
{
    Vec x;
    PetscInt _start[size], _end[size]; // Arrays of integers: start and end of the global DOFs for each process

    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, nDOFs);
    VecSetFromOptions(x);
    VecGetOwnershipRange(x, &_start[rank], &_end[rank]);

    VecDestroy(&x);

    // Broadcasting the start and end of the global DOFs for each process
    for (int i = 0; i < size; i++)
    {
        MPI_Bcast(&_start[i], 1, MPI_INT, i, PETSC_COMM_WORLD); // i is the root process
        MPI_Bcast(&_end[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
    }

    PetscInt &start = _start[rank];
    PetscInt &end = _end[rank];
}