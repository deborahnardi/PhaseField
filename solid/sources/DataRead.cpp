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

        if (result.size() > 4)
        {
            double gC = std::stod(result[4]);
            materials[i]->setGriffithCriterion(gC);
            double l0PF = std::stod(result[5]);
            materials[i]->setL0(l0PF);
        }
    }

    /*---------------------------------------------------------------------------------------
                                    PHYSICAL NAMES
    ----------------------------------------------------------------------------------------*/

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
            if (physicalEntities[i].dimension > elemDim)
                elemDim = physicalEntities[i].dimension;
        }
    }

    /*---------------------------------------------------------------------------------------
                                        NODES AND ELEMENTS
    ----------------------------------------------------------------------------------------*/

    while (line != "*NODES")
        std::getline(file, line);

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

        for (auto node : connectivity)
            node->setIsDiscritized();

        double value = physicalEntities[physicalEntity].value;
        int elemDim = physicalEntities[physicalEntity].dimension;

        switch (physicalEntities[physicalEntity].elementType)
        {
        case TRUSS_ELEMENT:
            numElNodes = 2;
            material = materials[physicalEntities[physicalEntity].material];
            elements.push_back(new Truss(elements.size(), elemDim, connectivity, material, physicalEntity, value));
            break;
        case SOLID_ELEMENT:
            numElNodes = 3; // 3-node triangle considered so far
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
            if (dof->getDOFType() != D)
            {
                dof->setIndex(globalDOFs.size());
                globalDOFs.push_back(dof);
            }
        }

    numElements = elements.size();
    nDOFs = globalDOFs.size(); // Only displacement DOFs are considered
    for (auto n : nodes)
        if (n->getIsDiscritized())
            discritizedNodes.push_back(n);

    for (auto n : nodes)
    {
        DOF *damageDOF = n->getDOFs()[2];
        damageDOF->setIndex(n->getIndex());
    }

    // Print damage DOFs
    for (auto n : nodes)
    {
        DOF *damageDOF = n->getDOFs()[2];
        std::cout << damageDOF->getIndex() << std::endl;
    }

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
    PetscPrintf(PETSC_COMM_WORLD, "Number of DOFs: %d\n", nDOFs);

    findNeighbours();
    decomposeElements(rhs, solution);
    matrixPreAllocation(IIIstart, IIIend);
    createPETScVariables(matrix, rhs, solution, nDOFs, true);
}

void FEM::findNeighbours()
{
    nodeNeighbours.resize(numNodes);

    /*
        std::vector<std::set<int>> nodeNeighbours declared in the private section of the class FEM;
        set is a data structure in c++ that is used to store unique elements;
        Insert the node itself and its neighbours, the insert method does not allow repeated elements.
    */

    for (auto elem : elements)
        for (auto node : elem->getElemConnectivity())
            for (auto node2 : elem->getElemConnectivity())
                nodeNeighbours[node->getIndex()].insert(node2->getIndex());
}

PetscErrorCode FEM::decomposeElements(Vec &b, Vec &x)
{
    MPI_Barrier(PETSC_COMM_WORLD); // Synchronizes all processes with PETSc communicator
    PetscPrintf(PETSC_COMM_WORLD, "Decomposing elements...\n");

    // ----------------------------------------------------------------
    // PARTIONING DOMAIN ELEMENTS
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, elements.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    // PARTIONING BOUNDARY ELEMENTS (FOR NEUMANN CONDITIONS)
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, bdElements.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IIstart, &IIend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    // PARTIONING GLOBAL DOFs
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, nDOFs);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IIIstart, &IIIend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    // PARTIONING PHASE FIELD DOFs
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, numNodes);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IstartPF, &IendPF);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);

    return ierr;
}

PetscErrorCode FEM::matrixPreAllocation(PetscInt start, PetscInt end)
{
    int rankLocalDOFs = end - start; // Number of nodes in the local partition

    d_nnz = new PetscInt[rankLocalDOFs]();
    o_nnz = new PetscInt[rankLocalDOFs]();

    for (auto node1 : nodes)
        for (auto dof1 : node1->getDOFs()) // Rows of the matrix
            if (dof1->getDOFType() != D)
                if (dof1->getIndex() >= IIIstart && dof1->getIndex() < IIIend)
                    for (auto node2 : nodeNeighbours[node1->getIndex()]) // Columns of the matrix
                        for (auto dof2 : nodes[node2]->getDOFs())
                            if (dof2->getDOFType() != D)
                                if (dof2->getIndex() >= IIIstart && dof2->getIndex() < IIIend)
                                    d_nnz[dof1->getIndex() - IIIstart]++; // - IIIstart to get the local index
                                else
                                    o_nnz[dof1->getIndex() - IIIstart]++;

    return ierr;
}

PetscErrorCode FEM::createPETScVariables(Mat &A, Vec &b, Vec &x, int mSize, bool showInfo) // mSize stands for matrix size, mSize = DOFs = rows = cols
{
    PetscLogDouble bytes;

    (size == 1)
        ? ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, mSize, mSize, NULL, d_nnz, &A)
        : ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, mSize, mSize, NULL, d_nnz, NULL, o_nnz, &A);
    CHKERRQ(ierr);

    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, mSize);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);

    if (showInfo && rank == 0)
    {
        ierr = PetscMemoryGetCurrentUsage(&bytes);
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes / (1024 * 1024));
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "Matrix and vectors created...\n");
        CHKERRQ(ierr);
    }

    delete[] d_nnz;
    delete[] o_nnz;

    return ierr;
}