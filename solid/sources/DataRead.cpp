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
        int j = std::stoi(result[1]) - 1;
        physicalEntities[j].dimension = std::stoi(result[0]);
        physicalEntities[j].indexType = std::stoi(result[1]);
        physicalEntities[j].name = result[2];

        if (result.size() > 3)
        {
            physicalEntities[j].material = std::stoi(result[3]) - 1;
            physicalEntities[j].value = std::stod(result[4]);
            physicalEntities[j].elementType = static_cast<ElementType>(std::stoi(result[5]));
            if (physicalEntities[j].dimension > elemDim)
                elemDim = physicalEntities[j].dimension;
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

        double value = physicalEntities[physicalEntity].value;
        int elemDim = physicalEntities[physicalEntity].dimension;

        switch (physicalEntities[physicalEntity].elementType)
        {
        case TRUSS_ELEMENT:
            numElNodes = 2;
            material = materials[physicalEntities[physicalEntity].material];
            elements.push_back(new Truss(elements.size(), elemDim, connectivity, material, physicalEntity, value, params));
            break;
        case SOLID_ELEMENT:
            numElNodes = 3; // 3-node triangle considered so far
            material = materials[physicalEntities[physicalEntity].material];
            elements.push_back(new Solid2D(elements.size(), elemDim, connectivity, material, physicalEntity, params));
            break;
        default:
            material = materials[physicalEntities[physicalEntity].material];
            bdElements.push_back(new BoundaryElement(bdElements.size(), elemDim, connectivity, material, physicalEntity, params));
            break;
        }
    }

    // ********** SETTING DOFS **********
    int idx = 0;
    for (auto n : nodes)
        if (n->getIsDiscritized())
        {
            discritizedNodes.push_back(n);
            n->setIndex(idx);
            idx++;
        }

    for (auto n : discritizedNodes)
        for (auto dof : n->getDOFs())
        {
            if (dof->getDOFType() != D)
            {
                dof->setIndex(globalDOFs.size());
                globalDOFs.push_back(dof);
            }
        }

    numNodes = discritizedNodes.size();
    numElements = elements.size();
    nDOFs = globalDOFs.size(); // Only displacement DOFs are considered

    for (auto n : discritizedNodes)
    {
        DOF *damageDOF = n->getDOFs()[2];
        damageDOF->setIndex(n->getIndex());
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
        BoundaryType bdType = static_cast<BoundaryType>(std::stoi(result[1])); // 0: Dirichlet; 1: Neumann, 2: Damage
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

    while (line != "*TRACTION BOUNDARY")
        std::getline(file, line);
    std::getline(file, line);
    int nTractionBCs = std::stoi(line);

    for (int i = 0; i < nTractionBCs; i++)
    {
        std::getline(file, line);
        std::vector<std::string> result = split(line, ' ');
        std::string name = result[1];
        int index = -1;
        for (auto pe : physicalEntities)
            if (pe.name == name)
                index = pe.indexType - 1;

        for (auto be : bdElements)
            if (be->getPhysicalEntity() == index)
                tractionBd.push_back(be);
    }

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
    int numNodesAux = nodes.size();
    nodeNeighbours.resize(numNodesAux);

    /*
        std::vector<std::set<int>> nodeNeighbours declared in the private section of the class FEM;
        set is a data structure in c++ that is used to store unique elements;
        Insert the node itself and its neighbours, the insert method does not allow repeated elements.

        - nodeNeighbours contains the node itself and its neighbours, the insert method does not allow repeated elements, the elements are stored in ascending order and if i < j, then the node is not inserted in the set of neighbours of node j (due to symmetry of the matrix);

        - CRSNodeNeighbours is a vector that stores the neighbours of each node in a CRS format, i.e., the neighbours of node 0 are stored in the first positions of the vector, the neighbours of node 1 are stored in the next positions, and so on;

        - n2nCSRTotal stores the cumulative length of the neighbours of each node starting from 0;

        - nodeElementMapping stores the elements that each node belongs to;
    */

    for (auto elem : elements)
        for (auto node : elem->getElemConnectivity())
            for (auto node2 : elem->getElemConnectivity())
                if (node2->getIndex() >= node->getIndex())
                    nodeNeighbours[node->getIndex()].insert(node2->getIndex()); // Equivalent to n2n_symm_mat_total Ayrton

    for (int i = 0; i < nodeNeighbours.size(); i++)
    {
        std::set<int> neighbours = nodeNeighbours[i];
        for (auto n : neighbours)
            CRSNodeNeighbours.push_back(n); // Equivalent to n2n_symm_total Ayrton
    }

    n2nCSRTotal.push_back(0);
    for (int i = 0; i < nodeNeighbours.size(); i++)
        n2nCSRTotal.push_back(n2nCSRTotal[i] + nodeNeighbours[i].size()); // Equivalent to n2n_csr_total Ayrton

    // Print CumulativeNodeNeighbours
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
        for (int i = 0; i < n2nCSRTotal.size(); i++)
            std::cout << n2nCSRTotal[i] << " ";
        std::cout << std::endl;
    }

    n2e.resize(numNodesAux);
    for (auto elem : elements)
        for (auto node : elem->getElemConnectivity())
            n2e[node->getIndex()].insert(elem->getIndex()); // Equivalent to n2e Ayrton

    MPI_Barrier(PETSC_COMM_WORLD);
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
    // PARTIONING PHASE FIELD DOFs (damage DOFs)
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

    // ----------------------------------------------------------------
    // PARTIONING NODES
    // size = 0;
    // rank = 0;
    int numNodesByRank = numNodes / size; // Decomposing the nodes among the processors
    int remainderNodes = numNodes % size; // Remainder of the division

    std::vector<int> numNodesForEachRank; // Array that stores the number of nodes for each processor
    for (int i = 0; i < size; i++)
        numNodesForEachRank.push_back(numNodesByRank); // Filling the array with the number of nodes for each processor

    for (int i = 0; i < remainderNodes; i++)
        numNodesForEachRank[i]++; // Equivalent to numNodesForEachRank Ayrton

    int start = 0, end = 0;
    for (int i = 0; i <= rank; i++)
    {
        start = end;
        end += numNodesForEachRank[i];
    }

    PetscInt n = end - start;

    // Create vector with the number of elements specific to the rank
    PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
    PetscCall(VecSetSizes(x, n, numNodes)); // Use end-start for local size, not total size
    PetscCall(VecSetFromOptions(x));

    int Istart, Iend;
    PetscCall(VecGetOwnershipRange(x, &Istart, &Iend));
    PetscCall(VecDestroy(&x));

    nodesForEachRankCSR.resize(size + 1, 0); // size + 1 because it starts at 0; 0 is the value assigned to all the elements of the vector
    for (int i = 1; i <= size; i++)
        nodesForEachRankCSR[i] = nodesForEachRankCSR[i - 1] + numNodesForEachRank[i - 1]; // Equivalent to n2n_csr Ayrton

    /*
        FROM THIS POINT ON, THE nodeNeighbours and nodeElementMapping VECTORS ARE SETTLED FOR THE LOCAL PARTITION
    */

    localRankNodeNeighbours.resize(n);
    n2eLocal.resize(n);

    MPI_Barrier(PETSC_COMM_WORLD);

    for (int i = start; i < end; i++)
    {
        localRankNodeNeighbours[i - start] = nodeNeighbours[i]; // Equivalent to n2n_symm_mat Ayrton
        n2eLocal[i - start] = n2e[i];
    }

    // Flattening localRankNodeNeighbours
    std::vector<int> localRankNodeNeighboursFlattened;
    for (int i = 0; i < localRankNodeNeighbours.size(); i++)
    {
        std::set<int> neighbours = localRankNodeNeighbours[i];
        for (auto n : neighbours)
            localRankNodeNeighboursFlattened.push_back(n); // Equivalent to n2n_symm Ayrton
    }
    // Print localRankNodeNeighbours
    if (rank == 0)
    {
        for (int i = 0; i < localRankNodeNeighbours.size(); i++)
        {
            std::set<int> neighbours = localRankNodeNeighbours[i];
            std::cout << "Node " << i << " neighbours: ";
            for (auto n : neighbours)
                std::cout << n << " ";
            std::cout << std::endl;
        }
        std::cout << "--------------------------------" << std::endl;
    }

    std::vector<int> n2nCSRLocal;
    n2nCSRLocal.push_back(0);
    for (int i = 0; i < localRankNodeNeighbours.size(); i++)
        n2nCSRLocal.push_back(n2nCSRLocal[i] + localRankNodeNeighbours[i].size()); // Equivalent to n2n_csr Ayrton

    // Print n2nCSRLocal
    if (rank == 0)
    {
        std::cout << "n2n_CSR:" << std::endl;
        for (int i = 0; i < n2nCSRLocal.size(); i++)
            std::cout << n2nCSRLocal[i] << " ";
        std::cout << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }

    /*
        EXPLAINING THE FOLLOWING CODE BLOCK:
        - n2nDRank is a vector that stores the number of neighbours of each node that are in other ranks;
        Ex.:
            - Let's suppose that there are 8 nodes (0 to 7) and 3 ranks;
            - The nodes 0, 1, and 2 are in rank 0;
            - The nodes 3, 4, and 5 are in rank 1;
            - The nodes 6 and 7 are in rank 2;

            CONSIDERING rank = 0:

            - localRankNodeNeighbours = [(0, 1, 4), (1, 2, 4, 5, 6), (2, 3, 6)], i.e, the neighbours of node 0, 1 and 2, respectively;
            - numNodesForEachRank = [3, 3, 2], i.e, the number of nodes for each rank;
            - nodesForEachRankCompressed = [0, 3, 6, 8], i.e, the number of nodes for each rank compressed (0-2 for rank 1, 3-5 for rank 2, 6-7 for rank 3);
            - n2nDRank for rank 0 = [1, 3, 2] -> (0, 1, 4): only node 4 is in rank 1; (1, 2, 4, 5, 6): 3 nodes are in other ranks (4, 5, 6); (2, 3, 6): only node 6 is in rank 2;

    */

    std::vector<int> n2nDRank(numNodesForEachRank[rank], 0); // For ex., if rank = 1, numNodesForEachRank[1] = 3, then n2nDrank.size() = 3
    for (int i = 0; i < numNodesForEachRank[rank]; i++)
    {
        std::set<int> neighbours = localRankNodeNeighbours[i];
        int greaterThanRankCounter = 0;
        for (auto n : neighbours)
            if (n > n2nCSRLocal[rank + 1] - 1) // -1 because the vector starts at 0
                greaterThanRankCounter++;
        n2nDRank[i] = greaterThanRankCounter;
    }

    // Print n2nDRank
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
        for (int i = 0; i < n2nDRank.size(); i++)
            std::cout << n2nDRank[i] << " ";
        std::cout << std::endl;
        std::cout << "--------------------------------" << std::endl;
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    /*
        EXPALINING THE FOLLOWING CODE BLOCK:

        FIRST LOOP:
        - nodesForEachRankCSR[rank] = [0, 3, 6, 8], i.e, the number of nodes for each rank compressed (0-2 for rank 1, 3-5 for rank 2, 6-7 for rank 3);
            Thus, node1 goes from 0 to 2 for rank 0, from 3 to 5 for rank 1, and from 6 to 7 for rank 2;

        SECOND LOOP:
        - jj goes over the neighbours of EACH NODE, i.e, in n2nCSRLocal we have [0, 3, 6, 8], so the neighbours of node 0 are stored in the positions 0 to 3, the neighbours of node 1 are stored in the positions 3 to 6, and so on;
            Thus, node2 for node 0 goes from 0 to 3 -> node2 assumes the values 0, 1, 4 (neighbours of node 0);
            node2 for node 1 goes from 3 to 6 -> node2 assumes the values 1, 2, 4, 5, 6 (neighbours of node 1);
            node2 for node 2 goes from 6 to 8 -> node2 assumes the values 2, 3, 6 (neighbours of node 2);

        eSameList is an array of arrays
    */

    std::vector<std::vector<int>> eSameList(localRankNodeNeighboursFlattened.size());
    for (int iNode1 = 0; iNode1 < n2nCSRLocal.size() - 1; iNode1++) // This loops iterates over the nodes of the local partition (size + 1)
    {
        int node1 = nodesForEachRankCSR[rank] + iNode1; // nodesForEachRankCSR[rank] is the first node of the rank

        for (int jj = n2nCSRLocal[iNode1]; jj < n2nCSRLocal[iNode1 + 1]; jj++) // This loop iterates over the neighbours of each node
        {

            int node2 = localRankNodeNeighboursFlattened[jj]; // The neighbour of the node

            std::vector<int> elemsFromNode1(n2e[node1].begin(), n2e[node1].end()); // Elements that node1 belongs to
            std::vector<int> elemsFromNode2(n2e[node2].begin(), n2e[node2].end()); // Elements that node2 belongs to

            // Couting how many elements are in common between node1 and node2
            int numSharedElems = 0;
            for (int i = 0; i < elemsFromNode1.size(); i++)
            {
                int elem1 = elemsFromNode1[i];
                for (int j = 0; j < elemsFromNode2.size(); j++)
                {
                    int elem2 = elemsFromNode2[j];
                    if (elem1 == elem2)
                        numSharedElems++;
                }
            }

            /*
                3 infos are needed to create the eSharedList array:
                - The element index;
                - The local index of node1;
                - The local index of node2;
            */

            std::vector<int> eSharedList(3 * numSharedElems, 0);
            int counter = -1;
            for (int i = 0; i < elemsFromNode1.size(); i++)
            {
                int elem1 = elemsFromNode1[i];
                for (int j = 0; j < elemsFromNode2.size(); j++)
                {
                    int elem2 = elemsFromNode2[j];
                    if (elem1 == elem2)
                    {
                        counter++;

                        // Get the local position of node1 and node2 for the element (from the connectivity)
                        std::vector<Node *> elemConnectivity = elements[elem1]->getElemConnectivity();
                        int localNode1 = elemConnectivity[0]->getIndex() == node1   ? 0
                                         : elemConnectivity[1]->getIndex() == node1 ? 1
                                         : elemConnectivity[2]->getIndex() == node1 ? 2
                                                                                    : -1;

                        int localNode2 = elemConnectivity[0]->getIndex() == node2   ? 0
                                         : elemConnectivity[1]->getIndex() == node2 ? 1
                                         : elemConnectivity[2]->getIndex() == node2 ? 2
                                                                                    : -1;

                        eSharedList[3 * counter] = elem1;
                        eSharedList[3 * counter + 1] = localNode1;
                        eSharedList[3 * counter + 2] = localNode2;
                    }
                }
            }

            eSameList[jj] = eSharedList;
        }
    }

    // Print eSameList
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
        for (int i = 0; i < eSameList.size(); i++)
        {
            std::vector<int> eSharedList = eSameList[i];
            std::cout << "Tuple " << i << " shared elements: ";
            for (int j = 0; j < eSharedList.size(); j++)
                std::cout << eSharedList[j] << " ";
            std::cout << std::endl;
        }
        std::cout << "--------------------------------" << std::endl;
    }

    // TOTAL NUMBER OF NONZERO (UPPER TRIANGULAR PART ONLY) ONLY FOR THE LOCAL PARTITION

    int totalNnz = 0;
    int nDOF = 2;
    for (int i = 0; i < numNodesForEachRank[rank]; i++)
    {
        std::set<int> neighbours = localRankNodeNeighbours[i];
        int numFriends = neighbours.size();
        totalNnz += nDOF * nDOF * numFriends - (nDOF * nDOF - (nDOF * (nDOF + 1)) / 2.0);
    }

    if (rank == 0)
        std::cout << "Total number of nonzeros for the local partition: " << totalNnz << std::endl;

    MPI_Barrier(PETSC_COMM_WORLD);

    // PREALLOCATION OF THE MPI MATRIX
    int m = numNodesForEachRank[rank] * nDOF; // Number of local lines at the rank. Ex.: 3 nodes at the local partition * 2 DOFs = 6 lines
    int nn = numNodes * nDOF;                 // Number of local columns at the rank. Ex.: 8 nodes * 2 DOFs = 16 columns

    d_nnz = new PetscInt[m]();
    o_nnz = new PetscInt[m]();

    for (int ii = 0; ii < numNodesForEachRank[rank]; ii++)
        for (int iDOF = 0; iDOF < nDOF; iDOF++)
        {
            d_nnz[nDOF * ii + iDOF] = nDOF * localRankNodeNeighbours[ii].size() - n2nDRank[ii] * nDOF - iDOF; // (n2nDRank[ii] * nDOF - iDOF) removes the DOFs that do not belong to the local partition
            o_nnz[nDOF * ii + iDOF] = n2nDRank[ii] * nDOF;
        }

    std::vector<int> val(totalNnz, 0);
    int kkn2n = 0, kk = 0;
    for (int iNode1 = 0; iNode1 < n2nCSRLocal.size() - 1; iNode1++)
    {
        int n1 = nodesForEachRankCSR[rank] + iNode1;
        std::vector<int> friendNodes(localRankNodeNeighbours[iNode1].begin(), localRankNodeNeighbours[iNode1].end());
        int numFriends = friendNodes.size();
        int iFriendCount = 0;

        for (auto n2 : friendNodes)
        {
            std::vector<int> elems = eSameList[kkn2n];

            /*  COMPUTE HERE THE Kglobal COMPONENTS ASSOCIATED TO THE INFLUENCE OF NODE n2 (COLUMN) ON NODE n1 (LINE)
                IF n1 == n2, ONLY 3 DIFFERENT COMPONENTS ARE COMPUTED
                CONSIDERING ONLY A 2D ANALYSIS, THOSE COMPONENTS MUST BE PLACED AT THE FOLLOWING POSITIONS ON VECTOR val:
             */

            if (n1 != n2)
            { // idof = 0,       jdof = 0
                int p1 = kk;
                // idof = 0,       jdof = 1
                int p2 = kk + 1;
                // idof = 1,       jdof = 0
                int p3 = kk + numFriends * nDOF - 1;
                // idof = 1,       jdof = 1
                int p4 = kk + numFriends * nDOF + 0;

                val[p1] += 1;
                val[p2] += 1;
                val[p3] += 1;
                val[p4] += 1;
            }
            else
            {
                // idof = 0,       jdof = 0
                int p1 = kk;
                // idof = 0,       jdof = 1
                int p2 = kk + 1;
                // idof = 1,       jdof = 1
                int p3 = kk + numFriends * nDOF;

                val[p1] += 1;
                val[p2] += 1;
                val[p3] += 1;
            }

            iFriendCount++;
            kkn2n++;
            kk += nDOF;
        }

        kk += nDOF * numFriends - (nDOF * nDOF - (nDOF * (nDOF + 1)) / 2.0);
    }

    // Print val
    MPI_Barrier(PETSC_COMM_WORLD);
    if (rank == 0)
    {
        for (int i = 0; i < val.size(); i++)
            std::cout << val[i] << " ";
        std::cout << std::endl;
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 1);

    return ierr;
}

PetscErrorCode FEM::matrixPreAllocation(PetscInt start, PetscInt end)
{
    // int rankLocalDOFs = end - start; // Number of nodes in the local partition

    // d_nnz = new PetscInt[rankLocalDOFs]();
    // o_nnz = new PetscInt[rankLocalDOFs]();

    // for (auto node1 : discritizedNodes)
    //     for (auto dof1 : node1->getDOFs()) // Rows of the matrix
    //         if (dof1->getDOFType() != D)
    //             if (dof1->getIndex() >= IIIstart && dof1->getIndex() < IIIend)
    //                 for (auto node2 : nodeNeighbours[node1->getIndex()]) // Columns of the matrix
    //                     for (auto dof2 : discritizedNodes[node2]->getDOFs())
    //                         if (dof2->getDOFType() != D)
    //                             if (dof2->getIndex() >= IIIstart && dof2->getIndex() < IIIend)
    //                                 d_nnz[dof1->getIndex() - IIIstart]++; // - IIIstart to get the local index
    //                             else
    //                                 o_nnz[dof1->getIndex() - IIIstart]++;

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

    PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE)); // or before solving the system
    PetscCall(MatSetOption(A, MAT_SPD, PETSC_TRUE));       // Positive definite matrix

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

    if (params->getCalculateReactionForces())
    {
        ierr = VecDuplicate(b, &nodalForces);
        CHKERRQ(ierr);
    }

    return ierr;
}