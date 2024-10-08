#include "../headers/FEM.h"

FEM::FEM() {}
FEM::FEM(const std::string _name)
    : name(_name)
{
    setResultsPath();
    deleteResults(true);
    createResultsPath();
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
};
FEM::~FEM() {}

/*----------------------------------------------------------------------------------
                            DATA INPUT METHODS
------------------------------------------------------------------------------------*/

void FEM::createResultsPath()
{
    if (rank == 0)
    {
        std::string command = "mkdir -p " + resultsPath + "results/hdf5";
        command += " && mkdir -p examples/";
        system(command.c_str());
    }
}

void FEM::deleteResults(bool deleteFiles)
{
    if (rank == 0)
    {
        if (deleteFiles)
        {
            std::cout << "Deleting files\n";
            std::string command = "rm -r " + resultsPath + "*";
            std::cout << command << "\n";
            system(command.c_str());
        }
    }
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

/*----------------------------------------------------------------------------------
                Assembling and solving problem PETSc
----------------------------------------------------------------------------------
*/
PetscErrorCode FEM::solveFEMProblem()
{
    findNeighbours();
    assembleProblem();
    solveLinearSystem(matrix, rhs, solution);
    if (rank == 0)
        showResults();

    return 0;
}

PetscErrorCode FEM::assembleProblem()
{
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Assembling problem...\n");

    decomposeElements(rhs, solution);
    matrixPreAllocation();
    createPETScVariables(matrix, rhs, solution, nDOFs, true);

    ierr = MatZeroEntries(matrix);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(rhs);
    CHKERRQ(ierr);
    ierr = VecZeroEntries(solution);
    CHKERRQ(ierr);

    for (int Ii = Istart; Ii < Iend; Ii++)
        elements[Ii]->getContribution(matrix, rhs);

    for (int Ii = IIstart; Ii < IIend; Ii++) // Neumann boundary conditions
        bdElements[Ii]->getContribution(rhs);

    // Assemble the matrix and the right-hand side vector
    ierr = VecAssemblyBegin(rhs);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);

    ierr = MatZeroRowsColumns(matrix, numDirichletDOFs, dirichletBC, 1., solution, rhs); // Apply Dirichlet boundary conditions
    CHKERRQ(ierr);

    // if (showMatrix && rank == 0) // Print the global stiffness matrix on the terminal
    // {
    //     ierr = PetscPrintf(PETSC_COMM_WORLD, " --- GLOBAL STIFFNESS MATRIX: ----\n");
    //     CHKERRQ(ierr);
    //     ierr = MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
    //     CHKERRQ(ierr);
    //     // printGlobalMatrix(matrix);
    // }

    return ierr;
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
    ierr = VecSetSizes(x, PETSC_DECIDE, globalDOFs.size());
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(x, &IIIstart, &IIIend);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    // ----------------------------------------------------------------
    return ierr;
}

PetscErrorCode FEM::matrixPreAllocation()
{
    int rankLocalDOFs = IIIend - IIIstart; // Number of nodes in the local partition

    d_nnz = new PetscInt[rankLocalDOFs]();
    o_nnz = new PetscInt[rankLocalDOFs]();

    for (auto node1 : nodes)
        for (auto dof1 : node1->getDOFs()) // Rows of the matrix
            if (dof1->getIndex() >= IIIstart && dof1->getIndex() < IIIend)
                for (auto node2 : nodeNeighbours[node1->getIndex()]) // Columns of the matrix
                    for (auto dof2 : nodes[node2]->getDOFs())
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

    return ierr;
}

PetscErrorCode FEM::printGlobalMatrix(Mat &A)
{
    PetscInt i, j, rows, cols;
    PetscScalar value;
    const int width = 10; // Columns width

    ierr = MatGetSize(A, &rows, &cols);
    CHKERRQ(ierr);

    std::cout << "Global stiffness matrix:" << std::endl;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            ierr = MatGetValue(A, i, j, &value);
            CHKERRQ(ierr);
            // std::cout << std::setw(width) << std::fixed << std::setprecision(0) << value << " ";
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    return ierr;
}

PetscErrorCode FEM::solveLinearSystem(Mat &A, Vec &b, Vec &x)
{
    KSP ksp;
    PC pc;
    PetscInt its;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    switch (size)
    {
    case 1:
        ierr = KSPGetPC(ksp, &pc);
        CHKERRQ(ierr);
        ierr = PCSetType(pc, PCJACOBI);
        CHKERRQ(ierr);
        break;
    }

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &its); // Gets the number of iterations
    CHKERRQ(ierr);

    ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // Prints the Krylov subspace method information
    CHKERRQ(ierr);

    ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); // Prints the solution vector
    CHKERRQ(ierr);

    /*----------------------------------------------------------------------------------
                                        POST-PROCESSING
    ------------------------------------------------------------------------------------
    */
    // Set the solution to the final coordinates of the nodes
    finalDisplacements = new double[globalDOFs.size()];
    int rankLocalDOFs = IIIend - IIIstart;
    double *localDisplacements;
    localDisplacements = new double[rankLocalDOFs];

    for (int i = IIIstart; i < IIIend; i++)
    {
        DOF *dof = globalDOFs[i];
        PetscScalar retrieved;
        PetscInt index = dof->getIndex();
        VecGetValues(x, 1, &index, &retrieved);

        localDisplacements[i - IIIstart] = retrieved;
    }

    // Transfering data between processes
    int *numEachProcess = new int[size];
    // Collecting the number of local DOFs from all processes
    MPI_Allgather(&rankLocalDOFs, 1, MPI_INT, numEachProcess, 1, MPI_INT, PETSC_COMM_WORLD);

    int *globalBuffer = new int[size];
    globalBuffer[0] = 0;
    for (int i = 1; i < size; i++)
        globalBuffer[i] = globalBuffer[i - 1] + numEachProcess[i - 1];

    // Collecting the displacements from all processes to the root process
    MPI_Gatherv(localDisplacements, rankLocalDOFs, MPI_DOUBLE, finalDisplacements, numEachProcess, globalBuffer, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    // Broadcasting the final displacements to all processes
    MPI_Bcast(finalDisplacements, nDOFs, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

    for (auto node : nodes)
    {
        double valueX = (!node->getDOF(0)->isDirichlet()) ? finalDisplacements[node->getDOF(0)->getIndex()] : node->getDOF(0)->getDirichletValue();
        double valueY = (!node->getDOF(1)->isDirichlet()) ? finalDisplacements[node->getDOF(1)->getIndex()] : node->getDOF(1)->getDirichletValue();

        node->setFinalDisplacement({valueX, valueY, 0.});
    }

    for (int i = 0; i < size; i++)
    {
        if (i == rank)
        {
            std::cout << "Final displacements: " << std::endl;
            for (auto node : nodes)
            {
                std::cout << node->getFinalDisplacement()[0] << std::endl;
                std::cout << node->getFinalDisplacement()[1] << std::endl;
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
    /*----------------------------------------------------------------------------------
                                        CLEANING UP
    ------------------------------------------------------------------------------------
    */

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);

    return ierr;
}