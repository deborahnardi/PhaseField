#include "../headers/PETScExs.h"

PETScExs::PETScExs() {}
PETScExs::~PETScExs() {}

PetscErrorCode PETScExs::PETScSequentialTest()
{
    /*
    A preconditioner is a matrix that approximates the inverse of the matrix that we need to solve.
    */
    Vec x, b, u;                     // Approx solution, rhs, exact solution
    Mat A;                           // Linear system matrix
    KSP ksp;                         // Defines the Krylov subspace method
    PC pc;                           // Preconditioner context
    PetscReal norm;                  // Norm of solution error
    PetscErrorCode ierr;             // Error code
    PetscInt i, n = 10, col[3], its; // Number of rows and columns, its -> Number of iterations
    PetscMPIInt size;                // Number of processes
    PetscScalar value[3];            // Values of the matrix

    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); // Get number of processes
    CHKERRQ(ierr);

    if (size != 1)
        SETERRQ(PETSC_COMM_WORLD, 1, "ATTENTION: This is a uniprocessor example only!");

    ierr = PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL); // This is a PETSc function that allows the user to set the number of rows and columns of the matrix via the command line -> THIS IS NOT NECESSARY IF NO COMMAND LINE ARGUMENTS ARE PASSED
    CHKERRQ(ierr);

    /*------------------------------------------------------------------
        Compute the matrix and righ-hand-side vector that define
        the linear system, Ax = b.
    --------------------------------------------------------------------
    */

    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x, "Solution"); // Sets the name of the vector x
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, n); // Sets the size of the vector x
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); // Good practice
    CHKERRQ(ierr);
    ierr = VecDuplicate(x, &b); // Creates a vector b that is a duplicate of x
    CHKERRQ(ierr);
    ierr = VecDuplicate(x, &u); // Creates a vector u that is a duplicate of x

    /*
        Create matrix. When using MatCreate(), the matrix format can be specified at runtime.

        Perfomance tuning note: For problems of substantial size, preallocation of matrix memory is crucial for attaining good performance.
    */

    ierr = MatCreate(PETSC_COMM_WORLD, &A); // Creates a matrix A
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n); // Sets the size of the matrix A
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); // Good practice
    CHKERRQ(ierr);
    ierr = MatSetUp(A); // Pre allocation of matrix memory
    CHKERRQ(ierr);

    /*
        Assemble matrix
    */

    value[0] = -1.0;
    value[1] = 2.0;
    value[2] = -1.0;

    for (int i = 1; i < n - 1; i++)
    {
        col[0] = i - 1;                                              // Left column of position i
        col[1] = i;                                                  // Central column of position i
        col[2] = i + 1;                                              // Right column of position i
        ierr = MatSetValues(A, 1, &i, 3, col, value, INSERT_VALUES); // Arguments: matrix, number of rows, line index, number of columns, column index, values to be inserted into the matrix, INSERT_VALUES
        CHKERRQ(ierr);
    }

    i = n - 1;
    col[0] = n - 2;
    col[1] = n - 1;

    ierr = MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES); // It has to be &i because i is a pointer to the memory address of line i
    CHKERRQ(ierr);

    i = 0;
    col[0] = 0;
    col[1] = 1;
    value[0] = 2.0;
    value[1] = -1.0;

    ierr = MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES);
    CHKERRQ(ierr);

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    /*
        Set exact solution; then compute right-hand side vector.
    */

    ierr = VecSet(u, 1.0); // All elements of the vector u are set to 1.0
    CHKERRQ(ierr);
    ierr = MatMult(A, u, b); // Multiplies the matrix A by the vector u and stores the result in the vector b
    CHKERRQ(ierr);

    /*
        Create the linear solver and set various options
    */

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); // Creates a Krylov subspace method
    CHKERRQ(ierr);

    ierr = KSPSetOperators(ksp, A, A); // ksp, A: matrix, A: preconditioner (preconditioner is the same as the matrix)
    CHKERRQ(ierr);

    /*
        The following lines are optional. They are used to set the type of the Krylov subspace method and the type of the preconditioner.
        If these lines are not used, PETSc will use the default values, which is GMRS (Generalized Minimal Residual) for the Krylov subspace method and Jacobi for the preconditioner.
    */

    ierr = KSPGetPC(ksp, &pc); // Gets the preconditioner context
    CHKERRQ(ierr);

    ierr = PCSetType(pc, PCJACOBI); // Sets the type of the preconditioner
    CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp); // This override the previous specifications
    CHKERRQ(ierr);

    /*
    --------------------------------------------------------------------
                            Solve the linear system
    --------------------------------------------------------------------
    */

    ierr = KSPSolve(ksp, b, x); // b is the right-hand side vector, x is the solution vector
    CHKERRQ(ierr);

    ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // Prints the Krylov subspace method information
    CHKERRQ(ierr);

    ierr = VecView(x, PETSC_VIEWER_STDOUT_WORLD); // Prints the solution vector
    CHKERRQ(ierr);

    /*
        Check the solution and clean up
    */

    ierr = VecAXPY(x, -1.0, u); // x = x - u, this defines the error, where u is the exact solution
    CHKERRQ(ierr);

    ierr = VecNorm(x, NORM_2, &norm); // Computes the norm of the vector x
    CHKERRQ(ierr);

    ierr = KSPGetIterationNumber(ksp, &its); // Gets the number of iterations
    CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g, Iterations %D\n", (double)norm, its); // Prints the norm of the error and the number of iterations

    /*
        Free work space. All PETSc objects should be destroyed when they are no longer needed.
    */

    ierr = VecDestroy(&x); // PETSc expects the address of the pointer to the object
    CHKERRQ(ierr);
    ierr = VecDestroy(&u);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
}

PetscErrorCode PETScExs::PETScParallelTest()
{
    Vec x, b, u;    // Approx solution, rhs, exact solution
    Mat A;          // Linear system matrix
    KSP ksp;        // Defines the Krylov subspace method
    PetscReal norm; // Norm of solution error
    PetscInt i, j, Ii, J, Istart, Iend, m = 8, n = 7, its;
    PetscErrorCode ierr;
    PetscBool flg;
    PetscScalar v;

    /*-------------------------------------------------------------------------------------
        Create the matrix and right-hand-side vector that define the linear system, Ax = b.
    ---------------------------------------------------------------------------------------*/

    ierr = MatCreate(PETSC_COMM_WORLD, &A);
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m * n, m * n);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    /*
        MatMPIAIJSetPreallocation(Mat mat, PetscInt d_nz, const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[]);

        mat: matrix;
        d_nz: number of non-zero (NZ) elements per row in the diagonal portion of the matrix; you can use PETSC_DECIDE to have PETSc determine the number of non-zeros;
        d_nnz: array containing the number of non-zeros in the various rows of the diagonal portion of the matrix; if you don't know the number of non-zeros, you can use NULL;
        o_nz: number of non-zero elements per row in the off-diagonal portion of the matrix; you can use PETSC_DECIDE to have PETSc determine the number of non-zeros;
        o_nnz: array containing the number of non-zeros in the various rows of the off-diagonal portion of the matrix; if you don't know the number of non-zeros, you can use NULL.
    */
    ierr = MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL); // Parallel version of the matrix
    CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A, 5, NULL); // Sequential version of the matrix
    CHKERRQ(ierr);
    ierr = MatSeqSBAIJSetPreallocation(A, 1, 5, NULL); // Symmetric block version of the matrix (Sequential)
    CHKERRQ(ierr);
    ierr = MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL); // Symmetric block version of the matrix (Parallel)
    CHKERRQ(ierr);
    /*
        MatMPISELLSetPreallocation: MPISELL means MPI Sparse Element-wise Linked List. This format is used for matrices that are not symmetric and have a general sparsity pattern.
    */
    ierr = MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
    CHKERRQ(ierr);
    ierr = MatSeqSELLSetPreallocation(A, 5, NULL);
    CHKERRQ(ierr);

    /*
        Determine wich rows of the matrix are locally owned.
    */
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);

    /*
        Set matrix elements for the 2-D, 5-point stencil in parallel.
    */

    for (Ii = Istart; Ii < Iend; Ii++)
    {
        v = -1.0;
        i = Ii / n;
        j = Ii - i * n;
        if (i > 0)
        {
            J = Ii - n;
            ierr =
                MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
            CHKERRQ(ierr);
        }
        if (i < m - 1)
        {
            J = Ii + n;
            ierr =
                MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
            CHKERRQ(ierr);
        }
        if (j > 0)
        {
            J = Ii - 1;
            ierr =
                MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
            CHKERRQ(ierr);
        }
        if (j < n - 1)
        {
            J = Ii + 1;
            ierr =
                MatSetValues(A, 1, &Ii, 1, &J, &v, ADD_VALUES);
            CHKERRQ(ierr);
        }
        v = 4.0;
        ierr = MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES);
        CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
    CHKERRQ(ierr);

    /*----------------------------------------------------------------
                        CREATE PARALLEL VECTORS
    ------------------------------------------------------------------
    */

    ierr = VecCreate(PETSC_COMM_WORLD, &u);
    CHKERRQ(ierr);
    ierr = VecSetSizes(u, PETSC_DECIDE, m * n);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);
    CHKERRQ(ierr);
    ierr = VecDuplicate(u, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(u, &b);
    CHKERRQ(ierr);

    /*--------------------------------------------------------------
                        Set exact solution
    ----------------------------------------------------------------
    */

    ierr = VecSet(u, 1.0);
    CHKERRQ(ierr);
    ierr = MatMult(A, u, b);

    /*
        View the exact solution vector
    */
    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRQ(ierr);

    /*---------------------------------------------------------------
            Create the linear solver and set various options
    -----------------------------------------------------------------
    */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1.e-2 / ((m + 1) * (n + 1)), 1e-50, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);

    /*---------------------------------------------------------------
                        Solve the linear system
    -----------------------------------------------------------------
    */
    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);

    /*----------------------------------------------------------------
                    Check the solution and clean up
    ------------------------------------------------------------------
    */
    ierr = VecAXPY(x, -1.0, u); // x = x - u to get the error
    CHKERRQ(ierr);
    ierr = VecNorm(x, NORM_2, &norm);
    CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp, &its);
    CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g iterations%D\n", double(norm), its);
    CHKERRQ(ierr);

    /*
        Free work space. All PETSc objects should be destroyed when they are no longer needed.
    */
    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&u);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
}