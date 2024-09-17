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