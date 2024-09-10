#include <petsc.h> // Include PETSc library

/*
AppCtx stands for Application Context. It is a class that holds the application context, which is a set of data that is used by the application to perform its tasks. It is a way to encapsulate the data that is used by the application, so that it can be passed around easily and accessed by different parts of the application.
*/
class AppCtx
{
public:
    AppCtx();
    ~AppCtx();
};

// static means that the function is only visible in the file where it is defined. It is a way to make the function private to the file.
static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx &options)
{
    PetscInt sol = 0, def = 0;
}