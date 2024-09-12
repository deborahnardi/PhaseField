#include <petsc.h> // Include PETSc library

#include "../../enumclass.hpp"

/*
AppCtx stands for Application Context. It is a class that holds the application context, which is a set of data that is used by the application to perform its tasks. It is a way to encapsulate the data that is used by the application, so that it can be passed around easily and accessed by different parts of the application.
*/
class AppCtx
{
private:
    std::string dmType;
    DeformType deform;
    PetscBag bag;
    PetscBool useNearNullSpace;

public:
    AppCtx();
    ~AppCtx();

    const char *deformTypes[NUM_DEFORM_TYPES + 1] = {"none", "shear", "step", "unknown"};
    const char *solutionTypes[NUM_SOLUTION_TYPES + 1] = {"vlap_quad", "elas_quad", "vlap_trig", "elas_trig", "elas_axial_disp", "elas_uniform_strain", "elas_ge", "mass_quad", "unknown"};
};

// static means that the function is only visible in the file where it is defined. It is a way to make the function private to the file.
// The linker will not link this function with other files, so it is not visible outside of this file.

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx &options)
{
    PetscInt sol = 0, def = 0;
}