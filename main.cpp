static char help[] = "Solves the phase-field problem using Finete Element Method "
                     "Start date: July 22nd, 2024 "
                     "Written by: Deborah Nardi"
                     "Advisor: Prof. Edson Denner Leonel"
                     "Co-advisor: Prof. Ayrton Ferreira"
                     "São Carlos School of Engineering - University of São Paulo";

#include "mesh_interface/headers/Geometry.h"
#include "solid/headers/FEM.h"
#include "solid/headers/Quadrature.h"
#include "solid/headers/ShapeFunction.h"
#include "petsc/headers/PETScExs.h"
#include "solid/headers/DenseEigen.h"

#include <omp.h>

// class MyClass
// {
// public:
//   void processData(int numThreads)
//   {
//     // Define o número de threads para OpenMP
//     omp_set_num_threads(numThreads);

// // Simula um loop de processamento de dados
// #pragma omp parallel for
//     for (int i = 0; i < 10; ++i)
//     {
//       int thread_id = omp_get_thread_num();
//       std::cout << "Thread " << thread_id << " is processing item " << i << std::endl;
//     }
//   }
// };

int main(int argc, char **argv)
{

  PetscInitialize(&argc, &argv, (char *)0, help); // Starts main program invoking PETSc

  //================================ TESTING OPENMP ================================

  // std::cout << "OpenMP is initialized (MAIN) with " << omp_get_max_threads() << " threads." << std::endl;

  //   int rank, size;
  //   MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Obtém o rank do processo MPI
  //   MPI_Comm_size(MPI_COMM_WORLD, &size); // Obtém o número total de processos MPI

  // // Definir número de threads pelo OpenMP
  // #pragma omp parallel
  //   {
  //     int thread_id = omp_get_thread_num();
  //     int num_threads = omp_get_num_threads();

  // // Usamos critical para evitar misturas na saída
  // #pragma omp critical
  //     {
  //       std::cout << "MPI Rank " << rank << " - Hello from thread "
  //                 << thread_id << " out of " << num_threads << " threads."
  //                 << std::endl;
  //     }
  //   }
  // ======================================================================================

  // DISCOMMENT THE ABOVE LINE BEFORE RUNNING THE NON PETSC EXAMPLES
  // #include "examples/pointerAndReference.hpp"
  // #include "examples/Ex01Inclusions.hpp"
  // #include "examples/square.hpp"
  //    #include "examples/squareEllipse.hpp"
  //           #include "examples/Ex02NumericalIntegration.hpp"
  //   #include "examples/Ex03Truss.hpp"
  //      #include "examples/Ex04Truss.hpp"

  // ================ Phase Field Examples =================
  // #include "examples/phaseField1D.hpp"
  //       #include "examples/phaseField2D-01.hpp"
  //      #include "examples/phaseField2D-02.hpp"
#include "examples/phaseField2D-03.hpp"
  //        #include "examples/phaseField2D-04.hpp"
  // #include "examples/phaseField2D-04copy.hpp"
  //           #include "examples/FEMHardWay.hpp"

  PetscFinalize(); // Finalize main program
}