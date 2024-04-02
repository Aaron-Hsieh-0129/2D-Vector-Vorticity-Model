#include "Iteration.hpp"
#include <petsc.h>

vvm model;
int main(int argc, char **argv) {
    Init::Init1d(model);
    Init::Init2d(model);
    Output::printInit(model);
    Output::create_all_directory();

    PetscInitialize(&argc, &argv, NULL, NULL);
    #if defined(POISSONTEST)
        vvm::PoissonSolver PoissonSolver;
        PoissonSolver.cal_w(model);
        PoissonSolver.cal_u(model);
    #else
        Iteration::TimeMarching(model);
    #endif
    PetscFinalize();
    return 0;
}
