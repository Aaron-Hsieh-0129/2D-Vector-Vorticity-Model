#include "Iteration.hpp"
#include <petsc.h>

vvm model;
int main(int argc, char **argv) {
    Init::Init1d(model);
    Init::Init2d(model);
    Output::printInit(model);
    Output::create_all_directory();

    PetscInitialize(&argc, &argv, NULL, NULL);
    Iteration::TimeMarching(model);
    PetscFinalize();
    return 0;
}
