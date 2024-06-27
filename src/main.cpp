#include "Declare.hpp"
#include <omp.h>
#if defined(PETSC)
    #include <petscsys.h>
    #include <petsc.h>
#endif

// Config(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, std::string outputpath, int OUTPUTSTEP
//        double Kx, double Kz, double TIMETS, double POISSONPARAMU, double POISSONPARAMW, double GRAVITY, double Cp, double Cv, double Rd, double Lv
//        double P0, double PSURF, double ADDFORCINGTIME)

int main(int argc, char **argv) {
    #if defined(PETSC)
        PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    #endif

    Config_VVM config(4., 200., 200., 100000, 20000, 40000., 10000, "/data/Aaron/TMIF/0619_test_Bubble/", 1, 
                    200., 200., 0.01, 1E-9, 0., 1E-20, 9.80665, 1003.5, 716.5, 287., 2.5E6, 
                    1E5, 96500., 10., 1);
    vvm model(config);
    
    #if defined(LOADFROMPREVIOUSFILE)
        vvm::Init::LoadFromPreviousFile(model);
    #elif defined(LOAD2DINIT)
        vvm::Init::Load2DInit(model);
    #else
        vvm::Init::Init1d(model);
        vvm::Init::Init2d(model);
    #endif

    vvm::Output::printInit(model);
    vvm::Output::create_all_directory(model);

    #if defined(POISSONTEST)
        vvm::PoissonSolver PoissonSolver;
        PoissonSolver.cal_w(model);
        PoissonSolver.cal_u(model);
    #else
        vvm::Iteration::TimeMarching(model);
    #endif

    #if defined(PETSC)
        PetscCall(PetscFinalize());
    #endif
    return 0;
}
