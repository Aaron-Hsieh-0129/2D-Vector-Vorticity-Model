#include "Declare.hpp"
#include <petsc.h>

// Config(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, std::string OUTPUTPATH, int OUTPUTSTEP
//        double Kx, double Kz, double TIMETS, double POISSONPARAMU, double POISSONPARAMW, double GRAVITY, double Cp, double Cv, double Rd, double Lv
//        double P0, double PSURF, double ADDFORCINGTIME)
Config_VVM config(1., 200., 200., 60000, 20000, 1500000., 10000, "/data/Aaron/TMIF/0527/vvm/", 10, 
              200, 200, 0.01, 1E-9, 0., 9.80665, 1003.5, 716.5, 287., 2.5E6, 
              1E5, 96500., 1200.);
vvm model(config);
int main(int argc, char **argv) {
    vvm::Init init;
    #if defined(LOADFROMPREVIOUSFILE)
        Init::LoadFromPreviousFile(model);
    #else
        vvm::Init::Init1d(model);
        vvm::Init::Init2d(model);
    #endif

    vvm::Output::printInit(model);
    vvm::Output::create_all_directory(model);

    PetscInitialize(&argc, &argv, NULL, NULL);
    #if defined(POISSONTEST)
        vvm::PoissonSolver PoissonSolver;
        PoissonSolver.cal_w(model);
        PoissonSolver.cal_u(model);
    #else
        vvm::Iteration::TimeMarching(model);
    #endif
    PetscFinalize();
    return 0;
}
