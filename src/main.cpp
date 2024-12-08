#include "Declare.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(PETSC)
    #include <petscsys.h>
    #include <petsc.h>
#endif
#include <fstream>

// Config(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, std::string outputpath, int OUTPUTSTEP
//        double Kx, double Kz, double TIMETS, double POISSONPARAMU, double POISSONPARAMW, double GRAVITY, double Cp, double Cv, double Rd, double Lv
//        double P0, double PSURF, double ADDFORCINGTIME)

int main(int argc, char **argv) {
    #if defined(PETSC)
        PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    #endif

    #ifdef _OPENMP
    omp_set_num_threads(8);
    Eigen::setNbThreads(8);
    #endif

    Config_VVM config(3., 200., 200., 100000, 20000, 90000., 10000, "/data/Aaron/2DVVM/Bubble_shear_eva/", 50, 
                    70., 70., 0.01, 1E-22, 9.80665, 1003.5, 716.5, 287., 2.5E6, 
                    1E5, 96500., -1., 2);
    vvm model(config);
    
    #if defined(LOADFROMPREVIOUSFILE)
        vvm::Init::LoadFromPreviousFile(model);
    #elif defined(LOAD2DINIT)
        vvm::Init::Load2DInit(model);
    #else
        vvm::Init::Init1d(model);
        vvm::Init::Init2d(model);
    #endif

    #ifndef PETSC
        vvm::PoissonSolver::InitPoissonMatrix(model);
    #endif

    vvm::Output::printInit(model);
    vvm::Output::create_all_directory(model);

    // // This initialization is for NGAC3F coupling comparison
    // std::ifstream inputFile("/data/Aaron/TMIF/th0.txt");
    // for (int k = 1; k < model.nz-1; k++) {
    //     for (int i = 1; i < model.nx-1; i++) {
    //         inputFile >> model.th[i][k];
    //         model.thm[i][k] = model.th[i][k];
    //     }
    // }
    

    #if defined(POISSONTEST)
        vvm::PoissonSolver::cal_w(model);
        vvm::PoissonSolver::cal_u(model);
    #else
        vvm::Iteration::TimeMarching(model);
    #endif

    #if defined(PETSC)
        PetscCall(PetscFinalize());
    #endif
    return 0;
}
