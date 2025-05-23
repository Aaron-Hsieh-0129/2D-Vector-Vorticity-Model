#include "Declare.hpp"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(PETSC)
    #include <petscsys.h>
    #include <petsc.h>
#endif
#include "ReadConfig.hpp"
#include "Timer.hpp"
#if defined(GPU_POISSON)
    #include <mpi.h>
#endif

// Config(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, std::string outputpath, int OUTPUTSTEP
//        double Kx, double Kz, double TIMETS, double POISSONPARAMU, double POISSONPARAMW, double GRAVITY, double Cp, double Cv, double Rd, double Lv
//        double P0, double PSURF, double ADDFORCINGTIME)

#if defined(P3_MICROPHY)
extern "C" {
    void __microphy_p3_MOD_p3_init(
        char* lookup_file_dir, int* nCat, bool* trplMomI, bool* liqfrac,
        char* model, int* stat, bool* abort_on_err, bool* dowr, size_t lookup_file_dir_len, size_t model_name_len
    );
}
#endif

int main(int argc, char **argv) {
    #if defined(GPU_POISSON)
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    #endif

    Timer timer_all;
    timer_all.reset();
    #if defined(PETSC)
        PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    #endif

    #ifdef _OPENMP
    omp_set_num_threads(8);
    Eigen::setNbThreads(8);
    #endif

    std::map<std::string, std::string> configs = vvm_read_config("../vvm_config.txt");
    std::string vvmoutputpath = configs["VVM_OUTPUTPATH"] + "/";
    double vvm_dt = std::stod(configs["VVM_DT"]);
    double vvm_dx = std::stod(configs["VVM_DX"]);
    double vvm_dz = std::stod(configs["VVM_DZ"]);
    int vvm_xrange = std::stoi(configs["VVM_XRANGE"]);
    int vvm_zrange = std::stoi(configs["VVM_ZRANGE"]);
    double vvm_timeend = std::stod(configs["VVM_TIMEEND"]);
    int vvm_outputstep = std::stoi(configs["VVM_OUTPUTSTEP"]);
    int vvm_case = std::stoi(configs["VVM_CASE"]);
    double vvm_moisture_nudge_time = std::stod(configs["VVM_MOISTURE_NUDGE_TIME"]);
    int vvm_year = std::stoi(configs["VVM_YEAR"]);
    int vvm_month = std::stoi(configs["VVM_MONTH"]);
    int vvm_day = std::stoi(configs["VVM_DAY"]);
    double vvm_hour = std::stod(configs["VVM_HOUR"]);
    double vvm_minute = std::stod(configs["VVM_MINUTE"]);
    double vvm_sec = std::stod(configs["VVM_SECOND"]);
    double vvm_lon = std::stod(configs["VVM_LON"]);
    double vvm_lat = std::stod(configs["VVM_LAT"]);

    Config_VVM config(vvm_dt, vvm_dx, vvm_dz, vvm_xrange, vvm_zrange, vvm_timeend, 10000, vvmoutputpath, vvm_outputstep, 
                    10., 10., 0.01, 1E-22, 9.80665, 1003.5, 716.5, 287., 2.5E6, 
                    1E5, 96500., 10, vvm_case, vvm_moisture_nudge_time, 
                    vvm_year, vvm_month, vvm_day, vvm_hour, vvm_minute, vvm_sec, 
                    vvm_lon, vvm_lat);
    
    vvm model(config);
    
    #if defined(LOADFROMPREVIOUSFILE)
        vvm::Init::LoadFromPreviousFile(model);
    #elif defined(LOAD2DINIT)
        vvm::Init::Load2DInit(model);
    #else
        vvm::Init::Init1d(model);
        vvm::Init::Init2d(model);
    #endif

    #if defined(GPU_POISSON)
        vvm::PoissonSolver::InitAMGX(model);
    #endif
    vvm::PoissonSolver::InitPoissonMatrix(model);

    vvm::Output::create_all_directory(model);
    vvm::Output::printInit(model);

    // p3 initialization
    #if defined(P3_MICROPHY)
    __microphy_p3_MOD_p3_init(
        vvm::P3::lookup_file_dir, &vvm::P3::nCat, &vvm::P3::trplMomI, &vvm::P3::liqfrac,
        vvm::P3::model_name, &vvm::P3::stat, &vvm::P3::abort_on_err, &vvm::P3::dowr, strlen(vvm::P3::lookup_file_dir), strlen(vvm::P3::model_name)
    );

    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) { 
            model.dz_all[i][k] = model.dz;
            model.w_all[i][k] = 0.;
            model.pb_all[i][k] = model.pb[k];
            model.zi_all[i][k] = 0.;
            model.ssat_all[i][k] = 0.;
        }
    }
    #endif

    // // This initialization is for NGAC3F coupling comparison
    // std::ifstream inputFile("/data/Aaron/TMIF/th0.txt");
    // for (int k = 1; k < model.nz-1; k++) {
    //     for (int i = 1; i < model.nx-1; i++) {
    //         inputFile >> model.th[i][k];
    //         model.thm[i][k] = model.th[i][k];
    //     }
    // }

    // Copy grads ctl file to the output directory
    std::string source = "../scripts/vvm.ctl";
    std::string destination = model.outputpath + "nc/vvm.ctl";

    // Construct the command
    std::string command = "cp " + source + " " + destination;

    // Execute the command
    system(command.c_str());
    

    #if defined(POISSONTEST)
        vvm::PoissonSolver::cal_w(model);
        vvm::PoissonSolver::cal_u(model);
    #else
        vvm::Iteration::TimeMarching(model);
    #endif

    #if defined(PETSC)
        PetscCall(PetscFinalize());
    #endif
    #if defined(GPU_POISSON)
        vvm::PoissonSolver::CleanupAMGX(model);
        std::cout << "Total time: " << timer_all.elapsed() << std::endl;
        MPI_Finalize();
    #endif
    return 0;
}
