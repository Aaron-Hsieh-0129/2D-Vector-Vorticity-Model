#include "Config.hpp"
#include <string>
#include <iostream>
#ifndef PETSC
    // #include "../include/Eigen/Sparse"
    #include <Eigen/Sparse>
#endif
#if defined(GPU_POISSON)
    #include <amgx_c.h>
    #include <cuda_runtime.h>
    #include <mpi.h>
    // #include <cusparse.h>
#endif

#if defined(RTERRTMGP)
#include <boost/algorithm/string.hpp>
#include <chrono>
#include <iomanip>

#include "../external/rte-rrtmgp-cpp/include_test/Status.h"
#include "../external/rte-rrtmgp-cpp/include_test/Netcdf_interface.h"
#include "../external/rte-rrtmgp-cpp/include/Array.h"
#include "../external/rte-rrtmgp-cpp/include/Aerosol_optics.h"
#include "../external/rte-rrtmgp-cpp/include_test/Radiation_solver.h"
#include "../external/rte-rrtmgp-cpp/include/types.h"
#endif

class Config_VVM {
public:
    Config_VVM(double dt, double dx, double dz, double dz1, int nz, int XRANGE, double TIMEEND, int TIMEROUTPUTSIZE, 
           std::string outputpath, int OUTPUTSTEP, double Kx, double Kz, double tolerance,
           double GRAVITY, double Cp, double Cv, double Rd, double Lv, double P0, double PSURF, double addforcingtime, int CASE, double mositure_nudge_time, 
           int year, int month, int day, double hour, double minute, double second, double lon, double lat)
        : dt(dt), dx(dx), dz(dz), dz1(dz1), nz(nz+2), XRANGE(XRANGE+2*dx), TIMEEND(TIMEEND), TIMEROUTPUTSIZE(TIMEROUTPUTSIZE), 
          outputpath(outputpath), OUTPUTSTEP(OUTPUTSTEP), Kx(Kx), Kz(Kz),
          tolerance(tolerance), GRAVITY(GRAVITY), Cp(Cp), Cv(Cv), Rd(Rd), Lv(Lv), P0(P0), PSURF(PSURF), addforcingtime(addforcingtime), CASE(CASE), mositure_nudge_time(mositure_nudge_time), 
          year(year), month(month), day(day), hour(hour), minute(minute), second(second), lon(lon), lat(lat) {}
    ~Config_VVM() {}

    double dt;              ///< Time step for vvm [s].
    double dx;              ///< Grid size in x-direction [m].
    double dz;              ///< z stretch coefficient1
    double dz1;              ///< z stretch coefficient2
    int nz;                 ///< z layer number
    int XRANGE;             ///< Domain size of the model in x-direction [m].
    double TIMEEND;         ///< End time of the simulation [s].
    int TIMEROUTPUTSIZE;    ///< The size of the timer output.
    std::string outputpath; ///< The path for the output file. It should be a directory, such as "/data/vvm/".
    int OUTPUTSTEP;         ///< The output interval for the output file.
    double Kx;              ///< The eddy diffusion coefficient in x-direction [m^2/s], this is activated when DIFFUSION flag is turned on. If the flag is not turned on, the coeffcient will be calculated through the turbulent closure.
    double Kz;              ///< The eddy diffusion coefficient in z-direction [m^2/s], this is activated when DIFFUSION flag is turned on. If the flag is not turned on, the coeffcient will be calculated through the turbulent closure.
    double tolerance;       ///< The tolerance for the Poisson Solver.
    double GRAVITY;         ///< The gravity acceleration [m/s^2]. It's 9.80665 m/s^2 for default.
    double Cp;              ///< The specific heat capacity at constant pressure [J/kg/K]. It's 1003.5 J/kg/K for default.
    double Cv;              ///< The specific heat capacity at constant volume [J/kg/K]. It's 716.5 J/kg/K for default.
    double Rd;              ///< The gas constant for dry air [J/kg/K]. It's 287 J/kg/K for default.
    double Lv;              ///< The latent heat of vaporization [J/kg]. It's 2.5E6 J/kg for default.
    double P0;              ///< The reference pressure [Pa]. It's 1E5 Pa for default.
    double PSURF;           ///< The surface pressure [Pa]. It's 96500 Pa for default.
    double addforcingtime;  ///< The time for adding the perturbation. The perturbation is used to break the symmetry of the model.
    int CASE;               ///< The case number for the model. It's used to specify the initial condition and the forcing. If CASE is 0, the initial condition is equal to mean state. If CASE is 1, the initial condition is equal to mean state plus a warm bubble. 
    double mositure_nudge_time; ///< The time for nudging the moisture. It's used for nudging the moisture to the mean state. It's used for nudging the moisture to the mean state.
    int year;
    int month;
    int day;
    double hour;
    double minute;
    double second;
    double lon;
    double lat;
};


class vvm {
public:
    /**
     * vvm constructor.
     * Used to initialize the model.
     */
    vvm(const Config_VVM& config)
        : rdx(1.0 / config.dx), r2dx(rdx / 2.0), rdz(1.0 / config.dz), 
          r2dz(rdz / 2.0), rdx2(rdx * rdx),
          rdz2(rdz * rdz), nx(config.XRANGE/config.dx), 
          dt(config.dt), d2t(2.0 * config.dt), 
          dx(config.dx), dz(config.dz), dz1(config.dz1), nz(config.nz),
          XRANGE(config.XRANGE), TIMEEND(config.TIMEEND),
          TIMEROUTPUTSIZE(config.TIMEROUTPUTSIZE), outputpath(config.outputpath), OUTPUTSTEP(config.OUTPUTSTEP), Kx(config.Kx), Kz(config.Kz),
          tolerance(config.tolerance),
          GRAVITY(config.GRAVITY),
          Cp(config.Cp), Cv(config.Cv),
          Rd(config.Rd), Lv(config.Lv),
          P0(config.P0), PSURF(config.PSURF), addforcingtime(config.addforcingtime), CASE(config.CASE), moisture_nudge_time(config.mositure_nudge_time), 
          year(config.year), month(config.month), day(config.day), hour(config.hour), minute(config.minute), second(config.second), 
          lon(config.lon), lat(config.lat)
    {
        allocateMemory();
    }
    ~vvm() {
        printf("Free vvm\n");
        deallocateMemory();
    }

    // public members
    std::string outputpath = "";                  ///< From Config_VVM given by users.


    // ========= HELPER METHODS =========
    // 1D Allocation
    template<typename T>
    void create_variable(T*& array, size_t size) {
        array = new T[size]();
    }

    // 2D Allocation
    template<typename T>
    void create_variable(T**& array, T*& cont_memory, size_t dim1, size_t dim2) {
        cont_memory = new T[dim1 * dim2]();
        array = new T*[dim1]();
        for (size_t i = 0; i < dim1; ++i) {
            array[i] = &cont_memory[i * dim2];
        }
    }

    // 3D Allocation
    template<typename T>
    void create_variable(T***& array, T*& cont_memory, size_t dim1, size_t dim2, size_t dim3) {
        cont_memory = new T[dim1 * dim2 * dim3]();
        array = new T**[dim1]();
        for (size_t i = 0; i < dim1; ++i) {
            array[i] = new T*[dim2]();
            for (size_t j = 0; j < dim2; ++j) {
                array[i][j] = &cont_memory[i * (dim2 * dim3) + j * dim3];
            }
        }
    }

    // 1D Deallocation
    template<typename T>
    void destroy_variable(T*& array) {
        delete[] array;
        array = nullptr; // Good practice to prevent dangling pointers
    }

    // 2D Deallocation
    template<typename T>
    void destroy_variable(T**& array, T*& cont_memory) {
        if (array != nullptr) {
            delete[] cont_memory;
            delete[] array;
            array = nullptr;
            cont_memory = nullptr;
        }
    }

    // 3D Deallocation
    template<typename T>
    void destroy_variable(T***& array, T*& cont_memory, size_t dim1) {
        if (array != nullptr) {
            delete[] cont_memory;
            for (size_t i = 0; i < dim1; ++i) {
                delete[] array[i];
            }
            delete[] array;
            array = nullptr;
            cont_memory = nullptr;
        }
    }
    
    // Boundary Process => BoundaryProcess.cpp
    // **********************************************************************
    /**
     * A member function that process the boundary of the 1D array where the varibles are at the center of the grid
     * @param var an one dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess1D_center(double *var, int nz);

    /**
     * A member function that process the boundary of the 2D array where the varibles are at the center of the grid
     * @param var an two dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess2D_center(double **var, int nx, int nz);

    /**
     * A member function that process the boundary of the 2D array where the varibles are at the southwestern side of the grid
     * @param var an two dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess2D_westdown(double **var, int nx, int nz);

    static void BoundaryProcess2D_all(vvm &model);
    // **********************************************************************


    // Advection Scheme => Advection.cpp
    // *********************************************************************************
    /**
     * A member function that do advection process to the vorticity (zeta) field.
     * @param model the vvm object which is used to advect the vorticity field and put into it.
     */
    static void Advection_zeta(vvm &model);

    /**
     * A member function that do advection process to the vorticity (zeta) field.
     * @param previous an two dimensional array that the timestep is the previous one such as zetam, thm.
     * @param now an two dimensional array that the timestep is now such as zeta, th.
     * @param future an two dimensional array that the timestep is the future one such as zetap, thp.
     * @param model the vvm object.
     */
    // static void Advection_thermo(double **previous, double **now, double **future, vvm &model);
    static void Advection_thermo(double **past, double **now, double **future, double ***dvar, vvm &model);

    #if defined(WATER) && defined(KESSLER_MICROPHY)
        static void Advection_qrVT(vvm &model);
    #endif
    // *********************************************************************************

    // Buoyancy term (including heating by microphysics and radiation) => Buoyancy.cpp
    // **********************************************************************
    /**
     * A member function that do buoyancy process to the vorticity (zeta) field and heating by microphysics and radiation.
     * @param model the vvm object.
     */
    static void Buoyancy(vvm &model);
    // **********************************************************************


    // Surface flux heating and moistening term => SurfaceFlux.cpp
    // **********************************************************************
    /**
     * A member function that do surface flux process to the vorticity (zeta) field and heating by microphysics and radiation.
     * @param model the vvm object.
     */
    #if defined(SFCFLX)
        static void SurfaceFlux(vvm &model);
    #endif
    // **********************************************************************

    double getScalar(std::string name) {
        if (name == "rdx") return rdx;
        else if (name == "rdz") return rdz;
        else if (name == "r2dx") return r2dx;
        else if (name == "r2dz") return r2dz;
        else if (name == "rdx2") return rdx2;
        else if (name == "rdz2") return rdz2;
        else if (name == "dt") return dt;
        else if (name == "d2t") return d2t;
        else if (name == "dx") return dx;
        else if (name == "dz") return dz;
        else if (name == "dz1") return dz1;
        else if (name == "nx") return nx;
        else if (name == "nz") return nz;
        else {
            std::cerr << "Unknown scalar name: " << name << std::endl;
            exit(1);
        }
    }


    // Poisson Solver => PoissonSolver.cpp
    // *********************************************************************************

    class PoissonSolver {
    public:
        #if defined(GPU_POISSON)
            inline static int gpu_id = 0;
            inline static cudaStream_t stream = nullptr;           // Static stream
            // inline static cusparseHandle_t cusparseHandle = nullptr; // Static cuSPARSE handle
            inline static bool initialized = false; // Flag to track initialization per process
            PoissonSolver() {
                // Set GPU for this instance using the static gpu_id
                if (cudaSetDevice(gpu_id) != cudaSuccess) {
                    int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    std::cerr << "Rank " << rank << ": Failed to set GPU " << gpu_id << " for PoissonSolver" << std::endl;
                    exit(1);
                }
            }
        #endif

        #ifndef PETSC
            static void InitPoissonMatrix(vvm &model);
            #if defined(GPU_POISSON)
                static void Initialize(); // Static initialization
                static void Finalize();   // Static cleanup
                static void InitAMGX(vvm &model);
                static void CleanupAMGX(vvm &model);
            #else
                Eigen::SparseMatrix<double> A;
                Eigen::SparseMatrix<double> G;
            #endif
        #endif
        #if defined(STREAMFUNCTION)
            static void calpsiuw(vvm &model);
        #else
            static void cal_w(vvm &, int p = 0, int i = 0, int j = 0);
            static void cal_u(vvm &model);
            static void pubarTop_pt(vvm &model);
        #endif
    };

    // *********************************************************************************

    // Diffusion and Time Filter => NumericalProcess.cpp
    // *********************************************************************************
    class NumericalProcess {
    public:
        static void interpolate(const std::vector<double>& known_heights,
                                const std::vector<std::vector<double>>& known_data_fields,
                                const std::vector<double>& new_heights,
                                std::vector<std::vector<double>>& interpolated_data_fields);

        #if defined(DIFFUSION_VVM)
            static void Diffusion(double **var_in, double **var_out, vvm &model);
            static void DiffusionAll(vvm &model);
        #endif
        #if defined(TIMERFILTER)
            static void TimeFilter(double **previous, double **now, double **future, vvm &model);
            static void timeFilterAll(vvm &model);
        #endif
        static void Nudge_qv(vvm &model);
        static void GravityWaveDampingExponential(vvm &model);
        static void NegativeValueProcess(double **var, int nx, int nz);
    };

    // *********************************************************************************

    #if defined(WATER)
    #if defined(KESSLER_MICROPHY)
    class MicroPhysics {
    public:
        static void condensation(vvm &model); 	// condensation of qc by qv
        static void autoconversion(vvm &model); 	// autoconversion of qc to qr
        static void accretion(vvm &model); 		// accretion of qc by qr
        static void evaporation(vvm &model); 	// evaporation of rain water
    };
    #endif

        static void AddForcing(vvm &model);
    #endif

    #if defined(P3_MICROPHY)
    class P3 {
    public:
        inline static char *lookup_file_dir = strdup("../external/P3-microphysics/lookup_tables");
        inline static int nCat = 1;
        inline static bool trplMomI = false; // 3-element array
        inline static bool liqfrac = false;
        inline static char *model_name = strdup("2DVVM");
        inline static int stat = 0;
        inline static bool abort_on_err = true;
        inline static bool dowr = true;

        inline static int n_diag_2d = 1;
        inline static int n_diag_3d = 1;

        inline static bool log_predictNc  = false; // 3-element array
        inline static double scpf_pfrac   = 0.;    // dummy variable (not used), set to 0
        inline static double scpf_resfact = 0.;    // dummy variable (not used), set to 0
        inline static double clbfact_dep  = 1.;    // calibration for deposition
        inline static double clbfact_sub  = 1.;    // calibration for sublimation
        inline static bool debug_on  = false;
        inline static bool scpf_on   = false;      // cloud fraction version not used
    };
    #endif

    // Variables for tropical forcing
    double** init_th_forcing = nullptr;
    bool status_for_adding_forcing = true;
    #if defined(TROPICALFORCING)
        double* Q1LS;
        double* Q2LS;
    #endif

    class Init {
    public:
        static void Init1d(vvm &model);
        static void Init2d(vvm &model);
        static void RandomPerturbation(vvm &model, int t, double min_range=-0.25, double max_range=0.25, double standard_deviation=1.);
        
        #if defined(LOADFILE)
            static void LoadFile(vvm &model);
        #elif defined(LOADFROMPREVIOUSFILE)
            static void LoadFromPreviousFile(vvm &model);
        #elif defined(LOAD2DINIT)
            static void Load2DInit(vvm &model);
        #endif
            
    private:
        static double GetTB(int i, vvm &model);
        static double GetTHRAD(int i, int k, vvm &model);
        static double GetTH(int i, int k, vvm &model);
        #if defined(WATER)
            static double GetQVB(int k, vvm &model);
        #endif
    };

    class Output {
    public:
        static void printInit(vvm &model);
        static void create_all_directory(vvm &model);
        static void create_directory(std::string path);
        #if defined(OUTPUTNC)
            static void output_nc(int step, vvm &model);
            static void output_time_nc(int step, vvm &model);
        #endif
        static void copy_files(const std::string &source_path, const std::string &destination_path);
        static void copy_source_project(vvm &model);
        static void grads_ctl_file(vvm &model);

        #if defined(OUTPUTTXT)
            static void output_zeta(int step, vvm &model);
            static void output_th(int step, vvm &model);
            static void output_u(int step, vvm &model);
            static void output_w(int step, vvm &model);
            #if defined(WATER)
                static void output_qv(int step, vvm &model);
                static void output_qc(int step, vvm &model);
                static void output_qr(int step, vvm &model);
                static void output_precip(int step, vvm &model);
            #endif
            static void outputalltxt(int step, vvm &model);
        #endif
        
    };


    class Iteration {
    public:
        static void pzeta_pt(vvm &model);
        static void pth_pt(vvm &model);
        #if defined(WATER)
            #if defined(KESSLER_MICROPHY)
                static void pqv_pt(vvm &model);
                static void pqc_pt(vvm &model);
                static void pqr_pt(vvm &model);
            #endif
            #if defined(P3_MICROPHY)
                static void pqmicrophy_pt(vvm &model);
            #endif
        #endif

        static void updateMean(vvm &model);
        static void TimeMarching(vvm &model);
        static void nextTimeStep(vvm &model);
    };

    class Turbulence {
    public:
        static void RKM_RKH(vvm &model);
        static void Mparam(vvm &model, double **var_now, double **var_future);
        static void Hparam(vvm &model, double **var_now, double **var_future);
        static void ubarTop(vvm &model);
    };

    #if defined(RTERRTMGP)
    class Radiation {
    public:
        static void solve_radiation(vvm &model);
        static bool is_leap_year(int year);
        static int day_of_year(int year, int month, int day);
        static double calculate_scaling_factor(int year, int month, int day);
        static double calculate_cos_zenith(int year, int month, int day, double hour, double minute, double second,
                                    double longitude, double latitude);

    };
    #endif

private:
    double rdx = 0;                              ///< 1/dx, calculated from Config_VVM given by users.
    double r2dx = 0;                             ///< 1 / (2dx), calculated from Config_VVM given by users.
    double rdz = 0;                              ///< 1 / dz, calculated from Config_VVM given by users.
    double r2dz = 0;                             ///< 1 / (2dz), calculated from Config_VVM given by users.
    double rdx2 = 0;                             ///< 1 / (dx^2), calculated from Config_VVM given by users.
    double rdz2 = 0;                             ///< 1 / (dz^2), calculated from Config_VVM given by users.
    int nx = 0;                                  ///< Number of grid points in x direction, calculated from Config_VVM given by users.
    int nz = 0;                                  ///< Number of grid points in z direction, calculated from Config_VVM given by users.
    double dt = 0;                               ///< From Config_VVM given by users.
    double d2t = 0;                              ///< From Config_VVM given by users.
    double dx = 0;                               ///< From Config_VVM given by users.
    double dz = 0;                               ///< From Config_VVM given by users.
    double dz1 = 0;                               ///< From Config_VVM given by users.
    int XRANGE = 0;                              ///< From Config_VVM given by users.
    int ZRANGE = 0;                              ///< From Config_VVM given by users.
    double TIMEEND = 0;                          ///< From Config_VVM given by users.
    int TIMEROUTPUTSIZE = 0;                     ///< From Config_VVM given by users.
    int OUTPUTSTEP = 0;                          ///< From Config_VVM given by users.
    double Kx = 0;                               ///< From Config_VVM given by users.
    double Kz = 0;                               ///< From Config_VVM given by users.
    double TIMETS = 0;                           ///< From Config_VVM given by users.
    double tolerance = 0;                        ///< From Config_VVM given by users.
    double GRAVITY = 0;                          ///< From Config_VVM given by users.
    double Cp = 0;                               ///< From Config_VVM given by users.
    double Cv = 0;                               ///< From Config_VVM given by users.
    double Rd = 0;                               ///< From Config_VVM given by users.
    double Lv = 0;                               ///< From Config_VVM given by users.
    double P0 = 0;                               ///< From Config_VVM given by users.
    double PSURF = 0;                            ///< From Config_VVM given by users.
    double addforcingtime = 0;                   ///< From Config_VVM given by users.
    int CASE = 0;                                ///< From Config_VVM given by users.
    int year = 2025;
    int month = 3;
    int day = 20;
    double hour = 12.;
    double minute = 0.;
    double second = 0.;
    double lon = 0.;
    double lat = 0.;
    int k_diff_start = 0;
    double CZ1 = 0.;
    double CZ2 = 0.;

    // 0D variables
    int step = 0;                            ///< The current time step.
    double ubarTopp = 0;                         ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.
    double ubarTop = 0;                          ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.
    double ubarTopm = 0;                         ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.
    double moisture_nudge_time = 0.;         ///< The time for nudging the moisture field.
    double dubarTop_advect[2] = {0,0};

    // 1D variables
    double *thb = nullptr;                              ///< Horizontal mean potential temperature profile.
    double *thb_init = nullptr;                         ///< Initial horizontal mean potential temperature profile.
    double *thbm = nullptr;                             ///< Horizontal mean potential temperature profile for previous step.
    double *thb_zeta = nullptr;                         ///< Horizontal mean potential temperature profile at grid upper edge.
    double *rhou = nullptr;                             ///< Horizontal mean density profile at grid center.
    double *rhow = nullptr;                             ///< Horizontal mean density profile at grid upper edge.
    double *pib = nullptr;                              ///< Horizontal mean non-dimensional height profile at grid center.
    double *pib_lev = nullptr;                              ///< Horizontal mean non-dimensional height profile at grid boundary.
    double *qvb = nullptr;                              ///< Horizontal mean water vapor profile at grid center.
    double *qvb0 = nullptr;                              ///< Horizontal mean water vapor profile at grid center.
    double *qvsb = nullptr;                             ///< Horizontal mean saturated water vapor profile at grid center.
    double *pb = nullptr;                               ///< Horizontal mean pressure profile at grid center.
    double *pb_lev = nullptr;                               ///< Horizontal mean pressure profile at grid boundary.
    double *xi = nullptr;                               ///< The velocity potential in x-direction at top boundary grid center.
    double *uxi = nullptr;
    double *thvb = nullptr;
    double *thvbm = nullptr;
    double *x = nullptr;
    double *z = nullptr;
    double *z_zeta = nullptr;
    double *lambda2 = nullptr;
    double *lambda2_zeta = nullptr;
    double *th_ground = nullptr;
    double *qvs_ground = nullptr;
    double *addflux = nullptr;
    double *heatflux = nullptr;
    double *waterflux = nullptr;
    double *nudge_tau = nullptr;
    double *RH = nullptr;
    double *dz_th = nullptr;
    double *dz_zeta = nullptr;
    double *flex_height_coef_th = nullptr; ///< Coefficient for flexible height adjustment, used in the model design part.
    double *flex_height_coef_zeta = nullptr; ///< Coefficient for flexible height adjustment, used in the model design part.
    double *flex_height_coef_th_mean = nullptr; ///< Coefficient for flexible height adjustment, used in the model design part.
    double *flex_height_coef_zeta_mean = nullptr; ///< Coefficient for flexible height adjustment, used in the model design part.
    double *ubar = nullptr;

    #if defined(GPU_POISSON)
        int *row_ptr_w = nullptr;
        int *col_idx_w = nullptr;
        double *values_w = nullptr;
        int *row_ptr_u = nullptr;
        int *col_idx_u = nullptr;
        double *values_u = nullptr;
        int nnz_w = 0;
        int nnz_u = 0;
        // AMGX members for w (A matrix)
        AMGX_config_handle cfg_w;
        AMGX_resources_handle rsc_w;
        AMGX_matrix_handle A;
        AMGX_vector_handle b_vec_w, x_vec_w;
        AMGX_solver_handle solver_w;

        // AMGX members for u (G matrix)
        AMGX_config_handle cfg_u;
        AMGX_resources_handle rsc_u;
        AMGX_matrix_handle G;
        AMGX_vector_handle h_vec_u, y_vec_u;
        AMGX_solver_handle solver_u;

        bool initialized = false;

        int *d_row_ptr_w, *d_col_idx_w;
        double *d_values_w;
        int *d_row_ptr_u, *d_col_idx_u;
        double *d_values_u;

        double *d_b_w, *d_x_w; // For cal_w
        double *d_b_u, *d_x_u; // For cal_u
    #endif

    // 2D variables
    double **zetap = nullptr;
    double **zeta = nullptr;
    double **zetam = nullptr;
    double **thp = nullptr;
    double **th = nullptr;
    double **thm = nullptr;
    double **u = nullptr;
    double **w = nullptr;
    double **RKM = nullptr;
    double **RKH = nullptr;
    double **U_w = nullptr;
    double **W_u = nullptr;
    #if defined(RTERRTMGP)
        double **T = nullptr;
        double **T_lev = nullptr;
        double **radiation_heating_rate = nullptr; // K/s
    #endif

    double *zetapcont = nullptr;
    double *zetacont = nullptr;
    double *zetamcont = nullptr;
    double *thpcont = nullptr;
    double *thcont = nullptr;
    double *thmcont = nullptr;
    double *ucont = nullptr;
    double *wcont = nullptr;
    double *init_th_forcingcont = nullptr;
    double *RKMcont = nullptr;
    double *RKHcont = nullptr;
    double *U_wcont = nullptr;
    double *W_ucont = nullptr;
    #if defined(RTERRTMGP)
        double *Tcont = nullptr;
        double *T_levcont = nullptr;
        double *radiation_heating_ratecont = nullptr;
    #endif
    
    
    #if defined(STREAMFUNCTION)
        double** psi;
    #endif

    #if defined(WATER)
        double **qvp = nullptr, **qv = nullptr, **qvm = nullptr;
        double **qcp = nullptr, **qc = nullptr, **qcm = nullptr;
        double **qrp = nullptr, **qr = nullptr, **qrm = nullptr;
        #if defined(KESSLER_MICROPHY)
            double **evaporation = nullptr;
            double **accretion = nullptr;
            double **autoconversion = nullptr;
            double **condensation = nullptr;
            double *precip = nullptr;
        #endif
        #if defined(P3_MICROPHY)
            double **ncp = nullptr, **nc = nullptr, **ncm = nullptr;
            double **nrp = nullptr, **nr = nullptr, **nrm = nullptr;
            double **qitotp = nullptr, **qitot = nullptr, **qitotm = nullptr;
            double **qirimp = nullptr, **qirim = nullptr, **qirimm = nullptr;
            double **qiliqp = nullptr, **qiliq = nullptr, **qiliqm = nullptr;
            double **nip = nullptr, **ni = nullptr, **nim = nullptr;
            double **birimp = nullptr, **birim = nullptr, **birimm = nullptr;
            double *precip_liq = nullptr, *precip_sol = nullptr, *precip = nullptr;
            double **diag_ze = nullptr, **diag_effc = nullptr, **diag_effi = nullptr;
            double **diag_vmi = nullptr, **diag_di = nullptr, **diag_rhoi = nullptr, **cldfrac = nullptr;
            double **diag_2d = nullptr, ***diag_3d = nullptr;
            double **dz_all = nullptr;
            double **w_all = nullptr;
            double **pb_all = nullptr;
            double **pb_lev_all = nullptr;
            double **zi_all = nullptr;
            double **ssat_all = nullptr;
        #endif


        double *qvpcont, *qvcont, *qvmcont;
        double *qcpcont, *qccont, *qcmcont;
        double *qrpcont, *qrcont, *qrmcont;
        #if defined(KESSLER_MICROPHY)
            double *evaporationcont = nullptr;
            double *accretioncont = nullptr;
            double *autoconversioncont = nullptr;
            double *condensationcont = nullptr;
        #endif

        #if defined(P3_MICROPHY)
            double *ncpcont = nullptr, *nccont = nullptr, *ncmcont = nullptr;
            double *nrpcont = nullptr, *nrcont = nullptr, *nrmcont = nullptr;
            double *qitotpcont = nullptr, *qitotcont = nullptr, *qitotmcont = nullptr;
            double *qirimpcont = nullptr, *qirimcont = nullptr, *qirimmcont = nullptr;
            double *qiliqpcont = nullptr, *qiliqcont = nullptr, *qiliqmcont = nullptr;
            double *nipcont = nullptr, *nicont = nullptr, *nimcont = nullptr;
            double *birimpcont = nullptr, *birimcont = nullptr, *birimmcont = nullptr;
            double *diag_zecont = nullptr, *diag_effccont = nullptr, *diag_efficont = nullptr;
            double *diag_vmicont = nullptr, *diag_dicont = nullptr, *diag_rhoicont = nullptr, *cldfraccont = nullptr;
            double *diag_2dcont = nullptr, *diag_3dcont = nullptr;
            double *dz_allcont = nullptr;
            double *w_allcont = nullptr;
            double *pb_allcont = nullptr;
            double *pb_lev_allcont = nullptr;
            double *zi_allcont = nullptr;
            double *ssat_allcont = nullptr;
        #endif
    #endif

    // #####################################################################################
    // Used for AB2. These variables are declared but not initialized if it's not AB2
    double ***dth_advect = nullptr;
    double ***dth_buoyancy = nullptr;
    double ***dzeta_advect = nullptr;
    
    double *dth_advectcont = nullptr;
    double *dth_buoyancycont = nullptr;
    double *dzeta_advectcont = nullptr;

    #if defined(WATER)
        double ***dqv_advect = nullptr;
        double ***dqc_advect = nullptr;
        double ***dqr_advect = nullptr;

        double *dqv_advectcont = nullptr;
        double *dqc_advectcont = nullptr;
        double *dqr_advectcont = nullptr;

        #if defined(KESSLER_MICROPHY)
            double ***dqr_VT = nullptr;
            double *dqr_VTcont = nullptr;
        #endif

        #if defined(P3_MICROPHY)
            double ***dnc_advect = nullptr;
            double ***dnr_advect = nullptr;
            double ***dni_advect = nullptr;
            double ***dqitot_advect = nullptr;
            double ***dqirim_advect = nullptr;
            double ***dqiliq_advect = nullptr;
            double ***dbirim_advect = nullptr;

            double *dnc_advectcont = nullptr;
            double *dnr_advectcont = nullptr;
            double *dni_advectcont = nullptr;
            double *dqitot_advectcont = nullptr;
            double *dqirim_advectcont = nullptr;
            double *dqiliq_advectcont = nullptr;
            double *dbirim_advectcont = nullptr;
        #endif
    #endif
    // #####################################################################################


    double *t_advection = nullptr;
    double *t_poisson = nullptr;
    double *t_diffusion = nullptr;
    double *t_microphysics = nullptr;
    double *t_all = nullptr;

    #if !defined(PETSC) && !defined(GPU_POISSON) 
        Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>((nx-2)*(nz-3), (nx-2)*(nz-3));
        Eigen::SparseMatrix<double> G = Eigen::SparseMatrix<double>(nx-2, nx-2);
    #endif

    static double getTHV(int i, int k, vvm &model);

    void deallocateMemory() {
        // Free the allocated memory
        destroy_variable(t_advection);
        destroy_variable(t_poisson);
        destroy_variable(t_diffusion);
        destroy_variable(t_microphysics);
        destroy_variable(t_all);

        destroy_variable(thb);
        destroy_variable(thbm);
        destroy_variable(thb_zeta);
        destroy_variable(thb_init);
        destroy_variable(rhou);
        destroy_variable(rhow);
        destroy_variable(pib);
        destroy_variable(pib_lev);
        destroy_variable(qvb);
        destroy_variable(qvb0);
        destroy_variable(qvsb);
        destroy_variable(pb);
        destroy_variable(pb_lev);
        destroy_variable(xi);
        destroy_variable(uxi);
        destroy_variable(thvb);
        destroy_variable(thvbm);
        destroy_variable(x);
        destroy_variable(z);
        destroy_variable(z_zeta);
        destroy_variable(lambda2);
        destroy_variable(lambda2_zeta);
        destroy_variable(th_ground);
        destroy_variable(qvs_ground);
        destroy_variable(addflux);
        destroy_variable(heatflux);
        destroy_variable(waterflux);
        destroy_variable(nudge_tau);
        destroy_variable(RH);
        destroy_variable(dz_th);
        destroy_variable(dz_zeta);
        destroy_variable(flex_height_coef_th);
        destroy_variable(flex_height_coef_zeta);
        destroy_variable(ubar);

        destroy_variable(zetap, zetapcont);
        destroy_variable(zeta, zetacont);
        destroy_variable(zetam, zetamcont);
        destroy_variable(thp, thpcont);
        destroy_variable(th, thcont);
        destroy_variable(thm, thmcont);
        destroy_variable(u, ucont);
        destroy_variable(w, wcont);
        destroy_variable(init_th_forcing, init_th_forcingcont);
        destroy_variable(RKM, RKMcont);
        destroy_variable(RKH, RKHcont);
        destroy_variable(U_w, U_wcont);
        destroy_variable(W_u, W_ucont);

        #if defined(RTERRTMGP)
            destroy_variable(T, Tcont);
            destroy_variable(radiation_heating_rate, radiation_heating_ratecont);
        #endif

        #if defined(STREAMFUNCTION)
            destroy_variable(psi, psicont);
        #endif

        #if defined(WATER)
            destroy_variable(qvp, qvpcont);
            destroy_variable(qv, qvcont);
            destroy_variable(qvm, qvmcont);
            destroy_variable(qcp, qcpcont);
            destroy_variable(qc, qccont);
            destroy_variable(qcm, qcmcont);
            destroy_variable(qrp, qrpcont);
            destroy_variable(qr, qrcont);
            destroy_variable(qrm, qrmcont);
            #if defined(KESSLER_MICROPHY)
                delete[] precip;

                destroy_variable(evaporation, evaporationcont);
                destroy_variable(accretion, accretioncont);
                destroy_variable(autoconversion, autoconversioncont);
                destroy_variable(condensation, condensationcont);
            #endif

            #if defined(P3_MICROPHY)
                
                destroy_variable(ncp, ncpcont);
                destroy_variable(nc, nccont);
                destroy_variable(ncm, ncmcont);
                destroy_variable(nrp, nrpcont);
                destroy_variable(nr, nrcont);
                destroy_variable(nrm, nrmcont);
                destroy_variable(qitotp, qitotpcont);
                destroy_variable(qitot, qitotcont);
                destroy_variable(qitotm, qitotmcont);
                destroy_variable(qirimp, qirimpcont);
                destroy_variable(qirim, qirimcont);
                destroy_variable(qirimm, qirimmcont);
                destroy_variable(qiliqp, qiliqpcont);
                destroy_variable(qiliq, qiliqcont);
                destroy_variable(qiliqm, qiliqmcont);
                destroy_variable(nip, nipcont);
                destroy_variable(ni, nicont);
                destroy_variable(nim, nimcont);
                destroy_variable(birimp, birimpcont);
                destroy_variable(birim, birimcont);
                destroy_variable(birimm, birimmcont);
                destroy_variable(diag_ze, diag_zecont);
                destroy_variable(diag_effc, diag_effccont);
                destroy_variable(diag_effi, diag_efficont);
                destroy_variable(diag_vmi, diag_vmicont);
                destroy_variable(diag_di, diag_dicont);
                destroy_variable(diag_rhoi, diag_rhoicont);
                destroy_variable(cldfrac, cldfraccont);
                destroy_variable(diag_2d, diag_2dcont);
                destroy_variable(dz_all, dz_allcont);
                destroy_variable(w_all, w_allcont);
                destroy_variable(pb_all, pb_allcont);
                destroy_variable(pb_lev_all, pb_lev_allcont);
                destroy_variable(zi_all, zi_allcont);
                destroy_variable(ssat_all, ssat_allcont);

                destroy_variable(precip_liq);
                destroy_variable(precip_sol);
                destroy_variable(precip);
            #endif
        #endif

        #if defined(AB2)
            destroy_variable(dth_advect, dth_advectcont, nx);
            destroy_variable(dth_buoyancy, dth_buoyancycont, nx);
            destroy_variable(dzeta_advect, dzeta_advectcont, nx);
            
            #if defined(WATER)
                destroy_variable(dqv_advect, dqv_advectcont, nx);
                destroy_variable(dqc_advect, dqc_advectcont, nx);
                destroy_variable(dqr_advect, dqr_advectcont, nx);
                #if defined(KESSLER_MICROPHY)
                    destroy_variable(dqr_VT, dqr_VTcont);
                #endif
                #if defined(P3_MICROPHY)
                    destroy_variable(dnc_advect, dnc_advectcont, nx);
                    destroy_variable(dnr_advect, dnr_advectcont, nx);
                    destroy_variable(dni_advect, dni_advectcont, nx);
                    destroy_variable(dqitot_advect, dqitot_advectcont, nx);
                    destroy_variable(dqirim_advect, dqirim_advectcont, nx);
                    destroy_variable(dqiliq_advect, dqiliq_advectcont, nx);
                    destroy_variable(dbirim_advect, dbirim_advectcont, nx);
                #endif
            #endif

        #endif
        #if defined(P3_MICROPHY)
            destroy_variable(diag_3d, diag_3dcont, nx);
        #endif

        #if defined(TROPICALFORCING)
            destroy_variable(Q1LS);
            destroy_variable(Q2LS);
        #endif
    }

    void allocateMemory() {
        // 1D arrays for timing
        create_variable(t_advection, TIMEROUTPUTSIZE);
        create_variable(t_poisson, TIMEROUTPUTSIZE);
        create_variable(t_diffusion, TIMEROUTPUTSIZE);
        create_variable(t_microphysics, TIMEROUTPUTSIZE);
        create_variable(t_all, TIMEROUTPUTSIZE);

        // 1D arrays
        create_variable(thb, nz);
        create_variable(thbm, nz);
        create_variable(thb_zeta, nz);
        create_variable(thb_init, nz);
        create_variable(rhou, nz);
        create_variable(rhow, nz);
        create_variable(pib, nz);
        create_variable(pib_lev, nz+1);
        create_variable(qvb, nz);
        create_variable(qvb0, nz);
        create_variable(qvsb, nz);
        create_variable(pb, nz);
        create_variable(pb_lev, nz+1);
        create_variable(xi, nx);
        create_variable(uxi, nx);
        create_variable(thvb, nz);
        create_variable(thvbm, nz);
        create_variable(x, nx);
        create_variable(z, nz);
        create_variable(z_zeta, nz);
        create_variable(lambda2, nz);
        create_variable(lambda2_zeta, nz);
        #if defined(TROPICALFORCING)
            create_variable(Q1LS, nz);
            create_variable(Q2LS, nz);
        #endif
        create_variable(th_ground, nx);
        create_variable(qvs_ground, nx);
        create_variable(addflux, nx);
        create_variable(heatflux, nx);
        create_variable(waterflux, nx);
        create_variable(nudge_tau, nz);
        create_variable(RH, nz);
        create_variable(dz_th, nz);
        create_variable(dz_zeta, nz);
        create_variable(flex_height_coef_th, nz);
        create_variable(flex_height_coef_zeta, nz);
        create_variable(flex_height_coef_th_mean, nz);
        create_variable(flex_height_coef_zeta_mean, nz);
        create_variable(ubar, nz);

        // 2D arrays
        create_variable(zetap, zetapcont, nx, nz);
        create_variable(zeta, zetacont, nx, nz);
        create_variable(zetam, zetamcont, nx, nz);
        create_variable(thp, thpcont, nx, nz);
        create_variable(th, thcont, nx, nz);
        create_variable(thm, thmcont, nx, nz);
        create_variable(u, ucont, nx, nz);
        create_variable(w, wcont, nx, nz);
        create_variable(init_th_forcing, init_th_forcingcont, nx, nz);
        create_variable(RKM, RKMcont, nx, nz);
        create_variable(RKH, RKHcont, nx, nz);
        create_variable(U_w, U_wcont, nx, nz);
        create_variable(W_u, W_ucont, nx, nz);

        #if defined(RTERRTMGP)
            create_variable(T, Tcont, nx, nz);
            create_variable(T_lev, T_levcont, nx, nz+1);
            create_variable(radiation_heating_rate, radiation_heating_ratecont, nx, nz);
        #endif

        #if defined(STREAMFUNCTION)
            create_variable(psi, psicont, nx, nz);
        #endif

        #if defined(WATER)
            #if defined(KESSLER_MICROPHY)
                create_variable(precip, nx);
                create_variable(evaporation, evaporationcont, nx, nz);
                create_variable(accretion, accretioncont, nx, nz);
                create_variable(autoconversion, autoconversioncont, nx, nz);
                create_variable(condensation, condensationcont, nx, nz);
            #endif

            create_variable(qvp, qvpcont, nx, nz);
            create_variable(qv, qvcont, nx, nz);
            create_variable(qvm, qvmcont, nx, nz);
            create_variable(qcp, qcpcont, nx, nz);
            create_variable(qc, qccont, nx, nz);
            create_variable(qcm, qcmcont, nx, nz);
            create_variable(qrp, qrpcont, nx, nz);
            create_variable(qr, qrcont, nx, nz);
            create_variable(qrm, qrmcont, nx, nz);

            #if defined(P3_MICROPHY)
                create_variable(ncp, ncpcont, nx, nz);
                create_variable(nc, nccont, nx, nz);
                create_variable(ncm, ncmcont, nx, nz);
                create_variable(nrp, nrpcont, nx, nz);
                create_variable(nr, nrcont, nx, nz);
                create_variable(nrm, nrmcont, nx, nz);
                create_variable(qitotp, qitotpcont, nx, nz);
                create_variable(qitot, qitotcont, nx, nz);
                create_variable(qitotm, qitotmcont, nx, nz);
                create_variable(qirimp, qirimpcont, nx, nz);
                create_variable(qirim, qirimcont, nx, nz);
                create_variable(qirimm, qirimmcont, nx, nz);
                create_variable(qiliqp, qiliqpcont, nx, nz);
                create_variable(qiliq, qiliqcont, nx, nz);
                create_variable(qiliqm, qiliqmcont, nx, nz);
                create_variable(nip, nipcont, nx, nz);
                create_variable(ni, nicont, nx, nz);
                create_variable(nim, nimcont, nx, nz);
                create_variable(birimp, birimpcont, nx, nz);
                create_variable(birim, birimcont, nx, nz);
                create_variable(birimm, birimmcont, nx, nz);
                create_variable(diag_ze, diag_zecont, nx, nz);
                create_variable(diag_effc, diag_effccont, nx, nz);
                create_variable(diag_effi, diag_efficont, nx, nz);
                create_variable(diag_vmi, diag_vmicont, nx, nz);
                create_variable(diag_di, diag_dicont, nx, nz);
                create_variable(diag_rhoi, diag_rhoicont, nx, nz);
                create_variable(cldfrac, cldfraccont, nx, nz);
                create_variable(diag_2d, diag_2dcont, nx, vvm::P3::n_diag_2d);
                create_variable(dz_all, dz_allcont, nx, nz);
                create_variable(w_all, w_allcont, nx, nz);
                create_variable(pb_all, pb_allcont, nx, nz);
                create_variable(pb_lev_all, pb_lev_allcont, nx, nz+1);
                create_variable(zi_all, zi_allcont, nx, nz);
                create_variable(ssat_all, ssat_allcont, nx, nz);

                create_variable(precip_liq, nx);
                create_variable(precip_sol, nx);
                create_variable(precip, nx);
            #endif
        #endif

        #if defined(AB2)
            create_variable(dth_advect, dth_advectcont, nx, nz, 2);
            create_variable(dth_buoyancy, dth_buoyancycont, nx, nz, 2);
            create_variable(dzeta_advect, dzeta_advectcont, nx, nz, 2);
            #if defined(WATER)
                create_variable(dqv_advect, dqv_advectcont, nx, nz, 2);
                create_variable(dqc_advect, dqc_advectcont, nx, nz, 2);
                create_variable(dqr_advect, dqr_advectcont, nx, nz, 2);
                #if defined(KESSLER_MICROPHY)
                    create_variable(dqr_VT, dqr_VTcont, nx, nz, 2);
                #endif

                #if defined(P3_MICROPHY)
                    create_variable(dnc_advect, dnc_advectcont, nx, nz, 2);
                    create_variable(dnr_advect, dnr_advectcont, nx, nz, 2);
                    create_variable(dni_advect, dni_advectcont, nx, nz, 2);
                    create_variable(dqitot_advect, dqitot_advectcont, nx, nz, 2);
                    create_variable(dqirim_advect, dqirim_advectcont, nx, nz, 2);
                    create_variable(dqiliq_advect, dqiliq_advectcont, nx, nz, 2);
                    create_variable(dbirim_advect, dbirim_advectcont, nx, nz, 2);
                #endif
            #endif

            #if defined(P3_MICROPHY)
                create_variable(diag_3d, diag_3dcont, nx, nz, vvm::P3::n_diag_3d);
            #endif
        #endif
    }
};

