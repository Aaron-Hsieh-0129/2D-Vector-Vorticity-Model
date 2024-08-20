#include "Config.hpp"
#include <string>
#ifndef PETSC
    // #include "../include/Eigen/Sparse"
    #include <Eigen/Sparse>
#endif

class Config_VVM {
public:
    Config_VVM(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, 
           std::string outputpath, int OUTPUTSTEP, double Kx, double Kz, double TIMETS, double tolerance,
           double GRAVITY, double Cp, double Cv, double Rd, double Lv, double P0, double PSURF, double addforcingtime, int CASE)
        : dt(dt), dx(dx), dz(dz), XRANGE(XRANGE+2*dx), ZRANGE(ZRANGE+2*dz), TIMEEND(TIMEEND), TIMEROUTPUTSIZE(TIMEROUTPUTSIZE), 
          outputpath(outputpath), OUTPUTSTEP(OUTPUTSTEP), Kx(Kx), Kz(Kz), TIMETS(TIMETS),
          tolerance(tolerance), GRAVITY(GRAVITY), Cp(Cp), Cv(Cv), Rd(Rd), Lv(Lv), P0(P0), PSURF(PSURF), addforcingtime(addforcingtime), CASE(CASE) {}
    ~Config_VVM() {}

    double dt;              ///< Time step for vvm [s].
    double dx;              ///< Grid size in x-direction [m].
    double dz;              ///< Grid size in z-direction [m]. It should be the same as dx.
    int XRANGE;             ///< Domain size of the model in x-direction [m].
    int ZRANGE;             ///< Domain size of the model in z-direction [m].
    double TIMEEND;         ///< End time of the simulation [s].
    int TIMEROUTPUTSIZE;    ///< The size of the timer output.
    std::string outputpath; ///< The path for the output file. It should be a directory, such as "/data/vvm/".
    int OUTPUTSTEP;         ///< The output interval for the output file.
    double Kx;              ///< The eddy diffusion coefficient in x-direction [m^2/s], this is activated when DIFFUSION flag is turned on. If the flag is not turned on, the coeffcient will be calculated through the turbulent closure.
    double Kz;              ///< The eddy diffusion coefficient in z-direction [m^2/s], this is activated when DIFFUSION flag is turned on. If the flag is not turned on, the coeffcient will be calculated through the turbulent closure.
    double TIMETS;          ///< The time filter coefficient [s] for Leapfrog. This is activated when TIMEFILTER flag is turned on. Only need to turn on when Leapfrog is used.
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
          nz(config.ZRANGE/config.dz), dt(config.dt), 
          d2t(2.0 * config.dt), dx(config.dx), dz(config.dz),
          XRANGE(config.XRANGE), ZRANGE(config.ZRANGE), TIMEEND(config.TIMEEND),
          TIMEROUTPUTSIZE(config.TIMEROUTPUTSIZE), outputpath(config.outputpath), OUTPUTSTEP(config.OUTPUTSTEP), Kx(config.Kx), Kz(config.Kz),
          TIMETS(config.TIMETS), tolerance(config.tolerance),
          GRAVITY(config.GRAVITY),
          Cp(config.Cp), Cv(config.Cv),
          Rd(config.Rd), Lv(config.Lv),
          P0(config.P0), PSURF(config.PSURF), addforcingtime(config.addforcingtime), CASE(config.CASE)
    {
        allocateMemory();
    }
    ~vvm() {
        printf("Free vvm\n");
        deallocateMemory();
    }

    void deallocateMemory() {
        // Free the allocated memory
        delete[] thb;
        delete[] thbm;
        delete[] thb_zeta;
        delete[] thb_init;
        delete[] rhou;
        delete[] rhow;
        delete[] pib;
        delete[] qvb;
        delete[] qvsb;
        delete[] pb;
        delete[] xi;
        delete[] uxi;
        delete[] thvb;
        delete[] thvbm;
        delete[] z;
        delete[] z_zeta;
        delete[] lambda2;

        deallocate2DContinuousArray(zetap, zetapcont);
        deallocate2DContinuousArray(zeta, zetacont);
        deallocate2DContinuousArray(zetam, zetamcont);
        deallocate2DContinuousArray(thp, thpcont);
        deallocate2DContinuousArray(th, thcont);
        deallocate2DContinuousArray(thm, thmcont);
        deallocate2DContinuousArray(u, ucont);
        deallocate2DContinuousArray(w, wcont);
        deallocate2DContinuousArray(init_th_forcing, init_th_forcingcont);
        deallocate2DContinuousArray(RKM, RKMcont);
        deallocate2DContinuousArray(RKH, RKHcont);
        deallocate2DContinuousArray(U_w, U_wcont);
        deallocate2DContinuousArray(W_u, W_ucont);

        #if defined(STREAMFUNCTION)
            deallocate2DContinuousArray(psi, psicont);
        #endif

        #if defined(WATER)
            delete[] precip;

            deallocate2DContinuousArray(qvp, qvpcont);
            deallocate2DContinuousArray(qv, qvcont);
            deallocate2DContinuousArray(qvm, qvmcont);
            deallocate2DContinuousArray(qcp, qcpcont);
            deallocate2DContinuousArray(qc, qccont);
            deallocate2DContinuousArray(qcm, qcmcont);
            deallocate2DContinuousArray(qrp, qrpcont);
            deallocate2DContinuousArray(qr, qrcont);
            deallocate2DContinuousArray(qrm, qrmcont);
            deallocate2DContinuousArray(evaporation, evaporationcont);
            deallocate2DContinuousArray(accretion, accretioncont);
            deallocate2DContinuousArray(autoconversion, autoconversioncont);
            deallocate2DContinuousArray(condensation, condensationcont);
        #endif

        #if defined(AB2)
            deallocate3DContinuousArray(dth_advect, dth_advectcont);
            deallocate3DContinuousArray(dth_buoyancy, dth_buoyancycont);
            deallocate3DContinuousArray(dzeta_advect, dzeta_advectcont);
            
            #if defined(WATER)
                deallocate3DContinuousArray(dqv_advect, dqv_advectcont);
                deallocate3DContinuousArray(dqc_advect, dqc_advectcont);
                deallocate3DContinuousArray(dqr_advect, dqr_advectcont);
                deallocate3DContinuousArray(dqr_VT, dqr_VTcont);
            #endif
        #endif

        delete[] t_advection;
        delete[] t_poisson;
        delete[] t_diffusion;
        delete[] t_microphysics;
        delete[] t_all;

        #if defined(TROPICALFORCING)
            delete[] Q1LS;
            delete[] Q2LS;
        #endif
    }

    void allocateMemory() {
        t_advection = new double[TIMEROUTPUTSIZE];
        t_poisson = new double[TIMEROUTPUTSIZE];
        t_diffusion = new double[TIMEROUTPUTSIZE];
        t_microphysics = new double[TIMEROUTPUTSIZE];
        t_all = new double[TIMEROUTPUTSIZE];

        // 1D arrays
        thb = new double[nz];
        thbm = new double[nz];
        thb_zeta = new double[nz];
        thb_init = new double[nz];
        rhou = new double[nz];
        rhow = new double[nz];
        pib = new double[nz];
        qvb = new double[nz];
        qvsb = new double[nz];
        pb = new double[nz];
        xi = new double[nx];
        uxi = new double[nx];
        thvb = new double[nz];
        thvbm = new double[nz];
        z = new double[nz];
        z_zeta = new double[nz];
        lambda2 = new double[nz];
        #if defined(TROPICALFORCING)
            Q1LS = new double[nz];
            Q2LS = new double[nz];
        #endif

        // 2D arrays
        zetap = allocate2DContinuousArray(nx, nz, zetapcont);
        zeta = allocate2DContinuousArray(nx, nz, zetacont);
        zetam = allocate2DContinuousArray(nx, nz, zetamcont);
        thp = allocate2DContinuousArray(nx, nz, thpcont);
        th = allocate2DContinuousArray(nx, nz, thcont);
        thm = allocate2DContinuousArray(nx, nz, thmcont);
        u = allocate2DContinuousArray(nx, nz, ucont);
        w = allocate2DContinuousArray(nx, nz, wcont);
        init_th_forcing = allocate2DContinuousArray(nx, nz, init_th_forcingcont);
        RKM = allocate2DContinuousArray(nx, nz, RKMcont);
        RKH = allocate2DContinuousArray(nx, nz, RKHcont);
        U_w = allocate2DContinuousArray(nx, nz, U_wcont);
        W_u = allocate2DContinuousArray(nx, nz, W_ucont);
        
        #if defined(STREAMFUNCTION)
            psi = allocate2DContinuousArray(nx, nz, psicont);
        #endif

        #if defined(WATER)
            precip = new double[nx];

            qvp = allocate2DContinuousArray(nx, nz, qvpcont);
            qv = allocate2DContinuousArray(nx, nz, qvcont);
            qvm = allocate2DContinuousArray(nx, nz, qvmcont);
            qcp = allocate2DContinuousArray(nx, nz, qcpcont);
            qc = allocate2DContinuousArray(nx, nz, qccont);
            qcm = allocate2DContinuousArray(nx, nz, qcmcont);
            qrp = allocate2DContinuousArray(nx, nz, qrpcont);
            qr = allocate2DContinuousArray(nx, nz, qrcont);
            qrm = allocate2DContinuousArray(nx, nz, qrmcont);
            evaporation = allocate2DContinuousArray(nx, nz, evaporationcont);
            accretion = allocate2DContinuousArray(nx, nz, accretioncont);
            autoconversion = allocate2DContinuousArray(nx, nz, autoconversioncont);
            condensation = allocate2DContinuousArray(nx, nz, condensationcont);
        #endif

        #if defined(AB2)
            dth_advect = allocate3DContinuousArray(nx, nz, 2, dth_advectcont);
            dth_buoyancy = allocate3DContinuousArray(nx, nz, 2, dth_buoyancycont);
            dzeta_advect = allocate3DContinuousArray(nx, nz, 2, dzeta_advectcont);
            #if defined(WATER)
                dqv_advect = allocate3DContinuousArray(nx, nz, 2, dqv_advectcont);
                dqc_advect = allocate3DContinuousArray(nx, nz, 2, dqc_advectcont);
                dqr_advect = allocate3DContinuousArray(nx, nz, 2, dqr_advectcont);
                dqr_VT = allocate3DContinuousArray(nx, nz, 2, dqr_VTcont);
            #endif
        #endif
    }

    static double** allocate2DContinuousArray(int rows, int cols, double*& contMemory) {
        double** array = new double*[rows];
        contMemory = new double[rows * cols]; // Allocate continuous memory block
        for (int i = 0; i < rows; ++i) {
            array[i] = &contMemory[i * cols]; // Point to segments within continuous block
        }
        return array;
    }

    static void deallocate2DContinuousArray(double** array, double* contMemory) {
        if (array != nullptr) {
            delete[] contMemory; // Deallocate the continuous block of memory
            delete[] array;      // Deallocate the array of pointers
        }
    }

    #if defined(AB2)
    double*** allocate3DContinuousArray(int dim1, int dim2, int dim3, double*& contMemory) {
        double*** array = new double**[dim1];
        contMemory = new double[dim1 * dim2 * dim3]; // Allocate continuous memory block
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double*[dim2];
            for (int j = 0; j < dim2; ++j) {
                array[i][j] = &contMemory[i * dim2 * dim3 + j * dim3]; // Point to segments within continuous block
            }
        }
        return array;
    }

    void deallocate3DContinuousArray(double*** array, double* contMemory) {
        if (array != nullptr) {
            delete[] contMemory; // Deallocate the continuous block of memory
            for (int i = 0; i < nx; ++i) {
                delete[] array[i]; // Deallocate the array of pointers
            }
            delete[] array;      // Deallocate the array of pointers
        }
    }
    #endif

    double rdx;                              ///< 1/dx, calculated from Config_VVM given by users.
    double r2dx;                             ///< 1 / (2dx), calculated from Config_VVM given by users.
    double rdz;                              ///< 1 / dz, calculated from Config_VVM given by users.
    double r2dz;                             ///< 1 / (2dz), calculated from Config_VVM given by users.
    double rdx2;                             ///< 1 / (dx^2), calculated from Config_VVM given by users.
    double rdz2;                             ///< 1 / (dz^2), calculated from Config_VVM given by users.
    int nx;                                  ///< Number of grid points in x direction, calculated from Config_VVM given by users.
    int nz;                                  ///< Number of grid points in z direction, calculated from Config_VVM given by users.
    double dt;                               ///< From Config_VVM given by users.
    double d2t;                              ///< From Config_VVM given by users.
    double dx;                               ///< From Config_VVM given by users.
    double dz;                               ///< From Config_VVM given by users.
    int XRANGE;                              ///< From Config_VVM given by users.
    int ZRANGE;                              ///< From Config_VVM given by users.
    double TIMEEND;                          ///< From Config_VVM given by users.
    int TIMEROUTPUTSIZE;                     ///< From Config_VVM given by users.
    std::string outputpath;                  ///< From Config_VVM given by users.
    int OUTPUTSTEP;                          ///< From Config_VVM given by users.
    double Kx;                               ///< From Config_VVM given by users.
    double Kz;                               ///< From Config_VVM given by users.
    double TIMETS;                           ///< From Config_VVM given by users.
    double tolerance;                        ///< From Config_VVM given by users.
    double GRAVITY;                          ///< From Config_VVM given by users.
    double Cp;                               ///< From Config_VVM given by users.
    double Cv;                               ///< From Config_VVM given by users.
    double Rd;                               ///< From Config_VVM given by users.
    double Lv;                               ///< From Config_VVM given by users.
    double P0;                               ///< From Config_VVM given by users.
    double PSURF;                            ///< From Config_VVM given by users.
    double addforcingtime;                   ///< From Config_VVM given by users.
    int CASE;                                ///< From Config_VVM given by users.
    double CRAD = 1. / 100.;                   ///< From Config_VVM given by users.

    // 0D variables
    int step = 0;                            ///< The current time step.
    double ubarTopp;                         ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.
    double ubarTop;                          ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.
    double ubarTopm;                         ///< The top boundary of the zonal wind for future time step. In the model design part, this is used to predict the mean top boundary of the zonal wind in the 9th governing equation.

    // 1D variables
    double *thb;                              ///< Horizontal mean potential temperature profile.
    double *thb_init;                         ///< Initial horizontal mean potential temperature profile.
    double *thbm;                             ///< Horizontal mean potential temperature profile for previous step.
    double *thb_zeta;                         ///< Horizontal mean potential temperature profile at grid upper edge.
    double *rhou;                             ///< Horizontal mean density profile at grid center.
    double *rhow;                             ///< Horizontal mean density profile at grid upper edge.
    double *pib;                              ///< Horizontal mean non-dimensional height profile at grid center.
    double *qvb;                              ///< Horizontal mean water vapor profile at grid center.
    double *qvsb;                             ///< Horizontal mean saturated water vapor profile at grid center.
    double *pb;                               ///< Horizontal mean pressure profile at grid center.
    double *xi;                               ///< The velocity potential in x-direction at top boundary grid center.
    double *uxi;
    double *thvb;
    double *thvbm;
    double *z;
    double *z_zeta;
    double *lambda2;

    // 2D variables
    double **zetap;
    double **zeta;
    double **zetam;
    double **thp;
    double **th;
    double **thm;
    double **u;
    double **w;
    double **RKM;
    double **RKH;
    double **U_w;
    double **W_u;


    double *zetapcont;
    double *zetacont;
    double *zetamcont;
    double *thpcont;
    double *thcont;
    double *thmcont;
    double *ucont;
    double *wcont;
    double *init_th_forcingcont;
    double *RKMcont;
    double *RKHcont;
    double *U_wcont;
    double *W_ucont;
    
    
    #if defined(STREAMFUNCTION)
        double** psi;
    #endif

    #if defined(WATER)
        double **qvp;
        double **qv;
        double **qvm;
        double **qcp;
        double **qc;
        double **qcm;
        double **qrp;
        double **qr;
        double **qrm;
        double **evaporation;
        double **accretion;
        double **autoconversion;
        double **condensation;
        double *precip;
        
        double *qvpcont;
        double *qvcont;
        double *qvmcont;
        double *qcpcont;
        double *qccont;
        double *qcmcont;
        double *qrpcont;
        double *qrcont;
        double *qrmcont;
        double *evaporationcont;
        double *accretioncont;
        double *autoconversioncont;
        double *condensationcont;
    #endif

    // #####################################################################################
    // Used for AB2. These variables are declared but not initialized if it's not AB2
    double ***dth_advect;
    double ***dth_buoyancy;
    double ***dzeta_advect;
    
    double *dth_advectcont;
    double *dth_buoyancycont;
    double *dzeta_advectcont;

    #if defined(WATER)
        double ***dqv_advect;
        double ***dqc_advect;
        double ***dqr_advect;
        double ***dqr_VT;

        double *dqv_advectcont;
        double *dqc_advectcont;
        double *dqr_advectcont;
        double *dqr_VTcont;
    #endif
    // #####################################################################################


    double *t_advection;
    double *t_poisson;
    double *t_diffusion;
    double *t_microphysics;
    double *t_all;

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
     * @param model the vvm object which is the model that will be used to do the diffusion (mainly the wind and the grid info).
     */
    // static void Advection_thermo(double **previous, double **now, double **future, vvm &model);
    static void Advection_thermo(double **past, double **now, double **future, double ***dvar, vvm &model);

    static void Advection_qrVT(vvm &model);
    // *********************************************************************************

    static void Bouyancy(vvm &model);


    // Poisson Solver => PoissonSolver.cpp
    // *********************************************************************************
    #ifndef PETSC
        Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>((nx-2)*(nz-3), (nx-2)*(nz-3));
        Eigen::SparseMatrix<double> G = Eigen::SparseMatrix<double>(nx-2, nx-2);
    #endif

    class PoissonSolver {
    public:
        #ifndef PETSC
            Eigen::SparseMatrix<double> A;
            Eigen::SparseMatrix<double> G;
            static void InitPoissonMatrix(vvm &model);
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
        static void Diffusion(double **var_in, double **var_out, vvm &model);
        static void DiffusionAll(vvm &model);
        static void TimeFilter(double **previous, double **now, double **future, vvm &model);
        static void timeFilterAll(vvm &model);
        static void Nudge_theta(vvm &model);
        static void Nudge_zeta(vvm &model);
    };

    // *********************************************************************************

    #if defined(WATER)
    class MicroPhysics {
    public:
        static void condensation(vvm &model); 	// condensation of qc by qv
        static void autoconversion(vvm &model); 	// autoconversion of qc to qr
        static void accretion(vvm &model); 		// accretion of qc by qr
        static void evaporation(vvm &model); 	// evaporation of rain water
        static void NegativeValueProcess(double **var, int nx, int nz);
    };

        static void AddForcing(vvm &model);
    #endif

    // Variables for tropical forcing
    double** init_th_forcing;
    bool status_for_adding_forcing = true;
    #if defined(TROPICALFORCING)
        double* Q1LS;
        double* Q2LS;
    #endif

    class Init {
    public:
        static void Init1d(vvm &model);
        static void Init2d(vvm &model);
        static void RandomPerturbation(vvm &, int seed);
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
            static double GetQVB(int k, int dz);
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
            static void pqv_pt(vvm &model);
            static void pqc_pt(vvm &model);
            static void pqr_pt(vvm &model);
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
    };

};

