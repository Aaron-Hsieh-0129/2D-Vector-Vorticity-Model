#include "Const.hpp"
#include <string>
#ifndef PETSC
    #include "../include/Eigen/Sparse"
#endif

class Config_VVM {
public:
    double dt;
    double dx;
    double dz;
    int XRANGE;
    int ZRANGE;
    double TIMEEND;
    int TIMEROUTPUTSIZE;
    std::string OUTPUTPATH; 
    int OUTPUTSTEP;
    double Kx;
    double Kz;
    double TIMETS;
    double POISSONPARAMU;
    double POISSONPARAMW;
    double GRAVITY;
    double Cp;
    double Cv;
    double Rd;
    double Lv;
    double P0;
    double PSURF;
    double ADDFORCINGTIME;

    Config_VVM(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, 
           std::string OUTPUTPATH, int OUTPUTSTEP, double Kx, double Kz, double TIMETS, double POISSONPARAMU, double POISSONPARAMW, 
           double GRAVITY, double Cp, double Cv, double Rd, double Lv, double P0, double PSURF, double ADDFORCINGTIME)
        : dt(dt), dx(dx), dz(dz), XRANGE(XRANGE+2*dx), ZRANGE(ZRANGE+2*dz), TIMEEND(TIMEEND), TIMEROUTPUTSIZE(TIMEROUTPUTSIZE), 
          OUTPUTPATH(OUTPUTPATH), OUTPUTSTEP(OUTPUTSTEP), Kx(Kx), Kz(Kz), TIMETS(TIMETS), POISSONPARAMU(POISSONPARAMU), POISSONPARAMW(POISSONPARAMW), 
          GRAVITY(GRAVITY), Cp(Cp), Cv(Cv), Rd(Rd), Lv(Lv), P0(P0), PSURF(PSURF), ADDFORCINGTIME(ADDFORCINGTIME) {}
    ~Config_VVM() {}
};


class vvm {
public:
    /**
     * vvm constructor.
     * Used to initialize the model.
     */
    vvm(const Config_VVM& config)
        : dt(config.dt), dx(config.dx), dz(config.dz), 
          XRANGE(config.XRANGE), ZRANGE(config.ZRANGE),
          TIMEEND(config.TIMEEND), TIMEROUTPUTSIZE(config.TIMEROUTPUTSIZE),
          OUTPUTPATH(config.OUTPUTPATH), OUTPUTSTEP(config.OUTPUTSTEP), 
          Kx(config.Kx), Kz(config.Kz), TIMETS(config.TIMETS),
          POISSONPARAMU(config.POISSONPARAMU), POISSONPARAMW(config.POISSONPARAMW),
          GRAVITY(config.GRAVITY), Cp(config.Cp), Cv(config.Cv), Rd(config.Rd), Lv(config.Lv),
          P0(config.P0), PSURF(config.PSURF), ADDFORCINGTIME(config.ADDFORCINGTIME),
          d2t(2.0 * config.dt),
          rdx(1.0 / config.dx), r2dx(rdx / 2.0),
          rdz(1.0 / config.dz), r2dz(rdz / 2.0),
          rdx2(rdx * rdx), rdz2(rdz * rdz), nx(config.XRANGE/config.dx), nz(config.ZRANGE/config.dz) {

        allocateMemory();
        initializaArrays();
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

        // for (int i = 0; i < nx; ++i) {
        //     delete[] zetap[i];
        //     delete[] zeta[i];
        //     delete[] zetam[i];
        //     delete[] thp[i];
        //     delete[] th[i];
        //     delete[] thm[i];
        //     delete[] u[i];
        //     delete[] w[i];
        //     delete[] um[i];
        //     delete[] wm[i];
        // }
        delete[] zetap;
        delete[] zeta;
        delete[] zetam;
        delete[] thp;
        delete[] th;
        delete[] thm;
        delete[] u;
        delete[] w;
        delete[] um;
        delete[] wm;

        delete[] zetapcont;
        delete[] zetacont;
        delete[] zetamcont;
        delete[] thpcont;
        delete[] thcont;
        delete[] thmcont;
        delete[] ucont;
        delete[] wcont;
        delete[] umcont;
        delete[] wmcont;


        #if defined(STREAMFUNCTION)
            for (int i = 0; i < nx; ++i) {
                delete[] psi[i];
            }
            delete[] psi;
            delete[] psicont;
        #endif

        #if defined(WATER)
            // for (int i = 0; i < nx; ++i) {
            //     delete[] qvp[i];
            //     delete[] qv[i];
            //     delete[] qvm[i];
            //     delete[] qcp[i];
            //     delete[] qc[i];
            //     delete[] qcm[i];
            //     delete[] qrp[i];
            //     delete[] qr[i];
            //     delete[] qrm[i];
            //     delete[] evaporation[i];
            //     delete[] accretion[i];
            //     delete[] autoconversion[i];
            // }
            delete[] qvp;
            delete[] qv;
            delete[] qvm;
            delete[] qcp;
            delete[] qc;
            delete[] qcm;
            delete[] qrp;
            delete[] qr;
            delete[] qrm;
            delete[] evaporation;
            delete[] accretion;
            delete[] autoconversion;
            delete[] precip;

            delete[] qvpcont;
            delete[] qvcont;
            delete[] qvmcont;
            delete[] qcpcont;
            delete[] qccont;
            delete[] qcmcont;
            delete[] qrpcont;
            delete[] qrcont;
            delete[] qrmcont;
            delete[] evaporationcont;
            delete[] accretioncont;
            delete[] autoconversioncont;
        #endif

        delete[] t_advection;
        delete[] t_poisson;
        delete[] t_diffusion;
        delete[] t_microphysics;
        delete[] t_all;

        #if defined(TROPICALFORCING)
            delete[] Q1LS;
            delete[] Q2LS;
            // for (int i = 0; i < nx; ++i) {
            //     delete[] init_th_forcing[i];
            // }
            delete[] init_th_forcing;
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
        #if defined(TROPICALFORCING)
            Q1LS = new double[nz];
            Q2LS = new double[nz];
            init_th_forcing = allocate2DArray(nx, nz);
        #endif

        // 2D arrays
        zetap = allocate2DArray(nx, nz);
        zeta = allocate2DArray(nx, nz);
        zetam = allocate2DArray(nx, nz);
        thp = allocate2DArray(nx, nz);
        th = allocate2DArray(nx, nz);
        thm = allocate2DArray(nx, nz);
        u = allocate2DArray(nx, nz);
        w = allocate2DArray(nx, nz);
        um = allocate2DArray(nx, nz);
        wm = allocate2DArray(nx, nz);

        zetapcont = new double[nx*nz];
        zetacont = new double[nx*nz];
        zetamcont = new double[nx*nz];
        thpcont = new double[nx*nz];
        thcont = new double[nx*nz];
        thmcont = new double[nx*nz];
        ucont = new double[nx*nz];
        wcont = new double[nx*nz];
        umcont = new double[nx*nz];
        wmcont = new double[nx*nz];

        #if defined(STREAMFUNCTION)
            psi = allocate2DArray(nx, nz);
            psicont = new double[nx*nz];
        #endif

        #if defined(WATER)
            qvp = allocate2DArray(nx, nz);
            qv = allocate2DArray(nx, nz);
            qvm = allocate2DArray(nx, nz);
            qcp = allocate2DArray(nx, nz);
            qc = allocate2DArray(nx, nz);
            qcm = allocate2DArray(nx, nz);
            qrp = allocate2DArray(nx, nz);
            qr = allocate2DArray(nx, nz);
            qrm = allocate2DArray(nx, nz);
            evaporation = allocate2DArray(nx, nz);
            accretion = allocate2DArray(nx, nz);
            autoconversion = allocate2DArray(nx, nz);

            qvpcont = new double[nx*nz];
            qvcont = new double[nx*nz];
            qvmcont = new double[nx*nz];
            qcpcont = new double[nx*nz];
            qccont = new double[nx*nz];
            qcmcont = new double[nx*nz];
            qrpcont = new double[nx*nz];
            qrcont = new double[nx*nz];
            qrmcont = new double[nx*nz];
            evaporationcont = new double[nx*nz];
            accretioncont = new double[nx*nz];
            autoconversioncont = new double[nx*nz];
        #endif
    }

    void initializaArrays() {
        for (int i = 0; i < nx; i++) {
            zetap[i] = &zetapcont[i*nz];
            zeta[i] = &zetacont[i*nz];
            zetam[i] = &zetamcont[i*nz];
            thp[i] = &thpcont[i*nz];
            th[i] = &thcont[i*nz];
            thm[i] = &thmcont[i*nz];
            u[i] = &ucont[i*nz];
            w[i] = &wcont[i*nz];
            um[i] = &umcont[i*nz];
            wm[i] = &wmcont[i*nz];
            #if defined(STREAMFUNCTION)
                psi[i] = &psicont[i*nz];
            #endif

            #if defined(WATER)
                qvp[i] = &qvpcont[i*nz];
                qv[i] = &qvcont[i*nz];
                qvm[i] = &qvmcont[i*nz];
                qcp[i] = &qcpcont[i*nz];
                qc[i] = &qccont[i*nz];
                qcm[i] = &qcmcont[i*nz];
                qrp[i] = &qrpcont[i*nz];
                qr[i] = &qrcont[i*nz];
                qrm[i] = &qrmcont[i*nz];
                evaporation[i] = &evaporationcont[i*nz];
                accretion[i] = &accretioncont[i*nz];
                autoconversion[i] = &autoconversioncont[i*nz];
            #endif
        }
        #if defined(WATER)
            precip = new double[nx];
        #endif
    }

    double** allocate2DArray(int dim1, int dim2) {
        double** array = new double*[dim1];
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double[dim2];
        }
        return array;
    }

    double rdx;                                               ///< 1 / dx
    double r2dx;                                              ///< 1 / (2dx)
    double rdz;                                               ///< 1 / dz
    double r2dz;                                              ///< 1 / (2dz)
    double rdx2;                                              ///< 1 / (dx^2)
    double rdz2;                                              ///< 1 / (dz^2)
    int nx;                                                   ///< Number of grid points in x direction
    int nz;                                                   ///< Number of grid points in z direction
    double dt;
    double d2t;                                               ///< Integration time step
    double dx;
    double dz;
    int XRANGE;
    int ZRANGE;
    double TIMEEND;
    int TIMEROUTPUTSIZE;
    std::string OUTPUTPATH;
    int OUTPUTSTEP;
    double Kx;
    double Kz;
    double TIMETS;
    double POISSONPARAMU;
    double POISSONPARAMW;
    double GRAVITY;
    double Cp;
    double Cv;
    double Rd;
    double Lv;
    double P0;
    double PSURF;
    double ADDFORCINGTIME;

    // 0D variables
    double ubarTopp, ubarTop, ubarTopm;

    // 1D variables
    double *thb;
    double *thbm;
    double *thb_zeta;
    double *rhou;
    double *rhow;
    double *pib;
    double *qvb;
    double *qvsb;
    double *pb;
    double *xi;
    double *uxi;
    double *thvb;
    double *thvbm;

    // 2D variables
    double **zetap;
    double **zeta;
    double **zetam;
    double **thp;
    double **th;
    double **thm;
    double **u;
    double **w;
    double **um;
    double **wm;

    double *zetapcont;
    double *zetacont;
    double *zetamcont;
    double *thpcont;
    double *thcont;
    double *thmcont;
    double *ucont;
    double *wcont;
    double *umcont;
    double *wmcont;
    
    
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
    #endif

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

    static void BoundaryProcess2D_all(vvm &);
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
    static void Advection_thermo(double **previous, double **now, double **future, vvm &model);

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
            static void InitPoissonMatrix(vvm &);
        #endif
        #if defined(STREAMFUNCTION)
            static void calpsiuw(vvm &);
        #else
            static void cal_w(vvm &);
            static void cal_u(vvm &);
            static void pubarTop_pt(vvm &);
        #endif
    };

    // *********************************************************************************

    // Diffusion and Time Filter => NumericalProcess.cpp
    // *********************************************************************************
    class NumericalProcess {
    public:
        static void Diffusion(double **var_in, double **var_out, vvm &model);
        static void DiffusionAll(vvm &);
        static void TimeFilter(double **previous, double **now, double **future, vvm &model);
        static void timeFilterAll(vvm &);
    };

    // *********************************************************************************

    #if defined(WATER)
    class MicroPhysics {
    public:
        static void condensation(vvm &); 	// condensation of qc by qv
        static void autoconversion(vvm &); 	// autoconversion of qc to qr
        static void accretion(vvm &); 		// accretion of qc by qr
        static void evaporation(vvm &); 	// evaporation of rain water
        static void NegativeValueProcess(double **var, int nx, int nz);
    };

        #if defined(TROPICALFORCING)
            static void AddForcing(vvm &model);
        #endif
    #endif

    // Variables for tropical forcing
    #if defined(TROPICALFORCING)
        double* Q1LS;
        double* Q2LS;
        double** init_th_forcing;
        bool status_for_adding_forcing = false;
    #endif

    class Init {
    public:
        static void Init1d(vvm &);
        static void Init2d(vvm &);
        #if defined(LOADFILE)
            static void LoadFile(vvm &);
        #elif defined(LOADFROMPREVIOUSFILE)
            static void LoadFromPreviousFile(vvm &);
        #endif
        #if defined(TROPICALFORCING)
            static void RandomPerturbation(vvm &, int);
        #endif
            
    private:
        static double GetTB(int i, vvm &);
        static double GetTHRAD(int i, int k, vvm &);
        static double GetTH(int, int, vvm &);
        #if defined(WATER)
            static double GetQVB(int, int dz);
        #endif
    };

    class Output {
    public:
        static void printInit(vvm &);
        static void create_all_directory(vvm &);
        #if defined(OUTPUTNC)
            static void output_nc(int, vvm &);
            static void output_time_nc(int, vvm &);
        #endif

        #if defined(OUTPUTTXT)
            static void output_zeta(int, vvm &);
            static void output_th(int, vvm &);
            static void output_u(int, vvm &);
            static void output_w(int, vvm &);
            #if defined(WATER)
                static void output_qv(int, vvm &);
                static void output_qc(int, vvm &);
                static void output_qr(int, vvm &);
                static void output_precip(int, vvm &);
            #endif
            static void outputalltxt(int, vvm &);
        #endif
        
    private:
        static void create_directory(std::string);
    };


    class Iteration {
    public:
        static void pzeta_pt(vvm &);
        static void pth_pt(vvm &);
        #if defined(WATER)
            static void pqv_pt(vvm &);
            static void pqc_pt(vvm &);
            static void pqr_pt(vvm &);
        #endif

        static void updateMean(vvm &);
        static void TimeMarching(vvm &);
        static void nextTimeStep(vvm &);
    };

};

