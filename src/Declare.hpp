#include "Config.hpp"
#include <string>
#ifndef PETSC
    // #include "../include/Eigen/Sparse"
    #include <Eigen/Sparse>
#endif
#if defined(GPU_POISSON)
    #include <amgx_c.h>
#endif

class Config_VVM {
public:
    Config_VVM(double dt, double dx, double dz, int XRANGE, int ZRANGE, double TIMEEND, int TIMEROUTPUTSIZE, 
           std::string outputpath, int OUTPUTSTEP, double Kx, double Kz, double TIMETS, double tolerance,
           double GRAVITY, double Cp, double Cv, double Rd, double Lv, double P0, double PSURF, double addforcingtime, int CASE, double mositure_nudge_time)
        : dt(dt), dx(dx), dz(dz), XRANGE(XRANGE+2*dx), ZRANGE(ZRANGE+2*dz), TIMEEND(TIMEEND), TIMEROUTPUTSIZE(TIMEROUTPUTSIZE), 
          outputpath(outputpath), OUTPUTSTEP(OUTPUTSTEP), Kx(Kx), Kz(Kz), TIMETS(TIMETS),
          tolerance(tolerance), GRAVITY(GRAVITY), Cp(Cp), Cv(Cv), Rd(Rd), Lv(Lv), P0(P0), PSURF(PSURF), addforcingtime(addforcingtime), CASE(CASE), mositure_nudge_time(mositure_nudge_time) {}
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
    double mositure_nudge_time; ///< The time for nudging the moisture. It's used for nudging the moisture to the mean state. It's used for nudging the moisture to the mean state.
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
          P0(config.P0), PSURF(config.PSURF), addforcingtime(config.addforcingtime), CASE(config.CASE), moisture_nudge_time(config.mositure_nudge_time)
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
        delete[] qvb0;
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
            

            deallocate2DContinuousArray(qvp, qvpcont);
            deallocate2DContinuousArray(qv, qvcont);
            deallocate2DContinuousArray(qvm, qvmcont);
            deallocate2DContinuousArray(qcp, qcpcont);
            deallocate2DContinuousArray(qc, qccont);
            deallocate2DContinuousArray(qcm, qcmcont);
            deallocate2DContinuousArray(qrp, qrpcont);
            deallocate2DContinuousArray(qr, qrcont);
            deallocate2DContinuousArray(qrm, qrmcont);
            #if defined(KESSLER_MICROPHY)
                delete[] precip;

                deallocate2DContinuousArray(evaporation, evaporationcont);
                deallocate2DContinuousArray(accretion, accretioncont);
                deallocate2DContinuousArray(autoconversion, autoconversioncont);
                deallocate2DContinuousArray(condensation, condensationcont);
            #endif

            #if defined(P3_MICROPHY)
                
                deallocate2DContinuousArray(ncp, ncpcont);
                deallocate2DContinuousArray(nc, nccont);
                deallocate2DContinuousArray(ncm, ncmcont);
                deallocate2DContinuousArray(nrp, nrpcont);
                deallocate2DContinuousArray(nr, nrcont);
                deallocate2DContinuousArray(nrm, nrmcont);
                deallocate2DContinuousArray(qitotp, qitotpcont);
                deallocate2DContinuousArray(qitot, qitotcont);
                deallocate2DContinuousArray(qitotm, qitotmcont);
                deallocate2DContinuousArray(qirimp, qirimpcont);
                deallocate2DContinuousArray(qirim, qirimcont);
                deallocate2DContinuousArray(qirimm, qirimmcont);
                deallocate2DContinuousArray(qiliqp, qiliqpcont);
                deallocate2DContinuousArray(qiliq, qiliqcont);
                deallocate2DContinuousArray(qiliqm, qiliqmcont);
                deallocate2DContinuousArray(nip, nipcont);
                deallocate2DContinuousArray(ni, nicont);
                deallocate2DContinuousArray(nim, nimcont);
                deallocate2DContinuousArray(birimp, birimpcont);
                deallocate2DContinuousArray(birim, birimcont);
                deallocate2DContinuousArray(birimm, birimmcont);
                deallocate2DContinuousArray(diag_ze, diag_zecont);
                deallocate2DContinuousArray(diag_effc, diag_effccont);
                deallocate2DContinuousArray(diag_effi, diag_efficont);
                deallocate2DContinuousArray(diag_vmi, diag_vmicont);
                deallocate2DContinuousArray(diag_di, diag_dicont);
                deallocate2DContinuousArray(diag_rhoi, diag_rhoicont);
                deallocate2DContinuousArray(cldfrac, cldfraccont);
                deallocate2DContinuousArray(diag_2d, diag_2dcont);
                deallocate2DContinuousArray(dz_all, dz_allcont);
                deallocate2DContinuousArray(w_all, w_allcont);
                deallocate2DContinuousArray(pb_all, pb_allcont);
                deallocate2DContinuousArray(zi_all, zi_allcont);
                deallocate2DContinuousArray(ssat_all, ssat_allcont);

                delete[] precip_liq;
                delete[] precip_sol;
                delete[] precip;
            #endif
        #endif

        #if defined(AB2)
            deallocate3DContinuousArray(dth_advect, dth_advectcont);
            deallocate3DContinuousArray(dth_buoyancy, dth_buoyancycont);
            deallocate3DContinuousArray(dzeta_advect, dzeta_advectcont);
            
            #if defined(WATER)
                deallocate3DContinuousArray(dqv_advect, dqv_advectcont);
                deallocate3DContinuousArray(dqc_advect, dqc_advectcont);
                deallocate3DContinuousArray(dqr_advect, dqr_advectcont);
                #if defined(KESSLER_MICROPHY)
                    deallocate3DContinuousArray(dqr_VT, dqr_VTcont);
                #endif
                #if defined(P3_MICROPHY)
                    
                    deallocate3DContinuousArray(dnc_advect, dnc_advectcont);
                    deallocate3DContinuousArray(dnr_advect, dnr_advectcont);
                    deallocate3DContinuousArray(dni_advect, dni_advectcont);
                    deallocate3DContinuousArray(dqitot_advect, dqitot_advectcont);
                    deallocate3DContinuousArray(dqirim_advect, dqirim_advectcont);
                    deallocate3DContinuousArray(dqiliq_advect, dqiliq_advectcont);
                    deallocate3DContinuousArray(dbirim_advect, dbirim_advectcont);
                #endif
            #endif

        #endif
        #if defined(P3_MICROPHY)
            deallocate3DContinuousArray(diag_3d, diag_3dcont);
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
        t_advection = new double[TIMEROUTPUTSIZE]();
        t_poisson = new double[TIMEROUTPUTSIZE]();
        t_diffusion = new double[TIMEROUTPUTSIZE]();
        t_microphysics = new double[TIMEROUTPUTSIZE]();
        t_all = new double[TIMEROUTPUTSIZE]();

        // 1D arrays
        thb = new double[nz]();
        thbm = new double[nz]();
        thb_zeta = new double[nz]();
        thb_init = new double[nz]();
        rhou = new double[nz]();
        rhow = new double[nz]();
        pib = new double[nz]();
        qvb = new double[nz]();
        qvb0 = new double[nz]();
        qvsb = new double[nz]();
        pb = new double[nz]();
        xi = new double[nx]();
        uxi = new double[nx]();
        thvb = new double[nz]();
        thvbm = new double[nz]();
        z = new double[nz]();
        z_zeta = new double[nz]();
        lambda2 = new double[nz]();
        #if defined(TROPICALFORCING)
            Q1LS = new double[nz]();
            Q2LS = new double[nz]();
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
            #if defined(KESSLER_MICROPHY)
                precip = new double[nx]();

                evaporation = allocate2DContinuousArray(nx, nz, evaporationcont);
                accretion = allocate2DContinuousArray(nx, nz, accretioncont);
                autoconversion = allocate2DContinuousArray(nx, nz, autoconversioncont);
                condensation = allocate2DContinuousArray(nx, nz, condensationcont);
            #endif

            qvp = allocate2DContinuousArray(nx, nz, qvpcont);
            qv = allocate2DContinuousArray(nx, nz, qvcont);
            qvm = allocate2DContinuousArray(nx, nz, qvmcont);
            qcp = allocate2DContinuousArray(nx, nz, qcpcont);
            qc = allocate2DContinuousArray(nx, nz, qccont);
            qcm = allocate2DContinuousArray(nx, nz, qcmcont);
            qrp = allocate2DContinuousArray(nx, nz, qrpcont);
            qr = allocate2DContinuousArray(nx, nz, qrcont);
            qrm = allocate2DContinuousArray(nx, nz, qrmcont);

            #if defined(P3_MICROPHY)
                ncp = allocate2DContinuousArray(nx, nz, ncpcont);
                nc = allocate2DContinuousArray(nx, nz, nccont);
                ncm = allocate2DContinuousArray(nx, nz, ncmcont);
                nrp = allocate2DContinuousArray(nx, nz, nrpcont);
                nr = allocate2DContinuousArray(nx, nz, nrcont);
                nrm = allocate2DContinuousArray(nx, nz, nrmcont);
                qitotp = allocate2DContinuousArray(nx, nz, qitotpcont);
                qitot = allocate2DContinuousArray(nx, nz, qitotcont);
                qitotm = allocate2DContinuousArray(nx, nz, qitotmcont);
                qirimp = allocate2DContinuousArray(nx, nz, qirimpcont);
                qirim = allocate2DContinuousArray(nx, nz, qirimcont);
                qirimm = allocate2DContinuousArray(nx, nz, qirimmcont);
                qiliqp = allocate2DContinuousArray(nx, nz, qiliqpcont);
                qiliq = allocate2DContinuousArray(nx, nz, qiliqcont);
                qiliqm = allocate2DContinuousArray(nx, nz, qiliqmcont);
                nip = allocate2DContinuousArray(nx, nz, nipcont);
                ni = allocate2DContinuousArray(nx, nz, nicont);
                nim = allocate2DContinuousArray(nx, nz, nimcont);
                birimp = allocate2DContinuousArray(nx, nz, birimpcont);
                birim = allocate2DContinuousArray(nx, nz, birimcont);
                birimm = allocate2DContinuousArray(nx, nz, birimmcont);
                diag_ze = allocate2DContinuousArray(nx, nz, diag_zecont);
                diag_effc = allocate2DContinuousArray(nx, nz, diag_effccont);
                diag_effi = allocate2DContinuousArray(nx, nz, diag_efficont);
                diag_vmi = allocate2DContinuousArray(nx, nz, diag_vmicont);
                diag_di = allocate2DContinuousArray(nx, nz, diag_dicont);
                diag_rhoi = allocate2DContinuousArray(nx, nz, diag_rhoicont);
                cldfrac = allocate2DContinuousArray(nx, nz, cldfraccont);
                diag_2d = allocate2DContinuousArray(nx, vvm::P3::n_diag_2d, diag_2dcont);
                dz_all = allocate2DContinuousArray(nx, nz, dz_allcont);
                w_all = allocate2DContinuousArray(nx, nz, w_allcont);
                pb_all = allocate2DContinuousArray(nx, nz, pb_allcont);
                zi_all = allocate2DContinuousArray(nx, nz, zi_allcont);
                ssat_all = allocate2DContinuousArray(nx, nz, ssat_allcont);

                precip_liq = new double[nx]();
                precip_sol = new double[nx]();
                precip = new double[nx]();
            #endif
        #endif

        #if defined(AB2)
            dth_advect = allocate3DContinuousArray(nx, nz, 2, dth_advectcont);
            dth_buoyancy = allocate3DContinuousArray(nx, nz, 2, dth_buoyancycont);
            dzeta_advect = allocate3DContinuousArray(nx, nz, 2, dzeta_advectcont);
            #if defined(WATER)
                dqv_advect = allocate3DContinuousArray(nx, nz, 2, dqv_advectcont);
                dqc_advect = allocate3DContinuousArray(nx, nz, 2, dqc_advectcont);
                dqr_advect = allocate3DContinuousArray(nx, nz, 2, dqr_advectcont);
                #if defined(KESSLER_MICROPHY)
                    dqr_VT = allocate3DContinuousArray(nx, nz, 2, dqr_VTcont);
                #endif

                #if defined(P3_MICROPHY)
                    dnc_advect = allocate3DContinuousArray(nx, nz, 2, dnc_advectcont);
                    dnr_advect = allocate3DContinuousArray(nx, nz, 2, dnr_advectcont);
                    dni_advect = allocate3DContinuousArray(nx, nz, 2, dni_advectcont);
                    dqitot_advect = allocate3DContinuousArray(nx, nz, 2, dqitot_advectcont);
                    dqirim_advect = allocate3DContinuousArray(nx, nz, 2, dqirim_advectcont);
                    dqiliq_advect = allocate3DContinuousArray(nx, nz, 2, dqiliq_advectcont);
                    dbirim_advect = allocate3DContinuousArray(nx, nz, 2, dbirim_advectcont);
                #endif
            #endif

            #if defined(P3_MICROPHY)
                diag_3d = allocate3DContinuousArray(nx, nz, vvm::P3::n_diag_3d, diag_3dcont);
            #endif
        #endif
    }

    static double** allocate2DContinuousArray(int rows, int cols, double*& contMemory) {
        double** array = new double*[rows]();
        contMemory = new double[rows * cols](); // Allocate continuous memory block
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
        double*** array = new double**[dim1]();
        contMemory = new double[dim1 * dim2 * dim3](); // Allocate continuous memory block
        for (int i = 0; i < dim1; ++i) {
            array[i] = new double*[dim2]();
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
    int XRANGE = 0;                              ///< From Config_VVM given by users.
    int ZRANGE = 0;                              ///< From Config_VVM given by users.
    double TIMEEND = 0;                          ///< From Config_VVM given by users.
    int TIMEROUTPUTSIZE = 0;                     ///< From Config_VVM given by users.
    std::string outputpath = "";                  ///< From Config_VVM given by users.
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
    double CRAD = 1. / 3600.;                   ///< From Config_VVM given by users.

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
    double *qvb = nullptr;                              ///< Horizontal mean water vapor profile at grid center.
    double *qvb0 = nullptr;                              ///< Horizontal mean water vapor profile at grid center.
    double *qvsb = nullptr;                             ///< Horizontal mean saturated water vapor profile at grid center.
    double *pb = nullptr;                               ///< Horizontal mean pressure profile at grid center.
    double *xi = nullptr;                               ///< The velocity potential in x-direction at top boundary grid center.
    double *uxi = nullptr;
    double *thvb = nullptr;
    double *thvbm = nullptr;
    double *z = nullptr;
    double *z_zeta = nullptr;
    double *lambda2 = nullptr;

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

    #if defined(WATER) && defined(KESSLER_MICROPHY)
        static void Advection_qrVT(vvm &model);
    #endif
    // *********************************************************************************

    static void Bouyancy(vvm &model);


    // Poisson Solver => PoissonSolver.cpp
    // *********************************************************************************
    #if !defined(PETSC) && !defined(GPU_POISSON) 
        Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>((nx-2)*(nz-3), (nx-2)*(nz-3));
        Eigen::SparseMatrix<double> G = Eigen::SparseMatrix<double>(nx-2, nx-2);
    #endif

    class PoissonSolver {
    public:
        #ifndef PETSC
            static void InitPoissonMatrix(vvm &model);
            #if defined(GPU_POISSON)
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
        #if defined(DIFFUSION_VVM)
            static void Diffusion(double **var_in, double **var_out, vvm &model);
            static void DiffusionAll(vvm &model);
        #endif
        #if defined(TIMERFILTER)
            static void TimeFilter(double **previous, double **now, double **future, vvm &model);
            static void timeFilterAll(vvm &model);
        #endif
        static void Nudge_theta(vvm &model);
        static void Nudge_zeta(vvm &model);
        static void Nudge_qv(vvm &model);
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
        static void NegativeValueProcess(double **var, int nx, int nz);
    };
    #endif

        static void AddForcing(vvm &model);
    #endif

    #if defined(P3_MICROPHY)
    class P3 {
    public:
        inline static char *lookup_file_dir = strdup("../lookup_tables");
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
    };

};

