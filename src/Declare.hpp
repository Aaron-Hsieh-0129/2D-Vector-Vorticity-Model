#include "Const.hpp"
#include "../include/Eigen/Sparse"


/**
 * @file Declare.hpp
 * @author Aaron Hsieh (b08209006@ntu.edu.tw)
 * @brief 
 * @version 0.1
 * @date 2024-03-13
 * 
 * @copyright Copyright (c) 2024
 * @brief This file is used to declare the class vvm.
 * 
*/


/**
 *  vvm class. This is class for the whole model. 
 *  The configuration for the model is put inside the class and some functions related to model proccessing such as advection, boundary process, etc.
 */
class vvm {
public:
    /**
     * vvm constructor.
     * Used to initialize the model.
     */
    vvm();
    double rdx = 1. / (double) dx;                                              ///< 1 / dx
    double r2dx = rdx / 2.;                                                     ///< 1 / (2dx)
    double rdz = 1. / (double) dz;                                              ///< 1 / dz
    double r2dz = rdz / 2.;                                                     ///< 1 / (2dz)
    double rdx2 = rdx * rdx;                                                    ///< 1 / (dx^2)
    double rdz2 = rdz * rdz;                                                    ///< 1 / (dz^2)
    int nx = NX;                                                                ///< Number of grid points in z direction 
    int nz = NZ;                                                                ///< Number of grid points in z direction
    double d2t = 2. * DT;                                                       ///< Integration time step

    // 0D variables
    double ubarTopp, ubarTop, ubarTopm;

    // 1D variables
    double thb[NZ], thbm[NZ], thb_zeta[NZ], rhou[NZ], rhow[NZ], pib[NZ], qvb[NZ], qvsb[NZ], pb[NZ];
    double xi[NX], uxi[NX];
    double thvb[NZ], thvbm[NZ];

    // 2D variables
    double zetap[NX][NZ], zeta[NX][NZ], zetam[NX][NZ];
    double thp[NX][NZ], th[NX][NZ], thm[NX][NZ];
    double u[NX][NZ], w[NX][NZ];
    double um[NX][NZ], wm[NX][NZ];
    #if defined(STREAMFUNCTION)
        double psi[NX][NZ];
    #endif

    #if defined(WATER)
        double qvp[NX][NZ], qv[NX][NZ], qvm[NX][NZ];
		double qcp[NX][NZ], qc[NX][NZ], qcm[NX][NZ];
		double qrp[NX][NZ], qr[NX][NZ], qrm[NX][NZ];
        double evaporation[NX][NZ], accretion[NX][NZ], autoconversion[NX][NZ];
        double qrAcc[NX];
        double precip[NX];
    #endif

    double t_advection[TIMEROUTPUTSIZE], t_poisson[TIMEROUTPUTSIZE], t_diffusion[TIMEROUTPUTSIZE], t_microphysics[TIMEROUTPUTSIZE];
    double t_all[TIMEROUTPUTSIZE];

// Boundary Process => BoundaryProcess.cpp
// **********************************************************************
    /**
     * A member function that process the boundary of the 1D array where the varibles are at the center of the grid
     * @param var an one dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess1D_center(double var[]);

    /**
     * A member function that process the boundary of the 2D array where the varibles are at the center of the grid
     * @param var an two dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess2D_center(double var[][NZ]);

    /**
     * A member function that process the boundary of the 2D array where the varibles are at the southwestern side of the grid
     * @param var an two dimensional array that should be put into boundary process.
     */
    static void BoundaryProcess2D_westdown(double var[][NZ]);
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
    static void Advection_thermo(double previous[][NZ], double now[][NZ], double future[][NZ], vvm &model);

    static void Advection_qrVT(vvm &model);
// *********************************************************************************

    static void Bouyancy(vvm &model);


// Poisson Solver => PoissonSolver.cpp
// *********************************************************************************
    
    Eigen::SparseMatrix<double> A;
    Eigen::SparseMatrix<double> G;
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
    /**
     * A member function that do advection process to the vorticity (zeta) field.
     * @param var_in an two dimensional array that should be put in to be diffused.
     * @param var_out an two dimensional array that the result will be put into.
     * @param model the vvm object which is the model that will be used to do the diffusion (mainly the grid information).
     */
    static void Diffusion(double var_in[][NZ], double var_out[][NZ], vvm &model);

    /**
     * A member function that do advection process to the vorticity (zeta) field.
     * @param previous an two dimensional array that the timestep is the previous one such as zetam, thm.
     * @param now an two dimensional array that the timestep is now such as zeta, th.
     * @param future an two dimensional array that the timestep is the future one such as zetap, thp.
     * @param model the vvm object which is the model that will be used to do the diffusion (mainly the grid information).
     */
    static void TimeFilter(double previous[][NZ], double now[][NZ], double future[][NZ], vvm &model);
// *********************************************************************************

    #if defined(WATER)
    class MicroPhysics {
    public:
        static void condensation(vvm &); 	// condensation of qc by qv
        static void autoconversion(vvm &); 	// autoconversion of qc to qr
        static void accretion(vvm &); 		// accretion of qc by qr
        static void evaporation(vvm &); 	// evaporation of rain water
        static void NegativeValueProcess(double var[][NZ]);
    };

        #if defined(TROPICALFORCING)
            static void AddForcing(vvm &model);
        #endif
    #endif

    // Variables for tropical forcing
    #if defined(TROPICALFORCING)
        double Q1LS[NZ], Q2LS[NZ];
        double init_th_forcing[NX][NZ];
        bool status_for_adding_forcing = false;
    #endif
};
