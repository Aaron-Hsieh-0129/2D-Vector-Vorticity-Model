// #define OUTPUTNC // output in .nc format. This can't work for openmp so turn it off if you want to use openmp
#define OUTPUTTXT

#define AB2 // More stable than Leapfrog (default)
// #define PETSC // Using PETSc to solve Poisson equation (default is Eigen solver)
// #define STREAMFUNCTION // Don't turn this on. Debugging!!!!!
#if defined(AB2)
    #define ALPHA (1.)
#else
    #define ALPHA (0.)
#endif
// #define DIFFUSION
// #define TIMEFILTER

// #define DRY
// #define RHO1
#define WATER
// #define TROPICALFORCING
// #define RADIATIONCOOLING
// #define LOADFILE
#if defined(LOADFILE)
    #define LOADINITPATH "/home/Aaron/TMIF_VVM_CSSWM/2DVVM/input/init.txt"
#endif

// #define LOADFROMPREVIOUSFILE
#if defined(LOADFROMPREVIOUSFILE)
    #define TIMENOW (1)
    #define LOADINITPATH "/home/Aaron/TMIF_VVM_CSSWM/2DVVM/input/init.txt"
    #define LOADPATH1 "/data/Aaron/TMIF/0613_test/nc/10000.nc"
    #define LOADPATH2 "/data/Aaron/TMIF/0613_test/nc/10005.nc"
#endif

// #define POISSONTEST

// ***********************************************************************************
// Documentation Part
/**
 * @file Config.hpp
 * @author Aaron Hsieh (b08209006@ntu.edu.tw)
 * @brief 
 * @version 0.1
 * @date 2024-03-13
 * 
 * @copyright Copyright (c) 2024
 * @brief Here are the parameters for the model that you can tune them to fit your needs.
 * 
*/

/*!
    \def OUTPUTNC
    The switch for output file in .nc datatype. Note that this can't work for openmp so turn it off if you want to use openmp.
    \def OUTPUTTXT
    The switch for output file in .txt datatype. This is the default output format.
    \def PETSC
    The switch for Poisson Solver PETSc. If you want to use Eigen solver (default), turn it off.
    \def AB2
    The switch for Adams-Bashforth Numerical Method. This method is more stable because the spatial discretization is third-order accurate. If you want to use Leapfrog, turn it off. 
    \def DIFFUSION
    The switch for the diffusion process.
    If the flag is turned off, the first-order turbulent closure is used.
    \def TIMEFILTER
    The switch for the time filter process.
    Turn this on if you want to use Leapfrog method. This will filter out the computational mode.
    \def DRY
    Switch for constant potential temperature (300K)
    \def RHO1
    Switch for constant density profile (1kg/m^3)
    \def WATER
    Switch for microphysics (qc, qv, qr) with warm rain scheme
    \def TROPICALFORCING
    Switch for tropical forcing test case (With Q1, Q2)
    \def ADDFORCINGTIME
    Time (s) for adding random perturbation into model
    \def LOADFILE
    Switch for Loading file for tropical forcing. Turn it on when you use tropical forcing test
    \def LOADFROMPREVIOUSFILE
    Switch for restarting the model. You can specify the file path (for two time steps) you want to load 
    \def STREAMFUNCTION
    Solving stream function (Krurger 1988) to get u, w rather than solving them directly (default, Jung 2008)
    \def POISSONTEST
    Switch for testing the Poisson matrix
*/
// ***********************************************************************************
