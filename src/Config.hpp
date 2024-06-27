// #define OUTPUTNC
#define OUTPUTTXT

// #define PETSC
#define AB2 // More stable than Leapfrog (default)
// #define STREAMFUNCTION // Don't turn this on. Debugging!!!!!
#if defined(AB2)
    #define ALPHA (1.)
#else
    #define ALPHA (0.)
#endif
#define DIFFUSION
#define TIMEFILTER

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
    \def dx
    The grid size in x direction (m)
    \def dz
    The grid size in z direction (m)
    \def XRANGE
    The domain size for this model in x direction (m)
    \def ZRANGE
    The domain size for this model in z direction (m)
    \def dt
    The integration time step (s)
    \def TIMEEND
    The total integration time (s)
    \def OUTPUTPATH
    The output file path
    \def OUTPUTSTEP
    The time step for output
    \def OUTPUTNC
    The switch for output file in .nc datatype
    \def OUTPUTTXT
    The switch for output file in .txt datatype
    \def PETSC
    The switch for Poisson Solver PETSc
    \def AB2
    The switch for Adams-Bashforth Numerical Method
    \def Diffusion
    The switch for the diffusion process
    \def Kx
    The diffusion coefficient in x direction
    \def Kz
    The diffusion coefficient in z direction
    \def POISSONPARAM
    When solving 2D Poisson equation, adding small values to main diagonal to make the matrix solvable.
    \def POISSONPARAMU
    When solving 1D Poisson equation, adding small values to main diagonal to make the matrix solvable.
    \def TIMETS
    The coefficient for removing the computational mode
    \def DRY
    Switch for constant potential temperature (300K)
    \def RHO1
    Switch for constant density profile (1kg/m3)
    \def WATER
    Switch for microphysics (qc, qv, qr) with warm rain scheme
    \def TROPICALFORCING
    Switch for tropical forcing test case (With Q1, Q2)
    \def ADDFORCINGTIME
    Time (s) for adding random perturbation into model
    \def LOADFILE
    Switch for Loading file for tropical forcing. Turn it on when you use tropical forcing test
    \def STREAMFUNCTION
    Solving stream function (Krurger 1988) to get u, w rather than solving them directly (Jung 2008)
    \def VTconst
    In the Kessler warm rain scheme, set the terminal velocity to a constant (here is 6) or not. If this is turned off the terminal velocity will be calculated according to an imperical function
    \def POISSONTEST
    Switch for testing the Poisson matrix
*/
// ***********************************************************************************
