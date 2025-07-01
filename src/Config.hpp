#define OUTPUTNC
#define OUTPUTNCSAMEFILE
// #define OUTPUTMICROPHYSICS
// #define OUTPUTTXT

#define AB2
// #define PETSC
// #define STREAMFUNCTION // Don't turn this on. Debugging!!!!!
// #define DIFFUSION_VVM
// #define TIMEFILTER

// #define DRY
// #define RHO1
#define WATER
// #define KESSLER_MICROPHY
#define P3_MICROPHY
// #define RTERRTMGP
// #define TROPICALFORCING
// #define RADIATIONCOOLING
#define LOADFILE
#if defined(LOADFILE)
    #define LOADINITPATH "/work/aaron900129/NextACC/2DVVM/input/bubble_init.txt"
#endif
#define GPU_POISSON

// #define LOADFROMPREVIOUSFILE
#if defined(LOADFROMPREVIOUSFILE)
    #define TIMENOW (1)
    #define LOADINITPATH "/home/Aaron/TMIF_VVM_CSSWM/2DVVM/input/init.txt"
    #define LOADPATH1 "/data/Aaron/TMIF/0613_test/nc/10000.nc"
    #define LOADPATH2 "/data/Aaron/TMIF/0613_test/nc/10005.nc"
#endif

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


/** 
 * \def OUTPUTNC
 * The flag for output file in .nc datatype. Note that this can't work for openmp so turn it off if you want to use openmp.
 */

/** 
 * \def OUTPUTTXT
 * The flag for output file in .txt datatype. This is the default output format.
 */

/** 
 * \def PETSC
 * The flag for Poisson Solver PETSc. If you want to use Eigen solver (default), turn it off.
 */

/** 
 * \def AB2
 * The flag for Adams-Bashforth Numerical Method. This method is more stable because the spatial discretization is third-order accurate. If you want to use Leapfrog, turn it off. 
 */

/** 
 * \def DIFFUSION
 * The flag for the diffusion process. If the flag is turned off, the first-order turbulent closure is used.
 */

/** 
 * \def TIMEFILTER
 * The flag for the time filter process. Turn this on if you want to use the Leapfrog method. This will filter out the computational mode.
 */

/** 
 * \def DRY
 * The flag for constant potential temperature (300K).
 */

/** 
 * \def RHO1
 * The flag for constant density profile (1kg/m^3).
 */

/** 
 * \def WATER
 * The flag for microphysics (qc, qv, qr) with warm rain scheme.
 */

/** 
 * \def TROPICALFORCING
 * The flag for tropical forcing test case (With Q1, Q2).
 */

/** 
 * \def ADDFORCINGTIME
 * Time (s) for adding random perturbation into the model.
 */

/** 
 * \def LOADFILE
 * The flag for Loading file for tropical forcing case. Turn it on when you use tropical forcing test.
 */

/** 
 * \def LOADFROMPREVIOUSFILE
 * The flag for restarting the model. You can specify the file path (for two time steps) you want to load.
 */

/** 
 * \def STREAMFUNCTION
 * Solving stream function (Kruger 1988) to get u, w rather than solving them directly (default, Jung 2008).
 */

/** 
 * \def POISSONTEST
 * The flag for testing the Poisson matrix.
 */

// ***********************************************************************************
