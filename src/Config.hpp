#define dx (200)
#define dz (200)
#define XRANGE (120000+2*dx)
#define ZRANGE (20000+2*dz)
#define DT (1.)
#define TIMEEND (10000.)

#define OUTPUTPATH "/data/Aaron/2DVVM/3rd/Water/"
#define OUTPUTSTEP (1)
#define OUTPUTNC

#define PETSC
// #define AB3
#define DIFFUSION
#define Kx (250.)
#define Kz (250.)
#define POISSONPARAM (1E-9)
#define POISSONPARAMU (1E-9)
#ifndef AB3
    #define TIMEFILTER
#endif
#define TIMETS (0.01)

// #define DRY
// #define RHO1
#define WATER
#define TROPICALFORCING
#define ADDFORCINGTIME (1200.)
#define RADIATIONCOOLING
#define LOADFILE
// #define STREAMFUNCTION

#define VTconst

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
 * @brief Here are the parameters for the model that you can tune to fit your needs.
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
    The switch for output file
    \def Kx
    The diffusion coefficient in x direction
    \def Kz
    The diffusion coefficient in z direction
    \def TIMETS
    The coefficient for removing the computational mode
    \def DRY
    Switch for constant potential temperature (300K)
    \def TROPICALFORCING
    Switch for tropical forcing test case (With Q1, Q2)
    \def LOADFILE
    Switch for Loading file for tropical forcing. Turn it on when you use tropical forcing test
    \def WATER
    Switch for microphysics (qc, qv, qr) with warm rain scheme
*/
// ***********************************************************************************
