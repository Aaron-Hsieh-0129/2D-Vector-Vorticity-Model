#include "Config.hpp"

#define GRAVITY (9.80665)
#define C_p (1003.5)
#define Rd (287.)
#define Cv (C_p - Rd)
#define PSURF (96500.)
#define P0 (100000.)
#define Lv (2500000.)
#define NX (XRANGE / dx)
#define NZ (ZRANGE / dz)

#if defined(AB3)
    #define ALPHA (1.)
#else
    #define ALPHA (0.)
#endif





// ***********************************************************************************
// Documentation Part
/**
 * @file Const.hpp
 * @author Aaron Hsieh (b08209006@ntu.edu.tw)
 * @brief 
 * @version 0.1
 * @date 2024-03-13
 * 
 * @copyright Copyright (c) 2024
 * @brief Here are some constants defined in this model. You don't need to change them.
 * 
*/

/*!
    \def GRAVITY
    The gravity constant (m/s^2)
    \def C_p
    The specific heat at constant pressure (J/kg/K)
    \def Rd
    The gas constant for dry air (J/kg/K)
    \def Cv
    The specific heat at constant volume (J/kg/K)
    \def PSURF
    The surface pressure (Pa) for the initial condition
    \def P0
    The reference pressure (Pa) for the initial condition
    \def Lv
    The latent heat of vaporization (J/kg)
    \def NX
    The number of grid points in x direction
    \def NZ
    The number of grid points in z direction
*/
// ***********************************************************************************
