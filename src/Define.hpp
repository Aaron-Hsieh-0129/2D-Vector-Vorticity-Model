#define WORKDIRECTORY "/Users/Aaron/2021-Summer-Research/vvm2d/src/"
#define gravity (9.80665)
#define C_p (1003.5)
#define Rd (287.)
#define Cv (C_p - Rd)
#define PSURF (96500.)
#define P0 (100000.)

#define dx (200)
#define dz (200)
#define XRANGE (120000+2*dx)
#define ZRANGE (20000+2*dz)
#define nx (XRANGE / dx)
#define nz (ZRANGE / dz)
#define dt (1.)
#define rdx (1. / (double) dx)
#define rdz (1. / (double) dz)
#define rdx2 (1. / ((double) dx * dx))
#define rdz2 (1. / ((double) dz * dz))
#define d2t (2. * dt)
#define TIMEEND (86400. * 3. + 1.)
#define Kx (250.)
#define Kz (250.)
#define TIMETS (0.01)
#define Lv (2500000.)

#define OUTPUTSTEP (10)
#define OUTPUTPATH "/data/Aaron/2DVVM/20km_full_qv/6Q16Q2_qvfull_diff250/"
// #define OUTPUTTXT
#define OUTPUTNC
// #define OUTPUTGRAPHMODE

// #define DRY
// #define RHO1
// #define ADVECTIONU
// #define ADVECTIONW
// #define NoBouyance
#define DIFFUSION
#define TIMEFILTER
#define WATER
#define TROPICALFORCING
#define RADIATIONCOOLING
// #define SHEAR
// #define CLOUDLESS
// #define HEATFLUX

#define LOADFILE
#define ADDFORCINGTIME (1200.)
// #define LINEARIZEDQV
#define LINEARIZEDTH