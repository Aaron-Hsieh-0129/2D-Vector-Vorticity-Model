#define WORKDIRECTORY "/Users/Aaron/2021-Summer-Research/vvm2d/src/"
#define gravity (9.80665)
#define C_p (1003.5)
#define Rd (287.)
#define Cv (C_p - Rd)
#define PSURF (96500.)
#define P0 (100000.)

#define dx (250)
#define dz (250)
#define XRANGE (150000)
#define ZRANGE (15000)
#define nx (XRANGE / dx)
#define nz (ZRANGE / dz)
#define dt (0.1)
#define rdx (1. / (double) dx)
#define rdz (1. / (double) dz)
#define rdx2 (1. / ((double) dx * dx))
#define rdz2 (1. / ((double) dz * dz))
#define d2t (2. * dt)
#define TIMEEND (2001.)
#define OUTPUTSTEP (250)
#define Kx (1000.)
#define Kz (1000.)
#define TIMETS (0.01)
#define Lv (2500000.)

#define OUTPUTPATH "../outputs/"
#define OUTPUTFILEMODE
#define OUTPUTNC
#define TXTOUTPUT
#define NCOUTPUT
// #define OUTPUTGRAPHMODE

// #define DRY
// #define RHO1
// #define ADVECTIONU
// #define ADVECTIONW
// #define NoBouyance
#define DIFFUSION
#define TIMEFILTER
#define WATER
// #define SHEAR
// #define CLOUDLESS
// #define HEATFLUX