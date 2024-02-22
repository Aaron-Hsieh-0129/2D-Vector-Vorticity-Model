#include "Init.hpp"
#include <netcdf>
// #include <vector>

#ifdef OUTPUTGRAPHMODE
    #include "include/python3.10/Python.h"
    #ifndef WITHOUT_NUMPY
    #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
    #include "include/numpy/arrayobject.h"
    #endif
#endif

class Output {
public:
    static void printInit(vvmArray &);
    static void output_zeta(int, vvmArray &);
    static void output_th(int, vvmArray &);
    static void output_u(int, vvmArray &);
    static void output_w(int, vvmArray &);
    static void output_qv(int, vvmArray &);
    static void output_qc(int, vvmArray &);
    static void output_qr(int, vvmArray &);
    static void output_nc(int, vvmArray &);
    static void create_all_directory();
private:
    static void create_directory(std::string);
};

#ifdef OUTPUTGRAPHMODE
    class Plot {
        public:
            static int *initNumpy(void);
            static void initPython(void);
            static void plotInit(vvmArray &);
            static void plot_zeta(int, vvmArray &);
            static void plot_th2d(int, vvmArray &);
            static void plot_u(int, vvmArray &);
            static void plot_w(int, vvmArray &);
            static void plot_qv(int, vvmArray &);
            static void plot_qc(int, vvmArray &);
            static void plot_qr(int, vvmArray &);
            static void plot_qc_qr(int, vvmArray &);
            static void plot_qv_qc(int, vvmArray &);
            static void plot_qc_qr_th_u_w(int, vvmArray &);
            static void plot_qr_th_u_w(int, vvmArray &);
    };
#endif
