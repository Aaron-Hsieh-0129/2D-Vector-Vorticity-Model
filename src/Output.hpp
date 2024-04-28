#include "Init.hpp"
#include <netcdf>

class Output {
public:
    static void printInit(vvm &);
    static void create_all_directory();
    #if defined(OUTPUTNC)
        static void output_nc(int, vvm &);
        static void output_time_nc(int, vvm &);
    #endif

    #if defined(OUTPUTTXT)
        static void output_zeta(int, vvm &);
        static void output_th(int, vvm &);
        static void output_u(int, vvm &);
        static void output_w(int, vvm &);
        #if defined(WATER)
            static void output_qv(int, vvm &);
            static void output_qc(int, vvm &);
            static void output_qr(int, vvm &);
            static void output_precip(int, vvm &);
            static void output_precipAcc(int, vvm &);
        #endif
    #endif
    
private:
    static void create_directory(std::string);
};
