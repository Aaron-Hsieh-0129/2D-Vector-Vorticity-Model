#include "Output.hpp"

class Iteration {
public:
    static void pzeta_pt(vvm &);
    static void pth_pt(vvm &);
    #if defined(WATER)
        static void pqv_pt(vvm &);
        static void pqc_pt(vvm &);
        static void pqr_pt(vvm &);
    #endif

    static void updateMean(vvm &);
    static void TimeMarching(vvm &);
};
