#include "Declare.hpp"

class Init {
public:
    static void Init1d(vvm &);
    static void Init2d(vvm &);
    #if defined(LOADFILE)
        static void LoadFile(vvm &);
    #elif defined(LOADFROMPREVIOUSFILE)
        static void LoadFromPreviousFile(vvm &);
    #endif
    #if defined(TROPICALFORCING)
        static void RandomPerturbation(vvm &, int);
    #endif
        
private:
    static double GetTB(int);
    static double GetTHRAD(int, int);
    static double GetTH(int, int);
    #if defined(WATER)
        static double GetQVB(int);
    #endif
};
