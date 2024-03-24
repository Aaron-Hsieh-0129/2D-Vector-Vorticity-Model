#include "Declare.hpp"

class Init {
public:
    static void Init1d(vvm &);
    static void Init2d(vvm &);
    static void LoadFile(vvm &);
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
