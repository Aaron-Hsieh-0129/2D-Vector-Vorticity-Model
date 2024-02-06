#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include "Declare.hpp"

class Init {
	public:
		static void Init1d(vvmArray &);
		static void Init2d(vvmArray &);
		static void LoadFile(vvmArray &);
		#if defined(TROPICALFORCING)
			static void RandomPerturbation(vvmArray &, int);
		#endif
			
	private:
		static double GetTB(int);
		static double GetTHRAD(int, int);
		static double GetTH(int, int);
		static double GetQVB(int);
		static void InitPoissonMatrix(vvmArray &);
};