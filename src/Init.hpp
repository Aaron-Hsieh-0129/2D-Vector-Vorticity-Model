#include <iostream>
#include <random>
#include "Const.hpp"
#include "Declare.hpp"

class Init {
	public:
			static void Init1d(vvmArray &);
			static void Init2d(vvmArray &);
			
	private:
			static double GetTB(int);
			static double GetTHRAD(int, int);
			static double GetTH(int, int);
			static double GetQVB(int);
};