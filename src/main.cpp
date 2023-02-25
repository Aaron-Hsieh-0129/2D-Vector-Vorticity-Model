#include "Iteration.hpp"
#ifdef OUTPUTGRAPHMODE
	#include "matplotlib-cpp-master/matplotlibcpp.h"
#endif 

vvmArray myArray;

int main(void) {
	Init::Init1d(myArray);
	Init::Init2d(myArray);
	Output::create_all_directory();
	Output::printInit(myArray);
	Iteration::LeapFrog(myArray);
	return 0;
}
