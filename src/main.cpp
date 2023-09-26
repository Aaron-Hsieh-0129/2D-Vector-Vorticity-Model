#include "Iteration.hpp"
#ifdef OUTPUTGRAPHMODE
	#include "matplotlib-cpp-master/matplotlibcpp.h"
#endif 

vvmArray model;

int main(void) {
	Init::Init1d(model);
	Init::Init2d(model);
	Output::create_all_directory();
	Output::printInit(model);
	Iteration::LeapFrog(model);
	return 0;
}
