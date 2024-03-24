#include "Init.hpp"
#include <netcdf>

class Output {
public:
    static void printInit(vvm &);
    static void output_nc(int, vvm &);
    static void create_all_directory();
private:
    static void create_directory(std::string);
};