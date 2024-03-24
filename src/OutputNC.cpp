#include "OutputNC.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <netcdf>

using namespace std;
using namespace netCDF;
void Output::printInit(vvm &model) {
    double z;
    std::cout << "z             thb        thb_zeta     rhou       rhow       qvb   	 RH      pib" << std::endl;
    for (int k = 0; k <= model.nz-1;k++){
        z = (double) (k - 0.5) * dz ;
        std::cout << std::fixed << std::setprecision(4) << z << "    " << model.thb[k] << "    " << model.thb_zeta[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "    " << model.qvb[k] * 1000 << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << std::endl;
    }
    std::fstream initout;
    string initName = OUTPUTPATH + (string) "init.txt";
    initout.open(initName, std::ios::out);
    for (int k = 0; k <= model.nz-1; k++) {
        z = (double) (k - 0.5) * dz ;
        initout << z << "    " << model.thb[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "   	 " << model.qvb[k] << "    " << model.qvsb[k] << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << std::endl;
    }
    return;
};

void Output::output_nc(int n, vvm &model) {
    string ncName = OUTPUTPATH + (string) "nc/" + std::to_string(n) + (string) ".nc";

    NcFile dataFile(ncName, NcFile::replace);

    // Create netCDF dimensions
    NcDim xDim = dataFile.addDim("x", model.nx);
    NcDim zDim = dataFile.addDim("z", model.nz);
    vector<NcDim> xzDim, zNcDim;
    xzDim.push_back(xDim);
    xzDim.push_back(zDim);
    zNcDim.push_back(zDim);

    NcVar thData = dataFile.addVar("th", ncDouble, xzDim);
    NcVar zetaData = dataFile.addVar("zeta", ncDouble, xzDim);
    NcVar uData = dataFile.addVar("u", ncDouble, xzDim);
    NcVar wData = dataFile.addVar("w", ncDouble, xzDim);
    
    thData.putVar(model.th);
    zetaData.putVar(model.zeta);
    uData.putVar(model.u);
    wData.putVar(model.w);

    #if defined(WATER)
        NcVar qcData = dataFile.addVar("qc", ncDouble, xzDim);
        NcVar qrData = dataFile.addVar("qr", ncDouble, xzDim);
        NcVar qvData = dataFile.addVar("qv", ncDouble, xzDim);
        NcVar qrAccData = dataFile.addVar("qrAcc", ncDouble, xDim);
        NcVar accretionData = dataFile.addVar("accretion", ncDouble, xzDim);
        NcVar autoconversionData = dataFile.addVar("autoconversion", ncDouble, xzDim);
        NcVar evaporationData = dataFile.addVar("evaporation", ncDouble, xzDim);
        qcData.putVar(model.qc);
        qrData.putVar(model.qr);
        qvData.putVar(model.qv);
        qrAccData.putVar(model.qrAcc);
        accretionData.putVar(model.accretion);
        autoconversionData.putVar(model.autoconversion);
        evaporationData.putVar(model.evaporation);
    #endif

    if (n == 0) {
        NcVar thbData = dataFile.addVar("tb", ncDouble, zNcDim);
        NcVar rhouData = dataFile.addVar("rhou", ncDouble, zNcDim);
        NcVar rhowData = dataFile.addVar("rhow", ncDouble, zNcDim);
        NcVar qvbData = dataFile.addVar("qvb", ncDouble, zNcDim);
        NcVar qvsbData = dataFile.addVar("qvsb", ncDouble, zNcDim);

        thbData.putVar(model.thb);
        rhouData.putVar(model.rhou);
        rhowData.putVar(model.rhow);
        qvbData.putVar(model.qvb);
        qvsbData.putVar(model.qvsb);
    }
}

void Output::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void Output::create_all_directory() {
    // data directory
    #ifdef OUTPUTNC
        create_directory(OUTPUTPATH + (string) "nc");
    #endif

    // plot directory
    create_directory(OUTPUTPATH + (string) "graphs");
}
