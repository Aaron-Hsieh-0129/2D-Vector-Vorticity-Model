#include "Output.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <netcdf>

using namespace std;
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

#if defined(OUTPUTNC)
using namespace netCDF;

void Output::output_nc(int n, vvm &model) {
    string ncName = OUTPUTPATH + (string) "nc/" + std::to_string(n) + (string) ".nc";

    NcFile dataFile(ncName, NcFile::replace);

    // Create netCDF dimensions
    NcDim zeroDim = dataFile.addDim("zero", 1);
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
    NcVar ubartopData = dataFile.addVar("ubarTop", ncDouble, zeroDim);
    
    thData.putVar(model.th);
    zetaData.putVar(model.zeta);
    uData.putVar(model.u);
    wData.putVar(model.w);
    double tmpubarTop[1];
    tmpubarTop[0] = model.ubarTopp;
    ubartopData.putVar(tmpubarTop);

    #if defined(WATER)
        NcVar qcData = dataFile.addVar("qc", ncDouble, xzDim);
        NcVar qrData = dataFile.addVar("qr", ncDouble, xzDim);
        NcVar qvData = dataFile.addVar("qv", ncDouble, xzDim);
        NcVar qrAccData = dataFile.addVar("qrAcc", ncDouble, xDim);
        NcVar precipData = dataFile.addVar("precip", ncDouble, xDim);
        NcVar accretionData = dataFile.addVar("accretion", ncDouble, xzDim);
        NcVar autoconversionData = dataFile.addVar("autoconversion", ncDouble, xzDim);
        NcVar evaporationData = dataFile.addVar("evaporation", ncDouble, xzDim);
        qcData.putVar(model.qc);
        qrData.putVar(model.qr);
        qvData.putVar(model.qv);
        qrAccData.putVar(model.qrAcc);
        precipData.putVar(model.precip);
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

void Output::output_time_nc(int n, vvm &model) {
    string ncName = OUTPUTPATH + (string) "timer/" + std::to_string(n) + (string) ".nc";

    NcFile dataFile(ncName, NcFile::replace);
    NcDim tDim = dataFile.addDim("steps", TIMEROUTPUTSIZE);

    NcVar advectionData = dataFile.addVar("advection", ncDouble, tDim);
    NcVar poissonData = dataFile.addVar("poisson", ncDouble, tDim);
    NcVar diffusionData = dataFile.addVar("diffusion", ncDouble, tDim);
    NcVar microphysicsData = dataFile.addVar("microphysics", ncDouble, tDim);
    NcVar allData = dataFile.addVar("all", ncDouble, tDim);

    advectionData.putVar(model.t_advection);
    poissonData.putVar(model.t_poisson);
    diffusionData.putVar(model.t_diffusion);
    microphysicsData.putVar(model.t_microphysics);
    allData.putVar(model.t_all);
    return;
}

#endif

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
        create_directory(OUTPUTPATH + (string) "timer");
    #endif

    #if defined(OUTPUTTXT)
        create_directory(OUTPUTPATH + (string) "txtoutputs");
        create_directory(OUTPUTPATH + (string) "txtoutputs/u");
        create_directory(OUTPUTPATH + (string) "txtoutputs/w");
        create_directory(OUTPUTPATH + (string) "txtoutputs/zeta");
        create_directory(OUTPUTPATH + (string) "txtoutputs/th");
        #if defined(WATER)
            create_directory(OUTPUTPATH + (string) "txtoutputs/qc");
            create_directory(OUTPUTPATH + (string) "txtoutputs/qr");
            create_directory(OUTPUTPATH + (string) "txtoutputs/qv");
            create_directory(OUTPUTPATH + (string) "txtoutputs/precip");
            create_directory(OUTPUTPATH + (string) "txtoutputs/precipAcc");
        #endif
    #endif

    // plot directory
    create_directory(OUTPUTPATH + (string) "graphs");
}


#if defined(OUTPUTTXT)
void Output::output_zeta(int n, vvm &model) {
    std::fstream foutzeta;
    string zetaName = OUTPUTPATH + (string) "txtoutputs/zeta/zeta_" + std::to_string(n) + (string) ".txt";
    foutzeta.open(zetaName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutzeta << model.zeta[i][k] << " ";
        }
    }
}

void Output::output_th(int n, vvm &model) {
    std::fstream foutth;
    string thName = OUTPUTPATH + (string) "txtoutputs/th/th_" + std::to_string(n) + (string) ".txt";
    foutth.open(thName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutth << model.th[i][k] << " ";
        }
    }
}

void Output::output_u(int n, vvm &model) {
    std::fstream foutu;
    string uName = OUTPUTPATH + (string) "txtoutputs/u/u_" + std::to_string(n) + (string) ".txt";
    foutu.open(uName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutu << model.u[i][k] << " ";
        }
    }
}

void Output::output_w(int n, vvm &model) {
    std::fstream foutw;
    string wName = OUTPUTPATH + (string) "txtoutputs/w/w_" + std::to_string(n) + (string) ".txt";
    foutw.open(wName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutw << model.w[i][k] << " ";
        }
    }
}

void Output::output_qv(int n, vvm &model) {
    std::fstream foutqv;
    string qvName = OUTPUTPATH + (string) "txtoutputs/qv/qv_" + std::to_string(n) + (string) ".txt";
    foutqv.open(qvName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqv << model.qv[i][k] << " ";
        }
    }
}

void Output::output_qc(int n, vvm &model) {
    std::fstream foutqc;
    string qcName = OUTPUTPATH + (string) "txtoutputs/qc/qc_" + std::to_string(n) + (string) ".txt";
    foutqc.open(qcName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqc << model.qc[i][k] << " ";
        }
    }
}

void Output::output_qr(int n, vvm &model) {
    std::fstream foutqr;
    string qrName = OUTPUTPATH + (string) "txtoutputs/qr/qr_" + std::to_string(n) + (string) ".txt";
    foutqr.open(qrName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqr << model.qr[i][k] << " ";
        }
    }
}

void Output::output_precip(int n, vvm &model) {
    std::fstream foutqr;
    string qrName = OUTPUTPATH + (string) "txtoutputs/precip/precip_" + std::to_string(n) + (string) ".txt";
    foutqr.open(qrName, std::ios::out);
    for (int i = 0; i < model.nx; i++) {
        foutqr << model.precip[i] << " ";
    }
}

void Output::output_precipAcc(int n, vvm &model) {
    std::fstream foutqr;
    string qrName = OUTPUTPATH + (string) "txtoutputs/precip/precipAll_" + std::to_string(n) + (string) ".txt";
    foutqr.open(qrName, std::ios::out);
    for (int i = 0; i < model.nx; i++) {
        foutqr << model.qrAcc[i] << " ";
    }
}
#endif
