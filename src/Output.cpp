#include "Declare.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <netcdf>

using namespace std;
void vvm::Output::printInit(vvm &model) {
    double z;
    std::cout << "z             thb        thb_zeta     rhou       rhow       qvb   	 RH      pib" << std::endl;
    for (int k = 0; k <= model.nz-1;k++){
        z = (double) (k - 0.5) * model.dz ;
        std::cout << std::fixed << std::setprecision(4) << z << "    " << model.thb[k] << "    " << model.thb_zeta[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "    " << model.qvb[k] * 1000 << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << std::endl;
    }
    std::fstream initout;
    string initName = model.outputpath + (string) "init.txt";
    initout.open(initName, std::ios::out);
    for (int k = 0; k <= model.nz-1; k++) {
        z = (double) (k - 0.5) * model.dz ;
        initout << z << "    " << model.thb[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "   	 " << model.qvb[k] << "    " << model.qvsb[k] << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << std::endl;
    }
    initout.close();
    return;
};

#if defined(OUTPUTNC)
using namespace netCDF;
void checkErr(int status, int line) {
    if (status != NC_NOERR) {
        cerr << "NetCDF error at line " << line << ": " << nc_strerror(status) << endl;
        exit(EXIT_FAILURE);
    }
}

void vvm::Output::output_nc(int n, vvm &model) {
    string ncName = model.outputpath + (string) "nc/" + std::to_string(n) + (string) ".nc";

    int ncid, x_dimid, z_dimid;
    int retval;

    int thid, zetaid, uid, wid, ubarTopid;

    if ((retval = nc_create(ncName.c_str(), NC_CLOBBER, &ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_def_dim(ncid, "x", model.nx, &x_dimid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_dim(ncid, "z", model.nz, &z_dimid))) checkErr(retval, __LINE__);

    int dimids[2] = {x_dimid, z_dimid};
    int dimx1d[1] = {x_dimid};
    int dimz1d[1] = {z_dimid};
    if ((retval = nc_def_var(ncid, "th", NC_DOUBLE, 2, dimids, &thid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "zeta", NC_DOUBLE, 2, dimids, &zetaid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "u", NC_DOUBLE, 2, dimids, &uid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "w", NC_DOUBLE, 2, dimids, &wid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "ubarTop", NC_DOUBLE, 0, nullptr, &ubarTopid))) checkErr(retval, __LINE__);

    #if defined(WATER)
        int qvid, qcid, qrid, precipid, accretionid, autoconversionid, evaporationid, condensationid;

        if ((retval = nc_def_var(ncid, "qv", NC_DOUBLE, 2, dimids, &qvid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "qc", NC_DOUBLE, 2, dimids, &qcid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "qr", NC_DOUBLE, 2, dimids, &qrid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "accretion", NC_DOUBLE, 2, dimids, &accretionid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "autoconversion", NC_DOUBLE, 2, dimids, &autoconversionid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "evaporation", NC_DOUBLE, 2, dimids, &evaporationid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "condensation", NC_DOUBLE, 2, dimids, &condensationid))) checkErr(retval, __LINE__);
        if ((retval = nc_def_var(ncid, "precip", NC_DOUBLE, 1, dimx1d, &precipid))) checkErr(retval, __LINE__);
    #endif

    if ((retval = nc_enddef(ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_put_var_double(ncid, thid, model.thcont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, zetaid, model.zetacont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, uid, model.ucont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, wid, model.wcont))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, ubarTopid, &model.ubarTopp))) checkErr(retval, __LINE__);

    #if defined(WATER)
        if ((retval = nc_put_var_double(ncid, qvid, model.qvcont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, qcid, model.qccont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, qrid, model.qrcont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, accretionid, model.accretioncont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, autoconversionid, model.autoconversioncont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, evaporationid, model.evaporationcont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, condensationid, model.condensationcont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_var_double(ncid, precipid, model.precip))) checkErr(retval, __LINE__);
    #endif

    if ((retval = nc_close(ncid))) checkErr(retval, __LINE__);
}

void vvm::Output::output_time_nc(int n, vvm &model) {
    string ncName = model.outputpath + (string) "timer/" + std::to_string(n) + (string) ".nc";

    int ncid, t_dimid;
    int retval;

    int advectionid, poissonid, diffusionid, microphysicsid, allid;

    if ((retval = nc_create(ncName.c_str(), NC_CLOBBER, &ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_def_dim(ncid, "x", model.nx, &t_dimid))) checkErr(retval, __LINE__);

    int dimt1d[1] = {t_dimid};

    if ((retval = nc_def_var(ncid, "advection", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "poisson", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "diffusion", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "microphysics", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "all", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);

    if ((retval = nc_put_var_double(ncid, advectionid, model.t_advection))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, poissonid, model.t_poisson))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, diffusionid, model.t_diffusion))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, microphysicsid, model.t_microphysics))) checkErr(retval, __LINE__);
    if ((retval = nc_put_var_double(ncid, allid, model.t_all))) checkErr(retval, __LINE__);

    if ((retval = nc_close(ncid))) checkErr(retval, __LINE__);
    return;
}

#endif

void vvm::Output::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void vvm::Output::create_all_directory(vvm &model) {
    // data directory
    #ifdef OUTPUTNC
        create_directory(model.outputpath + (string) "nc");
        create_directory(model.outputpath + (string) "timer");
    #endif

    #if defined(OUTPUTTXT)
        create_directory(model.outputpath + (string) "txtoutputs");
        create_directory(model.outputpath + (string) "txtoutputs/u");
        create_directory(model.outputpath + (string) "txtoutputs/w");
        create_directory(model.outputpath + (string) "txtoutputs/zeta");
        create_directory(model.outputpath + (string) "txtoutputs/th");
        #if defined(WATER)
            create_directory(model.outputpath + (string) "txtoutputs/qc");
            create_directory(model.outputpath + (string) "txtoutputs/qr");
            create_directory(model.outputpath + (string) "txtoutputs/qv");
            create_directory(model.outputpath + (string) "txtoutputs/precip");
            create_directory(model.outputpath + (string) "txtoutputs/precipAcc");
        #endif
    #endif

    // plot directory
    create_directory(model.outputpath + (string) "graphs");
}


#if defined(OUTPUTTXT)
void vvm::Output::output_zeta(int n, vvm &model) {
    std::fstream foutzeta;
    string zetaName = model.outputpath + (string) "txtoutputs/zeta/zeta_" + std::to_string(n) + (string) ".txt";
    foutzeta.open(zetaName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutzeta << model.zeta[i][k] << " ";
        }
    }
    foutzeta.close();
}

void vvm::Output::output_th(int n, vvm &model) {
    std::fstream foutth;
    string thName = model.outputpath + (string) "txtoutputs/th/th_" + std::to_string(n) + (string) ".txt";
    foutth.open(thName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutth << model.th[i][k] << " ";
        }
        foutth << std::endl;
    }
    foutth.close();
}

void vvm::Output::output_u(int n, vvm &model) {
    std::fstream foutu;
    string uName = model.outputpath + (string) "txtoutputs/u/u_" + std::to_string(n) + (string) ".txt";
    foutu.open(uName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutu << model.u[i][k] << " ";
        }
        foutu << std::endl;
    }
    foutu.close();
}

void vvm::Output::output_w(int n, vvm &model) {
    std::fstream foutw;
    string wName = model.outputpath + (string) "txtoutputs/w/w_" + std::to_string(n) + (string) ".txt";
    foutw.open(wName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutw << model.w[i][k] << " ";
        }
        foutw << std::endl;
    }
    foutw.close();
}

#if defined(WATER)
void vvm::Output::output_qv(int n, vvm &model) {
    std::fstream foutqv;
    string qvName = model.outputpath + (string) "txtoutputs/qv/qv_" + std::to_string(n) + (string) ".txt";
    foutqv.open(qvName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqv << model.qv[i][k] << " ";
        }
        foutqv << std::endl;
    }
    foutqv.close();
}

void vvm::Output::output_qc(int n, vvm &model) {
    std::fstream foutqc;
    string qcName = model.outputpath + (string) "txtoutputs/qc/qc_" + std::to_string(n) + (string) ".txt";
    foutqc.open(qcName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqc << model.qc[i][k] << " ";
        }
        foutqc << std::endl;
    }
    foutqc.close();
}

void vvm::Output::output_qr(int n, vvm &model) {
    std::fstream foutqr;
    string qrName = model.outputpath + (string) "txtoutputs/qr/qr_" + std::to_string(n) + (string) ".txt";
    foutqr.open(qrName, std::ios::out);
    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            foutqr << model.qr[i][k] << " ";
        }
        foutqr << std::endl;
    }
    foutqr.close();
}

void vvm::Output::output_precip(int n, vvm &model) {
    std::fstream foutqr;
    string qrName = model.outputpath + (string) "txtoutputs/precip/precip_" + std::to_string(n) + (string) ".txt";
    foutqr.open(qrName, std::ios::out);
    for (int i = 0; i < model.nx; i++) {
        foutqr << model.precip[i] << " ";
    }
    foutqr.close();
}
#endif

void vvm::Output::outputalltxt(int n, vvm &model) {
    vvm::Output::output_zeta(n, model);
    vvm::Output::output_th(n, model);
    vvm::Output::output_u(n, model);
    vvm::Output::output_w(n, model);
    #if defined(WATER)
        vvm::Output::output_qv(n, model);
        vvm::Output::output_qc(n, model);
        vvm::Output::output_qr(n, model);
        vvm::Output::output_precip(n, model);
    #endif
}
#endif
