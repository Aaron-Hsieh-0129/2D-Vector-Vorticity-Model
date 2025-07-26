#include "Declare.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <netcdf>
#include <filesystem>

using namespace std;
void vvm::Output::printInit(vvm &model) {
    std::cout << "z          z_zeta        dz_th      dz_zeta          Pbar          thb        thb_zeta     rhou       rhow       qvb   	 RH      pib" << std::endl;
    for (int k = 0; k <= model.nz-1;k++){
        std::cout << std::fixed << std::setprecision(4) << model.z[k] << "      " << model.z_zeta[k] << "    " << model.dz_th[k] << "    " << model.dz_zeta[k] << "    " << model.pb[k] << "     " << model.thb[k] << "    " << model.thb_zeta[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "    " << model.qvb[k] * 1000 << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << std::endl;
    }
    std::fstream initout;
    string initName = model.outputpath + (string) "init.txt";
    initout.open(initName, std::ios::out);
    for (int k = 1; k <= model.nz-2; k++) {
        initout << model.z[k] << "    " << model.z_zeta[k] << "    " << model.dz_th[k] << "    " << model.dz_zeta[k] << "    " << model.thb[k] << "    " << model.rhou[k] << "     " 
        << model.rhow[k] << "   	 " << model.qvb[k] << "    " << model.qvsb[k] << "    " << model.qvb[k] / model.qvsb[k] << "    "
        << model.pib[k] << "    " << model.pb[k] << std::endl;
    }
    initout.close();

    std::fstream initout2;
    string init2Name = model.outputpath + (string) "init2.txt";
    initout2.open(init2Name, std::ios::out);
    for (int k = 1; k <= model.nz-1; k++) {
        initout2 << model.z[k] << "      " << model.pb_lev[k] << "     " << model.pib_lev[k] << std::endl;
    }
    initout2.close();

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
#if !defined(OUTPUTNCSAMEFILE)
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
#else
void vvm::Output::output_nc(int n, vvm &model) {
    #define NC_ERR(e) { printf("NetCDF error: %s\n", nc_strerror(e)); exit(2); }
    double t = n * model.dt;

    int ncid, t_dimid, x_dimid, z_dimid;
    int t_varid, x_varid, z_varid, th_id, zeta_id, u_id, w_id, ubarTop_id;
    int waterflux_id, heatflux_id, momentumflux_id;
    #if defined(RTERRTMGP)
        int rad_net_heat_rate_id, rad_lw_heat_rate_id, rad_sw_heat_rate_id;
    #endif
    #if defined(WATER)
        int qvid, qcid, qrid, precipid;
        #if defined(KESSLER_MICROPHY)
        #if defined(OUTPUTMICROPHYSICS)
            int accretionid, autoconversionid, evaporationid, condensationid;
        #endif
        #endif
        #if defined(P3_MICROPHY)
            int qitotid;
        #endif
    #endif

    int retval;
    size_t t_index = 0;

    int file_num = (n / 1200000);
    std::string file_name = model.outputpath + "nc/" + std::to_string(file_num) + ".nc";

    if ((retval = nc_open(file_name.c_str(), NC_WRITE, &ncid)) == NC_NOERR) {
        // File exists, get the dimensions and variable IDs
        printf("File %s exists. Opening for appending.\n", file_name.c_str());

        // Get variable IDs
        if ((retval = nc_inq_varid(ncid, "t", &t_varid))) NC_ERR(retval);
        if ((retval = nc_inq_varid(ncid, "th", &th_id))) NC_ERR(retval);
        if ((retval = nc_inq_varid(ncid, "zeta", &zeta_id))) NC_ERR(retval);
        if ((retval = nc_inq_varid(ncid, "u", &u_id))) NC_ERR(retval);
        if ((retval = nc_inq_varid(ncid, "w", &w_id))) NC_ERR(retval);
        if ((retval = nc_inq_varid(ncid, "ubarTop", &ubarTop_id))) NC_ERR(retval);
        #if defined(RTERRTMGP)
            if ((retval = nc_inq_varid(ncid, "rad_net_heat_rate", &rad_net_heat_rate_id))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "rad_lw_heat_rate", &rad_lw_heat_rate_id))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "rad_sw_heat_rate", &rad_sw_heat_rate_id))) NC_ERR(retval);
        #endif
        #if defined(SFCFLX)
            if ((retval = nc_inq_varid(ncid, "heatflux_sfc", &heatflux_id))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "waterflux_sfc", &waterflux_id))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "momentumflux_sfc", &momentumflux_id))) NC_ERR(retval);
        #endif

        #if defined(WATER)
            if ((retval = nc_inq_varid(ncid, "qv", &qvid))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "qc", &qcid))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "qr", &qrid))) NC_ERR(retval);
            if ((retval = nc_inq_varid(ncid, "precip", &precipid))) NC_ERR(retval);
            #if defined(KESSLER_MICROPHY)
            #if defined(OUTPUTMICROPHYSICS)
                if ((retval = nc_inq_varid(ncid, "accretion", &accretionid))) NC_ERR(retval);
                if ((retval = nc_inq_varid(ncid, "autoconversion", &autoconversionid))) NC_ERR(retval);
                if ((retval = nc_inq_varid(ncid, "evaporation", &evaporationid))) NC_ERR(retval);
                if ((retval = nc_inq_varid(ncid, "condensation", &condensationid))) NC_ERR(retval);
            #endif
            #endif

            #if defined(P3_MICROPHY)
                // if ((retval = nc_inq_varid(ncid, "qnc", &qncid))) NC_ERR(retval);
                // if ((retval = nc_inq_varid(ncid, "qnr", &qnrid))) NC_ERR(retval);
                // if ((retval = nc_inq_varid(ncid, "qni", &qniid))) NC_ERR(retval);
                if ((retval = nc_inq_varid(ncid, "qitot", &qitotid))) NC_ERR(retval);
                // if ((retval = nc_inq_varid(ncid, "qirim", &qirimid))) NC_ERR(retval);
                // if ((retval = nc_inq_varid(ncid, "birim", &birimid))) NC_ERR(retval);
            #endif
        #endif

        // Get the current time index (size of the time dimension)
        size_t len;
        if ((retval = nc_inq_dimlen(ncid, t_varid, &len))) NC_ERR(retval);
        t_index = len;
    } 
    else {
        // File doesn't exist, create a new file
        printf("File %s does not exist. Creating new file.\n", file_name.c_str());

        if ((retval = nc_create(file_name.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid))) NC_ERR(retval);

        // Define global attributes
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen("CF-1.8"), "CF-1.8"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen("2DVVM Output"), "2DVVM Output"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "source", strlen("2DVVM Model"), "2DVVM Model"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen("Created by 2DVVM model"), "Created by 2DVVM model"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "comment", strlen("Output from 2DVVM atmospheric simulation"), "Output from 2DVVM atmospheric simulation"))) NC_ERR(retval);
        char creation_date[25];
        time_t now = time(NULL);
        strftime(creation_date, sizeof(creation_date), "%Y-%m-%dT%H:%M:%SZ", gmtime(&now));
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "creation_date", strlen(creation_date), creation_date))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "creator_name", strlen("Aaron Hsieh"), "Aaron Hsieh"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "creator_email", strlen("b08209006@ntu.edu.tw"), "b08209006@ntu.edu.tw"))) NC_ERR(retval);

        // Define dimensions
        if ((retval = nc_def_dim(ncid, "t", NC_UNLIMITED, &t_dimid))) NC_ERR(retval);
        if ((retval = nc_def_dim(ncid, "x", model.nx, &x_dimid))) NC_ERR(retval);
        if ((retval = nc_def_dim(ncid, "z", model.nz, &z_dimid))) NC_ERR(retval);

        // Define coordinate variables
        if ((retval = nc_def_var(ncid, "t", NC_DOUBLE, 1, &t_dimid, &t_varid))) NC_ERR(retval);
        if ((retval = nc_def_var(ncid, "x", NC_DOUBLE, 1, &x_dimid, &x_varid))) NC_ERR(retval);
        if ((retval = nc_def_var(ncid, "z", NC_DOUBLE, 1, &z_dimid, &z_varid))) NC_ERR(retval);

        // Define data variables
        int dimids[3] = {t_dimid, x_dimid, z_dimid};
        int dimx1d[2] = {t_dimid, x_dimid};
        size_t chunksizes_th[3] = {1, 50, 10};
        if ((retval = nc_def_var(ncid, "th", NC_DOUBLE, 3, dimids, &th_id))) NC_ERR(retval);
        if ((retval = nc_def_var_chunking(ncid, th_id, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
        if ((retval = nc_def_var(ncid, "zeta", NC_DOUBLE, 3, dimids, &zeta_id))) NC_ERR(retval);
        if ((retval = nc_def_var_chunking(ncid, zeta_id, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
        if ((retval = nc_def_var(ncid, "u", NC_DOUBLE, 3, dimids, &u_id))) NC_ERR(retval);
        if ((retval = nc_def_var_chunking(ncid, u_id, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
        if ((retval = nc_def_var(ncid, "w", NC_DOUBLE, 3, dimids, &w_id))) NC_ERR(retval);
        if ((retval = nc_def_var_chunking(ncid, w_id, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);

        if ((retval = nc_def_var(ncid, "ubarTop", NC_DOUBLE, 1, &t_dimid, &ubarTop_id))) NC_ERR(retval);
        #if defined(SFCFLX)
            if ((retval = nc_def_var(ncid, "heatflux_sfc", NC_DOUBLE, 2, dimx1d, &heatflux_id))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "waterflux_sfc", NC_DOUBLE, 2, dimx1d, &waterflux_id))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "momentumflux_sfc", NC_DOUBLE, 2, dimx1d, &momentumflux_id))) NC_ERR(retval);
        #endif
        #if defined(RTERRTMGP)
            if ((retval = nc_def_var(ncid, "rad_net_heat_rate", NC_DOUBLE, 3, dimids, &rad_net_heat_rate_id))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "rad_lw_heat_rate", NC_DOUBLE, 3, dimids, &rad_lw_heat_rate_id))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "rad_sw_heat_rate", NC_DOUBLE, 3, dimids, &rad_sw_heat_rate_id))) NC_ERR(retval);
        #endif
        #if defined(WATER)
            if ((retval = nc_def_var(ncid, "qv", NC_DOUBLE, 3, dimids, &qvid))) NC_ERR(retval);
            if ((retval = nc_def_var_chunking(ncid, qvid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "qc", NC_DOUBLE, 3, dimids, &qcid))) NC_ERR(retval);
            if ((retval = nc_def_var_chunking(ncid, qcid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "qr", NC_DOUBLE, 3, dimids, &qrid))) NC_ERR(retval);
            if ((retval = nc_def_var_chunking(ncid, qrid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
            if ((retval = nc_def_var(ncid, "precip", NC_DOUBLE, 2, dimx1d, &precipid))) NC_ERR(retval);
            #if defined(KESSLER_MICROPHY)
            #if defined(OUTPUTMICROPHYSICS)
                if ((retval = nc_def_var(ncid, "accretion", NC_DOUBLE, 3, dimids, &accretionid))) NC_ERR(retval);
                if ((retval = nc_def_var_chunking(ncid, accretionid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
                if ((retval = nc_def_var(ncid, "autoconversion", NC_DOUBLE, 3, dimids, &autoconversionid))) NC_ERR(retval);
                if ((retval = nc_def_var_chunking(ncid, autoconversionid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
                if ((retval = nc_def_var(ncid, "evaporation", NC_DOUBLE, 3, dimids, &evaporationid))) NC_ERR(retval);
                if ((retval = nc_def_var_chunking(ncid, evaporationid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
                if ((retval = nc_def_var(ncid, "condensation", NC_DOUBLE, 3, dimids, &condensationid))) NC_ERR(retval);
            #endif
            #endif
            #if defined(P3_MICROPHY)
                if ((retval = nc_def_var(ncid, "qitot", NC_DOUBLE, 3, dimids, &qitotid))) NC_ERR(retval);
                if ((retval = nc_def_var_chunking(ncid, qitotid, NC_CHUNKED, chunksizes_th))) NC_ERR(retval);
            #endif
        #endif

        // Define attributes for coordinate variables
        if ((retval = nc_put_att_text(ncid, t_varid, "standard_name", strlen("time"), "time"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, t_varid, "long_name", strlen("Time"), "Time"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, t_varid, "units", strlen("seconds since 1970-01-01 00:00:00"), "seconds since 1970-01-01 00:00:00"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, t_varid, "calendar", strlen("gregorian"), "gregorian"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, t_varid, "axis", strlen("T"), "T"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, x_varid, "standard_name", strlen("projection_x_coordinate"), "projection_x_coordinate"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, x_varid, "long_name", strlen("Horizontal Coordinate"), "Horizontal Coordinate"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, x_varid, "units", strlen("m"), "m"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, x_varid, "axis", strlen("X"), "X"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, z_varid, "standard_name", strlen("height_above_reference_ellipsoid"), "height_above_reference_ellipsoid"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, z_varid, "long_name", strlen("Height Above Surface"), "Height Above Surface"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, z_varid, "units", strlen("m"), "m"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, z_varid, "axis", strlen("Z"), "Z"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, z_varid, "positive", strlen("up"), "up"))) NC_ERR(retval);

        // Define attributes for data variables
        double fill_value = -9999.0;
        if ((retval = nc_put_att_text(ncid, th_id, "standard_name", strlen("air_potential_temperature"), "air_potential_temperature"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, th_id, "long_name", strlen("Potential Temperature"), "Potential Temperature"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, th_id, "units", strlen("K"), "K"))) NC_ERR(retval);
        if ((retval = nc_put_att_double(ncid, th_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, th_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, zeta_id, "standard_name", strlen("atmosphere_relative_vorticity"), "atmosphere_relative_vorticity"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, zeta_id, "long_name", strlen("Relative Vorticity"), "Relative Vorticity"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, zeta_id, "units", strlen("s-1"), "s-1"))) NC_ERR(retval);
        if ((retval = nc_put_att_double(ncid, zeta_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, zeta_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, u_id, "standard_name", strlen("eastward_wind"), "eastward_wind"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, u_id, "long_name", strlen("Zonal Wind Speed"), "Zonal Wind Speed"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, u_id, "units", strlen("m s-1"), "m s-1"))) NC_ERR(retval);
        if ((retval = nc_put_att_double(ncid, u_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, u_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, w_id, "standard_name", strlen("upward_air_velocity"), "upward_air_velocity"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, w_id, "long_name", strlen("Vertical Air Velocity"), "Vertical Air Velocity"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, w_id, "units", strlen("m s-1"), "m s-1"))) NC_ERR(retval);
        if ((retval = nc_put_att_double(ncid, w_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, w_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

        if ((retval = nc_put_att_text(ncid, ubarTop_id, "standard_name", strlen("eastward_wind"), "eastward_wind"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, ubarTop_id, "long_name", strlen("Mean Zonal Wind at Top Level"), "Mean Zonal Wind at Top Level"))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, ubarTop_id, "units", strlen("m s-1"), "m s-1"))) NC_ERR(retval);
        if ((retval = nc_put_att_double(ncid, ubarTop_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
        if ((retval = nc_put_att_text(ncid, ubarTop_id, "coordinates", strlen("t"), "t"))) NC_ERR(retval);

        #if defined(SFCFLX)
            if ((retval = nc_put_att_text(ncid, heatflux_id, "standard_name", strlen("surface_upward_sensible_heat_flux"), "surface_upward_sensible_heat_flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, heatflux_id, "long_name", strlen("Surface Sensible Heat Flux"), "Surface Sensible Heat Flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, heatflux_id, "units", strlen("W m-2"), "W m-2"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, heatflux_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, heatflux_id, "coordinates", strlen("t x"), "t x"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, waterflux_id, "standard_name", strlen("surface_upward_water_vapor_flux"), "surface_upward_water_vapor_flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, waterflux_id, "long_name", strlen("Surface Water Vapor Flux"), "Surface Water Vapor Flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, waterflux_id, "units", strlen("kg m-2 s-1"), "kg m-2 s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, waterflux_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, waterflux_id, "coordinates", strlen("t x"), "t x"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, momentumflux_id, "standard_name", strlen("surface_upward_momentum_flux"), "surface_upward_momentum_flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, momentumflux_id, "long_name", strlen("Surface Momentum Flux"), "Surface Momentum Flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, momentumflux_id, "units", strlen("kg m-2 s-1"), "kg m-2 s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, momentumflux_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, momentumflux_id, "coordinates", strlen("t x"), "t x"))) NC_ERR(retval);
        #endif

        #if defined(RTERRTMGP)
            if ((retval = nc_put_att_text(ncid, rad_net_heat_rate_id, "standard_name", strlen("tendency_of_air_temperature_due_to_radiative_heating"), "tendency_of_air_temperature_due_to_radiative_heating"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_net_heat_rate_id, "long_name", strlen("Net Radiative Heating Rate"), "Net Radiative Heating Rate"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_net_heat_rate_id, "units", strlen("K s-1"), "K s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, rad_net_heat_rate_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_net_heat_rate_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, rad_lw_heat_rate_id, "standard_name", strlen("tendency_of_air_temperature_due_to_radiative_heating"), "tendency_of_air_temperature_due_to_radiative_heating"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_lw_heat_rate_id, "long_name", strlen("Longwave Radiative Heating Rate"), "Longwave Radiative Heating Rate"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_lw_heat_rate_id, "units", strlen("K s-1"), "K s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, rad_lw_heat_rate_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_lw_heat_rate_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, rad_sw_heat_rate_id, "standard_name", strlen("tendency_of_air_temperature_due_to_radiative_heating"), "tendency_of_air_temperature_due_to_radiative_heating"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_sw_heat_rate_id, "long_name", strlen("Shortwave Radiative Heating Rate"), "Shortwave Radiative Heating Rate"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_sw_heat_rate_id, "units", strlen("K s-1"), "K s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, rad_sw_heat_rate_id, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, rad_sw_heat_rate_id, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);
        #endif
        #if defined(WATER)
            if ((retval = nc_put_att_text(ncid, qvid, "standard_name", strlen("specific_humidity"), "specific_humidity"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qvid, "long_name", strlen("Water Vapor Mixing Ratio"), "Water Vapor Mixing Ratio"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qvid, "units", strlen("kg kg-1"), "kg kg-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, qvid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qvid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, qcid, "standard_name", strlen("mass_fraction_of_cloud_liquid_water_in_air"), "mass_fraction_of_cloud_liquid_water_in_air"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qcid, "long_name", strlen("Cloud Liquid Water Mixing Ratio"), "Cloud Liquid Water Mixing Ratio"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qcid, "units", strlen("kg kg-1"), "kg kg-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, qcid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qcid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, qrid, "standard_name", strlen("mass_fraction_of_rain_water_in_air"), "mass_fraction_of_rain_water_in_air"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qrid, "long_name", strlen("Rain Water Mixing Ratio"), "Rain Water Mixing Ratio"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qrid, "units", strlen("kg kg-1"), "kg kg-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, qrid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, qrid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

            if ((retval = nc_put_att_text(ncid, precipid, "standard_name", strlen("precipitation_flux"), "precipitation_flux"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, precipid, "long_name", strlen("Precipitation Rate"), "Precipitation Rate"))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, precipid, "units", strlen("kg m-2 s-1"), "kg m-2 s-1"))) NC_ERR(retval);
            if ((retval = nc_put_att_double(ncid, precipid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
            if ((retval = nc_put_att_text(ncid, precipid, "coordinates", strlen("t x"), "t x"))) NC_ERR(retval);

            #if defined(KESSLER_MICROPHY)
            #if defined(OUTPUTMICROPHYSICS)
                if ((retval = nc_put_att_text(ncid, accretionid, "long_name", strlen("Accretion Rate of Cloud Water to Rain"), "Accretion Rate of Cloud Water to Rain"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, accretionid, "units", strlen("kg kg-1 s-1"), "kg kg-1 s-1"))) NC_ERR(retval);
                if ((retval = nc_put_att_double(ncid, accretionid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, accretionid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

                if ((retval = nc_put_att_text(ncid, autoconversionid, "long_name", strlen("Autoconversion Rate of Cloud Water to Rain"), "Autoconversion Rate of Cloud Water to Rain"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, autoconversionid, "units", strlen("kg kg-1 s-1"), "kg kg-1 s-1"))) NC_ERR(retval);
                if ((retval = nc_put_att_double(ncid, autoconversionid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, autoconversionid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

                if ((retval = nc_put_att_text(ncid, evaporationid, "long_name", strlen("Evaporation Rate of Rain Water"), "Evaporation Rate of Rain Water"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, evaporationid, "units", strlen("kg kg-1 s-1"), "kg kg-1 s-1"))) NC_ERR(retval);
                if ((retval = nc_put_att_double(ncid, evaporationid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, evaporationid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);

                if ((retval = nc_put_att_text(ncid, condensationid, "long_name", strlen("Condensation Rate of Water Vapor to Cloud Water"), "Condensation Rate of Water Vapor to Cloud Water"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, condensationid, "units", strlen("kg kg-1 s-1"), "kg kg-1 s-1"))) NC_ERR(retval);
                if ((retval = nc_put_att_double(ncid, condensationid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, condensationid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);
            #endif
            #endif
            #if defined(P3_MICROPHY)
                if ((retval = nc_put_att_text(ncid, qitotid, "standard_name", strlen("mass_fraction_of_cloud_ice_in_air"), "mass_fraction_of_cloud_ice_in_air"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, qitotid, "long_name", strlen("Total Ice Mixing Ratio"), "Total Ice Mixing Ratio"))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, qitotid, "units", strlen("kg kg-1"), "kg kg-1"))) NC_ERR(retval);
                if ((retval = nc_put_att_double(ncid, qitotid, "_FillValue", NC_DOUBLE, 1, &fill_value))) NC_ERR(retval);
                if ((retval = nc_put_att_text(ncid, qitotid, "coordinates", strlen("t x z"), "t x z"))) NC_ERR(retval);
            #endif
        #endif

        // End define mode
        if ((retval = nc_enddef(ncid))) NC_ERR(retval);
        
        // Write coordinate data (assuming uniform grid for simplicity)
        if ((retval = nc_put_var_double(ncid, x_varid, model.x))) NC_ERR(retval);
        if ((retval = nc_put_var_double(ncid, z_varid, model.z))) NC_ERR(retval);
    }

    if ((retval = nc_put_var1_double(ncid, t_varid, &t_index, &t))) NC_ERR(retval);

    // Write the 2D data for this time step
    size_t start[3] = {t_index, 0, 0}; // Starting point (time_index, x, y)
    size_t count[3] = {1, (size_t) model.nx, (size_t) model.nz}; // Write one time slice, all x and y values
    #if defined(KESSLER_MICROPHY) || defined(P3_MICROPHY)
        size_t start_precip[2] = {t_index, 0};
        size_t count_precip[2] = {1, (size_t) model.nx};
    #endif
    if ((retval = nc_put_vara_double(ncid, th_id, start, count, model.thcont))) NC_ERR(retval);
    if ((retval = nc_put_vara_double(ncid, zeta_id, start, count, model.zetacont))) NC_ERR(retval);
    if ((retval = nc_put_vara_double(ncid, u_id, start, count, model.ucont))) NC_ERR(retval);
    if ((retval = nc_put_vara_double(ncid, w_id, start, count, model.wcont))) NC_ERR(retval);
    #if defined(RTERRTMGP)
        if ((retval = nc_put_vara_double(ncid, rad_net_heat_rate_id, start, count, model.rad_net_heat_ratecont))) NC_ERR(retval);
        if ((retval = nc_put_vara_double(ncid, rad_lw_heat_rate_id, start, count, model.rad_lw_heat_ratecont))) NC_ERR(retval);
        if ((retval = nc_put_vara_double(ncid, rad_sw_heat_rate_id, start, count, model.rad_sw_heat_ratecont))) NC_ERR(retval);
    #endif

    if ((retval = nc_put_var1_double(ncid, ubarTop_id, &t_index, &model.ubarTopp))) checkErr(retval, __LINE__);
    #if defined(SFCFLX)
        if ((retval = nc_put_vara_double(ncid, heatflux_id, start_precip, count_precip, model.heatflux))) checkErr(retval, __LINE__);
        if ((retval = nc_put_vara_double(ncid, waterflux_id, start_precip, count_precip, model.waterflux))) checkErr(retval, __LINE__);
        if ((retval = nc_put_vara_double(ncid, momentumflux_id, start_precip, count_precip, model.momentumflux))) checkErr(retval, __LINE__);
    #endif

    #if defined(WATER)
        if ((retval = nc_put_vara_double(ncid, qvid, start, count, model.qvcont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_vara_double(ncid, qcid, start, count, model.qccont))) checkErr(retval, __LINE__);
        if ((retval = nc_put_vara_double(ncid, qrid, start, count, model.qrcont))) checkErr(retval, __LINE__);
        #if defined(KESSLER_MICROPHY)
           if ((retval = nc_put_vara_double(ncid, precipid, start_precip, count_precip, model.precip))) checkErr(retval, __LINE__);
        #elif defined(P3_MICROPHY)
           if ((retval = nc_put_vara_double(ncid, precipid, start_precip, count_precip, model.diag_2dcont))) checkErr(retval, __LINE__);
        #endif
        #if defined(KESSLER_MICROPHY)
        #if defined(OUTPUTMICROPHYSICS)
            if ((retval = nc_put_vara_double(ncid, accretionid, start, count, model.accretioncont))) checkErr(retval, __LINE__);
            if ((retval = nc_put_vara_double(ncid, autoconversionid, start, count, model.autoconversioncont))) checkErr(retval, __LINE__);
            if ((retval = nc_put_vara_double(ncid, evaporationid, start, count, model.evaporationcont))) checkErr(retval, __LINE__);
            if ((retval = nc_put_vara_double(ncid, condensationid, start, count, model.condensationcont))) checkErr(retval, __LINE__);
        #endif
        #endif

        #if defined(P3_MICROPHY)
            // if ((retval = nc_put_vara_double(ncid, qncid, start, count, model.nccont))) checkErr(retval, __LINE__);
            // if ((retval = nc_put_vara_double(ncid, qnrid, start, count, model.nrcont))) checkErr(retval, __LINE__);
            // if ((retval = nc_put_vara_double(ncid, qniid, start, count, model.nicont))) checkErr(retval, __LINE__);
            if ((retval = nc_put_vara_double(ncid, qitotid, start, count, model.qitotcont))) checkErr(retval, __LINE__);
            // if ((retval = nc_put_vara_double(ncid, qirimid, start, count, model.qirimcont))) checkErr(retval, __LINE__);
            // if ((retval = nc_put_vara_double(ncid, birimid, start, count, model.birimcont))) checkErr(retval, __LINE__);
        #endif
    #endif

    // Close the file
    if ((retval = nc_close(ncid))) NC_ERR(retval);
}
#endif

void vvm::Output::output_time_nc(int n, vvm &model) {
    string ncName = model.outputpath + (string) "timer/" + std::to_string(n) + (string) ".nc";

    int ncid, t_dimid;
    int retval;

    int advectionid, poissonid, diffusionid, microphysicsid, allid;

    if ((retval = nc_create(ncName.c_str(), NC_CLOBBER, &ncid))) checkErr(retval, __LINE__);

    if ((retval = nc_def_dim(ncid, "x", model.nx, &t_dimid))) checkErr(retval, __LINE__);

    int dimt1d[1] = {t_dimid};

    if ((retval = nc_def_var(ncid, "advection", NC_DOUBLE, 1, dimt1d, &advectionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "poisson", NC_DOUBLE, 1, dimt1d, &poissonid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "diffusion", NC_DOUBLE, 1, dimt1d, &diffusionid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "microphysics", NC_DOUBLE, 1, dimt1d, &microphysicsid))) checkErr(retval, __LINE__);
    if ((retval = nc_def_var(ncid, "all", NC_DOUBLE, 1, dimt1d, &allid))) checkErr(retval, __LINE__);

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


void vvm::Output::copy_files(const std::string &source_path, const std::string &destination_path) {
    string str = "cp -r " + source_path + " " + destination_path;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void vvm::Output::grads_ctl_file(vvm &model) {
    // Open output file
    std::ofstream outFile(model.outputpath + "/nc/vvm.ctl");
    if (!outFile.is_open()) {
        std::cerr << "Error opening ctl file!" << std::endl;
        return;
    }

    // Write the .ctl file content
    outFile << "DSET ^0.nc\n";
    outFile << "DTYPE netcdf\n";
    outFile << "OPTIONS template\n";
    outFile << "TITLE NetCDF Data for GrADS\n";
    outFile << "UNDEF -9999.0\n";
    outFile << "XDEF " << model.nx << " LINEAR 0 1\n";
    outFile << "YDEF " << 1 << " LINEAR 0.0 1.0\n";
    outFile << "ZDEF " << model.nz << " LEVELS ";
    for (int k = 0; k < model.nz; k++) outFile << static_cast<int> (model.z[k]) << ", ";
    outFile << "\n";

    int outnum = 6;
    #if defined(WATER)
        outnum += 4;
    #endif
    #if defined(RTERRTMGP)
        outnum += 3;
    #endif
    #if defined(SFCFLX)
        outnum += 3;
    #endif

    outFile << "TDEF " << (int) model.TIMEEND / (model.dt*model.OUTPUTSTEP) << " LINEAR 00:00Z01JAN2000 " << "1hr\n";
    outFile << "\n";
    outFile << "VARS " << outnum << "\n";
    outFile << "th=>th " << model.nz << " t,x,z theta\n";
    outFile << "u=>u " << model.nz << " t,x,z u\n";
    outFile << "w=>w " << model.nz << " t,x,z w\n";
    outFile << "zeta=>zeta " << model.nz << " t,x,z zeta\n";
    outFile << "qv=>qv " << model.nz << " t,x,z qv\n";
    #if defined(WATER)
        outFile << "qc=>qc " << model.nz << " t,x,z qc\n";
        outFile << "qr=>qr " << model.nz << " t,x,z qr\n";
        outFile << "qitot=>qitot " << model.nz << " t,x,z qitot\n";
        outFile << "precip=>precip 1 t,x precip\n";
    #endif
    #if defined(RTERRTMGP)
        outFile << "rad_net_heat_rate=>rad_net " << model.nz << " t,x,z rhr\n";
        outFile << "rad_lw_heat_rate=>rad_lw " << model.nz << " t,x,z rhr\n";
        outFile << "rad_lw_heat_rate=>rad_sw " << model.nz << " t,x,z rhr\n";
    #endif
    #if defined(SFCFLX)
        outFile << "heatflux_sfc=>htflx 1 t,x heat flux\n";
        outFile << "waterflux_sfc=>wtflx 1 t,x water flux\n";
        outFile << "momentumflux_sfc=>mtflx 1 t,x momentum flux\n";
    #endif
    outFile << "ubarTop=>ubarTop 1 t ubarTop\n";
    outFile << "ENDVARS\n";

    // Close the file
    outFile.close();
    return;
}


void vvm::Output::create_all_directory(vvm &model) {
    // data directory
    #ifdef OUTPUTNC
        create_directory(model.outputpath + (string) "nc");
        create_directory(model.outputpath + (string) "timer");
        create_directory(model.outputpath + (string) "run_files");
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
    // create_directory(model.outputpath + (string) "graphs");
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

void vvm::Output::copy_source_project(vvm &model) {
    // Copy source files to output directory
    vvm::Output::copy_files("../src", model.outputpath+"run_files/.");
    vvm::Output::copy_files("../input", model.outputpath+"run_files/.");
    vvm::Output::copy_files("../external", model.outputpath+"run_files/.");
    vvm::Output::copy_files("../vvm_config.txt", model.outputpath+"run_files/.");
    vvm::Output::copy_files("../CMakeLists.txt", model.outputpath+"run_files/.");
    vvm::Output::copy_files("../run.sh", model.outputpath+"run_files/.");
}
