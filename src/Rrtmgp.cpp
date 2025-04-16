#include "Declare.hpp"

#if defined(RTERRTMGP)
#include <omp.h>


bool vvm::Radiation::is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

int vvm::Radiation::day_of_year(int year, int month, int day) {
    const int months[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int doy = months[month - 1] + day;
    if (month > 2 && is_leap_year(year)) {
        doy += 1;
    }
    return doy;
}

double vvm::Radiation::calculate_scaling_factor(int year, int month, int day) {
    // This function can be used when you want to calculate the change of solar constant according to the distance between earth and sun.

    // Spencer (1971)
    int doy = day_of_year(year, month, day);

    // Calculate gamma = 2π * (doy - 1) / 365
    double gamma = 2 * M_PI * (doy - 1) / 365.0;

    // Spencer (1971) formula for scaling factor
    double scaling_factor = 1.000110 + 0.034221 * cos(gamma) + 0.001280 * sin(gamma) +
                            0.000719 * cos(2 * gamma) + 0.000077 * sin(2 * gamma);

    return scaling_factor;
}

double vvm::Radiation::calculate_cos_zenith(int year, int month, int day, double hour, double minute, double second,
                                            double longitude, double latitude) {
    // Calculate day of the year
    int doy = day_of_year(year, month, day);

    // Calculate solar declination
    double gamma = 2 * M_PI * (doy - 1) / 365.0;
    double delta = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
                   0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
                   0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma);

    // Calculate equation of time (in minutes)
    double eot = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
                           0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma));

    // Convert UTC time to hours
    double utc_hours = hour + minute / 60.0 + second / 3600.0;

    // Calculate local solar time (LST) in hours
    double lst = utc_hours + (longitude / 15.0) + (eot / 60.0);

    // Calculate hour angle in degrees
    double h = 15.0 * (lst - 12.0);

    // Convert angles to radians
    double phi_rad = latitude * M_PI / 180.0;
    double h_rad = h * M_PI / 180.0;

    // Calculate cosine of solar zenith angle
    double cos_zenith = sin(phi_rad) * sin(delta) + cos(phi_rad) * cos(delta) * cos(h_rad);

    return cos_zenith;
}


void vvm::Radiation::solve_radiation(vvm &model) {
    auto time_start = std::chrono::high_resolution_clock::now();
    Status::print_message("###### Starting RTE+RRTMGP solver ######");
    ////// FLOW CONTROL SWITCHES //////
    const bool switch_shortwave         = true;  // Enable computation of shortwave radiation.
    const bool switch_longwave          = true;  // Enable computation of longwave radiation.
    const bool switch_fluxes            = true;  // Enable computation of fluxes.
    const bool switch_cloud_optics      = true;  // Enable cloud optics.
    const bool switch_aerosol_optics    = false; // Enable aerosol optics.
    const bool switch_output_optical    = false; // Enable output of optical properties.
    const bool switch_output_bnd_fluxes = false; // Enable output of band fluxes.
    const bool switch_delta_cloud       = true;  // delta-scaling of cloud optical properties
    const bool switch_delta_aerosol     = false; // delta-scaling of aerosol optical properties
    const bool switch_col_dry           = false;

    ////// READ THE ATMOSPHERIC DATA //////
    const int n_col_x = model.nx-2;
    const int n_col_y = 1;
    const int n_col = 1;
    const int n_lay = model.nz-2;
    const int n_lev = n_lay+1;

    const double qvmin = 1e-6;
    const double qcmin = 1e-7;
    const double qimin = 1e-8;

    const double effcmin = 2.5;  // unit: 1e-6m
    const double effcmax = 21.5; // unit: 1e-6m
    const double effimin = 10;   // unit: 1e-6m
    const double effimax = 180;  // unit: 1e-6m
    // const double rel_val = 0.5*(effcmax+effcmin); // unit: 1e-6m
    // const double rei_val = 0.5*(effimax+effimin); // unit: 1e-6m

    const double sec_per_day = 86400.0; // Seconds per day
    // !     Calculated value (from constants above and input cpdair)
    // !     (grav) x (#sec/day) / (specific heat of dry air at const. p x 1.e2)
    // !     Here, cpdair is in units of J kg-1 K-1, and the constant (1.e2) 
    // !     converts mb to Pa when heatfac is multiplied by W m-2 mb-1.
    const double heat_factor = model.GRAVITY * sec_per_day / model.Cp; 


    ////// STORAGE ARRAYS FOR ALL COLUMNS //////
    std::vector<Array<Float,2>> p_lay_all(n_col_x);
    std::vector<Array<Float,2>> p_lev_all(n_col_x);
    
    // Longwave storage
    std::vector<Array<Float,2>> lw_flux_up_all(n_col_x);
    std::vector<Array<Float,2>> lw_flux_dn_all(n_col_x);
    std::vector<Array<Float,2>> lw_flux_net_all(n_col_x);
    std::vector<Array<Float,3>> lw_tau_all(n_col_x);
    std::vector<Array<Float,3>> lay_source_all(n_col_x);
    std::vector<Array<Float,3>> lev_source_inc_all(n_col_x);
    std::vector<Array<Float,3>> lev_source_dec_all(n_col_x);
    std::vector<Array<Float,2>> sfc_source_all(n_col_x);
    std::vector<Array<Float,3>> lw_bnd_flux_up_all(n_col_x);
    std::vector<Array<Float,3>> lw_bnd_flux_dn_all(n_col_x);
    std::vector<Array<Float,3>> lw_bnd_flux_net_all(n_col_x);

    // Shortwave storage
    std::vector<Array<Float,2>> sw_flux_up_all(n_col_x);
    std::vector<Array<Float,2>> sw_flux_dn_all(n_col_x);
    std::vector<Array<Float,2>> sw_flux_dn_dir_all(n_col_x);
    std::vector<Array<Float,2>> sw_flux_net_all(n_col_x);
    std::vector<Array<Float,3>> sw_tau_all(n_col_x);
    std::vector<Array<Float,3>> ssa_all(n_col_x);
    std::vector<Array<Float,3>> g_all(n_col_x);
    std::vector<Array<Float,2>> toa_source_all(n_col_x);
    std::vector<Array<Float,3>> sw_bnd_flux_up_all(n_col_x);
    std::vector<Array<Float,3>> sw_bnd_flux_dn_all(n_col_x);
    std::vector<Array<Float,3>> sw_bnd_flux_dn_dir_all(n_col_x);
    std::vector<Array<Float,3>> sw_bnd_flux_net_all(n_col_x);

    ////// INITIALIZE SOLVERS ONCE WITH DUMMY GAS_CONCS //////
    Gas_concs dummy_gas_concs;
    // Fully representative dummy gas concentrations (matching loop defaults)
    double dummy_h2o[n_lay]; // Will be updated per column, use min value as placeholder
    for (int k = 0; k < n_lay; k++) dummy_h2o[k] = std::max(model.qvb[k], 0.);
    Array<Float, 2> dummy_h2o_lay(dummy_h2o, {n_col, n_lay}, true);

    double g1 = 3.6478, g2 = 0.83209, g3 = 11.3515;
    double o3_min = 1e-13;
    double dummy_o3[n_lay];
    for (int k = 0; k < n_lay; k++) {
        double p_hpa = model.pb[k + 1] / 100.;
        dummy_o3[k] = std::max(o3_min, g1 * std::pow(p_hpa, g2) * std::exp(-p_hpa / g3) * 1e-6);
    }
    Array<Float, 2> dummy_o3_lay(dummy_o3, {n_col, n_lay}, true);

    dummy_gas_concs.set_vmr("h2o", dummy_h2o_lay); // Minimum water vapor as placeholder
    dummy_gas_concs.set_vmr("co2", 348e-6);        // ppm
    dummy_gas_concs.set_vmr("o3", dummy_o3_lay);   // Calculated ozone profile
    dummy_gas_concs.set_vmr("n2o", 306e-9);        // ppm
    dummy_gas_concs.set_vmr("ch4", 1650e-9);       // ppm
    dummy_gas_concs.set_vmr("o2", 0.2095);         // Volume mixing ratio
    dummy_gas_concs.set_vmr("n2", 0.7808);         // Volume mixing ratio


    Radiation_solver_longwave* rad_lw = nullptr;
    int n_bnd_lw = 0, n_gpt_lw = 0;
    if (switch_longwave) {
        rad_lw = new Radiation_solver_longwave(
            dummy_gas_concs,
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-lw-g256.nc",
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-lw.nc");
        n_bnd_lw = rad_lw->get_n_bnd();
        n_gpt_lw = rad_lw->get_n_gpt();
    }

    Radiation_solver_shortwave* rad_sw = nullptr;
    int n_bnd_sw = 0, n_gpt_sw = 0;
    if (switch_shortwave) {
        rad_sw = new Radiation_solver_shortwave(
            dummy_gas_concs, switch_cloud_optics, switch_aerosol_optics,
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-sw-g224.nc",
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-sw.nc",
            "../external/rte-rrtmgp-cpp/allsky/aerosol_optics.nc");
        n_bnd_sw = rad_sw->get_n_bnd();
        n_gpt_sw = rad_sw->get_n_gpt();
    }
    // Get solar constant to calculate tsi scaling factor: About 1360.8 W/m^2
    const Float tsi_ref = rad_sw->get_tsi();
    const double tsi = 340.; // total_solar_irradiance 
    // const double tsi = tsi_ref; // total_solar_irradiance 

    ////// COMMON DATA OUTSIDE LOOP //////
    // Initialize O3
    double p_hpa[n_lay];
    double o3[n_lay];
    double tmp_emis_sfc[n_lay];
    for (int k = 0; k < n_lay; k++) {
        p_hpa[k] = model.pb[k]/100.;
        o3[k] = std::max(o3_min, g1 * std::pow(p_hpa[k], g2) * std::exp(-p_hpa[k]/g3) * 1e-6);
        tmp_emis_sfc[k] = 0.98;
    }

    // Initialize mu0 [cos(Solar Zenith Angle)]
    double vvm_mu0 = calculate_cos_zenith(model.year, model.month, model.day, model.hour, model.minute, model.second, model.lon, model.lat);
    // double vvm_mu0 = cos(60.*M_PI/180);
    if (vvm_mu0 < 0) vvm_mu0 = 0.; // Set to 0 if the sun is below the sky.
    double tmp_mu0[n_col];
    for (int i = 0; i < n_col; i++) tmp_mu0[i] = vvm_mu0;

    // Initialize surface albedo related coefficient
    double *tmp_sfc_alb_dircont = nullptr;
    double **tmp_sfc_alb_dir = nullptr;
    double *tmp_sfc_alb_difcont = nullptr;
    double **tmp_sfc_alb_dif = nullptr;
    tmp_sfc_alb_dir = allocate2DContinuousArray(16, n_col, tmp_sfc_alb_dircont);
    tmp_sfc_alb_dif = allocate2DContinuousArray(16, n_col, tmp_sfc_alb_difcont);
    for (int bnd = 0; bnd < 16; bnd++) {
        for (int i = 0; i < n_col; i++) {
            tmp_sfc_alb_dif[bnd][i] = 0.06;
            tmp_sfc_alb_dir[bnd][i] = 0.06;
        }
    }

    // Precompute T and T_lev
    #pragma omp parallel for
    for (int i = 1; i <= n_col_x; i++) {
        for (int k = 1; k < model.nz - 1; k++) {
            model.T[i][k] = model.th[i][k] * model.pib[k];
        }
        for (int k = 0; k < model.nz - 1; k++) {
            model.T_lev[i][k] = 0.5 * (model.th[i][k] + model.th[i][k + 1]) * model.pib_lev[k];
        }
        model.T_lev[i][model.nz - 1] = model.th[i][model.nz - 2] * model.pib_lev[model.nz - 1];
        model.T_lev[i][model.nz] = model.T_lev[i][model.nz - 1];
    }
    model.BoundaryProcess2D_center(model.T, model.nx, model.nz);


    for (int i = 1; i <= n_col_x; i++) {
        for (int k = 0; k < model.nz; k++) {
            if (std::isnan(model.th[i][k]) || std::isnan(model.qv[i][k]) ||
                std::isnan(model.qc[i][k]) || std::isnan(model.qitot[i][k]) ||
                std::isnan(model.diag_effc[i][k]) || std::isnan(model.diag_effi[i][k]) ||
                std::isnan(model.rhou[k])) {
                std::cerr << "NaN in model input at i=" << i << ", k=" << k << std::endl;
                exit(1);
            }
            if (model.rhou[k] <= 0) {
                std::cerr << "Non-positive density at k=" << k << ": rhou=" << model.rhou[k] << std::endl;
                exit(1);
            }
        }
    }

    double *tmp_pb = model.pb+1;
    double *tmp_pb_lev = model.pb_lev+1;

#pragma omp parallel \
    shared(model, rad_lw, rad_sw, n_bnd_lw, n_gpt_lw, n_bnd_sw, n_gpt_sw, \
            p_lay_all, p_lev_all, lw_flux_up_all, lw_flux_dn_all, lw_flux_net_all, \
            lw_tau_all, lay_source_all, lev_source_inc_all, lev_source_dec_all, sfc_source_all, \
            lw_bnd_flux_up_all, lw_bnd_flux_dn_all, lw_bnd_flux_net_all, \
            sw_flux_up_all, sw_flux_dn_all, sw_flux_dn_dir_all, sw_flux_net_all, \
            sw_tau_all, ssa_all, g_all, toa_source_all, \
            sw_bnd_flux_up_all, sw_bnd_flux_dn_all, sw_bnd_flux_dn_dir_all, sw_bnd_flux_net_all, \
            tmp_emis_sfc, tmp_mu0, tmp_sfc_alb_dircont, tmp_sfc_alb_difcont, o3)
{
    int tid = omp_get_thread_num();

    ////// LOOP OVER X COLUMNS //////
    #pragma omp for
    for (int i = 1; i <= n_col_x; i++) {
        double *tmp_T = model.T[i]+1;
        double *tmp_T_lev = model.T_lev[i]+1;

        Array<Float,2> p_lay(tmp_pb, {n_col, n_lay}, true); // Start from k = 1 to k = model.nz-2
        Array<Float,2> t_lay(tmp_T, {n_col, n_lay}, true); // Start from k = 1 to k = model.nz-2
        Array<Float,2> p_lev(tmp_pb_lev, {n_col, n_lev}, true); // Start from k = 1 to k = model.nz-1
        Array<Float,2> t_lev(tmp_T_lev, {n_col, n_lev}, true); // Start from k = 1 to k = model.nz-1
        
        p_lay_all[i-1] = p_lay;
        p_lev_all[i-1] = p_lev;
        
        // Not needed in this case
        Array<Float,2> col_dry;
        if (switch_col_dry) {
            col_dry.set_dims({n_col, n_lay});
            // col_dry = std::move(input_nc.get_variable<Float>("col_dry", {n_lay, n_col_y, n_col_x}));
        }

        // Create container for the gas concentrations and read gases.
        double *tmp_qv = model.qv[i]+1;
        double filter_qv[n_lay];
        for (int k = 0; k < n_lay; k++) {
            filter_qv[k] = std::max(tmp_qv[k], qvmin);
        }

        Array<Float,2> h2o_lay(filter_qv, {n_col, n_lay}, true); // Start from k = 1 to k = model.nz-2
        Array<Float,2> o3_lay(o3, {n_col, n_lay}, true); // Start from k = 1 to k = model.nz-2
        Gas_concs gas_concs;

        gas_concs.set_vmr("h2o", h2o_lay);
        gas_concs.set_vmr("co2", 348e-6); // vmr (ppm)
        gas_concs.set_vmr("o3", o3_lay);
        gas_concs.set_vmr("n2o", 306e-9); // vmr (ppm)
        gas_concs.set_vmr("ch4", 1650e-9); // vmr (ppm)
        gas_concs.set_vmr("o2", 0.2095);
        gas_concs.set_vmr("n2", 0.7808);

        // Not need to input
        // gas_concs.set_vmr("co", 306e-9); // no input co
        // read_and_set_vmr("ccl4"   , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("cfc11"  , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("cfc12"  , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("cfc22"  , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("hfc143a", n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("hfc125" , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("hfc23"  , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("hfc32"  , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("hfc134a", n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("cf4"    , n_col_x, n_col_y, n_lay, input_nc, gas_concs);
        // read_and_set_vmr("no2"    , n_col_x, n_col_y, n_lay, input_nc, gas_concs);

        // Reinitialize kdist with updated Gas_concs
        // #pragma omp critical
        // if (switch_longwave) rad_lw->reinit_gas_optics(gas_concs);
        // #pragma omp critical
        // if (switch_shortwave) rad_sw->reinit_gas_optics(gas_concs);

        double *tmp_rho = model.rhou+1;
        double *tmp_qc = model.qc[i]+1;
        double *tmp_qi = model.qitot[i]+1;
        double *tmp_effc = model.diag_effc[i]+1;
        double *tmp_effi = model.diag_effi[i]+1;
        double filter_qc[n_lay];
        double filter_qi[n_lay];
        double filter_effc[n_lay];
        double filter_effi[n_lay];
        double tmp_lwp[n_lay];
        double tmp_iwp[n_lay];
        for (int k = 0; k < n_lay; k++) {
            filter_qc[k] = std::max(tmp_qc[k]*1e3, 0.); // unit: g/kg
            filter_qi[k] = std::max(tmp_qi[k]*1e3, 0.); // unit: g/kg
            tmp_lwp[k] = tmp_rho[k] * filter_qc[k] * model.dz;
            tmp_iwp[k] = tmp_rho[k] * filter_qi[k] * model.dz;
            filter_effc[k] = std::max(tmp_effc[k]*1e6, effcmin); // unit: 1e-6 m
            filter_effi[k] = std::max(tmp_effi[k]*1e6, effimin); // unit: 1e-6 m
            filter_effc[k] = std::min(filter_effc[k], effcmax); // unit: 1e-6 m
            filter_effi[k] = std::min(filter_effi[k], effimax); // unit: 1e-6 m
        }


        Array<Float,2> lwp;
        Array<Float,2> iwp;
        Array<Float,2> rel;
        Array<Float,2> rei;
        if (switch_cloud_optics) {
            lwp.set_dims({n_col, n_lay});
            for (int ic = 0; ic < n_col; ++ic) {
                for (int k = 0; k < n_lay; ++k) {
                    lwp({ic + 1, k + 1}) = tmp_lwp[k]; // Replicate tmp_lwp across all columns
                }
            }

            iwp.set_dims({n_col, n_lay});
            for (int ic = 0; ic < n_col; ++ic) {
                for (int k = 0; k < n_lay; ++k) {
                    iwp({ic + 1, k + 1}) = tmp_iwp[k]; // Replicate tmp_iwp across all columns
                }
            }

            rel.set_dims({n_col, n_lay});
            for (int ic = 0; ic < n_col; ++ic) {
                for (int k = 0; k < n_lay; ++k) {
                    rel({ic + 1, k + 1}) = filter_effc[k]; // Replicate tmp_iwp across all columns
                }
            }

            rei.set_dims({n_col, n_lay});
            for (int ic = 0; ic < n_col; ++ic) {
                for (int k = 0; k < n_lay; ++k) {
                    rei({ic + 1, k + 1}) = filter_effi[k]; // Replicate tmp_iwp across all columns
                }
            }
        }


        // No aerosol in this run
        Array<Float,2> rh;
        Aerosol_concs aerosol_concs;
        // if (switch_aerosol_optics)
        // {
        //     rh.set_dims({n_col, n_lay});
        //     rh = std::move(input_nc.get_variable<Float>("rh", {n_lay, n_col_y, n_col_x}));
        //
        //     read_and_set_aer("aermr01", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr02", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr03", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr04", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr05", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr06", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr07", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr08", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr09", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr10", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        //     read_and_set_aer("aermr11", n_col_x, n_col_y, n_lay, input_nc, aerosol_concs);
        // }


        ////// RUN THE LONGWAVE SOLVER //////
        if (switch_longwave) {
            // Initialize the solver.
            Array<Float,2> emis_sfc(tmp_emis_sfc, {n_col, n_lay}, true); // Start from k = 1 to k = model.nz-2
            Array<Float,1> t_sfc(tmp_T_lev, {n_col});

            // Create output arrays.
            Array<Float,3> lw_tau;
            Array<Float,3> lay_source;
            Array<Float,3> lev_source_inc;
            Array<Float,3> lev_source_dec;
            Array<Float,2> sfc_source;

            if (switch_output_optical) {
                lw_tau        .set_dims({n_col, n_lay, n_gpt_lw});
                lay_source    .set_dims({n_col, n_lay, n_gpt_lw});
                lev_source_inc.set_dims({n_col, n_lay, n_gpt_lw});
                lev_source_dec.set_dims({n_col, n_lay, n_gpt_lw});
                sfc_source    .set_dims({n_col, n_gpt_lw});
            }

            Array<Float,2> lw_flux_up, lw_flux_dn, lw_flux_net;
            if (switch_fluxes) {
                lw_flux_up .set_dims({n_col, n_lev});
                lw_flux_dn .set_dims({n_col, n_lev});
                lw_flux_net.set_dims({n_col, n_lev});
            }

            Array<Float,3> lw_bnd_flux_up, lw_bnd_flux_dn, lw_bnd_flux_net;
            if (switch_output_bnd_fluxes) {
                lw_bnd_flux_up .set_dims({n_col, n_lev, n_bnd_lw});
                lw_bnd_flux_dn .set_dims({n_col, n_lev, n_bnd_lw});
                lw_bnd_flux_net.set_dims({n_col, n_lev, n_bnd_lw});
            }


            // Solve the radiation.
        #pragma omp critical (longwave)
        {
            rad_lw->solve(
                    switch_fluxes,
                    switch_cloud_optics,
                    switch_output_optical,
                    switch_output_bnd_fluxes,
                    gas_concs,
                    p_lay, p_lev,
                    t_lay, t_lev,
                    col_dry,
                    t_sfc, emis_sfc,
                    lwp, iwp,
                    rel, rei,
                    lw_tau, lay_source, lev_source_inc, lev_source_dec, sfc_source,
                    lw_flux_up, lw_flux_dn, lw_flux_net,
                    lw_bnd_flux_up, lw_bnd_flux_dn, lw_bnd_flux_net
            );
        }

            // Store results
            if (switch_fluxes) {
                lw_flux_up_all[i-1] = lw_flux_up;
                lw_flux_dn_all[i-1] = lw_flux_dn;
                lw_flux_net_all[i-1] = lw_flux_net;
            }
            if (switch_output_optical) {
                lw_tau_all[i-1] = lw_tau;
                lay_source_all[i-1] = lay_source;
                lev_source_inc_all[i-1] = lev_source_inc;
                lev_source_dec_all[i-1] = lev_source_dec;
                sfc_source_all[i-1] = sfc_source;
            }
            if (switch_output_bnd_fluxes) {
                lw_bnd_flux_up_all[i-1] = lw_bnd_flux_up;
                lw_bnd_flux_dn_all[i-1] = lw_bnd_flux_dn;
                lw_bnd_flux_net_all[i-1] = lw_bnd_flux_net;
            }
        }


        ////// RUN THE SHORTWAVE SOLVER //////
        if (switch_shortwave) {
            // Initialize the solver.
            Array<Float,1> mu0(tmp_mu0, {n_col});
            Array<Float,2> sfc_alb_dir(tmp_sfc_alb_dircont, {n_bnd_sw, n_col}, true);
            Array<Float,2> sfc_alb_dif(tmp_sfc_alb_difcont, {n_bnd_sw, n_col}, true);

            Array<Float,1> tsi_scaling({n_col});
            tsi_scaling({1}) = tsi / tsi_ref;

            // Create output arrays.
            Array<Float,3> sw_tau;
            Array<Float,3> ssa;
            Array<Float,3> g;
            Array<Float,2> toa_source;

            if (switch_output_optical) {
                sw_tau    .set_dims({n_col, n_lay, n_gpt_sw});
                ssa       .set_dims({n_col, n_lay, n_gpt_sw});
                g         .set_dims({n_col, n_lay, n_gpt_sw});
                toa_source.set_dims({n_col, n_gpt_sw});
            }

            Array<Float,2> sw_flux_up;
            Array<Float,2> sw_flux_dn;
            Array<Float,2> sw_flux_dn_dir;
            Array<Float,2> sw_flux_net;

            if (switch_fluxes) {
                sw_flux_up    .set_dims({n_col, n_lev});
                sw_flux_dn    .set_dims({n_col, n_lev});
                sw_flux_dn_dir.set_dims({n_col, n_lev});
                sw_flux_net   .set_dims({n_col, n_lev});
            }

            Array<Float,3> sw_bnd_flux_up;
            Array<Float,3> sw_bnd_flux_dn;
            Array<Float,3> sw_bnd_flux_dn_dir;
            Array<Float,3> sw_bnd_flux_net;

            if (switch_output_bnd_fluxes) {
                sw_bnd_flux_up    .set_dims({n_col, n_lev, n_bnd_sw});
                sw_bnd_flux_dn    .set_dims({n_col, n_lev, n_bnd_sw});
                sw_bnd_flux_dn_dir.set_dims({n_col, n_lev, n_bnd_sw});
                sw_bnd_flux_net   .set_dims({n_col, n_lev, n_bnd_sw});
            }

        for (int k = 0; k < n_lay; k++) {
            if (std::isnan(tmp_lwp[k]) || std::isnan(tmp_iwp[k]) || 
            std::isnan(filter_effc[k]) || std::isnan(filter_effi[k])) {
                std::cerr << "NaN detected at i=" << i << ", k=" << k 
                        << ": lwp=" << tmp_lwp[k] << ", iwp=" << tmp_iwp[k]
                        << ", rel=" << filter_effc[k] << ", rei=" << filter_effi[k] << std::endl;
                exit(1);
            }
            if (filter_effc[k] < 0 || filter_effi[k] < 0) {
                std::cerr << "Zero or negative radius at i=" << i << ", k=" << k 
                        << ": rel=" << filter_effc[k] << ", rei=" << filter_effi[k] << std::endl;
                exit(1);
            }
        }

            // Solve the radiation.
        #pragma omp critical (shortwave)
        {
            rad_sw->solve(
                    switch_fluxes,
                    switch_cloud_optics,
                    switch_aerosol_optics,
                    switch_output_optical,
                    switch_output_bnd_fluxes,
                    switch_delta_cloud,
                    switch_delta_aerosol,
                    gas_concs,
                    p_lay, p_lev,
                    t_lay, t_lev,
                    col_dry,
                    sfc_alb_dir, sfc_alb_dif,
                    tsi_scaling, mu0,
                    lwp, iwp,
                    rel, rei,
                    rh,
                    aerosol_concs,
                    sw_tau, ssa, g,
                    toa_source,
                    sw_flux_up, sw_flux_dn,
                    sw_flux_dn_dir, sw_flux_net,
                    sw_bnd_flux_up, sw_bnd_flux_dn,
                    sw_bnd_flux_dn_dir, sw_bnd_flux_net
            );
        }
            for (int lev = 0; lev < n_lev; lev++) {
                if (std::isnan(sw_flux_up({1, lev+1}))) {
                    std::cerr << "NaN in sw_flux_up at i=" << i << ", lev=" << lev << std::endl;
                    exit(1);
                }
            }

            // Store results
            if (switch_fluxes) {
                sw_flux_up_all[i-1] = sw_flux_up;
                sw_flux_dn_all[i-1] = sw_flux_dn;
                sw_flux_dn_dir_all[i-1] = sw_flux_dn_dir;
                sw_flux_net_all[i-1] = sw_flux_net;
            }
            if (switch_output_optical) {
                sw_tau_all[i-1] = sw_tau;
                ssa_all[i-1] = ssa;
                g_all[i-1] = g;
                toa_source_all[i-1] = toa_source;
            }
            if (switch_output_bnd_fluxes) {
                sw_bnd_flux_up_all[i-1] = sw_bnd_flux_up;
                sw_bnd_flux_dn_all[i-1] = sw_bnd_flux_dn;
                sw_bnd_flux_dn_dir_all[i-1] = sw_bnd_flux_dn_dir;
                sw_bnd_flux_net_all[i-1] = sw_bnd_flux_net;
            }
        }
    }
}

    ////// COMPUTE HEATING RATE //////
    std::vector<Array<Float,2>> heating_rate_lw_all(n_col_x);
    std::vector<Array<Float,2>> heating_rate_sw_all(n_col_x);
    std::vector<Array<Float,2>> heating_rate_all(n_col_x);
    if (switch_fluxes && (switch_longwave || switch_shortwave)) {

        // Initialize heating rate arrays
        for (int i = 0; i < model.nx; ++i) {
            for (int k = 0; k < model.nz; ++k) {
                model.radiation_heating_rate[i][k] = 0.0;
            }
        }

        for (int i = 0; i < n_col_x; ++i) {
            int model_i = i + 1; // Map to model.radiation_heating_rate (i=0 -> model_i=1)
            Array<Float,2> heating_rate_lw({n_col, n_lay});
            Array<Float,2> heating_rate_sw({n_col, n_lay});
            Array<Float,2> heating_rate({n_col, n_lay});

            for (int k = 1; k <= n_lay; ++k) {
                int model_k = k;

                // Layer properties

                // Compute heating rates
                double dT_dt_lw = 0.0, dT_dt_sw = 0.0;

                if (switch_longwave) {
                    double fnet_lw_bot = lw_flux_net_all[i]({1, k});   // Interface k (bottom)
                    double fnet_lw_top = lw_flux_net_all[i]({1, k+1}); // Interface k+1 (top)
                    double p_bot = p_lev_all[i]({1, k});   // Interface k (bottom)
                    double p_top = p_lev_all[i]({1, k+1}); // Interface k+1 (top)
                    dT_dt_lw = heat_factor * (fnet_lw_top-fnet_lw_bot) / (p_bot-p_top); // Note the order
                    heating_rate_lw({1, k}) = dT_dt_lw;
                }

                if (switch_shortwave) {
                    double fnet_sw_bot = sw_flux_net_all[i]({1, k});   // Interface k (bottom)
                    double fnet_sw_top = sw_flux_net_all[i]({1, k+1}); // Interface k+1 (top)
                    double p_bot = p_lev_all[i]({1, k});   // Interface k (bottom)
                    double p_top = p_lev_all[i]({1, k+1}); // Interface k+1 (top)
                    dT_dt_sw = heat_factor * (fnet_sw_top-fnet_sw_bot) / (p_bot-p_top); // Note the order
                    heating_rate_sw({1, k}) = dT_dt_sw;
                }

                double dT_dt = dT_dt_lw + dT_dt_sw;
                if (k == n_lay) dT_dt = 0.; // top layer with 0 heating rate
                heating_rate({1, k}) = dT_dt;
                model.radiation_heating_rate[model_i][model_k] = dT_dt;

                // Debug checks
                if (std::isnan(dT_dt) || std::isinf(dT_dt)) {
                    std::cerr << "Invalid heating rate at i=" << model_i << ", k=" << model_k
                            << ": total=" << dT_dt * sec_per_day
                            << ", lw=" << dT_dt_lw * sec_per_day
                            << ", sw=" << dT_dt_sw * sec_per_day << std::endl;
                }
            }

            heating_rate_lw_all[i] = heating_rate_lw;
            heating_rate_sw_all[i] = heating_rate_sw;
            heating_rate_all[i] = heating_rate;
        }

        model.BoundaryProcess2D_center(model.radiation_heating_rate, model.nx, model.nz);
    }


    ////// NETCDF OUTPUT //////
    std::string output_nc_name = model.outputpath + "nc/Radiation_" + std::to_string(model.step) + ".nc";
    Status::print_message("Creating NetCDF file: " + output_nc_name);

    try {
        Netcdf_file output_nc(output_nc_name, Netcdf_mode::Create);

        // Define dimensions
        output_nc.add_dimension("x", n_col_x);
        output_nc.add_dimension("y", n_col_y);
        output_nc.add_dimension("lay", n_lay);
        output_nc.add_dimension("lev", n_lev);
        output_nc.add_dimension("pair", 2);

        if (switch_longwave) {
            output_nc.add_dimension("gpt_lw", n_gpt_lw);
            output_nc.add_dimension("band_lw", n_bnd_lw);
        }
        if (switch_shortwave) {
            output_nc.add_dimension("gpt_sw", n_gpt_sw);
            output_nc.add_dimension("band_sw", n_bnd_sw);
        }

        // Define variables
        auto nc_lay = output_nc.add_variable<Float>("p_lay", {"lay", "y", "x"});
        auto nc_lev = output_nc.add_variable<Float>("p_lev", {"lev", "y", "x"});
        auto nc_heating_rate = output_nc.add_variable<Float>("heating_rate", {"lay", "y", "x"});
        auto nc_lw_heating_rate = output_nc.add_variable<Float>("lw_heating_rate", {"lay", "y", "x"});
        auto nc_sw_heating_rate = output_nc.add_variable<Float>("sw_heating_rate", {"lay", "y", "x"});

        // Longwave variables
        std::map<std::string, Netcdf_variable<Float>*> lw_vars;
        if (switch_longwave) {
            auto nc_lw_band_lims_wvn = output_nc.add_variable<Float>("lw_band_lims_wvn", {"band_lw", "pair"});
            nc_lw_band_lims_wvn.insert(rad_lw->get_band_lims_wavenumber().v(), {0, 0}, {n_bnd_lw, 2});

            if (switch_fluxes) {
                lw_vars["lw_flux_up"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_flux_up", {"lev", "y", "x"}));
                lw_vars["lw_flux_dn"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_flux_dn", {"lev", "y", "x"}));
                lw_vars["lw_flux_net"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_flux_net", {"lev", "y", "x"}));
            }
            if (switch_output_optical) {
                lw_vars["lw_tau"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_tau", {"gpt_lw", "lay", "y", "x"}));
                lw_vars["lay_source"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lay_source", {"gpt_lw", "lay", "y", "x"}));
                lw_vars["lev_source_inc"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lev_source_inc", {"gpt_lw", "lay", "y", "x"}));
                lw_vars["lev_source_dec"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lev_source_dec", {"gpt_lw", "lay", "y", "x"}));
                lw_vars["sfc_source"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sfc_source", {"gpt_lw", "y", "x"}));
            }
            if (switch_output_bnd_fluxes) {
                lw_vars["lw_bnd_flux_up"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_bnd_flux_up", {"band_lw", "lev", "y", "x"}));
                lw_vars["lw_bnd_flux_dn"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_bnd_flux_dn", {"band_lw", "lev", "y", "x"}));
                lw_vars["lw_bnd_flux_net"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("lw_bnd_flux_net", {"band_lw", "lev", "y", "x"}));
            }
        }

        // Shortwave variables
        std::map<std::string, Netcdf_variable<Float>*> sw_vars;
        if (switch_shortwave) {
            auto nc_sw_band_lims_wvn = output_nc.add_variable<Float>("sw_band_lims_wvn", {"band_sw", "pair"});
            nc_sw_band_lims_wvn.insert(rad_sw->get_band_lims_wavenumber().v(), {0, 0}, {n_bnd_sw, 2});

            if (switch_fluxes) {
                sw_vars["sw_flux_up"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_flux_up", {"lev", "y", "x"}));
                sw_vars["sw_flux_dn"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_flux_dn", {"lev", "y", "x"}));
                sw_vars["sw_flux_dn_dir"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_flux_dn_dir", {"lev", "y", "x"}));
                sw_vars["sw_flux_net"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_flux_net", {"lev", "y", "x"}));
            }
            if (switch_output_optical) {
                sw_vars["sw_tau"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_tau", {"gpt_sw", "lay", "y", "x"}));
                sw_vars["ssa"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("ssa", {"gpt_sw", "lay", "y", "x"}));
                sw_vars["g"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("g", {"gpt_sw", "lay", "y", "x"}));
                sw_vars["toa_source"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("toa_source", {"gpt_sw", "y", "x"}));
            }
            if (switch_output_bnd_fluxes) {
                sw_vars["sw_bnd_flux_up"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_bnd_flux_up", {"band_sw", "lev", "y", "x"}));
                sw_vars["sw_bnd_flux_dn"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_bnd_flux_dn", {"band_sw", "lev", "y", "x"}));
                sw_vars["sw_bnd_flux_dn_dir"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_bnd_flux_dn_dir", {"band_sw", "lev", "y", "x"}));
                sw_vars["sw_bnd_flux_net"] = new Netcdf_variable<Float>(output_nc.add_variable<Float>("sw_bnd_flux_net", {"band_sw", "lev", "y", "x"}));
            }
        }

        // Write data with dimension checks
        for (int i = 0; i < n_col_x; i++) {
            nc_lay.insert(p_lay_all[i].v(), {0, 0, i}, {n_lay, n_col_y, 1});
            nc_lev.insert(p_lev_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});

            if (switch_fluxes && switch_longwave && switch_shortwave) {
                nc_heating_rate.insert(heating_rate_all[i].v(), {0, 0, i}, {n_lay, n_col_y, 1});
                nc_lw_heating_rate.insert(heating_rate_lw_all[i].v(), {0, 0, i}, {n_lay, n_col_y, 1});
                nc_sw_heating_rate.insert(heating_rate_sw_all[i].v(), {0, 0, i}, {n_lay, n_col_y, 1});
            }

            if (switch_longwave) {
                if (switch_fluxes) {
                    lw_vars["lw_flux_up"]->insert(lw_flux_up_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                    lw_vars["lw_flux_dn"]->insert(lw_flux_dn_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                    lw_vars["lw_flux_net"]->insert(lw_flux_net_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                }
                if (switch_output_optical) {
                    lw_vars["lw_tau"]->insert(lw_tau_all[i].v(), {0, 0, 0, i}, {n_gpt_lw, n_lay, n_col_y, 1});
                    lw_vars["lay_source"]->insert(lay_source_all[i].v(), {0, 0, 0, i}, {n_gpt_lw, n_lay, n_col_y, 1});
                    lw_vars["lev_source_inc"]->insert(lev_source_inc_all[i].v(), {0, 0, 0, i}, {n_gpt_lw, n_lay, n_col_y, 1});
                    lw_vars["lev_source_dec"]->insert(lev_source_dec_all[i].v(), {0, 0, 0, i}, {n_gpt_lw, n_lay, n_col_y, 1});
                    lw_vars["sfc_source"]->insert(sfc_source_all[i].v(), {0, 0, i}, {n_gpt_lw, n_col_y, 1});
                }
                if (switch_output_bnd_fluxes) {
                    lw_vars["lw_bnd_flux_up"]->insert(lw_bnd_flux_up_all[i].v(), {0, 0, 0, i}, {n_bnd_lw, n_lev, n_col_y, 1});
                    lw_vars["lw_bnd_flux_dn"]->insert(lw_bnd_flux_dn_all[i].v(), {0, 0, 0, i}, {n_bnd_lw, n_lev, n_col_y, 1});
                    lw_vars["lw_bnd_flux_net"]->insert(lw_bnd_flux_net_all[i].v(), {0, 0, 0, i}, {n_bnd_lw, n_lev, n_col_y, 1});
                }
            }

            if (switch_shortwave) {
                if (switch_fluxes) {
                    sw_vars["sw_flux_up"]->insert(sw_flux_up_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                    sw_vars["sw_flux_dn"]->insert(sw_flux_dn_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                    sw_vars["sw_flux_dn_dir"]->insert(sw_flux_dn_dir_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                    sw_vars["sw_flux_net"]->insert(sw_flux_net_all[i].v(), {0, 0, i}, {n_lev, n_col_y, 1});
                }
                if (switch_output_optical) {
                    sw_vars["sw_tau"]->insert(sw_tau_all[i].v(), {0, 0, 0, i}, {n_gpt_sw, n_lay, n_col_y, 1});
                    sw_vars["ssa"]->insert(ssa_all[i].v(), {0, 0, 0, i}, {n_gpt_sw, n_lay, n_col_y, 1});
                    sw_vars["g"]->insert(g_all[i].v(), {0, 0, 0, i}, {n_gpt_sw, n_lay, n_col_y, 1});
                    sw_vars["toa_source"]->insert(toa_source_all[i].v(), {0, 0, i}, {n_gpt_sw, n_col_y, 1});
                }
                if (switch_output_bnd_fluxes) {
                    sw_vars["sw_bnd_flux_up"]->insert(sw_bnd_flux_up_all[i].v(), {0, 0, 0, i}, {n_bnd_sw, n_lev, n_col_y, 1});
                    sw_vars["sw_bnd_flux_dn"]->insert(sw_bnd_flux_dn_all[i].v(), {0, 0, 0, i}, {n_bnd_sw, n_lev, n_col_y, 1});
                    sw_vars["sw_bnd_flux_dn_dir"]->insert(sw_bnd_flux_dn_dir_all[i].v(), {0, 0, 0, i}, {n_bnd_sw, n_lev, n_col_y, 1});
                    sw_vars["sw_bnd_flux_net"]->insert(sw_bnd_flux_net_all[i].v(), {0, 0, 0, i}, {n_bnd_sw, n_lev, n_col_y, 1});
                }
            }
        }
    } catch (const std::exception& e) {
        Status::print_message("NetCDF error: " + std::string(e.what()));
        throw;
    }

    // Cleanup
    if (switch_longwave) delete rad_lw;
    if (switch_shortwave) delete rad_sw;
    delete[] tmp_sfc_alb_dircont;
    delete[] tmp_sfc_alb_difcont;
    delete[] tmp_sfc_alb_dir;
    delete[] tmp_sfc_alb_dif;

    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<double, std::milli>(time_end-time_start).count();
    Status::print_message("Duration longwave solver: " + std::to_string(duration) + " (ms)");
    Status::print_message("###### Finished RTE+RRTMGP solver ######");
    return;
}

#endif
