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

Float vvm::Radiation::calculate_scaling_factor(int year, int month, int day) {
    // This function can be used when you want to calculate the change of solar constant according to the distance between earth and sun.

    // Spencer (1971)
    int doy = day_of_year(year, month, day);

    // Calculate gamma = 2π * (doy - 1) / 365
    Float gamma = 2 * M_PI * (doy - 1) / 365.0;

    // Spencer (1971) formula for scaling factor
    Float scaling_factor = 1.000110 + 0.034221 * cos(gamma) + 0.001280 * sin(gamma) +
                            0.000719 * cos(2 * gamma) + 0.000077 * sin(2 * gamma);

    return scaling_factor;
}

Float vvm::Radiation::calculate_cos_zenith(int year, int month, int day, double hour, double minute, double second,
                                            Float longitude, double latitude) {
    // Calculate day of the year
    int doy = day_of_year(year, month, day);

    // Calculate solar declination
    Float gamma = 2 * M_PI * (doy - 1) / 365.0;
    Float delta = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
                   0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
                   0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma);

    // Calculate equation of time (in minutes)
    Float eot = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -
                           0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma));

    // Convert UTC time to hours
    Float utc_hours = hour + minute / 60.0 + second / 3600.0;

    // Calculate local solar time (LST) in hours
    Float lst = utc_hours + (longitude / 15.0) + (eot / 60.0);

    // Calculate hour angle in degrees
    Float h = 15.0 * (lst - 12.0);

    // Convert angles to radians
    Float phi_rad = latitude * M_PI / 180.0;
    Float h_rad = h * M_PI / 180.0;

    // Calculate cosine of solar zenith angle
    Float cos_zenith = sin(phi_rad) * sin(delta) + cos(phi_rad) * cos(delta) * cos(h_rad);

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
    const int n_col = n_col_x*n_col_y;
    const int n_lay = model.nz-2;
    const int n_lev = n_lay+1;

    // const Float qvmin = 1e-6;
    // const Float qcmin = 1e-7;
    // const Float qimin = 1e-8;

    const Float effcmin = 2.5;  // unit: 1e-6m
    const Float effcmax = 21.5; // unit: 1e-6m
    const Float effimin = 10;   // unit: 1e-6m
    const Float effimax = 180;  // unit: 1e-6m
    // const Float rel_val = 0.5*(effcmax+effcmin); // unit: 1e-6m
    // const Float rei_val = 0.5*(effimax+effimin); // unit: 1e-6m

    const Float sec_per_day = 86400.0; // Seconds per day
    // !     Calculated value (from constants above and input cpdair)
    // !     (grav) x (#sec/day) / (specific heat of dry air at const. p x 1.e2)
    // !     Here, cpdair is in units of J kg-1 K-1, and the constant (1.e2) 
    // !     converts mb to Pa when heatfac is multiplied by W m-2 mb-1.
    const Float heat_factor = model.GRAVITY * sec_per_day / model.Cp; 
    

    Array<Float, 2> p_lay({n_col_x, n_lay});
    Array<Float, 2> t_lay({n_col_x, n_lay});
    Array<Float, 2> p_lev({n_col_x, n_lev});
    Array<Float, 2> t_lev({n_col_x, n_lev});
    Array<Float, 2> h2o_lay({n_col_x, n_lev});

    // Not needed in this case
    Array<Float,2> col_dry;
    if (switch_col_dry) {
        col_dry.set_dims({n_col_x, n_lay});
        // col_dry = std::move(input_nc.get_variable<Float>("col_dry", {n_lay, n_col_y, n_col_x}));
    }

    Gas_concs gas_concs;

    Float o3_min = 1e-13;
    Float g1 = 3.6478, g2 = 0.83209, g3 = 11.3515;
    Array<Float, 2> o3_lay({n_col_x, n_lay});
    for (int i = 1; i <= n_col_x; i++) {
        for (int k = 1; k <= n_lay; k++) {
            p_lay({i,k}) = model.pb[k];
            t_lay({i,k}) = model.th[i][k] * model.pib[k];
            o3_lay({i,k}) = std::max(o3_min, g1 * std::pow(model.pb[k]/100., g2) * std::exp(-model.pb[k]/100. / g3) * 1e-6);
            h2o_lay({i,k}) = model.qv[i][k];
        }
        for (int k = 1; k <= n_lev-1; k++) {
            p_lev({i,k}) = model.pb_lev[k];
            t_lev({i,k}) = 0.5 * (model.th[i][k] + model.th[i][k + 1]) * model.pib_lev[k];
        }
        t_lev({i,n_lev-1}) = model.th[i][n_lev-1] * model.pib_lev[n_lev];
        t_lev({i,n_lev}) = t_lev({i,n_lev-1});
        p_lev({i,n_lev}) = model.pb_lev[n_lev];
    }

    gas_concs.set_vmr("h2o", h2o_lay); // Minimum water vapor as placeholder
    gas_concs.set_vmr("co2", 348e-6);        // ppm
    gas_concs.set_vmr("o3", o3_lay);   // Calculated ozone profile
    gas_concs.set_vmr("n2o", 306e-9);        // ppm
    gas_concs.set_vmr("ch4", 1650e-9);       // ppm
    gas_concs.set_vmr("o2", 0.2095);         // Volume mixing ratio
    gas_concs.set_vmr("n2", 0.7808);         // Volume mixing ratio
    
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

    Array<Float,2> lwp({n_col_x, n_lay});
    Array<Float,2> iwp({n_col_x, n_lay});
    Array<Float,2> rel({n_col_x, n_lay});
    Array<Float,2> rei({n_col_x, n_lay});
    for (int i = 1; i <= n_col_x; i++) {
        for (int k = 1; k <= n_lay; k++) {
            lwp({i,k}) = std::max(model.qc[i][k] * 1e3, 0.) * model.rhou[k] * model.dz; // (g/kg) * (kg/m^3) * (m)
            iwp({i,k}) = std::max(model.qitot[i][k] * 1e3, 0.) * model.rhou[k] * model.dz; // (g/kg) * (kg/m^3) * (m)
            rel({i,k}) = std::min(std::max(model.diag_effc[i][k], effcmin), effcmax); // (m)
            rei({i,k}) = std::min(std::max(model.diag_effi[i][k], effimin), effimax); // (m)
        }
    }
    printf("a1\n");


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


    ////// CREATE THE OUTPUT FILE //////
    // Create the general dimensions and arrays.
    Status::print_message("Preparing NetCDF output file.");

     std::string output_nc_name = model.outputpath + "nc/Radiation_" + std::to_string(model.step) + ".nc";
     Netcdf_file output_nc(output_nc_name, Netcdf_mode::Create);
     output_nc.add_dimension("x", n_col_x);
     output_nc.add_dimension("y", n_col_y);
     output_nc.add_dimension("lay", n_lay);
     output_nc.add_dimension("lev", n_lev);
     output_nc.add_dimension("pair", 2);

     auto nc_lay = output_nc.add_variable<Float>("p_lay", {"lay", "y", "x"});
     auto nc_lev = output_nc.add_variable<Float>("p_lev", {"lev", "y", "x"});

     nc_lay.insert(p_lay.v(), {0, 0, 0});
     nc_lev.insert(p_lev.v(), {0, 0, 0});


    // Nan check
    for (int i = 1; i <= n_col_x; i++) {
        for (int k = 0; k < model.nz; k++) {
            if (std::isnan(model.th[i][k]) || std::isnan(model.qv[i][k]) ||
                std::isnan(model.qc[i][k]) || std::isnan(model.qitot[i][k]) ||
                std::isnan(model.diag_effc[i][k]) || std::isnan(model.diag_effi[i][k]) ||
                std::isnan(model.rhou[k])) {
                std::cout << "NaN in model input at i=" << i << ", k=" << k << std::endl;
                exit(1);
            }
            if (model.rhou[k] <= 0) {
                std::cout << "Non-positive density at k=" << k << ": rhou=" << model.rhou[k] << std::endl;
                exit(1);
            }
        }
    }

    ////// RUN THE LONGWAVE SOLVER //////
    Array<Float,2> net_heat_rate({n_col_x, n_lay});
    if (switch_longwave)
    {
        // Initialize the solver.
        Status::print_message("Initializing the longwave solver.");
        Radiation_solver_longwave rad_lw(
            gas_concs,
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-lw-g128.nc",
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-lw.nc"
        );

        // Read the boundary conditions.
        const int n_bnd_lw = rad_lw.get_n_bnd();
        const int n_gpt_lw = rad_lw.get_n_gpt();

        Array<Float,2> emis_sfc({n_bnd_lw, n_col_x}); // Start from k = 1 to k = model.nz-2
        Array<Float,1> t_sfc({n_col_x});

        for (int i = 1; i <= n_col_x; i++) {
            t_sfc({i}) = t_lev({i,1});
            for (int b = 1; b <= n_bnd_lw; b++) {
                emis_sfc({b,i}) = 0.98;
            }
        }
        

        // Create output arrays.
        Array<Float,3> lw_tau;
        Array<Float,3> lay_source;
        Array<Float,3> lev_source_inc;
        Array<Float,3> lev_source_dec;
        Array<Float,2> sfc_source;

        if (switch_output_optical)
        {
            lw_tau        .set_dims({n_col_x, n_lay, n_gpt_lw});
            lay_source    .set_dims({n_col_x, n_lay, n_gpt_lw});
            lev_source_inc.set_dims({n_col_x, n_lay, n_gpt_lw});
            lev_source_dec.set_dims({n_col_x, n_lay, n_gpt_lw});
            sfc_source    .set_dims({n_col_x, n_gpt_lw});
        }

        Array<Float,2> lw_flux_up;
        Array<Float,2> lw_flux_dn;
        Array<Float,2> lw_flux_net;

        if (switch_fluxes)
        {
            lw_flux_up .set_dims({n_col_x, n_lev});
            lw_flux_dn .set_dims({n_col_x, n_lev});
            lw_flux_net.set_dims({n_col_x, n_lev});
        }

        Array<Float,3> lw_bnd_flux_up;
        Array<Float,3> lw_bnd_flux_dn;
        Array<Float,3> lw_bnd_flux_net;

        if (switch_output_bnd_fluxes)
        {
            lw_bnd_flux_up .set_dims({n_col_x, n_lev, n_bnd_lw});
            lw_bnd_flux_dn .set_dims({n_col_x, n_lev, n_bnd_lw});
            lw_bnd_flux_net.set_dims({n_col_x, n_lev, n_bnd_lw});
        }


        // Solve the radiation.
        // Status::print_message("Solving the longwave radiation.");
        printf("Solving the longwave radiation.\n");

        auto time_start = std::chrono::high_resolution_clock::now();

        rad_lw.solve(
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


        // Store the output.
        Status::print_message("Storing the longwave output.");

        output_nc.add_dimension("gpt_lw", n_gpt_lw);
        output_nc.add_dimension("band_lw", n_bnd_lw);

        auto nc_lw_band_lims_wvn = output_nc.add_variable<Float>("lw_band_lims_wvn", {"band_lw", "pair"});
        nc_lw_band_lims_wvn.insert(rad_lw.get_band_lims_wavenumber().v(), {0, 0});

        if (switch_output_optical)
        {
            auto nc_lw_band_lims_gpt = output_nc.add_variable<int>("lw_band_lims_gpt", {"band_lw", "pair"});
            nc_lw_band_lims_gpt.insert(rad_lw.get_band_lims_gpoint().v(), {0, 0});

            auto nc_lw_tau = output_nc.add_variable<Float>("lw_tau", {"gpt_lw", "lay", "y", "x"});
            nc_lw_tau.insert(lw_tau.v(), {0, 0, 0, 0});

            auto nc_lay_source     = output_nc.add_variable<Float>("lay_source"    , {"gpt_lw", "lay", "y", "x"});
            auto nc_lev_source_inc = output_nc.add_variable<Float>("lev_source_inc", {"gpt_lw", "lay", "y", "x"});
            auto nc_lev_source_dec = output_nc.add_variable<Float>("lev_source_dec", {"gpt_lw", "lay", "y", "x"});

            auto nc_sfc_source = output_nc.add_variable<Float>("sfc_source", {"gpt_lw", "y", "x"});

            nc_lay_source.insert    (lay_source.v()    , {0, 0, 0, 0});
            nc_lev_source_inc.insert(lev_source_inc.v(), {0, 0, 0, 0});
            nc_lev_source_dec.insert(lev_source_dec.v(), {0, 0, 0, 0});

            nc_sfc_source.insert(sfc_source.v(), {0, 0, 0});
        }

        if (switch_fluxes)
        {
            auto nc_lw_flux_up  = output_nc.add_variable<Float>("lw_flux_up" , {"lev", "y", "x"});
            auto nc_lw_flux_dn  = output_nc.add_variable<Float>("lw_flux_dn" , {"lev", "y", "x"});
            auto nc_lw_flux_net = output_nc.add_variable<Float>("lw_flux_net", {"lev", "y", "x"});

            nc_lw_flux_up .insert(lw_flux_up .v(), {0, 0, 0});
            nc_lw_flux_dn .insert(lw_flux_dn .v(), {0, 0, 0});
            nc_lw_flux_net.insert(lw_flux_net.v(), {0, 0, 0});

            if (switch_output_bnd_fluxes)
            {
                auto nc_lw_bnd_flux_up  = output_nc.add_variable<Float>("lw_bnd_flux_up" , {"band_lw", "lev", "y", "x"});
                auto nc_lw_bnd_flux_dn  = output_nc.add_variable<Float>("lw_bnd_flux_dn" , {"band_lw", "lev", "y", "x"});
                auto nc_lw_bnd_flux_net = output_nc.add_variable<Float>("lw_bnd_flux_net", {"band_lw", "lev", "y", "x"});

                nc_lw_bnd_flux_up .insert(lw_bnd_flux_up .v(), {0, 0, 0, 0});
                nc_lw_bnd_flux_dn .insert(lw_bnd_flux_dn .v(), {0, 0, 0, 0});
                nc_lw_bnd_flux_net.insert(lw_bnd_flux_net.v(), {0, 0, 0, 0});
            }
            auto nc_lw_heat     = output_nc.add_variable<Float>("lw_heat_rate"    , {"lay", "y", "x"});
            Array<Float,2> lw_heat_rate({n_col_x, n_lay});
            for (int i = 1; i <= n_col_x; i++) {
                for (int k = 1; k <= n_lay; k++) {
                    lw_heat_rate({i,k}) = heat_factor * (lw_flux_net({i,k+1})-lw_flux_net({i,k})) / (p_lev({i,k})-p_lev({i,k+1}) );
                    net_heat_rate({i,k}) = lw_heat_rate({i,k});
                }
            }
        }
    }
    printf("test3\n");


    ////// RUN THE SHORTWAVE SOLVER //////
    if (switch_shortwave)
    {
        // Initialize the solver.
        Status::print_message("Initializing the shortwave solver.");

        Radiation_solver_shortwave rad_sw(
            gas_concs, switch_cloud_optics, switch_aerosol_optics,
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-gas-sw-g112.nc",
            "../external/rte-rrtmgp-cpp/rrtmgp-data/rrtmgp-clouds-sw.nc",
            "../external/rte-rrtmgp-cpp/allsky/aerosol_optics.nc"
        );

        // Read the boundary conditions.
        const int n_bnd_sw = rad_sw.get_n_bnd();
        const int n_gpt_sw = rad_sw.get_n_gpt();
        const Float tsi_ref = rad_sw.get_tsi();
        const Float tsi = 340.; // total_solar_irradiance 

        Array<Float,1> mu0({n_col_x});
        Array<Float,2> sfc_alb_dir({n_bnd_sw, n_col_x});
        Array<Float,2> sfc_alb_dif({n_bnd_sw, n_col_x});
        Array<Float,1> tsi_scaling({n_col_x});
        for (int i = 1; i <= n_col_x; i++) {
            mu0({i}) = calculate_cos_zenith(model.year, model.month, model.day, model.hour, model.minute, model.second, model.lon, model.lat);
            tsi_scaling({i}) = tsi / tsi_ref;
            for (int b = 1; b <= n_bnd_sw; b++) {
                sfc_alb_dir({b,i}) = 0.06;
                sfc_alb_dif({b,i}) = 0.06;
            }
        }

        // Create output arrays.
        Array<Float,3> sw_tau;
        Array<Float,3> ssa;
        Array<Float,3> g;
        Array<Float,2> toa_source;

        if (switch_output_optical)
        {
            sw_tau    .set_dims({n_col_x, n_lay, n_gpt_sw});
            ssa       .set_dims({n_col_x, n_lay, n_gpt_sw});
            g         .set_dims({n_col_x, n_lay, n_gpt_sw});
            toa_source.set_dims({n_col_x, n_gpt_sw});
        }
        Array<Float,2> sw_flux_up;
        Array<Float,2> sw_flux_dn;
        Array<Float,2> sw_flux_dn_dir;
        Array<Float,2> sw_flux_net;

        if (switch_fluxes)
        {
            sw_flux_up    .set_dims({n_col_x, n_lev});
            sw_flux_dn    .set_dims({n_col_x, n_lev});
            sw_flux_dn_dir.set_dims({n_col_x, n_lev});
            sw_flux_net   .set_dims({n_col_x, n_lev});
        }

        Array<Float,3> sw_bnd_flux_up;
        Array<Float,3> sw_bnd_flux_dn;
        Array<Float,3> sw_bnd_flux_dn_dir;
        Array<Float,3> sw_bnd_flux_net;

        if (switch_output_bnd_fluxes)
        {
            sw_bnd_flux_up    .set_dims({n_col_x, n_lev, n_bnd_sw});
            sw_bnd_flux_dn    .set_dims({n_col_x, n_lev, n_bnd_sw});
            sw_bnd_flux_dn_dir.set_dims({n_col_x, n_lev, n_bnd_sw});
            sw_bnd_flux_net   .set_dims({n_col_x, n_lev, n_bnd_sw});
        }

        // Solve the radiation.
        Status::print_message("Solving the shortwave radiation.");

        auto time_start = std::chrono::high_resolution_clock::now();

        rad_sw.solve(
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


        // Store the output.
        Status::print_message("Storing the shortwave output.");

        output_nc.add_dimension("gpt_sw", n_gpt_sw);
        output_nc.add_dimension("band_sw", n_bnd_sw);

        auto nc_sw_band_lims_wvn = output_nc.add_variable<Float>("sw_band_lims_wvn", {"band_sw", "pair"});
        nc_sw_band_lims_wvn.insert(rad_sw.get_band_lims_wavenumber().v(), {0, 0});

        if (switch_output_optical)
        {
            auto nc_sw_band_lims_gpt = output_nc.add_variable<int>("sw_band_lims_gpt", {"band_sw", "pair"});
            nc_sw_band_lims_gpt.insert(rad_sw.get_band_lims_gpoint().v(), {0, 0});

            auto nc_sw_tau = output_nc.add_variable<Float>("sw_tau", {"gpt_sw", "lay", "y", "x"});
            auto nc_ssa    = output_nc.add_variable<Float>("ssa"   , {"gpt_sw", "lay", "y", "x"});
            auto nc_g      = output_nc.add_variable<Float>("g"     , {"gpt_sw", "lay", "y", "x"});

            nc_sw_tau.insert(sw_tau.v(), {0, 0, 0, 0});
            nc_ssa   .insert(ssa   .v(), {0, 0, 0, 0});
            nc_g     .insert(g     .v(), {0, 0, 0, 0});

            auto nc_toa_source = output_nc.add_variable<Float>("toa_source", {"gpt_sw", "y", "x"});
            nc_toa_source.insert(toa_source.v(), {0, 0, 0});
        }

        if (switch_fluxes)
        {
            auto nc_sw_flux_up     = output_nc.add_variable<Float>("sw_flux_up"    , {"lev", "y", "x"});
            auto nc_sw_flux_dn     = output_nc.add_variable<Float>("sw_flux_dn"    , {"lev", "y", "x"});
            auto nc_sw_flux_dn_dir = output_nc.add_variable<Float>("sw_flux_dn_dir", {"lev", "y", "x"});
            auto nc_sw_flux_net    = output_nc.add_variable<Float>("sw_flux_net"   , {"lev", "y", "x"});

            nc_sw_flux_up    .insert(sw_flux_up    .v(), {0, 0, 0});
            nc_sw_flux_dn    .insert(sw_flux_dn    .v(), {0, 0, 0});
            nc_sw_flux_dn_dir.insert(sw_flux_dn_dir.v(), {0, 0, 0});
            nc_sw_flux_net   .insert(sw_flux_net   .v(), {0, 0, 0});

            if (switch_output_bnd_fluxes)
            {
                auto nc_sw_bnd_flux_up     = output_nc.add_variable<Float>("sw_bnd_flux_up"    , {"band_sw", "lev", "y", "x"});
                auto nc_sw_bnd_flux_dn     = output_nc.add_variable<Float>("sw_bnd_flux_dn"    , {"band_sw", "lev", "y", "x"});
                auto nc_sw_bnd_flux_dn_dir = output_nc.add_variable<Float>("sw_bnd_flux_dn_dir", {"band_sw", "lev", "y", "x"});
                auto nc_sw_bnd_flux_net    = output_nc.add_variable<Float>("sw_bnd_flux_net"   , {"band_sw", "lev", "y", "x"});

                nc_sw_bnd_flux_up    .insert(sw_bnd_flux_up    .v(), {0, 0, 0, 0});
                nc_sw_bnd_flux_dn    .insert(sw_bnd_flux_dn    .v(), {0, 0, 0, 0});
                nc_sw_bnd_flux_dn_dir.insert(sw_bnd_flux_dn_dir.v(), {0, 0, 0, 0});
                nc_sw_bnd_flux_net   .insert(sw_bnd_flux_net   .v(), {0, 0, 0, 0});
            }


            auto nc_sw_heat     = output_nc.add_variable<Float>("sw_heat_rate"    , {"lay", "y", "x"});
            Array<Float,2> sw_heat_rate({n_col_x, n_lay});
            for (int i = 1; i <= n_col_x; i++) {
                for (int k = 1; k <= n_lay; k++) {
                    sw_heat_rate({i,k}) = heat_factor * (sw_flux_net({i,k+1})-sw_flux_net({i,k})) / (p_lev({i,k})-p_lev({i,k+1}) );
                    net_heat_rate({i,k}) += sw_heat_rate({i,k});
                }
            }
        }
    }

    for (int i = 1; i <= n_col_x; ++i) {
        for (int k = 1; k < n_lay; ++k) {
            model.radiation_heating_rate[i][k] = net_heat_rate({i,k});
        }
        model.radiation_heating_rate[i][n_lay] = net_heat_rate({i,n_lay}) = 0.;
    }
    model.BoundaryProcess2D_center(model.radiation_heating_rate, model.nx, model.nz);

    auto nc_net_heat     = output_nc.add_variable<Float>("net_heat_rate"    , {"lay", "y", "x"});
    nc_net_heat.insert(net_heat_rate.v(), {0, 0, 0});

    auto time_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<Float, std::milli>(time_end-time_start).count();
    Status::print_message("Duration solver: " + std::to_string(duration) + " (ms)");
    Status::print_message("###### Finished RTE+RRTMGP solver ######");
    return;
}

#endif
