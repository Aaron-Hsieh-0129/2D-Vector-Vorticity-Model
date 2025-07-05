#include "Declare.hpp"
#include <iostream>
#include <random>
#if defined(LOADFILE)
    #include <fstream>
#endif
#if defined(PETSC)
    #include <petsc.h>
#endif
#if defined(LOADFROMPREVIOUSFILE)
    #include <netcdf>
    using namespace netCDF;
#endif
#include <iomanip>

void generateAddlfxArray(double *array, int size, double variation = 0.1) {
    std::mt19937 rng(static_cast<unsigned>(time(nullptr))); // Random number generator
    std::uniform_real_distribution<double> dist(-variation, variation); // Variation range

    for (int i = 0; i < size; ++i) {
        array[i] = 1.0 + dist(rng); // Mean of 1.0 with random variation
    }
    return;
}

void vvm::Init::Init1d(vvm &model) {
    for (int i = 0; i < model.nx; i++) model.x[i] = (i-0.5) * model.dx;

    // Height stretch setting 
    double DOMAIN = 15000.;
    model.CZ2 = (model.dz-model.dz1) / (model.dz * (DOMAIN-model.dz));
    model.CZ1 = 1. - model.CZ2 * DOMAIN;
    
    for (int k = 0; k < model.nz; k++) {
        model.z[k] = (k-0.5) * model.dz;
        model.z_zeta[k] = (k-1) * model.dz;
    }

    for (int k = 0; k < model.nz; k++) {
        // Make the coefficient array of flexible height before the height is modified to new height
        model.flex_height_coef_th[k] = 1. / (model.CZ1 + 2 * model.CZ2 * model.z[k]);
        model.flex_height_coef_zeta[k] = 1. / (model.CZ1 + 2 * model.CZ2 * model.z_zeta[k]);

        model.z[k] = model.z[k] * (model.CZ1 + model.CZ2 * model.z[k]);
        model.z_zeta[k] = model.z_zeta[k] * (model.CZ1 + model.CZ2 * model.z_zeta[k]);
    }

    for (int k = 1; k < model.nz-1; k++) {
        model.dz_th[k] = model.z_zeta[k+1] - model.z_zeta[k];
        model.dz_zeta[k] = model.z[k] - model.z[k-1];
    }
    model.BoundaryProcess1D_center(model.dz_th,model.nz);
    model.BoundaryProcess1D_center(model.dz_zeta,model.nz);

    // Initialization for p3 microphysics
    #if defined(P3_MICROPHY) && defined(WATER)
        for (int k = 0; k < model.nz; k++) {
            for (int i = 0; i < model.nx; i++) { 
                model.dz_all[i][k] = model.dz_th[k];
                model.w_all[i][k] = 0.;
                model.pb_all[i][k] = model.pb[k];
                model.zi_all[i][k] = 0.;
                model.ssat_all[i][k] = 0.;
            }
        }
    #endif


    #if defined(LOADFILE)
        LoadFile(model);
    #else
        // init tb
        model.thb[1] = 300.;
        for (int k = 2; k <= model.nz-2; k++) {
            #ifdef DRY
                model.thb[k] = 300.;
            #else
                model.thb[k] = GetTB(k, model);
                model.thb_init[k] = model.thb[k];
            #endif
        }
        model.BoundaryProcess1D_center(model.thb, model.nz);
        model.BoundaryProcess1D_center(model.thb_init, model.nz);

        // init qvb, tvb
        for (int k = 1; k <= model.nz-2; k++) {
            #if defined(WATER)
                model.qvb[k] = GetQVB(k, model);
            #else
                model.qvb[k] = 0.;
            #endif
            model.thvb[k] = model.thb[k] * (1. + 0.608 * model.qvb[k]);
        }
        model.BoundaryProcess1D_center(model.qvb, model.nz);
        model.BoundaryProcess1D_center(model.qvb0, model.nz);
        model.BoundaryProcess1D_center(model.thvb, model.nz);

        // init pib
        double pisfc = pow((model.PSURF / model.P0), model.Rd / model.Cp);
        for (int k = 1; k <= model.nz-2; k++) {
            if (k == 1) model.pib[k] = pisfc - model.GRAVITY * 0.5 * model.dz_th[k] / (model.Cp * model.thvb[k]);
            else {
                double tvbavg = 0.5*(model.thvb[k] + model.thvb[k-1]);
                model.pib[k] = model.pib[k-1] - model.GRAVITY * model.dz_th[k] / (model.Cp * tvbavg);
            }
        }
        model.BoundaryProcess1D_center(model.pib, model.nz);
        
        for (int k = 1; k <= model.nz-2; k++) {
            model.pib_lev[k] = 0.5*(model.pib[k]+model.pib[k-1]);
        }
        // extrapolation
        model.pib_lev[1] = model.pib_lev[2] + (model.pib_lev[2]-model.pib_lev[3]);
        model.pib_lev[model.nz-1] = model.pib_lev[model.nz-2] - (model.pib_lev[model.nz-3]-model.pib_lev[model.nz-2]);
        model.BoundaryProcess1D_center(model.pib_lev, model.nz+1);

        // init tb_zeta, rhou
        for (int k = 1; k <= model.nz-2; k++) {
            #ifdef RHO1
                model.rhou[k] = 1.;
            #else
                model.rhou[k] = model.P0 * pow(model.pib[k], model.Cv/model.Rd) / (model.Rd * model.thvb[k]);
            #endif
        }
        model.BoundaryProcess1D_center(model.rhou, model.nz);

        // init tb_zeta, rhow
        for (int k = 2; k <= model.nz-1; k++) {
            model.thb_zeta[k] = 0.5 * (model.thb[k] + model.thb[k-1]);
            model.rhow[k] = 0.5 * (model.rhou[k] + model.rhou[k-1]);
        }    
        model.thb_zeta[1] = model.thb_zeta[2] - (model.thb_zeta[3] - model.thb_zeta[2]);
        model.rhow[1] = model.rhow[2] - (model.rhow[3] - model.rhow[2]);
        model.thb_zeta[0] = model.thb_zeta[1];
        model.rhow[0] = model.rhow[1];
        model.rhou[0] = model.rhow[0];

        // init pb, qvsb
        for (int k = 1; k <= model.nz-2; k++) {
            model.pb[k] = model.P0 * pow(model.pib[k], model.Cp / model.Rd);
            double Tc = model.thb[k] * model.pib[k];
            double es = 611.2 * std::exp(17.67 * (Tc-273.15) / (Tc-273.15+243.5));
            model.qvsb[k] = 0.622 * es/ (model.pb[k] - 0.378 * es);
        }
        model.BoundaryProcess1D_center(model.pb, model.nz);
        model.BoundaryProcess1D_center(model.qvsb, model.nz);

        // init pb_lev (nz: model.nz+1)
        for (int k = 1; k <= model.nz-2; k++) {
            model.pb_lev[k] = 0.5*(model.pb[k]+model.pb[k-1]);
        }
        model.pb_lev[1] = model.PSURF;
        model.pb_lev[model.nz-1] = model.pb_lev[model.nz-2] - (model.pb_lev[model.nz-3]-model.pb_lev[model.nz-2]); // extrapolation
        model.BoundaryProcess1D_center(model.pb_lev, model.nz+1);

        #if defined(WATER)
            for (int k = 1; k <= model.nz-2; k++) {
                model.qvb[k] = GetQVB(k, model);
            }
            model.BoundaryProcess1D_center(model.qvb, model.nz);
        #endif

        #if defined(RHO1)
            for (int k = 0; k < model.nz; k++) {
                model.rhou[k] = model.rhow[k] = 1.;
            }
        #endif
    #endif

    for (int k = 0; k < model.nz; k++) {
        model.thbm[k] = model.thb[k];
        model.qvb0[k] = model.qvb[k];
        #if defined(WATER)
            model.thvb[k] = model.thvbm[k] = model.thb[k] * (1 + 0.608 * model.qvb[k]);
        #else
            model.thvb[k] = model.thvbm[k] = model.thb[k];
        #endif
    }

    for (int k = 1; k < model.nz-1; k++) {
        model.lambda2[k] = 1. / (1. / pow(0.23 * std::sqrt(model.dx*model.dz), 2) + 1. / pow(0.4* 0.4 * model.z[k], 2));
        model.lambda2_zeta[k] = 1. / (1. / pow(0.23 * std::sqrt(model.dx*model.dz), 2) + 1. / pow(0.4* 0.4 * model.z_zeta[k], 2));
    }
    model.BoundaryProcess1D_center(model.lambda2, model.nz);
    model.BoundaryProcess1D_center(model.lambda2_zeta, model.nz);

    double tau_min = 60., tau_max = 1800.;
    int k_diff_start = 0.;
    for (int k = 1; k < model.nz-1; k++) {
        if (model.z[k] >= 15000) {
            k_diff_start = k;
            break;
        }
    }
    for (int k = 0; k <= model.nz-1; k++) {
        if (k >= k_diff_start) {
            model.nudge_tau[k] = tau_min * std::pow(tau_max/tau_min, ((model.z[model.nz-1]-model.z[k])/(model.z[model.nz-1]-model.z[model.nz-1-k_diff_start])));
            model.nudge_tau[k] = 1. / model.nudge_tau[k];
        }
        else model.nudge_tau[k] = 0.;
    }
    model.BoundaryProcess1D_center(model.nudge_tau, model.nz);

    for (int i = 0; i <= model.nx-1; i++) {
        model.th_ground[i] = 303.;
        double Tc = model.th_ground[i] * model.pib[1];
        double es = 611.2 * std::exp(17.67 * (Tc-273.15) / (Tc-273.15+243.5));
        model.qvs_ground[i] = 0.622 * es/ (model.P0 - 0.378 * es);
    }
    generateAddlfxArray(model.addflux, model.nx);

    // Initialization for p3 microphysics
    #if defined(P3_MICROPHY) && defined(WATER)
        for (int k = 0; k < model.nz; k++) {
            for (int i = 0; i < model.nx; i++) { 
                model.pb_all[i][k] = model.pb[k];
            }
        }
    #endif

    return;
}

void vvm::Init::Init2d(vvm &model) {
    #if defined(WATER)
        // init qv: where th != 0, qv = qvs
        for (int k = 0; k <= model.nz-1; k++) {
            for (int i = 0; i <= model.nx-1; i++) {
                model.qv[i][k] = model.qvm[i][k] = model.qvb[k];
                model.qvp[i][k] = 0.;
                model.qc[i][k] = model.qcp[i][k] = model.qcm[i][k] = 0.;
                model.qr[i][k] = model.qrp[i][k] = model.qrm[i][k] = 0.;
                model.dqv_advect[i][k][0] = model.dqv_advect[i][k][1] = 0;
                model.dqc_advect[i][k][0] = model.dqc_advect[i][k][1] = 0;
                model.dqr_advect[i][k][0] = model.dqr_advect[i][k][1] = 0;
                #if defined(P3_MICROPHY)
                    model.nc[i][k] = model.ncp[i][k] = model.ncm[i][k] = 0.;
                    model.nr[i][k] = model.nrp[i][k] = model.nrm[i][k] = 0.;
                    model.ni[i][k] = model.nip[i][k] = model.nim[i][k] = 0.;
                    model.qitot[i][k] = model.qitotp[i][k] = model.qitotm[i][k] = 0.;
                    model.qirim[i][k] = model.qirimp[i][k] = model.qirimm[i][k] = 0.;
                    model.qiliq[i][k] = model.qiliqp[i][k] = model.qiliqm[i][k] = 0.;
                    model.birim[i][k] = model.birimp[i][k] = model.birimm[i][k] = 0.;
                    model.dnc_advect[i][k][0] = model.dnc_advect[i][k][1] = 0;
                    model.dnr_advect[i][k][0] = model.dnr_advect[i][k][1] = 0;
                    model.dni_advect[i][k][0] = model.dni_advect[i][k][1] = 0;
                    model.dqitot_advect[i][k][0] = model.dqitot_advect[i][k][1] = 0;
                    model.dqirim_advect[i][k][0] = model.dqirim_advect[i][k][1] = 0;
                    model.dqiliq_advect[i][k][0] = model.dqiliq_advect[i][k][1] = 0;
                    model.dbirim_advect[i][k][0] = model.dbirim_advect[i][k][1] = 0;
                #endif
            }
            model.qvb0[k] = model.qvb[k];
        }
    #endif
	#if defined(TROPICALFORCING)
		// Generate random 2D Gaussian noise array within the specified range
		RandomPerturbation(model, 0);

        for (int i = 1; i <= model.nx-2; i++) {
            for (int k = 1; k <= model.nz-2; k++) {
				model.th[i][k] = model.thb[k] + model.init_th_forcing[i][k];
                model.thm[i][k] = model.th[i][k];

                model.qv[i][k] = model.qvb[k];
                model.qvm[i][k] = model.qv[i][k];
            }
        }
        model.BoundaryProcess2D_center(model.th, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.thm, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.qv, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.qvm, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.u, model.nx, model.nz);
    #else
        // init th
        for (int i = 1; i <= model.nx-2; i++) {
            for (int k = 1; k <= model.nz-2; k++) {
                if (model.CASE == 0) model.th[i][k] = model.thb[k];
                else if (model.CASE == 1) model.th[i][k] = model.thb[k] + GetTH(i, k, model);

                if (model.addforcingtime > 0) {
                    RandomPerturbation(model, 0);
                    model.th[i][k] += model.init_th_forcing[i][k];
                }
                
                model.thm[i][k] = model.th[i][k];

                model.u[i][k] = model.ubar[k];
                model.w[i][k] = 0.;
            }
        }
        model.BoundaryProcess2D_center(model.th, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.thm, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.u, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.w, model.nx, model.nz);

        for (int k = 0; k <= model.nz-1; k++) {
            for (int i = 0; i <= model.nx-1; i++) {
                #if defined(LOADFILE)
                    model.qv[i][k] = model.qvm[i][k] = model.qvb[k];
                    model.qvb0[k] = model.qvb[k];
                #else
                    model.qv[i][k] = model.qvm[i][k] = model.qvb[k] * 0.9;
                    model.qvb0[k] = model.qvb[k] * 0.9;
                #endif
            }
        }
	#endif

	// init zeta
	double pu_pz = 0., pw_px = 0.;
	for (int i = 1; i <= model.nx-2; i++) {
		for (int k = 1; k <= model.nz-2; k++) {
			pw_px = (model.w[i][k] - model.w[i-1][k]) * model.rdx;
			pu_pz = (model.u[i][k] - model.u[i][k-1]) * model.rdz;
			model.zeta[i][k] = (pw_px - pu_pz) / model.rhow[k];
			model.zetam[i][k] = model.zeta[i][k];
		}
	}
	model.BoundaryProcess2D_westdown(model.zeta, model.nx, model.nz);
	model.BoundaryProcess2D_westdown(model.zetam, model.nx, model.nz);

	// init ubar at top
	for (int i = 1; i < model.nx-1; i++) {
		model.ubarTopm += model.u[i][model.nz-2];
	}
	model.ubarTopm /= ((double) (model.nx - 2.));
	model.ubarTop = model.ubarTopm;

    for (int k = 0; k < model.nz; k++) {
        for (int i = 0; i < model.nx; i++) {
            model.zetap[i][k] = 0.;
            model.thp[i][k] = 0.;
        }
    }


    #if defined(P3_MICROPHY)
    for (int i = 0; i < model.nx; i++) {
        for (int k = 0; k <= model.nz; k++) {
            model.pb_lev_all[i][k] = model.pb_lev[k];
        }
    }
    #endif

	return;
}

double vvm::Init::GetTB(int k, vvm &model) {
    double z_top = 12000., T_top = 213., tb_top = 343.;
    if (model.z[k] <= z_top) return 300. + 43. * pow(model.z[k] / z_top, 1.25);
    else return tb_top * exp(model.GRAVITY * (model.z[k] - z_top) / (model.Cp * T_top));
}

double vvm::Init::GetTHRAD(int i, int k, vvm &model) {
    double XC = model.XRANGE / 2., XR = 4000.;
    double ZC = 2500., ZR = 2000.;
    double x = (i-0.5) * model.dx;
    double rad = sqrt(pow((x - XC) / XR, 2) + pow((model.z[k]- ZC) / ZR, 2));
    return rad;
}

double vvm::Init::GetTH(int i, int k, vvm &model) {
    double rad = GetTHRAD(i, k, model);
    double delta = 3.;
    if (rad <= 1) return 0.5 * delta * (cos(M_PI * rad) + 1);
    else return 0.;
}

#if defined(WATER)
double vvm::Init::GetQVB(int k, vvm &model) {
    if (model.z[k] <= 4000) return 0.0161 - 0.000003375 * model.z[k];
    else if (4000 < model.z[k] && model.z[k] <= 8000) return 0.0026 - 0.00000065 * (model.z[k] - 4000);
    else return 0.;
}
#endif

#if defined(LOADFILE)
void print_data(const std::vector<double>& heights,
                const std::vector<std::vector<double>>& data_fields) {
    // Check if data is empty
    if (heights.empty() || data_fields.empty()) {
        std::cout << "No data to print." << std::endl;
        return;
    }

    // Print header
    std::cout << std::setw(12) << "Height";
    for (size_t j = 0; j < data_fields.size(); ++j) {
        std::cout << std::setw(12) << "Field " + std::to_string(j + 1);
    }
    std::cout << std::endl;

    // Print data
    for (size_t i = 0; i < heights.size(); ++i) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << heights[i];
        for (size_t j = 0; j < data_fields.size(); ++j) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(2) << data_fields[j][i];
        }
        std::cout << std::endl;
    }
}

void vvm::Init::LoadFile(vvm &model) {
    std::string filename = "../input/bubble_init.txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::vector<std::vector<double>> temp_data;
    std::string line;
    size_t num_columns = 0;

    // Read each line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        double value;

        // Read all values in the line
        while (ss >> value) {
            row.push_back(value);
        }

        // Skip empty lines
        if (row.empty()) {
            continue;
        }

        // Verify consistent number of columns
        if (num_columns == 0) {
            num_columns = row.size();
            if (num_columns < 2) {
                throw std::runtime_error("File must have at least two columns (height and one data field): " + filename);
            }
        } else if (row.size() != num_columns) {
            throw std::runtime_error("Inconsistent number of columns in file: " + filename);
        }

        temp_data.push_back(row);
    }

    file.close();

    if (temp_data.empty()) {
        throw std::runtime_error("File is empty: " + filename);
    }

    // Ensure at least two points for interpolation/extrapolation
    if (temp_data.size() < 2) {
        throw std::runtime_error("File must contain at least two data points: " + filename);
    }

    std::vector<double> heights;
    std::vector<std::vector<double>> data_fields;

    // Split into heights and data fields
    heights.resize(temp_data.size());
    data_fields.resize(num_columns - 1); // One column for height, rest for data fields
    for (auto& field : data_fields) {
        field.resize(temp_data.size());
    }

    for (size_t i = 0; i < temp_data.size(); ++i) {
        heights[i] = temp_data[i][0]; // First column is height
        for (size_t j = 1; j < num_columns; ++j) {
            data_fields[j - 1][i] = temp_data[i][j]; // Other columns are data fields
        }
    }

    // Verify that heights is sorted in ascending order
    if (!std::is_sorted(heights.begin(), heights.end())) {
        throw std::runtime_error("Height data (first column) must be sorted in ascending order: " + filename);
    }

    // print_data(heights, data_fields);

    // Interpolate loaded input to the stretched coordinate
    std::vector<double> new_heights(model.z+1, model.z + model.nz-1);
    std::vector<std::vector<double>> interpolated_data_fields;
    vvm::NumericalProcess::interpolate(heights, data_fields, new_heights, interpolated_data_fields);

    // print_data(new_heights, interpolated_data_fields);

    // Assign interpolated data to model variables
    for (int k = 1; k < model.nz-1; k++) {
        model.rhou[k] = interpolated_data_fields[0][k-1];
        model.thb[k] = interpolated_data_fields[1][k-1];
        model.thb_init[k] = model.thb[k];
        model.pb[k] = interpolated_data_fields[2][k-1];
        model.pib[k] = interpolated_data_fields[3][k-1];
        model.qvb[k] = interpolated_data_fields[4][k-1];
        model.qvb0[k] = model.qvb[k];
        model.thvb[k] = model.thb[k] * (1. + 0.608 * model.qvb[k]);
        model.thvbm[k] = model.thvb[k];
        #if defined(TROPICALFORCING)
            model.Q1LS[k] = interpolated_data_fields[5][k-1] * 6.;
            model.Q2LS[k] = interpolated_data_fields[6][k-1] * 6.;
        #endif
        model.rhow[k] = interpolated_data_fields[7][k-1];
        model.RH[k] = interpolated_data_fields[8][k-1];
        model.ubar[k] = interpolated_data_fields[9][k-1];
    }
    model.BoundaryProcess1D_center(model.rhou, model.nz);
    model.BoundaryProcess1D_center(model.thb, model.nz);
    model.BoundaryProcess1D_center(model.thb_init, model.nz);
    model.BoundaryProcess1D_center(model.pb, model.nz);
    model.BoundaryProcess1D_center(model.pib, model.nz);
    model.BoundaryProcess1D_center(model.qvb, model.nz);
    model.BoundaryProcess1D_center(model.qvb0, model.nz);
    model.BoundaryProcess1D_center(model.thvb, model.nz);
    #if defined(TROPICALFORCING)
        model.BoundaryProcess1D_center(model.Q1LS, model.nz);
        model.BoundaryProcess1D_center(model.Q2LS, model.nz);
    #endif
    model.BoundaryProcess1D_center(model.rhow, model.nz);
    model.BoundaryProcess1D_center(model.RH, model.nz);
    model.BoundaryProcess1D_center(model.ubar, model.nz);
    model.rhow[model.nz-1] = model.rhou[model.nz-2];

    // self defined RH
    double z1 = 6000.;
    for (int k = 1; k <= model.nz-2; k++) {
        if (model.z[k] <= z1) model.RH[k] = 0.9;
        else if (z1 <=model.z[k] && model.z[k] <= model.z[model.nz-2]) {
            model.RH[k] = 0.9 - (0.9 / (model.z[model.nz-2] - z1)) * (model.z[k] - z1);
        }
        else model.RH[k] = 0.;
    }

    for (int k = 1; k < model.nz-1; k++) {
        // Give pibar by pbar rather than given input
        model.pib[k] = std::pow(model.pb[k]/100000., model.Rd/model.Cp);

        double Tc = model.thb[k] * model.pib[k];
        double es = 611.2 * std::exp(17.67 * (Tc-273.15) / (Tc-273.15+243.5));
        model.qvsb[k] = 0.622 * es/ (model.pb[k] - 0.378 * es);

        #if !defined(TROPICALFORCING)
            // Give qvb by RH rather than given input
            model.qvb[k] = model.RH[k] * model.qvsb[k];
        #endif
        model.thvb[k] = model.thb[k] * (1. + 0.608 * model.qvb[k]);
    }
    model.BoundaryProcess1D_center(model.pib, model.nz);
    model.BoundaryProcess1D_center(model.thvb, model.nz);
    model.BoundaryProcess1D_center(model.qvsb, model.nz);
    model.BoundaryProcess1D_center(model.qvb, model.nz);

    for (int k = 1; k <= model.nz-2; k++) {
        model.thb_zeta[k] = 0.5 * (model.thb[k] + model.thb[k-1]);
    }
	model.thb_zeta[1] = model.thb_zeta[2] - (model.thb_zeta[3] - model.thb_zeta[2]);
    model.BoundaryProcess1D_center(model.thb_zeta, model.nz);
    model.thb_zeta[model.nz-1] = model.thb[model.nz-2];


    for (int k = 1; k <= model.nz-2; k++) {
        model.pib_lev[k] = 0.5*(model.pib[k]+model.pib[k-1]);
    }
    // extrapolation
    model.pib_lev[1] = model.pib_lev[2] + (model.pib_lev[2]-model.pib_lev[3]);
    model.pib_lev[model.nz-1] = model.pib_lev[model.nz-2]; // Give a small value for top boundary
    model.BoundaryProcess1D_center(model.pib_lev, model.nz+1);


    // init pb_lev (nz: model.nz+1)
    for (int k = 1; k <= model.nz-2; k++) {
        model.pb_lev[k] = 0.5*(model.pb[k]+model.pb[k-1]);
    }
    // model.pb_lev[model.nz-1] = model.pb_lev[model.nz-2]; // Give a small value for model top
    model.pb_lev[1] = model.pb_lev[2] - (model.pb_lev[3]-model.pb_lev[2]);
    model.pb_lev[model.nz-1] = model.pb_lev[model.nz-2] - (model.pb_lev[model.nz-3]-model.pb_lev[model.nz-2]); // extrapolation
    model.PSURF = model.pb_lev[1];
    model.BoundaryProcess1D_center(model.pb_lev, model.nz+1);
    return;
}
#elif defined(LOADFROMPREVIOUSFILE)
void vvm::Init::LoadFromPreviousFile(vvm &model) {
    std::ifstream inputFile;

    inputFile.open(LOADINITPATH);
    std::string line;
    std::getline(inputFile, line);
    std::getline(inputFile, line); // Skip the zero level
    double ZZ, ZT, RHO, THBAR, PBAR, PIBAR, QVBAR, Q1LS, Q2LS, RHOZ;

    int i = 1;
    while (inputFile >> ZZ >> ZT >> RHO >> THBAR >> PBAR >> PIBAR >> QVBAR >> Q1LS >> Q2LS >> RHOZ) {
        model.thb[i] = THBAR;
        model.qvb[i] = QVBAR;
        model.pib[i] = PIBAR;
        model.pb[i] = PBAR;
        model.rhou[i] = RHOZ;
        model.rhow[i] = RHO;
        #if defined(TROPICALFORCING)
            model.Q1LS[i] = Q1LS * 6.;
            model.Q2LS[i] = Q2LS * 6.;
        #endif
        i++;
    }

    model.BoundaryProcess1D_center(model.pib, model.nz);
    model.BoundaryProcess1D_center(model.pb, model.nz);
    model.BoundaryProcess1D_center(model.rhou, model.nz);
    model.BoundaryProcess1D_center(model.rhow, model.nz);
    model.rhow[model.nz-1] = model.rhou[model.nz-2];

    std::string data_m = LOADPATH1;
    std::string data = LOADPATH2;
    NcFile df_m(data_m, NcFile::read);
    NcFile df(data, NcFile::read);

    auto thm_in = df_m.getVar("th");
    auto th_in = df.getVar("th");
    auto zetam_in = df_m.getVar("zeta");
    auto zeta_in = df.getVar("zeta");
    auto qvm_in = df_m.getVar("qv");
    auto qv_in = df.getVar("qv");
    auto qcm_in = df_m.getVar("qc");
    auto qc_in = df.getVar("qc");
    auto qrm_in = df_m.getVar("qr");
    auto qr_in = df.getVar("qr");
    auto u_in = df.getVar("u");
    auto w_in = df.getVar("w");
    auto ubarm_in = df_m.getVar("ubarTop");

    thm_in.getVar(model.thmcont);
    th_in.getVar(model.thcont);
    zetam_in.getVar(model.zetamcont);
    zeta_in.getVar(model.zetacont);
    qvm_in.getVar(model.qvmcont);
    qv_in.getVar(model.qvcont);
    qcm_in.getVar(model.qcmcont);
    qc_in.getVar(model.qccont);
    qrm_in.getVar(model.qrmcont);
    qr_in.getVar(model.qrcont);
    u_in.getVar(model.ucont);
    w_in.getVar(model.wcont);
    double tmp[1];
    ubarm_in.getVar(tmp);
    model.ubarTopm = tmp[0];

    if (model.CASE == 1) {
        for (int k = 1; k < model.nz-1; k++) {
            for (int i = 1; i < model.nx-1; i++) {
                model.thm[i][k] += vvm::Init::GetTH(i, k, model);
                model.th[i][k] += vvm::Init::GetTH(i, k, model);
            }
        }
    }
    vvm::BoundaryProcess2D_center(model.thm, model.nx, model.nz);
    vvm::BoundaryProcess2D_center(model.th, model.nx, model.nz);

    double tb = 0.;
    #if defined(WATER)
        double qvb = 0.;
    #endif
    for (int k = 1; k < model.nz-1; k++) {
        tb = 0.;
        #if defined(WATER)
            qvb = 0.;
        #endif
        
        for (int i = 1; i < model.nx-1; i++) {
            tb += model.th[i][k];
            #if defined(WATER)
                qvb += model.qv[i][k];
            #endif
        }
        model.thb[k] = tb / (double) (model.nx - 2.);
        #if defined(WATER)
            model.qvb[k] = qvb / (double) (model.nx - 2.);
            model.thvb[k] = model.thb[k] * (1 + 0.608 * model.qvb[k]);
        #else
            model.thvb[k] = model.thb[k];
        #endif
        model.qvsb[k] = (380. / model.pb[k]) * exp((17.27 * (model.thb[k] * model.pib[k] - 273.)) / (model.thb[k] * model.pib[k] - 36.));
    }
    model.BoundaryProcess1D_center(model.thb, model.nz);
    model.BoundaryProcess1D_center(model.thvb, model.nz);
    #if defined(WATER)
        model.BoundaryProcess1D_center(model.qvb, model.nz);
    #endif
    model.BoundaryProcess1D_center(model.qvsb, model.nz);

    for (int k = 1; k < model.nz-1; k++) {
        model.thb_zeta[k] = 0.5 * (model.thb[k-1] + model.thb[k]);
    }
    model.BoundaryProcess1D_center(model.thb_zeta, model.nz);
    return;
}

#endif

void vvm::Init::RandomPerturbation(vvm &model, int t, double min_range, double max_range, double standard_deviation) {
    std::mt19937 gen(t); // Mersenne Twister engine for random numbers
    gen.seed(t);
    std::normal_distribution<> distribution(0.0, 1.0); // Gaussian distribution with mean 0 and standard deviation 1

    // Parameters for the 2D Gaussian noise array
    double mean = 0.; // Mean of the Gaussian distribution

    for (int k = 1; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            if (model.z[k] < 200) {
                double random_noise = 0.;
                do {
                    random_noise = mean + standard_deviation * distribution(gen);
                } while (random_noise < min_range || random_noise > max_range);

                model.init_th_forcing[i][k] = random_noise;
            }
            else {
                model.init_th_forcing[i][k] = 0.;
            }
        }
    }
    model.BoundaryProcess2D_center(model.init_th_forcing, model.nx, model.nz);
}


