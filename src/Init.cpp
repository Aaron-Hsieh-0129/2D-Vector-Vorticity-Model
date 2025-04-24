#include "Declare.hpp"
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

void generateAddlfxArray(double *array, int size, double variation = 0.1) {
    std::mt19937 rng(static_cast<unsigned>(time(nullptr))); // Random number generator
    std::uniform_real_distribution<double> dist(-variation, variation); // Variation range

    for (int i = 0; i < size; ++i) {
        array[i] = 1.0 + dist(rng); // Mean of 1.0 with random variation
    }
    return;
}
void vvm::Init::Init1d(vvm &model) {
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
                model.qvb[k] = GetQVB(k, model.dz);
            #else
                model.qvb[k] = 0.;
            #endif
            model.thvb[k] = model.thb[k] * (1. + 0.61 * model.qvb[k]);
        }
        model.BoundaryProcess1D_center(model.qvb, model.nz);
        model.BoundaryProcess1D_center(model.qvb0, model.nz);
        model.BoundaryProcess1D_center(model.thvb, model.nz);

        // init pib
        double pisfc = pow((model.PSURF / model.P0), model.Rd / model.Cp);
        for (int k = 1; k <= model.nz-2; k++) {
            if (k == 1) model.pib[k] = pisfc - model.GRAVITY * 0.5 * model.dz / (model.Cp * model.thvb[k]);
            else {
                double tvbavg = 0.5*(model.thvb[k] + model.thvb[k-1]);
                model.pib[k] = model.pib[k-1] - model.GRAVITY * model.dz / (model.Cp * tvbavg);
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
            double Tbar = model.thb[k] * model.pib[k];
            model.qvsb[k] = (380. / model.pb[k]) * exp((17.27 * (Tbar - 273.)) / (Tbar - 36.));
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
                model.qvb[k] = GetQVB(k, model.dz);
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
            model.thvb[k] = model.thvbm[k] = model.thb[k] + 0.61 * model.qvb[k];
        #else
            model.thvb[k] = model.thvbm[k] = model.thb[k];
        #endif
    }

    for (int k = 1; k < model.nz-1; k++) {
        model.z[k] = (k-0.5) * model.dz;
        model.z_zeta[k] = (k-1) * model.dz;
        model.lambda2[k] = 1. / (1. / pow(0.23 * std::sqrt(model.dx*model.dz), 2) + 1. / pow(0.4* 0.4 * model.z[k], 2));
        model.lambda2_zeta[k] = 1. / (1. / pow(0.23 * std::sqrt(model.dx*model.dz), 2) + 1. / pow(0.4* 0.4 * model.z_zeta[k], 2));
    }
    model.BoundaryProcess1D_center(model.z, model.nz);
    model.BoundaryProcess1D_center(model.z_zeta, model.nz);
    model.BoundaryProcess1D_center(model.lambda2, model.nz);
    model.BoundaryProcess1D_center(model.lambda2_zeta, model.nz);


    for (int i = 0; i <= model.nx-1; i++) {
        model.th_ground[i] = 303.;
        double Tc = model.th_ground[i] * model.pib[1];
        double es = 611.2 * std::exp(17.67 * (Tc-273.15) / (Tc-273.15+243.5));
        model.qvs_ground[i] = 0.622 * es/ (model.P0 - 0.378 * es);
    }
    generateAddlfxArray(model.addflux, model.nx);
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
    #else
        // init th
        for (int i = 1; i <= model.nx-2; i++) {
            for (int k = 1; k <= model.nz-2; k++) {
                if (model.CASE == 0) model.th[i][k] = model.thb[k];
                else if (model.CASE == 1 || model.CASE == 2) model.th[i][k] = model.thb[k] + GetTH(i, k, model);

                if (model.addforcingtime > 0) {
                    RandomPerturbation(model, 0);
                    model.th[i][k] += model.init_th_forcing[i][k];
                }
                
                model.thm[i][k] = model.th[i][k];

                model.u[i][k] = 0.;
                model.w[i][k] = 0.;
            }
        }
        model.BoundaryProcess2D_center(model.th, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.thm, model.nx, model.nz);

        for (int k = 0; k <= model.nz-1; k++) {
            for (int i = 0; i <= model.nx-1; i++) {
                model.qv[i][k] = model.qvm[i][k] = model.qvb[k];
            }
            model.qvb0[k] = model.qvb[k];
        }

		// init u
		if (model.CASE == 2) {
            // From level to bubble center (umin -> 0), from bubble center to top (0 -> umax)
            double u_tmp = 0.;
            for (int k = 1; k <= model.nz-2; k++) {
                double z_u = (k-1) * model.dz;
                if (z_u <= 2000.) u_tmp = -8. + (8. / 2000.0) * z_u;

                for (int i = 1; i <= model.nx-2; i++) {
                    model.u[i][k] = u_tmp;
                }
            }
            model.BoundaryProcess2D_center(model.u, model.nx, model.nz);
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

    
    for (int i = 0; i < model.nx; i++) {
        for (int k = 0; k <= model.nz; k++) {
            model.pb_lev_all[i][k] = model.pb_lev[k];
        }
    }


	return;
}



double vvm::Init::GetTB(int k, vvm &model) {
    double z_top = 12000., T_top = 213., tb_top = 343.;
    double z_t = model.dz * (k - 0.5);
    if (z_t <= z_top) return 300. + 43. * pow(z_t / z_top, 1.25);
    else return tb_top * exp(model.GRAVITY * (z_t - z_top) / (model.Cp * T_top));
}

double vvm::Init::GetTHRAD(int i, int k, vvm &model) {
    double XC = model.XRANGE / 2., XR = 4000.;
    double ZC = 2500., ZR = 2000.;
    double x = (i-0.5) * model.dx, z = (k-0.5) * model.dz;
    double rad = sqrt(pow((x - XC) / XR, 2) + pow((z- ZC) / ZR, 2));
    return rad;
}

double vvm::Init::GetTH(int i, int k, vvm &model) {
    double rad = GetTHRAD(i, k, model);
    double delta = 3.;
    if (rad <= 1) return 0.5 * delta * (cos(M_PI * rad) + 1);
    else return 0.;
}

#if defined(WATER)
double vvm::Init::GetQVB(int k, int dz) {
    double z_t = (k - 0.5) * dz;
    if (z_t <= 4000) return 0.0161 - 0.000003375 * z_t;
    else if (4000 < z_t && z_t <= 8000) return 0.0026 - 0.00000065 * (z_t - 4000);
    else return 0.;
}
#endif

#if defined(LOADFILE)
void vvm::Init::LoadFile(vvm &model) {
    std::ifstream inputFile;

    inputFile.open(LOADINITPATH);
    std::string line;
    std::getline(inputFile, line);
    std::getline(inputFile, line); // Skip the zero level
    double ZZ, ZT, RHO, THBAR, PBAR, PIBAR, QVBAR, Q1LS, Q2LS, RHOZ;

    int i = 1;
    while (inputFile >> ZZ >> ZT >> RHO >> THBAR >> PBAR >> PIBAR >> QVBAR >> Q1LS >> Q2LS >> RHOZ) {
        model.thb[i] = THBAR;
        model.thb_init[i] = model.thb[i];
        model.qvb[i] = QVBAR;
        model.pib[i] = PIBAR;
        model.pb[i] = PBAR;
        model.rhou[i] = RHOZ;
        model.rhow[i] = RHO;
        #if defined(TROPICALFORCING)
            model.Q1LS[i] = Q1LS * 6.;
            model.Q2LS[i] = Q2LS * 6.;
        #endif

        model.thvb[i] = model.thb[i] * (1. + 0.61 * model.qvb[i]);
        model.qvsb[i] = (380. / model.pb[i]) * exp((17.27 * (model.thb[i] * model.pib[i] - 273.)) / (model.thb[i] * model.pib[i] - 36.));
        i++;
    }

    model.BoundaryProcess1D_center(model.thb, model.nz);
    model.BoundaryProcess1D_center(model.thb_init, model.nz);
    model.BoundaryProcess1D_center(model.qvb, model.nz);
    model.BoundaryProcess1D_center(model.pib, model.nz);
    model.BoundaryProcess1D_center(model.pb, model.nz);
    model.BoundaryProcess1D_center(model.rhou, model.nz);
    model.BoundaryProcess1D_center(model.rhow, model.nz);
    model.BoundaryProcess1D_center(model.thvb, model.nz);
    model.BoundaryProcess1D_center(model.qvsb, model.nz);
    model.rhow[model.nz-1] = model.rhou[model.nz-2];

    // Give pibar by pbar rather than given input
    for (int k = 0; k < model.nz-1; k++) model.pib[k] = std::pow(model.pb[k]/100000., model.Rd/model.Cp);

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
    model.pb_lev[1] = 100000;
    model.pb_lev[model.nz-1] = model.pb_lev[model.nz-2] - (model.pb_lev[model.nz-3]-model.pb_lev[model.nz-2]); // extrapolation
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

    if (model.CASE == 1 || model.CASE == 2) {
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
            model.thvb[k] = model.thb[k] * (1 + 0.61 * model.qvb[k]);
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

    double z = 0;

    for (int k = 1; k < model.nz-1; k++) {
        z = (k - 0.5) * model.dz;
        for (int i = 1; i < model.nx-1; i++) {
            if (z < 200) {
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


