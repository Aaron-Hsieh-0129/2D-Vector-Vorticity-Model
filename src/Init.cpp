#include "Init.hpp"
#include "Eigen/src/SparseCore/SparseMatrixBase.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <petsc.h>
#include <random>

void Init::Init1d(vvm &model) {
    #if defined(LOADFILE)
        LoadFile(model);
    #else
        // init tb
        model.thb[1] = 300.;
        for (int k = 2; k <= model.nz-2; k++) {
            #ifdef DRY
                model.thb[k] = 300.;
            #else
                model.thb[k] = GetTB(k);
            #endif
        }
        model.BoundaryProcess1D_center(model.thb);

        // init qvb, tvb
        for (int k = 1; k <= model.nz-2; k++) {
            #if defined(WATER)
                model.qvb[k] = GetQVB(k);
            #else
                model.qvb[k] = 0.;
            #endif
            model.tvb[k] = model.thb[k] * (1. + 0.61 * model.qvb[k]);
        }
        model.BoundaryProcess1D_center(model.qvb);
        model.BoundaryProcess1D_center(model.tvb);

        // init pib
        double pisfc = pow((PSURF / P0), Rd / C_p);
        for (int k = 1; k <= model.nz-2; k++) {
            if (k == 1) model.pib[k] = pisfc - GRAVITY * 0.5 * dz / (C_p * model.tvb[k]);
            else {
                double tvbavg = 0.5*(model.tvb[k] + model.tvb[k-1]);
                model.pib[k] = model.pib[k-1] - GRAVITY * dz / (C_p * tvbavg);
            }
        }
        model.BoundaryProcess1D_center(model.pib);

        // init tb_zeta, rhou
        for (int k = 1; k <= model.nz-2; k++) {
            #ifdef RHO1
                model.rhou[k] = 1.;
            #else
                model.rhou[k] = P0 * pow(model.pib[k], Cv/Rd) / (Rd * model.tvb[k]);
            #endif
        }
        model.BoundaryProcess1D_center(model.rhou);

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
            model.pb[k] = P0 * pow(model.pib[k], C_p / Rd);
            double Tbar = model.thb[k] * model.pib[k];
            model.qvsb[k] = (380. / model.pb[k]) * exp((17.27 * (Tbar - 273.)) / (Tbar - 36.));
        }
        model.BoundaryProcess1D_center(model.pb);
        model.BoundaryProcess1D_center(model.qvsb);

        #if defined(WATER)
            for (int k = 1; k <= model.nz-2; k++) {
                model.qvb[k] = GetQVB(k);
            }
            model.BoundaryProcess1D_center(model.qvb);
        #endif

        #if defined(RHO1)
            for (int k = 0; k < model.nz; k++) {
                model.rhou[k] = model.rhow[k] = 1.;
            }
        #endif
    #endif

    for (int k = 0; k < model.nz; k++) {
        model.thbm[k] = model.thb[k];
    }
    return;
}

void Init::Init2d(vvm &model) {
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
        model.BoundaryProcess2D_center(model.th);
        model.BoundaryProcess2D_center(model.thm);
        model.BoundaryProcess2D_center(model.qv);
        model.BoundaryProcess2D_center(model.qvm);
    #else
        // init th
        for (int i = 1; i <= model.nx-2; i++) {
            for (int k = 1; k <= model.nz-2; k++) {
                model.th[i][k] = GetTH(i, k);
                model.th[i][k] = model.thb[k] + GetTH(i, k);
                model.thm[i][k] = model.th[i][k];
            }
        }
        model.BoundaryProcess2D_center(model.th);
        model.BoundaryProcess2D_center(model.thm);

        #if defined(WATER)
            // init qv: where th != 0, qv = qvs
            for (int i = 1; i <= model.nx-2; i++) {
                for (int k = 1; k <= model.nz-2; k++) {
                    model.qv[i][k] = model.qvm[i][k] = model.qvb[k];
                }
            }
            model.BoundaryProcess2D_center(model.qv);
            model.BoundaryProcess2D_center(model.qvm);
        #endif

		// init u
		#if defined(SHEAR)
			// From level to bubble center (umin -> 0), from bubble center to top (0 -> umax)
			double umin = -10, umax = 10.;
			int bubble_center_idx = 2500/dz + 1;
			for (int i = 1; i <= model.nx-2; i++) {
				for (int k = 1; k <= model.nz-2; k++) {
					if (k <= bubble_center_idx) model.u[i][k] = umin - umin / (bubble_center_idx-1) * (k-1);
					else model.u[i][k] = umax / (model.nz-2 - (bubble_center_idx-1)) * (k-bubble_center_idx+1);
				}
			}
			model.BoundaryProcess2D_center(model.u);
		#endif
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
	model.BoundaryProcess2D_westdown(model.zeta);
	model.BoundaryProcess2D_westdown(model.zetam);

	// init ubar at top
	for (int i = 1; i < model.nx-1; i++) {
		model.ubarTopm += model.u[i][model.nz-2];
	}
	model.ubarTopm /= ((double) (model.nx - 2.));
	model.ubarTop = model.ubarTopm;

    #ifndef PETSC
	    // Assign values to the matrices that solve the Poisson equation for u and w
        vvm::PoissonSolver::InitPoissonMatrix(model);
    #endif
	return;
}

double Init::GetTB(int k) {
    double z_top = 12000., T_top = 213., tb_top = 343.;
    double z_t = dz * (k - 0.5);
    if (z_t <= z_top) return 300. + 43. * pow(z_t / z_top, 1.25);
    else return tb_top * exp(GRAVITY * (z_t - z_top) / (C_p * T_top));
}

double Init::GetTHRAD(int i, int k) {
    double XC = XRANGE / 2., XR = 4000.;
    double ZC = 2500., ZR = 2000.;
    double x = (i-0.5) * dx, z = (k-0.5) * dz;
    double rad = sqrt(pow((x - XC) / XR, 2) + pow((z- ZC) / ZR, 2));
    return rad;
}

double Init::GetTH(int i, int k) {
    double rad = GetTHRAD(i, k);
    double delta = 3.;
    if (rad <= 1) return 0.5 * delta * (cos(M_PI * rad) + 1);
    // if (k >= 10 && k <= 20) return 3.;
    else return 0.;
}

#if defined(WATER)
double Init::GetQVB(int k) {
    double z_t = (k - 0.5) * dz;
    if (z_t <= 4000) return 0.0161 - 0.000003375 * z_t;
    else if (4000 < z_t && z_t <= 8000) return 0.0026 - 0.00000065 * (z_t - 4000);
    else return 0.;
}
#endif

void Init::LoadFile(vvm &model) {
    std::ifstream inputFile;

    inputFile.open("../input/init.txt");
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

        model.tvb[i] = model.thb[i] * (1. + 0.61 * model.qvb[i]);
        model.qvsb[i] = (380. / model.pb[i]) * exp((17.27 * (model.thb[i] * model.pib[i] - 273.)) / (model.thb[i] * model.pib[i] - 36.));
        i++;
    }

    model.BoundaryProcess1D_center(model.thb);
    model.BoundaryProcess1D_center(model.qvb);
    model.BoundaryProcess1D_center(model.pib);
    model.BoundaryProcess1D_center(model.pb);
    model.BoundaryProcess1D_center(model.rhou);
    model.BoundaryProcess1D_center(model.rhow);
    model.BoundaryProcess1D_center(model.tvb);
    model.BoundaryProcess1D_center(model.qvsb);
    model.rhow[NZ-1] = model.rhou[NZ-2];

    for (int k = 1; k <= model.nz-2; k++) {
        model.thb_zeta[k] = 0.5 * (model.thb[k] + model.thb[k-1]);
    }
	model.thb_zeta[1] = model.thb_zeta[2] - (model.thb_zeta[3] - model.thb_zeta[2]);
    model.BoundaryProcess1D_center(model.thb_zeta);
    model.thb_zeta[NZ-1] = model.thb[NZ-2];
    return;
}

#if defined(TROPICALFORCING)
void Init::RandomPerturbation(vvm &model, int t) {
    std::mt19937 gen(t); // Mersenne Twister engine for random numbers
    std::normal_distribution<> distribution(0.0, 1.0); // Gaussian distribution with mean 0 and standard deviation 1

    // Parameters for the 2D Gaussian noise array
    double mean = 0.; // Mean of the Gaussian distribution
    double standard_deviation = 1.; // Standard deviation of the Gaussian distribution
    double min_range = -1.; // Minimum value of the generated noise
    double max_range = 1.; // Maximum value of the generated noise

    for (int i = 1; i < model.nx-1; i++) {
        for (int k = 1; k < model.nz-1; k++) {
            if (k <= model.nz / 15) {
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
    model.BoundaryProcess2D_center(model.init_th_forcing);
}
#endif


