#include "Init.hpp"

void Init::Init1d(vvmArray &model) {
	#if defined(LOADFILE)
		LoadFile(model);
	#else
		// init tb
		model.tb[1] = 300.;
		for (int k = 2; k <= nz-2; k++) {
			#ifdef DRY
				model.tb[k] = 300.;
			#else
				model.tb[k] = GetTB(k);
			#endif
		}
		model.tb[0] = model.tb[1];
		model.tb[nz-1] = model.tb[nz-2];

		// init qvb, tvb
		for (int k = 1; k <= nz-2; k++) {
			#if defined(WATER)
				model.qvb[k] = GetQVB(k);
			#else
				model.qvb[k] = 0.;
			#endif
			model.tvb[k] = model.tb[k] * (1. + 0.61 * model.qvb[k]);
		}
		model.qvb[0] = model.qvb[1];
		model.qvb[nz-1] = model.qvb[nz-2];
		model.tvb[0] = model.tvb[1];
		model.tvb[nz-1] = model.tvb[nz-2];

		// init pib
		double pisfc = pow((PSURF / P0), Rd / C_p);
		for (int k = 1; k <= nz-2; k++) {
			if (k == 1) model.pib[k] = pisfc - gravity * 0.5 * dz / (C_p * model.tvb[k]);
			else {
				double tvbavg = 0.5*(model.tvb[k] + model.tvb[k-1]);
				model.pib[k] = model.pib[k-1] - gravity * dz / (C_p * tvbavg);
			}
		}
		model.pib[0] = model.pib[1];
		model.pib[nz-1] = model.pib[nz-2];

		// init tb_zeta, rhou
		for (int k = 1; k <= nz-2; k++) {
			model.tb_zeta[k] = 0.5 * (model.tb[k-1] + model.tb[k]);
			#ifdef RHO1
				model.rhou[k] = 1.;
			#else
				model.rhou[k] = P0 * pow(model.pib[k], Cv/Rd) / (Rd * model.tvb[k]);
			#endif
		}
		model.rhou[0] = model.rhou[1];
		model.rhou[nz-1] = model.rhou[nz-2];

		// init tb_zeta, rhow
		for (int k = 1; k <= nz-1; k++) {
			model.tb_zeta[k] = 0.5 * (model.tb[k] + model.tb[k-1]);
			model.rhow[k] = 0.5 * (model.rhou[k] + model.rhou[k-1]);
		}
		model.tb_zeta[0] = model.tb_zeta[1];
		model.rhow[0] = model.rhow[1];

		// init pb, qvsb
		for (int k = 1; k <= nz-2; k++) {
			model.pb[k] = P0 * pow(model.pib[k], C_p / Rd);
			double Tbar = model.tb[k] * model.pib[k];
			model.qvsb[k] = (380. / model.pb[k]) * exp((17.27 * (Tbar - 273.)) / (Tbar - 36.));
		}
		model.pb[0] = model.pb[1];
		model.pb[nz-1] = model.pb[nz-2];
		model.qvsb[0] = model.qvsb[1];
		model.qvsb[nz-1] = model.qvsb[nz-2];

		// heat flux init
		#if defined(HEATFLUX)
			mt19937 mt(20210831);
			uniform_real_distribution<> distr(-1, 1);
			for (int i = 1; i <= nx-2; i++) {
				model.addflx[i] += distr(mt);
			}
		#endif
	#endif
	return;
}

void Init::Init2d(vvmArray &model) {
	#if defined(TROPICALFORCING)
		// Generate random 2D Gaussian noise array within the specified range
		RandomPerturbation(model, 0);

        for (int i = 1; i <= nx-2; i++) {
            for (int k = 1; k <= nz-2; k++) {
				#if defined(LINEARIZEDTH)
                	model.th[i][k] = model.init_th_forcing[i][k];
				#else
					model.th[i][k] = model.tb[k] + model.init_th_forcing[i][k];
				#endif
                model.thm[i][k] = model.th[i][k];

				#if defined(LINEARIZEDQV)
					model.qv[i][k] = 0.;
				#else
                	model.qv[i][k] = model.qvb[k];
				#endif
                model.qvm[i][k] = model.qv[i][k];

                #if defined(SHEAR)
                    if ((k-0.5) * dz <= 5000) {
                        model.u[i][k] = 0.004 * (k - 0.5) * dz - 10.5;
                    }
                    else {
                        model.u[i][k] = 0.001 * (k - 0.5) * dz + 5.5;
                    }
                #else
                    model.u[i][k] = 0.;
                    model.w[i][k] = 0.;
                #endif
            }
        }
        model.BoundaryProcess(model.th);
        model.BoundaryProcess(model.thm);
        model.BoundaryProcess(model.qv);
        model.BoundaryProcess(model.qvm);
        model.BoundaryProcess(model.u);
        model.BoundaryProcess(model.w);

    #else
		// init th
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				#if defined(LINEARIZEDTH)
					model.th[i][k] = GetTH(i, k);
				#else
					model.th[i][k] = model.th[k] + GetTH(i, k);
				#endif
				model.thm[i][k] = model.th[i][k];
			}
		}
		model.BoundaryProcess(model.th);
		model.BoundaryProcess(model.thm);

		// init qv: where th != 0, qv = qvs
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				// if (model.th[i][k] != 0) model.qv[i][k] = model.qvsb[k] - model.qvb[k];
				// else model.qv[i][k] = 0.;
				#if defined(LINEARIZEDQV)
					model.qv[i][k] = 0.;
				#else
					model.qv[i][k] = model.qvb[k];
				#endif
				model.qvm[i][k] = model.qv[i][k];
			}
		}
		model.BoundaryProcess(model.qv);
		model.BoundaryProcess(model.qvm);

		// init u
		#if defined(SHEAR)
			// From level to bubble center (umin -> 0), from bubble center to top (0 -> umax)
			double umin = -10, umax = 10.;
			int bubble_center_idx = 2500/dz + 1;
			for (int i = 1; i <= nx-2; i++) {
				for (int k = 1; k <= nz-2; k++) {
					if (k <= bubble_center_idx) model.u[i][k] = umin - umin / (bubble_center_idx-1) * (k-1);
					else model.u[i][k] = umax / (nz-2 - (bubble_center_idx-1)) * (k-bubble_center_idx+1);
				}
			}
			model.BoundaryProcess(model.u);
		#endif
	#endif

	// init zeta
	double pu_pz = 0., pw_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			pw_px = (model.w[i][k] - model.w[i-1][k]) * rdx;
			pu_pz = (model.u[i][k] - model.u[i][k-1]) * rdz;
			model.zeta[i][k] = pw_px - pu_pz;
			model.zetam[i][k] = model.zeta[i][k];
		}
	}
	model.BoundaryProcessZETA(model.zeta);
	model.BoundaryProcessZETA(model.zetam);
}


double Init::GetTB(int k) {
	double z_top = 12000., T_top = 213., tb_top = 343.;
	double z_t = dz * (k - 0.5);
	if (z_t <= z_top) return 300. + 43. * pow(z_t / z_top, 1.25);
	else return tb_top * exp(gravity * (z_t - z_top) / (C_p * T_top));
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
	else return 0.;
}

double Init::GetQVB(int k) {
	double z_t = (k - 0.5) * dz;
	if (z_t <= 4000) return 0.0161 - 0.000003375 * z_t;
	else if (4000 < z_t && z_t <= 8000) return 0.0026 - 0.00000065 * (z_t - 4000);
	else return 0.;
}

void Init::LoadFile(vvmArray &model) {
	std::ifstream inputFile;

	inputFile.open("../input/init.txt");
	std::string line;
	std::getline(inputFile, line);
	std::getline(inputFile, line); // Skip the zero level
	double ZZ, ZT, RHO, THBAR, PBAR, PIBAR, QVBAR, Q1LS, Q2LS, RHOZ;

	int i = 1;
	while (inputFile >> ZZ >> ZT >> RHO >> THBAR >> PBAR >> PIBAR >> QVBAR >> Q1LS >> Q2LS >> RHOZ) {
		model.tb[i] = THBAR;
		model.qvb[i] = QVBAR;
		model.pib[i] = PIBAR;
		model.pb[i] = PBAR;
		model.rhou[i] = RHOZ;
		model.rhow[i] = RHO;
		#if defined(TROPICALFORCING)
			model.Q1LS[i] = Q1LS * 6.;
			model.Q2LS[i] = Q2LS * 6.;
		#endif

		model.tvb[i] = model.tb[i] * (1. + 0.61 * model.qvb[i]);
		model.qvsb[i] = (380. / model.pb[i]) * exp((17.27 * (model.tb[i] * model.pib[i] - 273.)) / (model.tb[i] * model.pib[i] - 36.));
		i++;
	}

	model.tb[0] = model.tb[1];
	model.tb[nz-1] = model.tb[nz-2];
	model.qvb[0] = model.qvb[1];
	model.qvb[nz-1] = model.qvb[nz-2];
	model.pib[0] = model.pib[1];
	model.pib[nz-1] = model.pib[nz-2];
	model.pb[0] = model.pb[1];
	model.pb[nz-1] = model.pb[nz-2];
	model.rhou[0] = model.rhou[1];
	model.rhou[nz-1] = model.rhou[nz-2];
	model.rhow[0] = model.rhow[1];
	model.rhow[nz-1] = model.rhow[nz-2];
	model.tvb[0] = model.tvb[1];
	model.tvb[nz-1] = model.tvb[nz-2];
	model.qvsb[0] = model.qvsb[1];
	model.qvsb[nz-1] = model.qvsb[nz-2];

	for (int k = 1; k < nz-1; k++) {
		model.tb_zeta[k] = 0.5 * (model.tb[k-1] + model.tb[k]);
	}
	model.tb_zeta[0] = model.tb_zeta[1];
	model.tb_zeta[nz-1] = model.tb_zeta[nz-2];

	return;
}

#if defined(TROPICALFORCING)
void Init::RandomPerturbation(vvmArray &model, int t) {
    std::mt19937 gen(t); // Mersenne Twister engine for random numbers
    std::normal_distribution<> distribution(0.0, 1.0); // Gaussian distribution with mean 0 and standard deviation 1

    // Parameters for the 2D Gaussian noise array
    double mean = 0.; // Mean of the Gaussian distribution
    double standard_deviation = 1.; // Standard deviation of the Gaussian distribution
    double min_range = -1.; // Minimum value of the generated noise
    double max_range = 1.; // Maximum value of the generated noise

    for (int i = 1; i < nx-1; i++) {
        for (int k = 1; k < nz-1; k++) {
            if (k <= nz / 15) {
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
    model.BoundaryProcess(model.init_th_forcing);
}
#endif