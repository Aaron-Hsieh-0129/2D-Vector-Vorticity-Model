#include "Init.hpp"

void Init::Init1d(vvmArray & myArray) {
	// init tb
	myArray.tb[1] = 300.;
	for (int k = 2; k <= nz-2; k++) {
		#ifdef DRY
			myArray.tb[k] = 300.;
		#else
			myArray.tb[k] = GetTB(k);
		#endif
	}
	myArray.tb[0] = myArray.tb[1];
	myArray.tb[nz-1] = myArray.tb[nz-2];

	// init qvb, tvb
	for (int k = 1; k <= nz-2; k++) {
		#if defined(WATER)
			myArray.qvb[k] = GetQVB(k);
		#else
			myArray.qvb[k] = 0.;
		#endif
		myArray.tvb[k] = myArray.tb[k] * (1. + 0.61 * myArray.qvb[k]);
	}
	myArray.qvb[0] = myArray.qvb[1];
	myArray.qvb[nz-1] = myArray.qvb[nz-2];
	myArray.tvb[0] = myArray.tvb[1];
	myArray.tvb[nz-1] = myArray.tvb[nz-2];

	// init pib
	double pisfc = pow((PSURF / P0), Rd / C_p);
	for (int k = 1; k <= nz-2; k++) {
		if (k == 1) myArray.pib[k] = pisfc - gravity * 0.5 * dz / (C_p * myArray.tvb[k]);
		else {
			double tvbavg = 0.5*(myArray.tvb[k] + myArray.tvb[k-1]);
			myArray.pib[k] = myArray.pib[k-1] - gravity * dz / (C_p * tvbavg);
		}
	}
	myArray.pib[0] = myArray.pib[1];
	myArray.pib[nz-1] = myArray.pib[nz-2];

	// init tb_zeta, rhou
	for (int k = 1; k <= nz-2; k++) {
		myArray.tb_zeta[k] = 0.5 * (myArray.tvb[k-1] + myArray.tvb[k]);
		#ifdef RHO1
			myArray.rhou[k] = 1.;
		#else
			myArray.rhou[k] = P0 * pow(myArray.pib[k], Cv/Rd) / (Rd * myArray.tvb[k]);
		#endif
	}
	myArray.rhou[0] = myArray.rhou[1];
	myArray.rhou[nz-1] = myArray.rhou[nz-2];

	// init tb_zeta, rhow
	for (int k = 1; k <= nz-1; k++) {
		myArray.tb_zeta[k] = 0.5 * (myArray.tb[k] + myArray.tb[k-1]);
		myArray.rhow[k] = 0.5 * (myArray.rhou[k] + myArray.rhou[k-1]);
	}
	myArray.tb_zeta[0] = myArray.tb_zeta[1];
	myArray.rhow[0] = myArray.rhow[1];

	// init pb, qvsb
	for (int k = 1; k <= nz-2; k++) {
		myArray.pb[k] = P0 * pow(myArray.pib[k], C_p / Rd);
		double Tbar = myArray.tb[k] * myArray.pib[k];
		myArray.qvsb[k] = (380. / myArray.pb[k]) * exp((17.27 * (Tbar - 273.)) / (Tbar - 36.));
	}
	myArray.pb[0] = myArray.pb[1];
	myArray.pb[nz-1] = myArray.pb[nz-2];
	myArray.qvsb[0] = myArray.qvsb[1];
	myArray.qvsb[nz-1] = myArray.qvsb[nz-2];

	// heat flux init
	#if defined(HEATFLUX)
		mt19937 mt(20210831);
		uniform_real_distribution<> distr(-1, 1);
		for (int i = 1; i <= nx-2; i++) {
			myArray.addflx[i] += distr(mt);
		}
	#endif
	return;
}

void Init::Init2d(vvmArray & myArray) {
	// init th
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			myArray.th[i][k] = GetTH(i, k);
			myArray.thm[i][k] = myArray.th[i][k];
		}
	}
	myArray.BoundaryProcess(myArray.th);
	myArray.BoundaryProcess(myArray.thm);

	// init qv: where th != 0, qv = qvs
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			// if (myArray.th[i][k] != 0) myArray.qv[i][k] = myArray.qvsb[k] - myArray.qvb[k];
			// else myArray.qv[i][k] = 0.;
			myArray.qv[i][k] = 0.;
			myArray.qvm[i][k] = myArray.qv[i][k];
		}
	}
	myArray.BoundaryProcess(myArray.qv);
	myArray.BoundaryProcess(myArray.qvm);

	// init u
	#if defined(SHEAR)
		double umax = 10.;
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				if ((k-0.5) * dz <= 5000) {
					myArray.u[i][k] = 0.004 * (k - 0.5) * dz - 10.5;
				}
				else {
					myArray.u[i][k] = 0.001 * (k - 0.5) * dz + 5.5;
				}
			}
		}
		myArray.BoundaryProcess(myArray.u);

	#elif defined(ADVECTIONU)
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
                myArray.u[i][k] = 100.;
            }
        }
		myArray.BoundaryProcess(myArray.u);

    #elif defined(ADVECTIONW)
        for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
                myArray.w[i][k] = 100.;
            }
		}
		myArray.BoundaryProcess(myArray.w);
	#endif

	// init zeta
	double pu_pz = 0., pw_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			pw_px = (myArray.w[i][k] - myArray.w[i-1][k]) * rdx;
			pu_pz = (myArray.u[i][k] - myArray.u[i][k-1]) * rdz;
			myArray.zeta[i][k] = pw_px - pu_pz;
			myArray.zetam[i][k] = myArray.zeta[i][k];
		}
	}
	myArray.BoundaryProcessZETA(myArray.zeta);
	myArray.BoundaryProcessZETA(myArray.zetam);
}


double Init::GetTB(int k) {
	double z_top = 12000., T_top = 213., tb_top = 343.;
	double z_t = dz * (k - 0.5);
	if (z_t <= z_top) return 300. + 43. * pow(z_t / z_top, 1.25);
	else return tb_top * exp(gravity * (z_t - z_top) / (C_p * T_top));
}

double Init::GetTHRAD(int i, int k) {
	double XC = 74875., XR = 4000.;
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
