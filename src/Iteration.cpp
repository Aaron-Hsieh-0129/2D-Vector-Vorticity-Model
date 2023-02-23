#include "Iteration.hpp"

void Iteration::pzeta_pt(vvmArray& myArray) {
	double puzeta_px = 0., pwzeta_pz = 0., g_tb_pth_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puzeta_px = (0.25*(myArray.u[i+1][k] + myArray.u[i+1][k-1] + myArray.u[i][k] + myArray.u[i][k-1]) * 0.5*(myArray.zeta[i+1][k] + myArray.zeta[i][k]) - 
						 0.25*(myArray.u[i][k] + myArray.u[i][k-1] + myArray.u[i-1][k] + myArray.u[i-1][k-1]) * 0.5*(myArray.zeta[i][k] + myArray.zeta[i-1][k])) * rdx;
			pwzeta_pz = (0.25*(myArray.w[i][k+1] + myArray.w[i-1][k+1] + myArray.w[i][k] + myArray.w[i-1][k]) * 0.5*(myArray.zeta[i][k+1] + myArray.zeta[i][k]) - 
						 0.25*(myArray.w[i][k] + myArray.w[i-1][k] + myArray.w[i][k-1] + myArray.w[i-1][k-1]) * 0.5*(myArray.zeta[i][k] + myArray.zeta[i][k-1])) * rdz;

			// Advection test
			#if defined(ADVECTIONU)
				g_tb_pth_px = 0.;
			#elif defined(ADVECTIONW)
				g_tb_pth_px = 0.;
			#elif defined(NoBouyance)
				g_tb_pth_px = 0.;
			#else
				g_tb_pth_px = gravity / myArray.tb_zeta[k] * (0.5*(myArray.th[i][k] + myArray.th[i][k-1]) - 0.5*(myArray.th[i-1][k] + myArray.th[i-1][k-1])) * rdx;
			#endif

			// Add water 
			#if defined(WATER)
				double g_pqv_px = gravity * (0.5*(myArray.qv[i][k] + myArray.qv[i][k-1]) - 0.5*(myArray.qv[i-1][k] + myArray.qv[i-1][k-1])) * rdx;
				double gpqc_px = gravity * (0.5*(myArray.qc[i][k] + myArray.qc[i][k-1]) - 0.5*(myArray.qc[i-1][k] + myArray.qc[i-1][k-1])) * rdx;
				double gpqr_px = gravity * (0.5*(myArray.qr[i][k] + myArray.qr[i][k-1]) - 0.5*(myArray.qr[i-1][k] + myArray.qr[i-1][k-1])) * rdx;
				myArray.zetap[i][k] = myArray.zetam[i][k] + d2t * (g_tb_pth_px - puzeta_px - pwzeta_pz + 0.61 * g_pqv_px - gpqc_px - gpqr_px);
			#else
				myArray.zetap[i][k] = myArray.zetam[i][k] + d2t * (g_tb_pth_px - puzeta_px - pwzeta_pz);
			#endif

			// Add diffusion
			#ifdef DIFFUSION
				myArray.zetap[i][k] += d2t * Kx * rdx2 * (myArray.zetam[i+1][k] - 2. * myArray.zetam[i][k] + myArray.zetam[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.zetam[i][k+1] - 2. * myArray.zetam[i][k] + myArray.zetam[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		myArray.BoundaryProcessDouble(myArray.zetap);
	#elif defined(ADVECTIONW)
		myArray.BoundaryProcessDouble(myArray.zetap);
	#else
		myArray.BoundaryProcessZETA(myArray.zetap);
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.zeta[i][k] += TIMETS * (myArray.zetap[i][k] - 2 * myArray.zeta[i][k] + myArray.zetam[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pth_pt(vvmArray& myArray) {
	double puth_px = 0., prhowth_pz_rho = 0., wptb_pz = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puth_px = (myArray.u[i+1][k] * 0.5*(myArray.th[i+1][k] + myArray.th[i][k]) - myArray.u[i][k] * 0.5*(myArray.th[i][k] + myArray.th[i-1][k])) * rdx;
			prhowth_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.th[i][k+1] + myArray.th[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.th[i][k] + myArray.th[i][k-1])) * rdz / myArray.rhou[k];
			wptb_pz = 0.5*(myArray.w[i][k+1] + myArray.w[i][k]) * (0.5*(myArray.tb[k+1] + myArray.tb[k]) - 0.5*(myArray.tb[k] + myArray.tb[k-1])) * rdz;

			myArray.thp[i][k] = myArray.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho - wptb_pz);

			#ifdef DIFFUSION
				myArray.thp[i][k] += d2t * Kx * rdx2 * (myArray.thm[i+1][k] - 2. * myArray.thm[i][k] + myArray.thm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.thm[i][k+1] - 2. * myArray.thm[i][k] + myArray.thm[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		myArray.BoundaryProcessDouble(myArray.thp);
	#elif defined(ADVECTIONW)
		myArray.BoundaryProcessDouble(myArray.thp);
	#else
		myArray.BoundaryProcess(myArray.thp);
	#endif

	#if defined(HEATFLUX)
		for (int i = 1; i <= nx-2; i++) {
			if (i <= nx / 2 - 1) heatflux(myArray, i, 1, 1);
		}
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.th[i][k] += TIMETS * (myArray.thp[i][k] - 2 * myArray.th[i][k] + myArray.thm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::cal_w(vvmArray& myArray) {
	#if defined(ADVECTIONU)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { myArray.u[i][k] = 10.;}}
	#elif defined(ADVECTIONW)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { myArray.u[i][k] = 0.;}}
	#else
		// eigen
		// A: i = 1~NX-2, k = 2~NZ-2
		int NX = nx, NZ = nz;
		int k = 1;
		Eigen::VectorXd x((NX-2)*(NZ-3)), b((NX-2)*(NZ-3));
		Eigen::SparseMatrix<double> A((NX-2)*(NZ-3), (NX-2)*(NZ-3));
		std::vector<T> coeff;
		for (int idx = 1; idx < (NX-2)*(NZ-3); idx++) {
			// Height
			if (idx % (NX-2) == 1) k++;

			// D
			coeff.push_back(T(idx-1, idx-1, 4. -  (- 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k+1]
												   - 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k] 
												   + 2*(myArray.rhow[k+1] - 2*myArray.rhow[k] + myArray.rhow[k-1]) / (myArray.rhou[k+1] + myArray.rhou[k]))));

			// left/right: -1
			if ((idx-1) % (NX-2) != 0) coeff.push_back(T(idx-1, idx-2, -1.));
			if (idx % (NX-2) != 0) coeff.push_back(T(idx-1, idx, -1.));

			// Boundary
			if ((idx-1) % (NX-2) == 0) {
				coeff.push_back(T(idx-1, idx-1+(NX-3), -1.));
				coeff.push_back(T(idx-1+(NX-3), idx-1, -1.));
			}
		}
		coeff.push_back(T((NX-2)*(NZ-3)-1, (NX-2)*(NZ-3)-1, 4. - (- 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k+1]
															      - 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k] 
															      + 2*(myArray.rhow[k+1] - 2*myArray.rhow[k] + myArray.rhow[k-1]) / (myArray.rhou[k+1] + myArray.rhou[k]))));
		coeff.push_back(T((NX-2)*(NZ-3)-1, (NX-2)*(NZ-3)-1-1, -1.)); // left

		k = 1;
		for (int idx = 1; idx <= (NX-2)*(NZ-4); idx++) {
			// Height
			if (idx % (NX-2) == 1) k++;

			// E
			coeff.push_back(T(idx-1, idx+(NX-2)-1, -1. - 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k+1]));
			
			// F
			coeff.push_back(T(idx+(NX-2)-1, idx-1, -1. - 0.5*(myArray.rhou[k+1] - myArray.rhou[k]) / myArray.rhou[k]));
		}
		A.setFromTriplets(coeff.begin(), coeff.end());

		// b
		int count = 0;
		for (int k = 2; k <= NZ-2; k++) {
			for (int i = 1; i <= NX-2; i++) {
				b(count) = -(myArray.zetap[i+1][k] - myArray.zetap[i][k]) * dx;
				count++;
			}
		}

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		x = solver.solve(b);

		int cnt = 0;
		for (int k = 2; k <= nz-2; k++) {
			for (int i = 1; i <= nx-2; i++) {
				myArray.w[i][k] = x[cnt];
				cnt++;
			}
		}
		myArray.BoundaryProcessZETA(myArray.w);
	#endif
	return;
}

void Iteration::cal_u(vvmArray& myArray) {
	#if defined(ADVECTIONU)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { myArray.w[i][k] = 0.;}}
	#elif defined(ADVECTIONW)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { myArray.w[i][k] = 10.;}}
	#else
		const double ubar = 0.;
		// G
		int NX = nx;
		Eigen::VectorXd y(NX-2), h(NX-2);
		Eigen::SparseMatrix<double> G(NX-2, NX-2);
		std::vector<T> coeff_xi;

		for (int k = 1; k < NX-2; k++) {
			// D
			coeff_xi.push_back(T(k-1, k-1, 2.));
			coeff_xi.push_back(T(k, k-1, -1.));
			coeff_xi.push_back(T(k-1, k, -1.));
		}
		// Boundary
		coeff_xi.push_back(T(0, NX-3, -1.));
		coeff_xi.push_back(T(NX-3, 0, -1.));
		coeff_xi.push_back(T(NX-3, NX-3, -2.));
		G.setFromTriplets(coeff_xi.begin(), coeff_xi.end());

		// h
		for (int i = 1; i <= NX-2; i++) {
			h(i-1) = (-myArray.rhow[nz-2] * myArray.w[i][nz-2]) / myArray.rhou[nz-2] * dx;
		}

		// solve
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solve_xi;
		y = solve_xi.compute(G).solve(h);
		// xi
		for (int i = 1; i <= nx-2; i++) {
			myArray.xi[i] = y[i-1];
		}
		myArray.xi[0] = myArray.xi[nx-2];
		myArray.xi[nx-1] = myArray.xi[1];

		// uxi
		for (int i = 1; i <= nx-2; i++) {
			myArray.uxi[i] = (myArray.xi[i] - myArray.xi[i-1]) * rdx;
		}
		myArray.uxi[0] = myArray.uxi[nx-2];
		myArray.uxi[nx-1] = myArray.uxi[1];

		// u_top
		for (int i = 1; i <= nx-2; i++) {
			#if defined(SHEAR)
				myArray.u[i][nz-2] = myArray.uxi[i] + ubar + 10.;
			#else
				myArray.u[i][nz-2] = myArray.uxi[i] + ubar;
			#endif
		}
		myArray.u[0][nz-2] = myArray.u[nx-2][nz-2];
		myArray.u[nx-1][nz-2] = myArray.u[1][nz-2];

		// u
		double uper = 0., now = 0.;
		for (int i = 1; i <= nx-2; i++) {
			double area = 0.;
			for (int k = nz-3; k >= 1; k--) {
				uper = (0.5*(myArray.w[i][k+2] + myArray.w[i][k+1]) - 0.5*(myArray.w[i-1][k+2] + myArray.w[i-1][k+1])) * rdx - 0.5*(myArray.zetap[i][k+2] + myArray.zetap[i][k+1]);
				now = (0.5*(myArray.w[i][k+1] + myArray.w[i][k]) - 0.5*(myArray.w[i-1][k+1] + myArray.w[i-1][k])) * rdx - 0.5*(myArray.zetap[i][k+1] + myArray.zetap[i][k]);
				area += (uper + now) * 0.5 * (-dz);
				myArray.u[i][k] = area + myArray.u[i][nz-2];
			}
		}
		myArray.BoundaryProcess(myArray.u);
	#endif
	return;
}

void Iteration::pqv_pt(vvmArray & myArray) {
	double puqv_px = 0., prhowqv_pz_rho = 0., wpqvb_pz = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqv_px = (myArray.u[i+1][k] * 0.5*(myArray.qv[i+1][k] + myArray.qv[i][k]) - myArray.u[i][k] * 0.5*(myArray.qv[i][k] + myArray.qv[i-1][k])) * rdx;
			prhowqv_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.qv[i][k+1] + myArray.qv[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.qv[i][k] + myArray.qv[i][k-1])) * rdz / myArray.rhou[k];
			wpqvb_pz = 0.5*(myArray.w[i][k+1] + myArray.w[i][k]) * (0.5*(myArray.qvb[k+1] + myArray.qvb[k]) - 0.5*(myArray.qvb[k] + myArray.qvb[k-1])) * rdz;
			myArray.qvp[i][k] = myArray.qvm[i][k] + d2t * (-puqv_px - prhowqv_pz_rho - wpqvb_pz);

			// diffusion
			#ifdef DIFFUSION
				myArray.qvp[i][k] += d2t * Kx * rdx2 * (myArray.qvm[i+1][k] - 2. * myArray.qvm[i][k] + myArray.qvm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qvm[i][k+1] - 2. * myArray.qvm[i][k] + myArray.qvm[i][k-1]);
			#endif

			// negative qv process: (source)
			if (myArray.qvp[i][k] + myArray.qvb[k] < 0.) myArray.qvp[i][k] = -myArray.qvb[k];
		}
	}
	myArray.BoundaryProcess(myArray.qvp);

	// saturation process: if cloudless.
	#ifdef CLOUDLESS
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0);   // coefficient
				double pth = myArray.thp[i][k] + myArray.tb[k];
				double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));
				if (myArray.qvp[i][k] + myArray.qvb[k] > qvs) {
					myArray.qvp[i][k] = qvs - myArray.qvb[k];
				}
			}
		}
		myArray.BoundaryProcess(myArray.qvp);
	#endif


	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qv[i][k] += TIMETS * (myArray.qvp[i][k] - 2 * myArray.qv[i][k] + myArray.qvm[i][k]);
			}
		}
	#endif

	return;
}

void Iteration::pqc_pt(vvmArray & myArray) {
	double puqc_px = 0., prhowqc_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqc_px = (myArray.u[i+1][k] * 0.5*(myArray.qc[i+1][k] + myArray.qc[i][k]) - myArray.u[i][k] * 0.5*(myArray.qc[i][k] + myArray.qc[i-1][k])) * rdx;
			prhowqc_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.qc[i][k+1] + myArray.qc[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.qc[i][k] + myArray.qc[i][k-1])) * rdz / myArray.rhou[k];

			myArray.qcp[i][k] = myArray.qcm[i][k] + d2t * (-puqc_px - prhowqc_pz_rho);
			// negative qc process
			if (myArray.qcp[i][k] < 0.) myArray.qcp[i][k] = 0.;

			// saturation process: sink and source (qv <--> qc)
			condensation(myArray, i, k);

			#ifdef DIFFUSION
				myArray.qcp[i][k] += d2t * Kx * rdx2 * (myArray.qcm[i+1][k] - 2. * myArray.qcm[i][k] + myArray.qcm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qcm[i][k+1] - 2. * myArray.qcm[i][k] + myArray.qcm[i][k-1]);
			#endif
		}
	}
	myArray.BoundaryProcess(myArray.qcp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qc[i][k] += TIMETS * (myArray.qcp[i][k] - 2 * myArray.qc[i][k] + myArray.qcm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pqr_pt(vvmArray & myArray) {
	double puqr_px = 0., prhowVTqr_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqr_px = (myArray.u[i+1][k] * 0.5*(myArray.qr[i+1][k] + myArray.qr[i][k]) - myArray.u[i][k] * 0.5*(myArray.qr[i][k] + myArray.qr[i-1][k])) * rdx;
			// TODO: VT
			double VT = 6.;
			prhowVTqr_pz_rho = (myArray.rhow[k+1] * (myArray.w[i][k+1] - VT) * 0.5*(myArray.qr[i][k+1] + myArray.qr[i][k]) - 
							    myArray.rhow[k] * (myArray.w[i][k] - VT) * 0.5*(myArray.qr[i][k] + myArray.qr[i][k-1])) * rdz / myArray.rhou[k];
			myArray.qrp[i][k] = myArray.qrm[i][k] + d2t * (-puqr_px - prhowVTqr_pz_rho);
			// negative qr process
			if (myArray.qrp[i][k] < 0.) myArray.qrp[i][k] = 0.;
			
			#ifdef DIFFUSION
				myArray.qrp[i][k] += d2t * Kx * rdx2 * (myArray.qrm[i+1][k] - 2. * myArray.qrm[i][k] + myArray.qrm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qrm[i][k+1] - 2. * myArray.qrm[i][k] + myArray.qrm[i][k-1]);
			#endif

			autoconversion(myArray, i, k);
			accretion(myArray, i, k);
			evaporation(myArray, i, k);
		}
	}
	myArray.BoundaryProcess(myArray.qrp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qr[i][k] += TIMETS * (myArray.qrp[i][k] - 2 * myArray.qr[i][k] + myArray.qrm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::condensation(vvmArray & myArray, int i, int k) {
	double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0); 	// coefficient
	double pth = myArray.thp[i][k] + myArray.tb[k];
	double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));
	double phi = qvs * (17.27 * 237. * Lv) / (C_p * pow(pth * myArray.pib[k] - 36., 2));

	double C = (myArray.qvp[i][k] + myArray.qvb[k] - qvs) / (1 + phi); 

	// C should less than qc
	if (fabs(C) > myArray.qcp[i][k] && C < 0) C = -myArray.qcp[i][k];
	
	myArray.qvp[i][k] = myArray.qvp[i][k] - C;
	myArray.qcp[i][k] = myArray.qcp[i][k] + C;
	myArray.thp[i][k] = myArray.thp[i][k] + Lv / (C_p * myArray.pib[k]) * C;
}

// autoconversion of qc to qr
void Iteration::autoconversion(vvmArray & myArray, int i, int k) {
	double autort = 0.001, autotr = 0.001; // autocon rate [1/sec], autocon threshold [kg/kg]
	double qcplus = std::max(0., myArray.qcp[i][k]);
	double ar = autort * (qcplus - autotr);
	ar = std::max(0., ar);
	double arcrdt = std::min(ar * d2t, qcplus);
	myArray.qcp[i][k] = myArray.qcp[i][k] - arcrdt;
	myArray.qrp[i][k] = myArray.qrp[i][k] + arcrdt;
	return;
}

// accretion of qc by qr
void Iteration::accretion(vvmArray & myArray, int i, int k) {
	double accrrt = 2.2; // accretion rate [1/sec]
	double qcplus = std::max(0., myArray.qcp[i][k]);
	double qrplus = std::max(0., myArray.qrp[i][k]);

	double cr = myArray.rhou[k] * accrrt * qcplus * pow(qrplus, 0.875);
	double arcrdt = std::min(cr * d2t, qcplus);

	myArray.qcp[i][k] = myArray.qcp[i][k] - arcrdt;
	myArray.qrp[i][k] = myArray.qrp[i][k] + arcrdt;
	return;
}

// evaporation of rain water
void Iteration::evaporation(vvmArray & myArray, int i, int k) {
	double qrplus = std::max(0., myArray.qrp[i][k]);
	double qvplus = std::max(0., myArray.qvp[i][k] + myArray.qvb[k]);

	double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0);	 // coefficient
	double pth = myArray.thp[i][k] + myArray.tb[k];
	double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));	// Tetens equation

	double coef = 1.6 + 30.39 * pow((myArray.rhou[k] * qrplus), 0.2046);	// ventilation coef.
	double deficit = std::max((1. - qvplus / qvs), 0.);							// saturation dificit (RH < 100%)

	double er = coef * deficit * (pow(myArray.rhou[k] * qrplus, 0.525)) / 
				((2.03e4 + 9.584e6 / (myArray.pb[k] * qvs)) * myArray.rhou[k]);
	double erdt = std::min(qrplus, std::max(0., er * d2t));

	myArray.qrp[i][k] = myArray.qrp[i][k] - erdt;
	myArray.qvp[i][k] = myArray.qvp[i][k] + erdt;
	myArray.thp[i][k] = myArray.thp[i][k] - Lv * erdt / (C_p * myArray.pib[k]);
	return;
}

void Iteration::heatflux(vvmArray & myArray, int i, int k, int ishflux) {
	if (ishflux == 1 && k == 1) {
		double cdh = 7e-3;
		double tground = 303.;
		double tdif = std::max(tground - (myArray.thm[i][k] + myArray.tb[k]), 0.);
		double avgu = 0.5 * abs(myArray.u[i+1][k] + myArray.u[i][k]);
		avgu = std::max(avgu, 2.);
		double wnetc = 2. * sqrt(tdif);      // "convective velocity" adjustment
		double vel = sqrt(pow(avgu, 2) + pow(wnetc, 2));
		myArray.thp[i][k] = myArray.thp[i][k] + d2t * cdh * vel * myArray.addflx[i] * tdif * rdz;
	}
}

void Iteration::LeapFrog(vvmArray & myArray) {
	int n = 0;
	// double timenow = 0.;
	double temp = TIMEEND / dt;
	int nmax = (int) temp;
	while (n < nmax) {
		std::cout << n << std::endl;
		// output
		Output::create_all_directory();
		if (n % OUTPUTSTEP == 0) {
			#if defined(OUTPUTFILEMODE)
				Output::output_zeta(n, myArray);
				Output::output_th(n, myArray);
				Output::output_u(n, myArray);
				Output::output_w(n, myArray);
				Output::output_qv(n, myArray);
				Output::output_qc(n, myArray);
				Output::output_qr(n, myArray);
			#endif
			#if defined(OUTPUTGRAPHMODE)
				Plot::plot_zeta(n, myArray);
			#endif
			#if defined(OUTPUTNC)
				Output::output_nc(n, myArray);
			#endif
		}
		n++;
		// timenow = n * dt;

		// calculate
		pzeta_pt(myArray);
		pth_pt(myArray);
		cal_w(myArray);
		cal_u(myArray);
		pqv_pt(myArray);
		pqc_pt(myArray);
		pqr_pt(myArray);

		// next step
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.zetam[i][k] = myArray.zeta[i][k];
				myArray.zeta[i][k] = myArray.zetap[i][k];

				myArray.thm[i][k] = myArray.th[i][k];
				myArray.th[i][k] = myArray.thp[i][k];

				myArray.qvm[i][k] = myArray.qv[i][k];
				myArray.qv[i][k] = myArray.qvp[i][k];

				myArray.qcm[i][k] = myArray.qc[i][k];
				myArray.qc[i][k] = myArray.qcp[i][k];

				myArray.qrm[i][k] = myArray.qr[i][k];
				myArray.qr[i][k] = myArray.qrp[i][k];
			}
		}

	}
	return;

}


