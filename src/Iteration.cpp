#include "Iteration.hpp"

void Iteration::pzeta_pt(vvmArray &model) {
	double puzeta_px = 0., pwzeta_pz = 0., g_tb_pth_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puzeta_px = (0.25*(model.u[i+1][k] + model.u[i+1][k-1] + model.u[i][k] + model.u[i][k-1]) * 0.5*(model.zeta[i+1][k] + model.zeta[i][k]) - 
						 0.25*(model.u[i][k] + model.u[i][k-1] + model.u[i-1][k] + model.u[i-1][k-1]) * 0.5*(model.zeta[i][k] + model.zeta[i-1][k])) * rdx;
			pwzeta_pz = (0.25*(model.w[i][k+1] + model.w[i-1][k+1] + model.w[i][k] + model.w[i-1][k]) * 0.5*(model.zeta[i][k+1] + model.zeta[i][k]) - 
						 0.25*(model.w[i][k] + model.w[i-1][k] + model.w[i][k-1] + model.w[i-1][k-1]) * 0.5*(model.zeta[i][k] + model.zeta[i][k-1])) * rdz;

			// Advection test
			#if defined(ADVECTIONU)
				g_tb_pth_px = 0.;
			#elif defined(ADVECTIONW)
				g_tb_pth_px = 0.;
			#elif defined(NoBouyance)
				g_tb_pth_px = 0.;
			#else
				g_tb_pth_px = gravity / model.tb_zeta[k] * (0.5*(model.th[i][k] + model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * rdx;
			#endif

			// Add water 
			#if defined(WATER)
				double g_pqv_px = gravity * (0.5*(model.qv[i][k] + model.qv[i][k-1]) - 0.5*(model.qv[i-1][k] + model.qv[i-1][k-1])) * rdx;
				double gpqc_px = gravity * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * rdx;
				double gpqr_px = gravity * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * rdx;
				model.zetap[i][k] = model.zetam[i][k] + d2t * (g_tb_pth_px - puzeta_px - pwzeta_pz + 0.61 * g_pqv_px - gpqc_px - gpqr_px);
			#else
				model.zetap[i][k] = model.zetam[i][k] + d2t * (g_tb_pth_px - puzeta_px - pwzeta_pz);
			#endif

			// Add diffusion
			#ifdef DIFFUSION
				model.zetap[i][k] += d2t * Kx * rdx2 * (model.zetam[i+1][k] - 2. * model.zetam[i][k] + model.zetam[i-1][k]) + 
									 d2t * Kz * rdz2 * (model.zetam[i][k+1] - 2. * model.zetam[i][k] + model.zetam[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		model.BoundaryProcessDouble(model.zetap);
	#elif defined(ADVECTIONW)
		model.BoundaryProcessDouble(model.zetap);
	#else
		model.BoundaryProcessZETA(model.zetap);
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.zeta[i][k] += TIMETS * (model.zetap[i][k] - 2 * model.zeta[i][k] + model.zetam[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pth_pt(vvmArray &model) {
    double puth_px = 0., prhowth_pz_rho = 0., forcing = 0.;
	#if defined(LINEARIZEDTH)
		double wptb_pz = 0.;
	#endif
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puth_px = (model.u[i+1][k] * 0.5*(model.th[i+1][k] + model.th[i][k]) - model.u[i][k] * 0.5*(model.th[i][k] + model.th[i-1][k])) * rdx;
			prhowth_pz_rho = (model.rhow[k+1] * model.w[i][k+1] * 0.5*(model.th[i][k+1] + model.th[i][k]) - 
							  model.rhow[k] * model.w[i][k] * 0.5*(model.th[i][k] + model.th[i][k-1])) * rdz / model.rhou[k];

			#if defined(TROPICALFORCING)
                forcing = model.Q1LS[k];
            #else
                forcing = 0.;
            #endif

			#if defined(LINEARIZEDTH)
				wptb_pz = 0.5*(model.w[i][k+1] + model.w[i][k]) * (0.5*(model.tb[k+1] + model.tb[k]) - 0.5*(model.tb[k] + model.tb[k-1])) * rdz;
				#if defined(TROPICALFORCING)
					if (model.status_for_adding_forcing == true) {
						model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho - wptb_pz + forcing + model.init_th_forcing[i][k]);
					}
					else {
						model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho - wptb_pz + forcing);
					}

					#if defined(RADIATIONCOOLING) 
						model.thp[i][k] += d2t * (-2 / 86400);
					#endif
				#else
					model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho - wptb_pz + forcing);
				#endif
            #else
                #if defined(TROPICALFORCING)
					if (model.status_for_adding_forcing == true) {
						model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho + forcing + model.init_th_forcing[i][k]);
					}
					else {
						model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho + forcing);
					}

					#if defined(RADIATIONCOOLING) 
						model.thp[i][k] += d2t * (-2 / 86400);
					#endif
				#else
					model.thp[i][k] = model.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho + forcing);
				#endif
            #endif

			#ifdef DIFFUSION
				model.thp[i][k] += d2t * Kx * rdx2 * (model.thm[i+1][k] - 2. * model.thm[i][k] + model.thm[i-1][k]) + 
									 d2t * Kz * rdz2 * (model.thm[i][k+1] - 2. * model.thm[i][k] + model.thm[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		model.BoundaryProcessDouble(model.thp);
	#elif defined(ADVECTIONW)
		model.BoundaryProcessDouble(model.thp);
	#else
		model.BoundaryProcess(model.thp);
	#endif

	#if defined(HEATFLUX)
		for (int i = 1; i <= nx-2; i++) {
			if (i <= nx / 2 - 1) heatflux(model, i, 1, 1);
		}
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.th[i][k] += TIMETS * (model.thp[i][k] - 2 * model.th[i][k] + model.thm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::cal_w(vvmArray &model) {
	#if defined(ADVECTIONU)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { model.u[i][k] = 10.;}}
	#elif defined(ADVECTIONW)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { model.u[i][k] = 0.;}}
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
			coeff.push_back(T(idx-1, idx-1, 4. -  (- 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k+1]
												   - 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k] 
												   + 2*(model.rhow[k+1] - 2*model.rhow[k] + model.rhow[k-1]) / (model.rhou[k+1] + model.rhou[k]))));

			// left/right: -1
			if ((idx-1) % (NX-2) != 0) coeff.push_back(T(idx-1, idx-2, -1.));
			if (idx % (NX-2) != 0) coeff.push_back(T(idx-1, idx, -1.));

			// Boundary
			if ((idx-1) % (NX-2) == 0) {
				coeff.push_back(T(idx-1, idx-1+(NX-3), -1.));
				coeff.push_back(T(idx-1+(NX-3), idx-1, -1.));
			}
		}
		coeff.push_back(T((NX-2)*(NZ-3)-1, (NX-2)*(NZ-3)-1, 4. - (- 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k+1]
															      - 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k] 
															      + 2*(model.rhow[k+1] - 2*model.rhow[k] + model.rhow[k-1]) / (model.rhou[k+1] + model.rhou[k]))));
		coeff.push_back(T((NX-2)*(NZ-3)-1, (NX-2)*(NZ-3)-1-1, -1.)); // left

		k = 1;
		for (int idx = 1; idx <= (NX-2)*(NZ-4); idx++) {
			// Height
			if (idx % (NX-2) == 1) k++;

			// E
			coeff.push_back(T(idx-1, idx+(NX-2)-1, -1. - 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k+1]));
			
			// F
			coeff.push_back(T(idx+(NX-2)-1, idx-1, -1. - 0.5*(model.rhou[k+1] - model.rhou[k]) / model.rhou[k]));
		}
		A.setFromTriplets(coeff.begin(), coeff.end());

		// b
		int count = 0;
		for (int k = 2; k <= NZ-2; k++) {
			for (int i = 1; i <= NX-2; i++) {
				b(count) = -(model.zetap[i+1][k] - model.zetap[i][k]) * dx;
				count++;
			}
		}

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		x = solver.solve(b);

		int cnt = 0;
		for (int k = 2; k <= nz-2; k++) {
			for (int i = 1; i <= nx-2; i++) {
				model.w[i][k] = x[cnt];
				cnt++;
			}
		}
		model.BoundaryProcessZETA(model.w);
	#endif
	return;
}

void Iteration::cal_u(vvmArray &model) {
	#if defined(ADVECTIONU)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { model.w[i][k] = 0.;}}
	#elif defined(ADVECTIONW)
		for (int i = 0; i <= nx-1; i++) { for (int k = 0; k <= nz-1; k++) { model.w[i][k] = 10.;}}
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
			h(i-1) = (-model.rhow[nz-2] * model.w[i][nz-2]) / model.rhou[nz-2] * dx;
		}

		// solve
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solve_xi;
		y = solve_xi.compute(G).solve(h);
		// xi
		for (int i = 1; i <= nx-2; i++) {
			model.xi[i] = y[i-1];
		}
		model.xi[0] = model.xi[nx-2];
		model.xi[nx-1] = model.xi[1];

		// uxi
		for (int i = 1; i <= nx-2; i++) {
			model.uxi[i] = (model.xi[i] - model.xi[i-1]) * rdx;
		}
		model.uxi[0] = model.uxi[nx-2];
		model.uxi[nx-1] = model.uxi[1];

		// u_top
		for (int i = 1; i <= nx-2; i++) {
			#if defined(SHEAR)
			// You should modify the last number which means the fixed u at the top layer
				model.u[i][nz-2] = model.uxi[i] + ubar + 10.;
			#else
				model.u[i][nz-2] = model.uxi[i] + ubar;
			#endif
		}
		model.u[0][nz-2] = model.u[nx-2][nz-2];
		model.u[nx-1][nz-2] = model.u[1][nz-2];

		// u
		double uper = 0., now = 0.;
		for (int i = 1; i <= nx-2; i++) {
			double area = 0.;
			for (int k = nz-3; k >= 1; k--) {
				uper = (0.5*(model.w[i][k+2] + model.w[i][k+1]) - 0.5*(model.w[i-1][k+2] + model.w[i-1][k+1])) * rdx - 0.5*(model.zetap[i][k+2] + model.zetap[i][k+1]);
				now = (0.5*(model.w[i][k+1] + model.w[i][k]) - 0.5*(model.w[i-1][k+1] + model.w[i-1][k])) * rdx - 0.5*(model.zetap[i][k+1] + model.zetap[i][k]);
				area += (uper + now) * 0.5 * (-dz);
				model.u[i][k] = area + model.u[i][nz-2];
			}
		}
		model.BoundaryProcess(model.u);
	#endif
	return;
}

void Iteration::pqv_pt(vvmArray &model) {
	double puqv_px = 0., prhowqv_pz_rho = 0., forcing = 0.;
	#if defined(LINEARIZEDQV)
		double wpqvb_pz = 0.;
	#endif
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqv_px = (model.u[i+1][k] * 0.5*(model.qv[i+1][k] + model.qv[i][k]) - model.u[i][k] * 0.5*(model.qv[i][k] + model.qv[i-1][k])) * rdx;
			prhowqv_pz_rho = (model.rhow[k+1] * model.w[i][k+1] * 0.5*(model.qv[i][k+1] + model.qv[i][k]) - 
							  model.rhow[k] * model.w[i][k] * 0.5*(model.qv[i][k] + model.qv[i][k-1])) * rdz / model.rhou[k];
			
			#if defined(TROPICALFORCING)
                forcing = model.Q2LS[k];
            #else
                forcing = 0.;
            #endif

            #if defined(LINEARIZEDQV)
				wpqvb_pz = 0.5*(model.w[i][k+1] + model.w[i][k]) * (0.5*(model.qvb[k+1] + model.qvb[k]) - 0.5*(model.qvb[k] + model.qvb[k-1])) * rdz;
                model.qvp[i][k] = model.qvm[i][k] + d2t * (-puqv_px - prhowqv_pz_rho - wpqvb_pz + forcing);
            #else
                model.qvp[i][k] = model.qvm[i][k] + d2t * (-puqv_px - prhowqv_pz_rho + forcing);
            #endif

			// diffusion
			#ifdef DIFFUSION
				model.qvp[i][k] += d2t * Kx * rdx2 * (model.qvm[i+1][k] - 2. * model.qvm[i][k] + model.qvm[i-1][k]) + 
								   d2t * Kz * rdz2 * (model.qvm[i][k+1] - 2. * model.qvm[i][k] + model.qvm[i][k-1]);
			#endif

			// negative qv process: (source)
			#if defined(LINEARIZEDQV)
				if (model.qvp[i][k] + model.qvb[k] < 0.) model.qvp[i][k] = -model.qvb[k];
			#else
				if (model.qvp[i][k] < 0.) model.qvp[i][k] = 0.;
			#endif
		}
	}
	model.BoundaryProcess(model.qvp);

	// saturation process: if cloudless.
	#ifdef CLOUDLESS
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				double pc = 380. / (pow(model.pib[k], C_p / Rd) * P0);   // coefficient
				double pth = model.thp[i][k] + model.tb[k];
				double qvs = pc * exp(17.27 * (model.pib[k] * pth - 273.) / (model.pib[k] * pth - 36.));
				if (model.qvp[i][k] + model.qvb[k] > qvs) {
					model.qvp[i][k] = qvs - model.qvb[k];
				}
			}
		}
		model.BoundaryProcess(model.qvp);
	#endif


	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.qv[i][k] += TIMETS * (model.qvp[i][k] - 2 * model.qv[i][k] + model.qvm[i][k]);
			}
		}
	#endif

	return;
}

void Iteration::pqc_pt(vvmArray &model) {
	double puqc_px = 0., prhowqc_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqc_px = (model.u[i+1][k] * 0.5*(model.qc[i+1][k] + model.qc[i][k]) - model.u[i][k] * 0.5*(model.qc[i][k] + model.qc[i-1][k])) * rdx;
			prhowqc_pz_rho = (model.rhow[k+1] * model.w[i][k+1] * 0.5*(model.qc[i][k+1] + model.qc[i][k]) - 
							  model.rhow[k] * model.w[i][k] * 0.5*(model.qc[i][k] + model.qc[i][k-1])) * rdz / model.rhou[k];

			model.qcp[i][k] = model.qcm[i][k] + d2t * (-puqc_px - prhowqc_pz_rho);
			// negative qc process
			if (model.qcp[i][k] < 0.) model.qcp[i][k] = 0.;

			// saturation process: sink and source (qv <--> qc)
			condensation(model, i, k);

			#ifdef DIFFUSION
				model.qcp[i][k] += d2t * Kx * rdx2 * (model.qcm[i+1][k] - 2. * model.qcm[i][k] + model.qcm[i-1][k]) + 
									 d2t * Kz * rdz2 * (model.qcm[i][k+1] - 2. * model.qcm[i][k] + model.qcm[i][k-1]);
			#endif
		}
	}
	model.BoundaryProcess(model.qcp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.qc[i][k] += TIMETS * (model.qcp[i][k] - 2 * model.qc[i][k] + model.qcm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pqr_pt(vvmArray &model) {
	double puqr_px = 0., prhowVTqr_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = nz-2; k >= 1; k--) {
			if (k == 2) model.qr[i][1] = model.qr[i][2];
			if (k == 1) {
				model.qr[i][1] = 0.;
				model.qr[i][0] = 0.;
			}

			puqr_px = (model.u[i+1][k] * 0.5*(model.qr[i+1][k] + model.qr[i][k]) - model.u[i][k] * 0.5*(model.qr[i][k] + model.qr[i-1][k])) * rdx;
			// TODO: VT
			double VT = 6.;
			prhowVTqr_pz_rho = (model.rhow[k+1] * (model.w[i][k+1] - VT) * 0.5*(model.qr[i][k+1] + model.qr[i][k]) - 
							    model.rhow[k] * (model.w[i][k] - VT) * 0.5*(model.qr[i][k] + model.qr[i][k-1])) * rdz / model.rhou[k];
			
			if (k == 1) model.qrp[i][k] = d2t * (-prhowVTqr_pz_rho);
			else model.qrp[i][k] = model.qrm[i][k] + d2t * (-puqr_px - prhowVTqr_pz_rho);

			// negative qr process
			if (model.qrp[i][k] < 0.) model.qrp[i][k] = 0.;
			
			#ifdef DIFFUSION
				model.qrp[i][k] += d2t * Kx * rdx2 * (model.qrm[i+1][k] - 2. * model.qrm[i][k] + model.qrm[i-1][k]) + 
								   d2t * Kz * rdz2 * (model.qrm[i][k+1] - 2. * model.qrm[i][k] + model.qrm[i][k-1]);
			#endif

			if (k >= 2) {
				autoconversion(model, i, k);
				accretion(model, i, k);
				evaporation(model, i, k);
			}

			model.qrAcc[i] += model.qrp[i][1];
		}
	}
	model.BoundaryProcess(model.qrp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.qr[i][k] += TIMETS * (model.qrp[i][k] - 2 * model.qr[i][k] + model.qrm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::condensation(vvmArray &model, int i, int k) {
	double pc = 380. / (pow(model.pib[k], C_p / Rd) * P0); 	// coefficient
	#if defined(LINEARIZEDTH)
		double pth = model.thp[i][k] + model.tb[k];
	#else
		double pth = model.thp[i][k];
	#endif
	double qvs = pc * exp(17.27 * (model.pib[k] * pth - 273.) / (model.pib[k] * pth - 36.));
	double phi = qvs * (17.27 * 237. * Lv) / (C_p * pow(pth * model.pib[k] - 36., 2));

	#if defined(LINEARIZEDQV)
		double C = (model.qvp[i][k] + model.qvb[k] - qvs) / (1 + phi); 
	#else
		double C = (model.qvp[i][k] - qvs) / (1 + phi); 
	#endif

	// C should less than qc (C can be sink for qc and source for qv, so it should not excess qc)
	if (fabs(C) > model.qcp[i][k] && C < 0) C = -model.qcp[i][k];
	
	model.qvp[i][k] = model.qvp[i][k] - C;
	model.qcp[i][k] = model.qcp[i][k] + C;
	model.thp[i][k] = model.thp[i][k] + Lv / (C_p * model.pib[k]) * C;
}

// autoconversion of qc to qr
void Iteration::autoconversion(vvmArray & model, int i, int k) {
	double autort = 0.001, autotr = 0.001; // autocon rate [1/sec], autocon threshold [kg/kg]
	double qcplus = std::max(0., model.qcp[i][k]);
	double ar = autort * (qcplus - autotr);
	ar = std::max(0., ar);
	double arcrdt = std::min(ar * d2t, qcplus);
	model.qcp[i][k] = model.qcp[i][k] - arcrdt;
	model.qrp[i][k] = model.qrp[i][k] + arcrdt;
	return;
}

// accretion of qc by qr
void Iteration::accretion(vvmArray &model, int i, int k) {
	double accrrt = 2.2; // accretion rate [1/sec]
	double qcplus = std::max(0., model.qcp[i][k]);
	double qrplus = std::max(0., model.qrp[i][k]);

	double cr = model.rhou[k] * accrrt * qcplus * pow(qrplus, 0.875);
	double arcrdt = std::min(cr * d2t, qcplus);

	model.qcp[i][k] = model.qcp[i][k] - arcrdt;
	model.qrp[i][k] = model.qrp[i][k] + arcrdt;
	return;
}

// evaporation of rain water
void Iteration::evaporation(vvmArray &model, int i, int k) {
	double qrplus = std::max(0., model.qrp[i][k]);
	double qvplus = std::max(0., model.qvp[i][k]);

	double pc = 380. / (pow(model.pib[k], C_p / Rd) * P0);	 // coefficient
	#if defined(LINEARIZEDTH)
		double pth = model.thp[i][k] + model.tb[k];
	#else
		double pth = model.thp[i][k];
	#endif
	double qvs = pc * exp(17.27 * (model.pib[k] * pth - 273.) / (model.pib[k] * pth - 36.));	// Tetens equation

	double coef = 1.6 + 30.39 * pow((model.rhou[k] * qrplus), 0.2046);	// ventilation coef.
	double deficit = std::max((1. - qvplus / qvs), 0.);					// saturation dificit (RH < 100%)

	double er = coef * deficit * (pow(model.rhou[k] * qrplus, 0.525)) / 
				((2.03e4 + 9.584e6 / (model.pb[k] * qvs)) * model.rhou[k]);
	double erdt = std::min(qrplus, std::max(0., er * d2t));

	model.qrp[i][k] = model.qrp[i][k] - erdt;
	model.qvp[i][k] = model.qvp[i][k] + erdt;
	model.thp[i][k] = model.thp[i][k] - Lv * erdt / (C_p * model.pib[k]);
	return;
}

void Iteration::heatflux(vvmArray &model, int i, int k, int ishflux) {
	if (ishflux == 1 && k == 1) {
		double cdh = 7e-3;
		double tground = 303.;
		double tdif = std::max(tground - (model.thm[i][k] + model.tb[k]), 0.);
		double avgu = 0.5 * abs(model.u[i+1][k] + model.u[i][k]);
		avgu = std::max(avgu, 2.);
		double wnetc = 2. * sqrt(tdif);      // "convective velocity" adjustment
		double vel = sqrt(pow(avgu, 2) + pow(wnetc, 2));
		model.thp[i][k] = model.thp[i][k] + d2t * cdh * vel * model.addflx[i] * tdif * rdz;
	}
}

void Iteration::updateMean(vvmArray &model) {
	double tb = 0.;
	for (int k = 1; k < nz-1; k++) {
		tb = 0.;
		for (int i = 1; i < nx-1; i++) {
			tb += model.thp[i][k];
		}
		model.tb[k] = tb / (nx-2);
	}
	model.tb[0] = model.tb[1];
	model.tb[nz-1] = model.tb[nz-2];

	for (int k = 1; k < nz-1; k++) {
		model.tb_zeta[k] = 0.5 * (model.tb[k-1] + model.tb[k]);
	}
	model.tb_zeta[0] = model.tb_zeta[1];
	model.tb_zeta[nz-1] = model.tb_zeta[nz-2];
	
}

void Iteration::LeapFrog(vvmArray &model) {
	int n = 0;
	double temp = TIMEEND / dt;
	int nmax = (int) temp;
	while (n < nmax) {
		std::cout << n << std::endl;
		// output
		if (n % OUTPUTSTEP == 0) {
			#if defined(OUTPUTTXT)
				Output::output_zeta(n, model);
				Output::output_th(n, model);
				Output::output_u(n, model);
				Output::output_w(n, model);
				Output::output_qv(n, model);
				Output::output_qc(n, model);
				Output::output_qr(n, model);
			#endif
			#if defined(OUTPUTGRAPHMODE)
				Plot::plot_zeta(n, model);
			#endif
			#if defined(OUTPUTNC)
				Output::output_nc(n, model);
			#endif
		}
		n++;
		
        #if defined(TROPICALFORCING)
            if (n * dt <= ADDFORCINGTIME) model.status_for_adding_forcing = true;
            else model.status_for_adding_forcing = false;

            // Generate new random th perturbation for tropical forcing case
            if (model.status_for_adding_forcing == true) {
                Init::RandomPerturbation(model, n);
            }
        #endif

		// calculate
		pzeta_pt(model);
		pth_pt(model);
		cal_w(model);
		cal_u(model);
		#if defined(WATER)
			pqv_pt(model);
			pqc_pt(model);
			pqr_pt(model);
		#endif
		cal_w(model);
		cal_u(model);

		#ifndef LINEARIZEDTH
			updateMean(model);
		#endif 

		// next step
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				model.zetam[i][k] = model.zeta[i][k];
				model.zeta[i][k] = model.zetap[i][k];

				model.thm[i][k] = model.th[i][k];
				model.th[i][k] = model.thp[i][k];

				#if defined(WATER)
					model.qvm[i][k] = model.qv[i][k];
					model.qv[i][k] = model.qvp[i][k];

					model.qcm[i][k] = model.qc[i][k];
					model.qc[i][k] = model.qcp[i][k];

					model.qrm[i][k] = model.qr[i][k];
					model.qr[i][k] = model.qrp[i][k];
				#endif
			}
		}

	}
	return;

}


