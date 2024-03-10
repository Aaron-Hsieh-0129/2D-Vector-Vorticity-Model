#include "Iteration.hpp"
#include "Eigen/src/IterativeLinearSolvers/BiCGSTAB.h"

// Thomas method from wikipedia
void solve(double* a, double* b, double* c, double* d, int n) {
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
    return;
}

// #if defined(STREAMFUNCTION)
// void Iteration::pzeta_pt(vvmArray &model) {
// 	double g_tbrho_pth_px = 0.;
//     double jacobian = 0.;
// 	for (int i = 1; i <= model.nx-2; i++) {
// 		for (int k = 1; k <= model.nz-2; k++) {
//             jacobian = ((model.zeta[i+1][k] - model.zeta[i-1][k]) * (model.psi[i][k+1] - model.psi[i][k-1]) - 
//                         (model.zeta[i][k+1] - model.zeta[i][k-1]) * (model.psi[i+1][k] - model.psi[i-1][k])) * 0.25*model.rdx2 / model.rhow[k];

// 			// Advection test
// 			#if defined(ADVECTIONU) || defined(ADVECTIONW) || defined(NoBouyance)
// 				g_tbrho_pth_px = 0.;
// 			#else
// 				#if defined(LINEARIZEDTH)
// 					g_tbrho_pth_px = gravity / model.tb_zeta[k] / model.rhow[k] * (0.5*(model.th[i][k] + model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * rdx;
// 				// g_tbrho_pth_px = gravity / model.tb_zeta[k] * (0.5*(model.th[i][k] + model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * rdx;
// 				// g_tbrho_pth_px = gravity / 300. / model.rhow[k] * (0.5*(model.th[i][k] + model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * rdx;
// 				#else
// 					g_tbrho_pth_px = gravity / model.tb_zeta[k] / model.rhow[k] * (0.5*(model.th[i][k]+model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * model.rdx;
// 				#endif
// 			#endif

// 			// Add water 
// 			#if defined(WATER)
// 				double g_pqv_px = gravity / model.rhow[k] * (0.5*(model.qv[i][k] + model.qv[i][k-1]) - 0.5*(model.qv[i-1][k] + model.qv[i-1][k-1])) * rdx;
// 				double gpqc_px = gravity / model.rhow[k] * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * rdx;
// 				double gpqr_px = gravity / model.rhow[k] * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * rdx;
// 				model.zetap[i][k] = model.zetam[i][k] + d2t * (g_tbrho_pth_px - puzeta_px - prhowzeta_pz_rho + 0.61 * g_pqv_px - gpqc_px - gpqr_px);
// 			#else
// 				model.zetap[i][k] = model.zetam[i][k] + d2t * (g_tbrho_pth_px + jacobian);
// 			#endif

// 			// Add diffusion
// 			#ifdef DIFFUSION
// 				model.zetap[i][k] += d2t * Kx * model.rdx2 * (model.zetam[i+1][k] - 2. * model.zetam[i][k] + model.zetam[i-1][k]) + 
// 									 d2t * Kz * model.rdz2 * (model.zetam[i][k+1] - 2. * model.zetam[i][k] + model.zetam[i][k-1]);
// 			#endif
// 		}
// 	}
// 	#if defined(ADVECTIONU) || defined(ADVECTIONW)
// 		model.BoundaryProcessDouble(model.zetap);
// 	#else
// 		model.BoundaryProcessZETA(model.zetap);
// 	#endif

// 	// Time filter
// 	#ifdef TIMEFILTER
// 		for (int i = 0; i <= model.nx-1; i++) {
// 			for (int k = 0; k <= model.nz-1; k++) {
// 				model.zeta[i][k] += TIMETS * (model.zetap[i][k] - 2 * model.zeta[i][k] + model.zetam[i][k]);
// 			}
// 		}
// 	#endif
// 	return;
// }
// #else
void Iteration::pzeta_pt(vvmArray &model) {
    double upzeta_px = 0., zetapu_px = 0., wpzeta_pz = 0., zetapw_pz = 0., wzeta_rhoprho_pz = 0., g_tbrho_pth_px = 0.;
    for (int i = 1; i <= model.nx-2; i++) {
        for (int k = 2; k <= model.nz-2; k++) {
            upzeta_px = 0.5*(model.u[i][k]+model.u[i][k-1]) * (model.zeta[i+1][k]-model.zeta[i-1][k]) * model.r2dx;
            wpzeta_pz = 0.5*(model.w[i][k]+model.w[i-1][k]) * (model.zeta[i][k+1]-model.zeta[i][k-1]) * model.r2dz;

            model.zetap[i][k] = model.zetam[i][k] + d2t * (-upzeta_px - wpzeta_pz);

            #if defined(FLUXFORM)
                zetapu_px = model.zeta[i][k] * (0.5*(model.u[i+1][k]+model.u[i+1][k-1]) - 0.5*(model.u[i-1][k]+model.u[i-1][k-1])) * model.r2dx;
                zetapw_pz = model.zeta[i][k] * (0.5*(model.w[i][k+1]+model.w[i-1][k+1]) - 0.5*(model.w[i][k-1]+model.w[i-1][k-1])) * model.r2dz;
                wzeta_rhoprho_pz = 0.5*(model.w[i][k]+model.w[i-1][k]) * model.zeta[i][k] * (model.rhou[k]-model.rhou[k-1]) * model.rdz / model.rhow[k];
                model.zetap[i][k] += d2t * (-zetapu_px - zetapw_pz - wzeta_rhoprho_pz);
            #endif

           // Advection test
            #if defined(ADVECTIONU) || defined(ADVECTIONW) || defined(NoBouyance)
                g_tbrho_pth_px = 0.;
            #else
                g_tbrho_pth_px = gravity / model.tb_zeta[k] / model.rhow[k] * (0.5*(model.th[i][k]+model.th[i][k-1]) - 0.5*(model.th[i-1][k] + model.th[i-1][k-1])) * model.rdx;
            #endif
            model.zetap[i][k] += d2t * (g_tbrho_pth_px);

            #if defined(WATER)
                double gpqv_px = gravity / model.rhow[k] * (0.5*(model.qv[i][k] + model.qv[i][k-1]) - 0.5*(model.qv[i-1][k] + model.qv[i-1][k-1])) * model.rdx;
                double gpqc_px = gravity / model.rhow[k] * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * model.rdx;
                double gpqr_px = gravity / model.rhow[k] * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * model.rdx;
                model.zetap[i][k] += d2t * (0.61 * gpqv_px - gpqc_px - gpqr_px);
            #endif

            // Add diffusion
            #ifdef DIFFUSION
                model.zetap[i][k] += d2t * Kx * model.rdx2 * (model.zetam[i+1][k] - 2. * model.zetam[i][k] + model.zetam[i-1][k]) + 
                                        d2t * Kz * model.rdz2 * (model.zetam[i][k+1] - 2. * model.zetam[i][k] + model.zetam[i][k-1]);
            #endif
        }
    }
    #if defined(ADVECTIONU) || defined(ADVECTIONW)
        model.BoundaryProcessDouble(model.zetap);
    #else
        model.BoundaryProcessZETA(model.zetap);
    #endif

    // Time filter
    #ifdef TIMEFILTER
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.zeta[i][k] += TIMETS * (model.zetap[i][k] - 2 * model.zeta[i][k] + model.zetam[i][k]);
            }
        }
    #endif
    return;
}
// #endif

void Iteration::pth_pt(vvmArray &model) {
    double upth_px = 0., wpth_pz = 0.;
    #if defined(FLUXFORM)
        double thpu_px = 0., thpw_pz = 0.;
        #if defined(LINEARIZEDTH)
            double wptb_pz = 0.;
        #else
            double wth_rhoprho_pz = 0.;
        #endif
    #endif

    for (int i = 1; i <= model.nx-2; i++) {
		for (int k = 1; k <= model.nz-2; k++) {
            upth_px = 0.5*(model.u[i+1][k]+model.u[i][k]) * (model.th[i+1][k] - model.th[i-1][k]) * model.r2dx;
            wpth_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * (model.th[i][k+1] - model.th[i][k-1]) * model.r2dz;

            model.thp[i][k] = model.thm[i][k] + d2t * (-upth_px - wpth_pz);

            #if defined(FLUXFORM)
                thpu_px = model.th[i][k] * (model.u[i+1][k] - model.u[i][k]) * model.rdx;
                thpw_pz = model.th[i][k] * (model.w[i][k+1] - model.w[i][k]) * model.rdz;
                // model.thp[i][k] += d2t * (-thpu_px - thpw_pz);
                #if defined(LINEARIZEDTH)
                    wptb_pz = 0.5*(model.w[i][k+1] + model.w[i][k]) * (0.5*(model.tb[k+1] + model.tb[k]) - 0.5*(model.tb[k] + model.tb[k-1])) * rdz;
                    model.thp[i][k] += d2t * (-wptb_pz);
                #else
                    wth_rhoprho_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * model.th[i][k] / model.rhou[k] * (model.rhow[k+1] - model.rhow[k]) * model.rdz;
                    // model.thp[i][k] += d2t * (-wth_rhoprho_pz);
                #endif
            #endif

            #if defined(TROPICALFORCING)
                model.thp[i][k] += d2t * model.Q1LS[k];
                if (model.status_for_adding_forcing == true) model.thp[i][k] += d2t * model.init_th_forcing[i][k];
                
                #if defined(RADIATIONCOOLING) 
                    model.thp[i][k] += d2t * (-2. / 86400.);
                #endif
            #endif

			#ifdef DIFFUSION
				model.thp[i][k] += d2t * Kx * model.rdx2 * (model.thm[i+1][k] - 2. * model.thm[i][k] + model.thm[i-1][k]) + 
								   d2t * Kz * model.rdz2 * (model.thm[i][k+1] - 2. * model.thm[i][k] + model.thm[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU) || defined(ADVECTIONW)
		model.BoundaryProcessDouble(model.thp);
	#else
		model.BoundaryProcess(model.thp);
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= model.nx-1; i++) {
			for (int k = 0; k <= model.nz-1; k++) {
				model.th[i][k] += TIMETS * (model.thp[i][k] - 2 * model.th[i][k] + model.thm[i][k]);
			}
		}
	#endif
	return;
}

#if defined(STREAMFUNCTION)
void Iteration::calpsiuw(vvmArray &model) {
    // eigen solver
    Eigen::VectorXd x((model.nx-2)*(model.nz-3)), b((model.nx-2)*(model.nz-3));

    // b
    int count = 0;
    for (int k = 2; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            b(count) = -model.rhow[k] * model.zetap[i][k] * dx * dx;
            count++;
        }
    }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.setTolerance(1e-16);
    solver.compute(model.A);
    x = solver.solve(b);

    int cnt = 0;
    for (int k = 2; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            model.psi[i][k] = x[cnt];
            cnt++;
        }
    }
    for (int k = 2; k <= model.nz-2; k++) {
        model.psi[0][k] = model.psi[model.nx-2][k];
        model.psi[model.nx-1][k] = model.psi[1][k];
    }
    for (int i = 0; i <= model.nx-1; i++) {
        model.psi[i][0] = model.psi[i][1] = model.psi[i][model.nz-1] = 0.; 
    }

    for (int i = 1; i <= model.nx-2; i++) {
        for (int k = 1; k <= model.nz-2; k++) {
            model.u[i][k] = -(model.psi[i][k+1]-model.psi[i][k]) * model.rdz / model.rhou[k];
            model.w[i][k] = (model.psi[i+1][k]-model.psi[i][k]) * model.rdx / model.rhow[k];
        }
    }
	for (int i = 0; i < model.nz; i++) {
        model.w[i][model.nz-1] = model.w[i][0] = model.w[i][1] = 0.;
    }
    model.BoundaryProcess(model.u);
    model.BoundaryProcess(model.w);
    
    return;
}
#else
void Iteration::pubarTop_pt(vvmArray &model) {
    double uwUp = 0., uwDown = 0.;
    for (int i = 1; i < model.nx-1; i++) {
        uwUp += 0.5*(model.u[i][model.nz-1]+model.u[i][model.nz-2]) * 0.5*(model.w[i][model.nz-1]+model.w[i-1][model.nz-1]);
        uwDown += 0.5*(model.u[i][model.nz-2]+model.u[i][model.nz-3]) * 0.5*(model.w[i][model.nz-2]+model.w[i-1][model.nz-2]);
    }
    uwUp /= ((double) (model.nx - 2.)); // Cast the denominator to double
    uwDown /= ((double) (model.nx - 2.)); // Cast the denominator to double
    double puwb_pz = (uwUp - uwDown) * model.rdz;
    double uwb_rhoprho_pz = 0.5*(uwUp+uwDown) / model.rhou[model.nz-2] * (model.rhow[model.nz-2] - model.rhow[model.nz-3]) * model.rdz;
    model.ubarTopp = model.ubarTopm + d2t * (-puwb_pz - uwb_rhoprho_pz);
    return;
}

void Iteration::cal_w(vvmArray &model) {
    #if defined(ADVECTIONU)
        for (int i = 0; i <= model.nx-1; i++) { for (int k = 0; k <= model.nz-1; k++) { model.u[i][k] = 10.;}}
    #elif defined(ADVECTIONW)
        for (int i = 0; i <= model.nx-1; i++) { for (int k = 0; k <= model.nz-1; k++) { model.u[i][k] = 0.;}}
    #else
        // eigen solver
        Eigen::VectorXd x((model.nx-2)*(model.nz-3)), b((model.nx-2)*(model.nz-3));
        
        // b
        int count = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                b(count) = -model.rhow[k]*(model.zetap[i+1][k] - model.zetap[i][k]) * dx;
                count++;
            }
        }

        // Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        solver.setTolerance(1e-16);
        solver.compute(model.A);
        x = solver.solve(b);

        int cnt = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                model.w[i][k] = x[cnt];
                cnt++;
            }
        }

        // make sure w at surface and top = 0
        // for (int i = 0; i < model.nx; i++) model.w[i][1] = 0.;
        model.BoundaryProcessZETA(model.w);
    #endif
    return;
}

void Iteration::cal_u(vvmArray &model) {
    #if defined(ADVECTIONU)
        for (int i = 0; i <= model.nx-1; i++) { for (int k = 0; k <= model.nz-1; k++) { model.w[i][k] = 0.;}}
    #elif defined(ADVECTIONW)
        for (int i = 0; i <= model.nx-1; i++) { for (int k = 0; k <= model.nz-1; k++) { model.w[i][k] = 10.;}}
    #else
        /*
        Eigen::VectorXd y(model.nx-2), h(model.nx-2);
        
        // h
        for (int i = 1; i <= model.nx-2; i++) {
            h(i-1) = (0.-model.w[i][model.nz-2]) * dx + 0.5*(0.+model.w[i][model.nz-2]) / model.rhou[model.nz-2] * (model.rhow[model.nz-1] - model.rhow[model.nz-2]) * dx;
        }

        // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver_xi;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_xi;
        solver_xi.setTolerance(1e-16);
        y = solver_xi.compute(model.G).solve(h);

        // xi
        for (int i = 1; i <= model.nx-2; i++) {
            model.xi[i] = y[i-1];
        }
        model.xi[0] = model.xi[model.nx-2];
        model.xi[model.nx-1] = model.xi[1];

        */
        double a[model.nx-2];
        double b[model.nx-2];
        double c[model.nx-2];
        double d[model.nx-2];
        int i_k = 0;
        for (int i = 0; i < model.nx-2; i++) {
            i_k = i + 1;
            if (i == 0) a[i] = 0.;
            if (i == model.nx-3) c[i] = 0.;
            a[i] = -1.;
            b[i] = 2.;
            c[i] = -1.;
            d[i] = (0.-model.w[i_k][model.nz-2]) * dx + 0.5*(0.+model.w[i_k][model.nz-2]) / model.rhou[model.nz-2] * (model.rhow[model.nz-1] - model.rhow[model.nz-2]) * dx;
        }
        solve(a, b, c, d, model.nx-2);
        
        for (int i = 1; i <= model.nx-2; i++) {
            model.xi[i] = d[i-1];
        }
        model.xi[0] = model.xi[model.nx-2];
        model.xi[model.nx-1] = model.xi[1];


        // std::cout << "================" << std::endl;
        // std::cout << model.G * y << std::endl;

        // for (int i = 0; i < model.nx; i++) {
        // 	std::cout << model.xi[i] << " ";
        // }
        // std::cout << std::endl;

        // uxi
        for (int i = 1; i <= model.nx-2; i++) {
            model.uxi[i] = (model.xi[i] - model.xi[i-1]) * model.rdx;
        }
        model.uxi[0] = model.uxi[model.nx-2];
        model.uxi[model.nx-1] = model.uxi[1];

        // u_top
        for (int i = 1; i <= model.nx-2; i++) {
            #if defined(SHEAR)
            // You should modify the last number which means the fixed u at the top layer
                model.u[i][model.nz-2] = model.uxi[i] + model.ubarTopp + 10.;
            #else
                model.u[i][model.nz-2] = model.uxi[i] + model.ubarTopp;
                // model.u[i][model.nz-2] = model.uxi[i];
            #endif
        }
        model.u[0][model.nz-2] = model.u[model.nx-2][model.nz-2];
        model.u[model.nx-1][model.nz-2] = model.u[1][model.nz-2];

        // u
        for (int i = 1; i <= model.nx-2; i++) {
            double area = 0.;
            for (int k = model.nz-3; k >= 1; k--) {
                area += ((0.5*(model.w[i][k+2] + model.w[i][k+1]) - 0.5*(model.w[i-1][k+2] + model.w[i-1][k+1])) * model.rdx - 0.5*(model.zetap[i][k+2] + model.zetap[i][k+1])*model.rhou[k+1]) * -dz;
                model.u[i][k] = area + model.u[i][model.nz-2];
            }
        }
        model.BoundaryProcess(model.u);
    #endif
    return;
}
#endif

#if defined(WATER)
void Iteration::pqv_pt(vvmArray &model) {
    double upqv_px = 0., wpqv_pz = 0.;
    #if defined(FLUXFORM)
        double qvpu_px = 0., qvpw_pz = 0.;
        #if defined(LINEARIZEDQV)
            double wpqvb_pz = 0.;
        #else
            double wqv_rhoprho_pz = 0.;
        #endif
    #endif
    
    for (int i = 1; i <= model.nx-2; i++) {
        for (int k = 1; k <= model.nz-2; k++) {
            upqv_px = 0.5*(model.u[i+1][k]+model.u[i][k]) * (model.qv[i+1][k] - model.qv[i-1][k]) * model.r2dx;
            wpqv_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * (model.qv[i][k+1] - model.qv[i][k-1]) * model.r2dz;

            model.qvp[i][k] = model.qvm[i][k] + d2t * (-upqv_px - wpqv_pz);

            #if defined(FLUXFORM)
                qvpu_px = model.qv[i][k] * (model.u[i+1][k] - model.u[i][k]) * model.rdx;
                qvpw_pz = model.qv[i][k] * (model.w[i][k+1] - model.w[i][k]) * model.rdz;
                model.qvp[i][k] += d2t * (-qvpu_px - qvpw_pz);
                #if defined(LINEARIZEDQV)
                    wpqvb_pz = 0.5*(model.w[i][k+1] + model.w[i][k]) * (model.qvb[k+1] - model.qvb[k-1]) * r2dz;
                    model.qvp[i][k] += d2t * (-wpqvb_pz);
                #else
                    wqv_rhoprho_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * model.qv[i][k] / model.rhou[k] * (model.rhow[k+1] - model.rhow[k]) * model.rdz;
                    model.qvp[i][k] += d2t * (-wqv_rhoprho_pz);
                #endif
            #endif

            #if defined(TROPICALFORCING)
                model.qvp[i][k] += d2t * model.Q2LS[k];
            #endif

            // diffusion
            #ifdef DIFFUSION
                model.qvp[i][k] += d2t * Kx * model.rdx2 * (model.qvm[i+1][k] - 2. * model.qvm[i][k] + model.qvm[i-1][k]) + 
                                   d2t * Kz * model.rdz2 * (model.qvm[i][k+1] - 2. * model.qvm[i][k] + model.qvm[i][k-1]);
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
        for (int i = 1; i <= model.nx-2; i++) {
            for (int k = 1; k <= model.nz-2; k++) {
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
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.qv[i][k] += TIMETS * (model.qvp[i][k] - 2 * model.qv[i][k] + model.qvm[i][k]);
            }
        }
    #endif

    return;
}

void Iteration::pqc_pt(vvmArray &model) {
    double upqc_px = 0., wpqc_pz = 0.;
    #if defined(FLUXFORM)
        double qcpu_px = 0., qcpw_pz = 0., wqc_rhoprho_pz = 0.;;
    #endif

    for (int i = 1; i <= model.nx-2; i++) {
        for (int k = 1; k <= model.nz-2; k++) {
            upqc_px = 0.5*(model.u[i+1][k]+model.u[i][k]) * (model.qc[i+1][k] - model.qc[i-1][k]) * model.r2dx;
            wpqc_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * (model.qc[i][k+1] - model.qc[i][k-1]) * model.r2dz;

            model.qcp[i][k] = model.qcm[i][k] + d2t * (-upqc_px - wpqc_pz);

            #if defined(FLUXFORM)
                qcpu_px = model.qc[i][k] * (model.u[i+1][k] - model.u[i][k]) * model.rdx;
                qcpw_pz = model.qc[i][k] * (model.w[i][k+1] - model.w[i][k]) * model.rdz;
                wqc_rhoprho_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * model.qc[i][k] / model.rhou[k] * (model.rhow[k+1] - model.rhow[k]) * model.rdz;
                model.qcp[i][k] += d2t * (-qcpu_px - qcpw_pz - wqc_rhoprho_pz);
            #endif

            // negative qc process
            if (model.qcp[i][k] < 0.) model.qcp[i][k] = 0.;

            // saturation process: sink and source (qv <--> qc)
            condensation(model, i, k);

            #ifdef DIFFUSION
                model.qcp[i][k] += d2t * Kx * model.rdx2 * (model.qcm[i+1][k] - 2. * model.qcm[i][k] + model.qcm[i-1][k]) + 
                                   d2t * Kz * model.rdz2 * (model.qcm[i][k+1] - 2. * model.qcm[i][k] + model.qcm[i][k-1]);
            #endif
        }
    }
    model.BoundaryProcess(model.qcp);

    // Time filter
    #ifdef TIMEFILTER
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.qc[i][k] += TIMETS * (model.qcp[i][k] - 2 * model.qc[i][k] + model.qcm[i][k]);
            }
        }
    #endif
    return;
}

void Iteration::pqr_pt(vvmArray &model) {
    double upqr_px = 0., wpqr_pz = 0., VTpqr_pz = 0.;
    #if defined(FLUXFORM)
        double qrpu_px = 0., qrpw_pz = 0., wqr_rhoprho_pz = 0.;
        double qrpVT_pz = 0., VTqr_rhoprho_pz = 0.;
        double VT_up = 0., VT_down = 0.;
    #endif

    double VT = 0.;
    double VT0 = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            VT = 36.34 * pow(model.rhou[k]*model.qr[i][k], 0.1346) * pow(model.rhou[k]/model.rhow[1], -0.5);
            upqr_px = 0.5*(model.u[i+1][k]+model.u[i][k]) * (model.qr[i+1][k] - model.qr[i-1][k]) * model.r2dx;
            wpqr_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * (model.qr[i][k+1] - model.qr[i][k-1]) * model.r2dz;
            VTpqr_pz = VT * (model.qr[i][k+1] - model.qr[i][k-1]) * model.r2dz;

            model.qcp[i][k] = model.qcm[i][k] + d2t * (-upqr_px - wpqr_pz + VTpqr_pz);

            #if defined(FLUXFORM)
                VT_up = 36.34 * pow(model.rhou[k+1]*model.qr[i][k+1], 0.1346) * pow(model.rhou[k+1]/model.rhow[1], -0.5);
                VT_down = 36.34 * pow(model.rhou[k-1]*model.qr[i][k-1], 0.1346) * pow(model.rhou[k-1]/model.rhow[1], -0.5);

                qrpu_px = model.qr[i][k] * (model.u[i+1][k] - model.u[i][k]) * model.rdx;
                qrpw_pz = model.qr[i][k] * (model.w[i][k+1] - model.w[i][k]) * model.rdz;
                wqr_rhoprho_pz = 0.5*(model.w[i][k+1]+model.w[i][k]) * model.qr[i][k] / model.rhou[k] * (model.rhow[k+1] - model.rhow[k]) * model.rdz;
                qrpVT_pz = model.qr[i][k] * (VT_up - VT_down) * model.r2dz;
                VTqr_rhoprho_pz = VT * model.qr[i][k] / model.rhou[k] * (model.rhow[k+1] - model.rhow[k]) * model.rdz;

                model.qrp[i][k] += d2t * (-qrpu_px - qrpw_pz - wqr_rhoprho_pz + qrpVT_pz + VTqr_rhoprho_pz);
            #endif

            if (k == 1) {
                VT0 = 36.34 * pow(model.rhow[1]*model.qr[i][1], 0.1346);
                model.qrAcc[i] += (-model.rhou[1] * (0.5*(model.w[i][2]+model.w[i][0]) - VT0) * model.qr[i][1]);
            }

            #ifdef DIFFUSION
                model.qrp[i][k] += d2t * Kx * model.rdx2 * (model.qrm[i+1][k] - 2. * model.qrm[i][k] + model.qrm[i-1][k]) + 
                                    d2t * Kz * model.rdz2 * (model.qrm[i][k+1] - 2. * model.qrm[i][k] + model.qrm[i][k-1]);
            #endif

            autoconversion(model, i, k);
            accretion(model, i, k);
            evaporation(model, i, k);
        
            // negative qr process
            if (model.qrp[i][k] < 0.) model.qrp[i][k] = 0.;
        }
    }
    model.BoundaryProcess(model.qrp);

    // Time filter
    #ifdef TIMEFILTER
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
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
    model.autoconversion[i][k] = arcrdt;
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
    model.accretion[i][k] = arcrdt;
    return;
}

// evaporation of rain water
void Iteration::evaporation(vvmArray &model, int i, int k) {
    double qrplus = std::max(0., model.qrp[i][k]);
    #if defined(LINEARIZEDQV)
        double qvplus = std::max(0., model.qvp[i][k] + model.qvb[k]);
    #else
        double qvplus = std::max(0., model.qvp[i][k]);
    #endif

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

    if (erdt < 0.) {
        std::cout << "Evaporation Wrong" << std::endl;
        return;
    }

    model.qrp[i][k] = model.qrp[i][k] - erdt;
    model.qvp[i][k] = model.qvp[i][k] + erdt;
    model.thp[i][k] = model.thp[i][k] - Lv * erdt / (C_p * model.pib[k]);
    model.evaporation[i][k] = erdt;
    return;
}
#endif


void Iteration::updateMean(vvmArray &model) {
    double tb = 0.;
    for (int k = 1; k < model.nz-1; k++) {
        tb = 0.;
        for (int i = 1; i < model.nx-1; i++) {
            tb += model.thp[i][k];
        }
        model.tb[k] = tb / (double) (model.nx - 2.);
    }
    model.tb[0] = model.tb[1];
    model.tb[model.nz-1] = model.tb[model.nz-2];

    for (int k = 1; k < model.nz-1; k++) {
        model.tb_zeta[k] = 0.5 * (model.tb[k-1] + model.tb[k]);
    }
    model.tb_zeta[0] = model.tb_zeta[1];
    model.tb_zeta[model.nz-1] = model.tb_zeta[model.nz-2];
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
        #if defined(WATER)
            pqv_pt(model);
            pqc_pt(model);
            pqr_pt(model);
        #endif
        #if defined(STREAMFUNCTION)
            calpsiuw(model);
        #else
            pubarTop_pt(model);
            cal_w(model);
            cal_u(model);
        #endif

        #ifndef LINEARIZEDTH
            updateMean(model);
        #endif 

        // next step
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.zetam[i][k] = model.zeta[i][k];
                model.zeta[i][k] = model.zetap[i][k];

                model.thm[i][k] = model.th[i][k];
                model.th[i][k] = model.thp[i][k];

                model.ubarTopm = model.ubarTop;
                model.ubarTop = model.ubarTopp;

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