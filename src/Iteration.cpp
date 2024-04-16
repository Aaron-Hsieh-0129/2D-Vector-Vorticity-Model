#include "Iteration.hpp"
#include <iomanip>

vvm::PoissonSolver poissonSolver;
#if defined(WATER)
    vvm::MicroPhysics microphysics;
#endif

void Iteration::pzeta_pt(vvm &model) {
    model.Advection_zeta(model);
    model.Bouyancy(model);
    model.BoundaryProcess2D_westdown(model.zetap);
    return;
}

void Iteration::pth_pt(vvm &model) {
    model.Advection_thermo(model.thm, model.th, model.thp, model);
    model.BoundaryProcess2D_center(model.thp);
    return;
}

void Iteration::updateMean(vvm &model) {
    #if defined(AB3)
        for (int k = 0; k < model.nz; k++) {
            model.thbm[k] = model.thb[k];
            model.thvbm[k] = model.thvb[k];
        }
    #endif
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
            tb += model.thp[i][k];
            #if defined(WATER)
                qvb += model.qvp[i][k];
            #endif
        }
        model.thb[k] = tb / (double) (model.nx - 2.);
        #if defined(WATER)
            model.qvb[k] = qvb / (double) (model.nx - 2.);
            model.thvb[k] = model.thb[k] + 0.61 * model.qvb[k];
        #else
            model.thvb[k] = model.thb[k];
        #endif
    }
    model.BoundaryProcess1D_center(model.thb);
    model.BoundaryProcess1D_center(model.thvb);
    #if defined(WATER)
        model.BoundaryProcess1D_center(model.qvb);
    #endif

    for (int k = 1; k < model.nz-1; k++) {
        model.thb_zeta[k] = 0.5 * (model.thb[k-1] + model.thb[k]);
    }
    model.thb_zeta[0] = model.thb_zeta[1];
    model.thb_zeta[model.nz-1] = model.thb_zeta[model.nz-2];
    return;
}

#if defined(WATER)
void Iteration::pqv_pt(vvm &model) {
    model.Advection_thermo(model.qvm, model.qv, model.qvp, model);
    microphysics.NegativeValueProcess(model.qvp);
    model.BoundaryProcess2D_center(model.qvp);

    return;
}

void Iteration::pqc_pt(vvm &model) {
    model.Advection_thermo(model.qcm, model.qc, model.qcp, model);
    microphysics.NegativeValueProcess(model.qcp);
    model.BoundaryProcess2D_center(model.qcp);
    return;
}

void Iteration::pqr_pt(vvm &model) {
    model.Advection_thermo(model.qrm, model.qr, model.qrp, model);
    model.Advection_qrVT(model);
    microphysics.NegativeValueProcess(model.qrp);

    for (int i = 1; i <= model.nx-2; i++) {
        double VT = 36.34 * pow(model.rhou[1]*model.qr[i][1], 0.1346) * pow(model.rhou[1]/model.rhow[1], -0.5);
        double rain = -model.rhou[1] * (0.5*(model.w[i][2]+0.) - VT) * model.qr[i][1];
        if (rain < 0) {
            model.precip[i] = 0.;
        }
        else {
            model.qrAcc[i] += rain;
            model.precip[i] = rain;
        }
    }
            
    model.BoundaryProcess2D_center(model.qrp);
    return;
}
#endif

void Iteration::TimeMarching(vvm &model) {
    int n = 0;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;
    #ifndef PETSC
        poissonSolver.InitPoissonMatrix(model);
    #endif
    while (n < nmax) {
        std::cout << n << std::endl;
        // output
        if (n % OUTPUTSTEP == 0 || n == TIMEEND-1) {
            #if defined(OUTPUTNC)
                Output::output_nc(n, model);
            #endif

            #if defined(OUTPUTTXT)
                Output::output_zeta(n, model);
                Output::output_th(n, model);
                Output::output_u(n, model);
                Output::output_w(n, model);
                #if defined(WATER)
                    Output::output_qv(n, model);
                    Output::output_qc(n, model);
                    Output::output_qr(n, model);
                    Output::output_precip(n, model);
                    Output::output_precipAcc(n, model);
                #endif
            #endif
        }
        n++;
        
        #if defined(TROPICALFORCING)
            if (n * DT <= ADDFORCINGTIME) model.status_for_adding_forcing = true;
            else model.status_for_adding_forcing = false;

            // Generate new random th perturbation for tropical forcing case
            if (model.status_for_adding_forcing == true) {
                Init::RandomPerturbation(model, n);
            }
        #endif

        #if defined(AB3)
        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.um[i][k] = model.u[i][k];
                model.wm[i][k] = model.w[i][k];
            }
        }
        #endif

        pzeta_pt(model);
        pth_pt(model);
        #if defined(WATER)
            pqv_pt(model);
            pqc_pt(model);
            pqr_pt(model);
            #if defined(TROPICALFORCING)
                model.AddForcing(model);
            #endif
        #endif
        #if defined(STREAMFUNCTION)
            poissonSolver.calpsiuw(model);
        #else
            poissonSolver.pubarTop_pt(model);
            poissonSolver.cal_w(model);
            poissonSolver.cal_u(model);
        #endif
        

        #if defined(DIFFUSION)
            model.Diffusion(model.zetam, model.zetap, model);
            model.Diffusion(model.thm, model.thp, model);
            #if defined(WATER)
                model.Diffusion(model.qvm, model.qvp, model);
                model.Diffusion(model.qcm, model.qcp, model);
                model.Diffusion(model.qrm, model.qrp, model);
            #endif
        #endif
        model.BoundaryProcess2D_westdown(model.zetap);
        model.BoundaryProcess2D_center(model.thp);
        model.BoundaryProcess2D_westdown(model.w);
        model.BoundaryProcess2D_center(model.u);
        #if defined(WATER)
            model.BoundaryProcess2D_center(model.qvp);
            model.BoundaryProcess2D_center(model.qcp);
            model.BoundaryProcess2D_center(model.qrp);
        #endif

        #if defined(WATER)
            microphysics.autoconversion(model);
            microphysics.accretion(model);
            microphysics.evaporation(model);
            microphysics.condensation(model); // saturation adjustment

            // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
            microphysics.NegativeValueProcess(model.qvp);
            microphysics.NegativeValueProcess(model.qcp);
            microphysics.NegativeValueProcess(model.qrp);
        #endif

        model.BoundaryProcess2D_westdown(model.zetap);
        model.BoundaryProcess2D_center(model.thp);
        model.BoundaryProcess2D_westdown(model.w);
        model.BoundaryProcess2D_center(model.u);
        #if defined(WATER)
            model.BoundaryProcess2D_center(model.qvp);
            model.BoundaryProcess2D_center(model.qcp);
            model.BoundaryProcess2D_center(model.qrp);
        #endif

        #ifndef AB3
            #if defined(TIMEFILTER)
                model.TimeFilter(model.zetam, model.zeta, model.zetap, model);
                model.TimeFilter(model.thm, model.th, model.thp, model);
                #if defined(WATER)
                    model.TimeFilter(model.qvm, model.qv, model.qvp, model);
                    model.TimeFilter(model.qcm, model.qc, model.qcp, model);
                    model.TimeFilter(model.qrm, model.qr, model.qrp, model);
                #endif
            #endif
        #endif


        updateMean(model);

        // time marching
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
