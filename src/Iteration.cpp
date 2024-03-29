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
        for (int k = 0; k < model.nz; k++) model.thbm[k] = model.thb[k];
    #endif
    double tb = 0.;
    for (int k = 1; k < model.nz-1; k++) {
        tb = 0.;
        
        for (int i = 1; i < model.nx-1; i++) {
            tb += model.thp[i][k];
        }
        model.thb[k] = tb / (double) (model.nx - 2.);
    }
    model.thb[0] = model.thb[1];
    model.thb[model.nz-1] = model.thb[model.nz-2];

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

    double VT_u = 6., VT_d = 6.;
    for (int i = 1; i <= model.nx-2; i++) {
        double rain = -model.rhow[1] * (model.w[i][2] - 0.5*(VT_u + VT_d)) * model.qr[i][1];
        if (rain < 0) continue;
        model.qrAcc[i] += rain;
        model.precip[i] = rain;
    }
            
    model.BoundaryProcess2D_center(model.qrp);
    return;
}
#endif

void Iteration::TimeMarching(vvm &model) {
    int n = 0;
    double temp = TIMEEND / DT;
    int nmax = (int) temp;
    while (n < nmax) {
        std::cout << n << std::endl;
        // output
        if (n % OUTPUTSTEP == 0 || n == TIMEEND-1) {
            #if defined(OUTPUTNC)
                Output::output_nc(n, model);
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

        for (int i = 0; i <= model.nx-1; i++) {
            for (int k = 0; k <= model.nz-1; k++) {
                model.um[i][k] = model.u[i][k];
                model.wm[i][k] = model.w[i][k];
            }
        }

        pzeta_pt(model);
        pth_pt(model);
        #if defined(STREAMFUNCTION)
            calpsiuw(model);
        #else
            poissonSolver.pubarTop_pt(model);
            poissonSolver.cal_w(model);
            poissonSolver.cal_u(model);
        #endif
        
        #if defined(WATER)
            pqv_pt(model);
            pqc_pt(model);
            pqr_pt(model);
            #if defined(TROPICALFORCING)
                model.AddForcing(model);
            #endif
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
