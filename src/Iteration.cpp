#include "Declare.hpp"
#include <cmath>
#include "Timer.hpp"

void vvm::Iteration::pzeta_pt(vvm &model) {
    model.Advection_zeta(model);
    model.Bouyancy(model);
    model.BoundaryProcess2D_westdown(model.zetap, model.nx, model.nz);
    return;
}

void vvm::Iteration::pth_pt(vvm &model) {
    model.Advection_thermo(model.thm, model.th, model.thp, model);
    model.BoundaryProcess2D_center(model.thp, model.nx, model.nz);
    return;
}

void vvm::Iteration::updateMean(vvm &model) {
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
    model.BoundaryProcess1D_center(model.thb, model.nz);
    model.BoundaryProcess1D_center(model.thvb, model.nz);
    #if defined(WATER)
        model.BoundaryProcess1D_center(model.qvb, model.nz);
    #endif

    for (int k = 1; k < model.nz-1; k++) {
        model.thb_zeta[k] = 0.5 * (model.thb[k-1] + model.thb[k]);
    }
    model.thb_zeta[0] = model.thb_zeta[1];
    model.thb_zeta[model.nz-1] = model.thb_zeta[model.nz-2];
    return;
}

#if defined(WATER)
void vvm::Iteration::pqv_pt(vvm &model) {
    model.Advection_thermo(model.qvm, model.qv, model.qvp, model);
    vvm::MicroPhysics::NegativeValueProcess(model.qvp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qvp, model.nx, model.nz);

    return;
}

void vvm::Iteration::pqc_pt(vvm &model) {
    model.Advection_thermo(model.qcm, model.qc, model.qcp, model);
    vvm::MicroPhysics::NegativeValueProcess(model.qcp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qcp, model.nx, model.nz);
    return;
}

void vvm::Iteration::pqr_pt(vvm &model) {
    model.Advection_thermo(model.qrm, model.qr, model.qrp, model);
    model.Advection_qrVT(model);
    vvm::MicroPhysics::NegativeValueProcess(model.qrp, model.nx, model.nz);

    for (int i = 1; i <= model.nx-2; i++) {
        double VT = 1E-2 * (3634 * std::pow(1E-3*model.rhou[1] * model.qr[i][1], 0.1346) * std::pow(model.rhou[1]/model.rhow[1], -0.5));
        double rain = -model.rhou[1] * (0.5*(model.w[i][2]+0.) - VT) * model.qr[i][1];
        if (rain < 0) model.precip[i] = 0.;
        else model.precip[i] = rain;
    }
            
    model.BoundaryProcess2D_center(model.qrp, model.nx, model.nz);
    return;
}
#endif

void vvm::Iteration::nextTimeStep(vvm &model) {
    updateMean(model);
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



void vvm::Iteration::TimeMarching(vvm &model) {
    Timer timer;
    Timer time_all;

    int n = 0;
    double temp = model.TIMEEND / model.dt;
    int nmax = (int) temp;

    #if defined(LOADFROMPREVIOUSFILE)
        n = TIMENOW;
        std::cout << "timenow: " << n << std::endl;
    #endif

    #ifndef PETSC
        poissonSolver.InitPoissonMatrix(model);
    #endif
    while (n < nmax) {
        time_all.reset();
        std::cout << n << std::endl;
        // output
        if (n % model.OUTPUTSTEP == 0 || n == model.TIMEEND-1 || n == model.TIMEEND-2 || n == 550001) {
            #if defined(OUTPUTNC)
                vvm::Output::output_nc(n, model);
            #endif
            #if defined(OUTPUTTXT)
                vvm::Output::outputalltxt(n, model);
            #endif
        }
        n++;

        if (n % model.TIMEROUTPUTSIZE == 0) {
            vvm::Output::output_time_nc(n, model);
        }

        #if defined(TROPICALFORCING)
            if (n * model.dt <= model.ADDFORCINGTIME) model.status_for_adding_forcing = true;
            else model.status_for_adding_forcing = false;

            // Generate new random th perturbation for tropical forcing case
            if (model.status_for_adding_forcing == true) {
                vvm::Init::RandomPerturbation(model, n);
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

        timer.reset();
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
        model.t_advection[(n-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();

        timer.reset();
        #if defined(STREAMFUNCTION)
            vvm::PoissonSolver::calpsiuw(model);
        #else
            vvm::PoissonSolver::pubarTop_pt(model);
            vvm::PoissonSolver::cal_w(model);
            vvm::PoissonSolver::cal_u(model);
        #endif
        model.t_poisson[(n-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();
        
        timer.reset();
        #if defined(DIFFUSION)
            vvm::NumericalProcess::DiffusionAll(model);
        #endif
        model.t_diffusion[(n-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();


        vvm::BoundaryProcess2D_all(model);

        timer.reset();
        #if defined(WATER)
            vvm::MicroPhysics::autoconversion(model);
            vvm::MicroPhysics::accretion(model);
            vvm::MicroPhysics::evaporation(model);
            vvm::MicroPhysics::condensation(model); // saturation adjustment

            // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
            vvm::MicroPhysics::NegativeValueProcess(model.qvp, model.nx, model.nz);
            vvm::MicroPhysics::NegativeValueProcess(model.qcp, model.nx, model.nz);
            vvm::MicroPhysics::NegativeValueProcess(model.qrp, model.nx, model.nz);
        #endif
        model.t_microphysics[(n-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();

        vvm::BoundaryProcess2D_all(model);

        #if defined(TIMEFILTER)
            vvm::NumericalProcess::timeFilterAll(model);
        #endif

        vvm::Iteration::nextTimeStep(model);
        model.t_all[(n-1)%model.TIMEROUTPUTSIZE] = time_all.elapsed();
    }
    return;
}
