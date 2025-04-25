#include "Declare.hpp"
#include "Timer.hpp"

#if defined(P3_MICROPHY)
extern "C" {
    void __microphy_p3_MOD_p3_main(
        double* qc, double* nc, double* qr, double* nr, 
        double* th_old, double* th, double* qv_old, double* qv, 
        double* dt, double* qitot, double* qirim, double* qiliq, 
        double* ni, double* birim, double* zi, double* ssat, 
        double* w, double* p, double* dz, int* itt, 
        double* precip_liq, double* precip_sol, 
        int* one, int* ncols, int* one2, int* nz, int* nCat, 
        double* diag_ze, double* diag_effc, double* diag_effi, 
        double* diag_vmi, double* diag_di, double* diag_rhoi, 
        int* n_diag_2d, double* diag_2d, int* n_diag_3d, double* diag_3d, 
        bool* log_predictNc, char* model, 
        double* clbfact_dep, double* clbfact_sub, 
        bool* debug_on, bool* scpf_on, double* scpf_pfrac, 
        double* scpf_resfact, double* cldfrac, 
        bool* trplMomI, bool* liqfrac,
        double* , double*, double*, double*, double*,
        double* , double*, double*, double*, double*,
        double* , double*, double*, double*, double*, size_t model_name_len
    );
}
#endif

void vvm::Iteration::pzeta_pt(vvm &model) {
    model.Advection_zeta(model);
    model.Bouyancy(model);
    model.BoundaryProcess2D_westdown(model.zetap, model.nx, model.nz);
    return;
}

void vvm::Iteration::pth_pt(vvm &model) {
    model.Advection_thermo(model.thm, model.th, model.thp, model.dth_advect, model);
    model.BoundaryProcess2D_center(model.thp, model.nx, model.nz);
    return;
}

void vvm::Iteration::updateMean(vvm &model) {
    #if defined(AB2)
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
            model.thvb[k] = model.thb[k] * (1 + 0.608 * model.qvb[k]);
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
#if defined(KESSLER_MICROPHY)
void vvm::Iteration::pqv_pt(vvm &model) {
    model.Advection_thermo(model.qvm, model.qv, model.qvp, model.dqv_advect, model);
    #if defined(KESSLER_MICROPHY)
    vvm::MicroPhysics::NegativeValueProcess(model.qvp, model.nx, model.nz);
    #endif
    model.BoundaryProcess2D_center(model.qvp, model.nx, model.nz);

    return;
}

void vvm::Iteration::pqc_pt(vvm &model) {
    model.Advection_thermo(model.qcm, model.qc, model.qcp, model.dqc_advect, model);
    #if defined(KESSLER_MICROPHY)
    vvm::MicroPhysics::NegativeValueProcess(model.qcp, model.nx, model.nz);
    #endif
    model.BoundaryProcess2D_center(model.qcp, model.nx, model.nz);
    return;
}

void vvm::Iteration::pqr_pt(vvm &model) {
    model.Advection_thermo(model.qrm, model.qr, model.qrp, model.dqr_advect, model);
    model.Advection_qrVT(model);
    #if defined(KESSLER_MICROPHY)
    vvm::MicroPhysics::NegativeValueProcess(model.qrp, model.nx, model.nz);
    for (int i = 1; i <= model.nx-2; i++) {
        double VT = 1E-2 * (3634 * std::pow(1E-3*model.rhou[1] * model.qr[i][1], 0.1346) * std::pow(model.rhou[1]/model.rhow[1], -0.5));
        double rain = -model.rhou[1] * (0.5*(model.w[i][2]+0.) - VT) * model.qr[i][1];
        if (rain < 0) model.precip[i] = 0.;
        else model.precip[i] = rain;
    }
    #endif
            
    model.BoundaryProcess2D_center(model.qrp, model.nx, model.nz);
    return;
}
#endif

#if defined(P3_MICROPHY)
void vvm::Iteration::pqmicrophy_pt(vvm &model) {
    model.Advection_thermo(model.qvm, model.qv, model.qvp, model.dqv_advect, model);
    model.Advection_thermo(model.qcm, model.qc, model.qcp, model.dqc_advect, model);
    model.Advection_thermo(model.qrm, model.qr, model.qrp, model.dqr_advect, model);
    model.Advection_thermo(model.ncm, model.nc, model.ncp, model.dnc_advect, model);
    model.Advection_thermo(model.nrm, model.nr, model.nrp, model.dnr_advect, model);
    model.Advection_thermo(model.nim, model.ni, model.nip, model.dni_advect, model);
    model.Advection_thermo(model.qitotm, model.qitot, model.qitotp, model.dqitot_advect, model);
    model.Advection_thermo(model.qirimm, model.qirim, model.qirimp, model.dqirim_advect, model);
    /* model.Advection_thermo(model.qiliqm, model.qiliq, model.qiliqp, model.dqiliq_advect, model); */
    model.Advection_thermo(model.birimm, model.birim, model.birimp, model.dbirim_advect, model);
    model.BoundaryProcess2D_center(model.qvp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qcp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qrp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.ncp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.nrp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.nip, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qitotp, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.qirimp, model.nx, model.nz);
    /* model.BoundaryProcess2D_center(model.qiliqp, model.nx, model.nz); */
    model.BoundaryProcess2D_center(model.birimp, model.nx, model.nz);
}
#endif
#endif

void vvm::Iteration::nextTimeStep(vvm &model) {
    updateMean(model);
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
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

                #if defined(P3_MICROPHY)
                    model.ncm[i][k] = model.nc[i][k];
                    model.nc[i][k] = model.ncp[i][k];

                    model.nrm[i][k] = model.nr[i][k];
                    model.nr[i][k] = model.nrp[i][k];

                    model.nim[i][k] = model.ni[i][k];
                    model.ni[i][k] = model.nip[i][k];

                    model.qitotm[i][k] = model.qitot[i][k];
                    model.qitot[i][k] = model.qitotp[i][k];

                    model.qirimm[i][k] = model.qirim[i][k];
                    model.qirim[i][k] = model.qirimp[i][k];

                    /* model.qiliqm[i][k] = model.qiliq[i][k]; */
                    /* model.qiliq[i][k] = model.qiliqp[i][k]; */

                    model.birimm[i][k] = model.birim[i][k];
                    model.birim[i][k] = model.birimp[i][k];
                #endif
            #endif
        }
    }
}



void vvm::Iteration::TimeMarching(vvm &model) {
    Timer timer;
    Timer time_all;

    model.step = 0;
    double temp = model.TIMEEND / model.dt;
    int nmax = (int) temp;

    #if defined(LOADFROMPREVIOUSFILE)
        model.step = TIMENOW;
        std::cout << "timenow: " << model.step << std::endl;
    #endif


    // Radiation scheme call for step 0
    #if defined(RTERRTMGP)
        vvm::Radiation::solve_radiation(model);
    #endif

    while (model.step < nmax) {
        time_all.reset();
        std::cout << model.step << std::endl;
        // output
        if (model.step % model.OUTPUTSTEP == 0 || model.step == model.TIMEEND-1 || model.step == model.TIMEEND-2 || model.step == 550001) {
            #if defined(_OPENMP)
            #pragma omp critical
            {
            #endif
                #if defined(OUTPUTNC)
                    vvm::Output::output_nc(model.step, model);
                #endif
                #if defined(OUTPUTTXT)
                    vvm::Output::outputalltxt(model.step, model);
                #endif
            #if defined(_OPENMP)
            }
            #endif
        }
        model.step++;


        // if (model.step % model.TIMEROUTPUTSIZE == 0) {
        //     #if defined(_OPENMP)
        //     #pragma omp critical
        //     {
        //     #endif
        //         #if defined(OUTPUTNC)
        //             vvm::Output::output_time_nc(model.step, model);
        //         #endif
        //     #if defined(_OPENMP)
        //     }
        //     #endif
        // }

        timer.reset();
        pzeta_pt(model);
        pth_pt(model);
        #if defined(WATER)
            #if defined(KESSLER_MICROPHY)
                pqv_pt(model);
                pqc_pt(model);
                pqr_pt(model);
            #endif

            #if defined(P3_MICROPHY)
                pqmicrophy_pt(model);
            #endif

            if (model.step * model.dt <= model.addforcingtime) model.status_for_adding_forcing = true;
            else model.status_for_adding_forcing = false;

            // Generate new random th perturbation for tropical forcing case
            if (model.status_for_adding_forcing == true) {
                vvm::Init::RandomPerturbation(model, model.step+0, -0.001, 0.001, 1.);
            }
            model.AddForcing(model);
        #endif
        vvm::BoundaryProcess2D_all(model);
        model.t_advection[(model.step-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();

        timer.reset();
        #if defined(STREAMFUNCTION)
            vvm::PoissonSolver::calpsiuw(model);
        #else
            vvm::PoissonSolver::pubarTop_pt(model);
            vvm::PoissonSolver::cal_w(model);
            vvm::PoissonSolver::cal_u(model);
        #endif
        vvm::BoundaryProcess2D_all(model);
        model.t_poisson[(model.step-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();


        timer.reset();
        #if defined(WATER)
            #if defined(KESSLER_MICROPHY)
                vvm::MicroPhysics::autoconversion(model);
                vvm::MicroPhysics::accretion(model);
                vvm::MicroPhysics::evaporation(model);
                vvm::MicroPhysics::condensation(model); // saturation adjustment
            #endif
            // It is supposed to not have negative values. But due to numerical process, it might produce some teeny-tiny values.
            vvm::NumericalProcess::NegativeValueProcess(model.qvp, model.nx, model.nz);
            vvm::NumericalProcess::NegativeValueProcess(model.qcp, model.nx, model.nz);
            vvm::NumericalProcess::NegativeValueProcess(model.qrp, model.nx, model.nz);

            #if defined(P3_MICROPHY)
            vvm::NumericalProcess::NegativeValueProcess(model.ncp, model.nx, model.nz);
            vvm::NumericalProcess::NegativeValueProcess(model.nrp, model.nx, model.nz);
            vvm::NumericalProcess::NegativeValueProcess(model.qitotp, model.nx, model.nz);
            vvm::NumericalProcess::NegativeValueProcess(model.nip, model.nx, model.nz);


            for (int k = 0; k < model.nz; k++) {
                for (int i = 0; i < model.nx; i++) {
                    model.qiliqp[i][k] = 0.;
                }
            }
            int one = 1;

            double *qcp1d = nullptr;
            double *ncp1d = nullptr;
            double *qrp1d = nullptr;
            double *nrp1d = nullptr;
            double *th1d = nullptr;
            double *thp1d = nullptr;
            double *qv1d = nullptr;
            double *qvp1d = nullptr;
            double *qitotp1d = nullptr;
            double *qirimp1d = nullptr;
            double *qiliqp1d = nullptr;
            double *nip1d = nullptr;
            double *birimp1d = nullptr;
            double *zi_all1d = nullptr;
            double *ssat_all1d = nullptr;
            double *w_all1d = nullptr;
            double *pb_all1d = nullptr;
            double *dz_all1d = nullptr;
            double *precip_liq1d = nullptr;
            double *precip_sol1d = nullptr;
            double *diag_ze1d = nullptr;
            double *diag_effc1d = nullptr;
            double *diag_effi1d = nullptr;
            double *diag_vmi1d = nullptr;
            double *diag_di1d = nullptr;
            double *diag_rhoi1d = nullptr;
            double *diag_2d1d = nullptr;
            double *diag_3d1d = nullptr;
            double *cldfrac1d = nullptr;
            for (int i = 0; i < model.nx; i++) {
                qcp1d = model.qcp[i]; 
                ncp1d = model.ncp[i]; 
                qrp1d = model.qrp[i]; 
                nrp1d = model.nrp[i];
                th1d = model.th[i]; 
                thp1d = model.thp[i]; 
                qv1d = model.qv[i]; 
                qvp1d = model.qvp[i];
                qitotp1d = model.qitotp[i]; 
                qirimp1d = model.qirimp[i]; 
                qiliqp1d = model.qiliqp[i]; 
                nip1d = model.nip[i];
                birimp1d = model.birimp[i]; 
                zi_all1d = model.zi_all[i]; 
                ssat_all1d = model.ssat_all[i]; 
                w_all1d = model.w[i]; 
                pb_all1d = model.pb_all[i];
                dz_all1d = model.dz_all[i];
                precip_liq1d = &model.precip_liq[i];
                precip_sol1d = &model.precip_sol[i];
                diag_ze1d = model.diag_ze[i];
                diag_effc1d = model.diag_effc[i];
                diag_effi1d = model.diag_effi[i];
                diag_vmi1d = model.diag_vmi[i];
                diag_di1d = model.diag_di[i];
                diag_rhoi1d = model.diag_rhoi[i];
                diag_2d1d = model.diag_2d[i];
                diag_3d1d = model.diag_3d[i][0];
                cldfrac1d = model.cldfrac[i];

                __microphy_p3_MOD_p3_main(
                    qcp1d, ncp1d, qrp1d, nrp1d, 
                    th1d, thp1d, qv1d, qvp1d, &model.dt,
                    qitotp1d, qirimp1d, qiliqp1d, nip1d,
                    birimp1d, zi_all1d, ssat_all1d, w_all1d, pb_all1d,
                    dz_all1d, &model.step, precip_liq1d, precip_sol1d, &one, &one, &one, &model.nz, 
                    &vvm::P3::nCat, diag_ze1d, diag_effc1d, diag_effi1d,
                    diag_vmi1d, diag_di1d, diag_rhoi1d, 
                    &vvm::P3::n_diag_2d, diag_2d1d, &vvm::P3::n_diag_3d, diag_3d1d,
                    &vvm::P3::log_predictNc, vvm::P3::model_name, &vvm::P3::clbfact_dep, 
                    &vvm::P3::clbfact_sub, &vvm::P3::debug_on, &vvm::P3::scpf_on, 
                    &vvm::P3::scpf_pfrac, &vvm::P3::scpf_resfact, cldfrac1d, 
                    &vvm::P3::trplMomI, &vvm::P3::liqfrac, 
                    nullptr, nullptr, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr,
                    nullptr, nullptr, nullptr, nullptr, nullptr, strlen(vvm::P3::model_name)

                );
            }
            #endif
        #endif
        vvm::BoundaryProcess2D_all(model);
        model.t_microphysics[(model.step-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();


        // model.SurfaceFlux(model);

        timer.reset();
        updateMean(model);
        #if defined(DIFFUSION_VVM)
            vvm::NumericalProcess::DiffusionAll(model);
        #else
            vvm::Turbulence::RKM_RKH(model);
        #endif
        // Nudging process to damp the gravity wave
        // vvm::NumericalProcess::Nudge_theta(model);
        // if (model.CASE != 2) vvm::NumericalProcess::Nudge_zeta(model);
        // vvm::NumericalProcess::Nudge_qv(model);
        vvm::BoundaryProcess2D_all(model);
        model.t_diffusion[(model.step-1)%model.TIMEROUTPUTSIZE] = timer.elapsed();

        #if defined(TIMEFILTER)
            vvm::NumericalProcess::timeFilterAll(model);
        #endif

        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        #if defined(RTERRTMGP)
            if (model.step % 40 == 0) {
                vvm::Radiation::solve_radiation(model);
            }
        #endif

        vvm::Iteration::nextTimeStep(model);
        model.t_all[(model.step-1)%model.TIMEROUTPUTSIZE] = time_all.elapsed();
    }
    return;
}
