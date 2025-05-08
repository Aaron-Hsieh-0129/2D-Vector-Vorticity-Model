#include "Declare.hpp"

#if defined(TIMEFILTER)
void vvm::NumericalProcess::TimeFilter(double **previous, double **now, double **future, vvm &model) {
    for (int i = 0; i <= model.nx-1; i++) {
        for (int k = 0; k <= model.nz-1; k++) {
            now[i][k] += model.TIMETS * (future[i][k] - 2.*now[i][k] + previous[i][k]);
        }
    }
    return;
}

void vvm::NumericalProcess::timeFilterAll(vvm &model) {
    #if defined(TIMEFILTER)
        TimeFilter(model.zetam, model.zeta, model.zetap, model);
        TimeFilter(model.thm, model.th, model.thp, model);
        #if defined(WATER)
            TimeFilter(model.qvm, model.qv, model.qvp, model);
            TimeFilter(model.qcm, model.qc, model.qcp, model);
            TimeFilter(model.qrm, model.qr, model.qrp, model);
        #endif
    #endif
}
#endif

#if defined(DIFFUSION_VVM)
void vvm::NumericalProcess::Diffusion(double **var_in, double **var_out, vvm &model) {
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            var_out[i][k] += model.d2t * model.Kx * model.rdx2 * (var_in[i+1][k] - 2. * var_in[i][k] + var_in[i-1][k]) + 
                             model.d2t * model.Kz * model.rdz2 * (var_in[i][k+1] - 2. * var_in[i][k] + var_in[i][k-1]);
        }
    }
    return;
}


void vvm::NumericalProcess::DiffusionAll(vvm &model) {
    Diffusion(model.zetam, model.zetap, model);
    Diffusion(model.thm, model.thp, model);
    #if defined(WATER)
        Diffusion(model.qvm, model.qvp, model);
        Diffusion(model.qcm, model.qcp, model);
        Diffusion(model.qrm, model.qrp, model);
        #if defined(P3_MICROPHY)
            Diffusion(model.ncm, model.ncp, model);
            Diffusion(model.nrm, model.nrp, model);
            Diffusion(model.nim, model.nip, model);
            Diffusion(model.qitotm, model.qitotp, model);
            Diffusion(model.qirimm, model.qirimp, model);
            Diffusion(model.birimm, model.birimp, model);
        #endif
    #endif
}
#endif


void vvm::NumericalProcess::GravityWaveDampingExponential(vvm &model) {
    int n_damp = 17000./model.dz + 1;

    for (int k = n_damp; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            model.thp[i][k] -= model.dt * model.nudge_tau[k] * (model.thp[i][k] - model.thb_init[k]);
            model.zetap[i][k] -= model.dt * model.nudge_tau[k] * (model.zetap[i][k]);
            #if defined(WATER)
                model.qvp[i][k] -= model.dt * model.nudge_tau[k] * (model.qvp[i][k] - model.qvb0[k]);
                model.qcp[i][k] -= model.dt * model.nudge_tau[k] * model.qcp[i][k];
                model.qrp[i][k] -= model.dt * model.nudge_tau[k] * model.qrp[i][k];
                #if defined(P3_MICROPHY)
                    model.ncp[i][k] -= model.dt * model.nudge_tau[k] * model.ncp[i][k];
                    model.nrp[i][k] -= model.dt * model.nudge_tau[k] * model.nrp[i][k];
                    model.nip[i][k] -= model.dt * model.nudge_tau[k] * model.nip[i][k];
                    model.qitotp[i][k] -= model.dt * model.nudge_tau[k] * model.qitotp[i][k];
                    model.qirimp[i][k] -= model.dt * model.nudge_tau[k] * model.qirimp[i][k];
                    model.birimp[i][k] -= model.dt * model.nudge_tau[k] * model.birimp[i][k];
                #endif
            #endif
        }
    }
}


void vvm::NumericalProcess::Nudge_qv(vvm &model) {
    if (model.moisture_nudge_time == 0) return;

    for (int i = 0; i <= model.nx-1; i++) {
        for (int k = 0; k <= model.nz-1; k++) {
            model.qvp[i][k] = model.qvp[i][k] + model.dt * (model.qvb0[k] - model.qv[i][k]) / model.moisture_nudge_time;
        }
    }
    return;
}

void vvm::NumericalProcess::NegativeValueProcess(double **var, int nx, int nz) {
    double positive = 0.;
    double negative = 0.;
    for (int k = 1; k <= nz-2; k++) {
        for (int i = 1; i <= nx-2; i++) {
            if (var[i][k] >= 0.) positive += var[i][k];
            else {
                negative += var[i][k];
                var[i][k] = 0.;
            }
        }
    }

    if (positive == 0. || std::abs(negative) > positive) return;

    double correctionRatio = 1. - std::abs(negative/positive);
    for (int k = 1; k <= nz-2; k++) {
        for (int i = 1; i <= nx-2; i++) {
            if (var[i][k] > 0) var[i][k] = var[i][k] * correctionRatio;
        }
    }
    return;
}
