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

#if defined(DIFFUSION)
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

void vvm::NumericalProcess::Nudge_theta(vvm &model) {
    // Nudge the top layers (larger than 15 km)
    int k_start = 12000. / model.dz + 1;
    double CGR = 0.;
    for (int k = k_start; k < model.nz-1; k++) {
        CGR = model.CRAD * (model.z[k] - model.z[(model.nz-1)-k_start]) / (model.z[model.nz-1] - model.z[(model.nz-1)-k_start]);
        for (int i = 1; i < model.nx-1; i++) {
            model.thp[i][k] -= model.dt * CGR * (model.thp[i][k] - model.thb_init[k]);
        }
    }
    return;
}

void vvm::NumericalProcess::Nudge_zeta(vvm &model) {
    // Nudge the top layers (larger than 15 km)
    int k_start = 12000. / model.dz + 1;
    double CGR = 0.;
    for (int k = k_start; k < model.nz-1; k++) {
        /*CGR = model.CRAD * (model.z_zeta[k] - model.z_zeta[(model.nz-1)-k_start]) / (model.z_zeta[model.nz-1] - model.z_zeta[(model.nz-1)-k_start]);*/
        CGR = (1./360.) * (model.z_zeta[k] - model.z_zeta[(model.nz-1)-k_start]) / (model.z_zeta[model.nz-1] - model.z_zeta[(model.nz-1)-k_start]);
        for (int i = 1; i < model.nx-1; i++) {
            model.zetap[i][k] -= model.dt * CGR * (model.zetap[i][k]);
        }
    }
    return;
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

