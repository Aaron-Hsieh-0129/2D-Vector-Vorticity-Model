#include "Declare.hpp"

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

void vvm::NumericalProcess::DiffusionAll(vvm &model) {
    Diffusion(model.zetam, model.zetap, model);
    Diffusion(model.thm, model.thp, model);
    #if defined(WATER)
        // Diffusion(model.qvm, model.qvp, model);
        #pragma omp parallel for collapse(2)
        for (int k = 1; k < model.nz-1; k++) {
            for (int i = 1; i < model.nx-1; i++) {
                model.qvp[i][k] += model.d2t * 1. * model.rdx2 * (model.qvm[i+1][k] - 2. * model.qvm[i][k] + model.qvm[i-1][k]) + 
                                   model.d2t * 1. * model.rdz2 * (model.qvm[i][k+1] - 2. * model.qvm[i][k] + model.qvm[i][k-1]);
            }
        }
        Diffusion(model.qcm, model.qcp, model);
        Diffusion(model.qrm, model.qrp, model);
    #endif
}