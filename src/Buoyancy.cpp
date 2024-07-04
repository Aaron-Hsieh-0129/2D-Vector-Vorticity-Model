#include "Declare.hpp"

double getTHV(int i, int k, vvm &model) {
    #if defined(WATER)
        return model.th[i][k] + 0.608 * model.qv[i][k];
    #else
        return model.th[i][k];
    #endif
}

void vvm::Bouyancy(vvm &model) {
    double g_rhothvbpthv_px = 0.;
    #if defined(WATER)
        double g_rhopqc_px = 0., g_rhopqr_px = 0.;
    #endif
    
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 2; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            g_rhothvbpthv_px = model.GRAVITY / model.rhow[k] * 0.5*((getTHV(i, k, model) - getTHV(i-1, k, model))/model.thvb[k] + (getTHV(i, k-1, model) - getTHV(i-1, k-1, model))/model.thvb[k-1]) * model.rdx;

            #if defined(AB2)
                model.dth_buoyancy[i][k][(model.step+1)%2] = g_rhothvbpthv_px;
                if (model.step == 0) model.dth_buoyancy[i][k][0] = model.dth_buoyancy[i][k][1];
            #endif

            #if defined(WATER)
                g_rhopqc_px = model.GRAVITY / model.rhow[k] * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * model.rdx;
                g_rhopqr_px = model.GRAVITY / model.rhow[k] * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * model.rdx;
                #if defined(AB2)
                    model.dth_buoyancy[i][k][(model.step+1)%2] += -g_rhopqc_px - g_rhopqr_px;
                    if (model.step == 0) model.dth_buoyancy[i][k][0] = model.dth_buoyancy[i][k][1];
                #endif
            #endif

            #if defined(AB2)
                model.zetap[i][k] += (1.5*model.dt*model.dth_buoyancy[i][k][(model.step+1)%2]) - (0.5*model.dt*model.dth_buoyancy[i][k][model.step%2]);
            #else
                model.zetap[i][k] += model.d2t * g_rhothvbpthv_px;
                #if defined(WATER)
                    model.zetap[i][k] += model.d2t * (-g_rhopqc_px - g_rhopqr_px);
                #endif
            #endif
        }
    }
    return;
}
