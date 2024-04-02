#include "Declare.hpp"

double getTHV(int i, int k, vvm &model) {
    #if defined(WATER)
        return model.th[i][k] + 0.61 * model.qv[i][k];
    #else
        return model.th[i][k];
    #endif
}

#if defined(AB3)
double getTHVm(int i, int k, vvm &model) {
    #if defined(WATER)
        return model.th[i][k] + 0.61 * model.qv[i][k];
    #else
        return model.th[i][k];
    #endif
}
#endif

void vvm::Bouyancy(vvm &model) {
    double g_rhothbpth_px = 0.;
    double g_rhothbpth_px_m = 0.;
    double thv = 0.;
    for (int k = 2; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            // g_rhothbpth_px = GRAVITY / model.rhow[k] * 0.5*((model.th[i][k] - model.th[i-1][k])/model.thb[k] + (model.th[i][k-1] - model.th[i-1][k-1])/model.thb[k-1]) * model.rdx;
            g_rhothbpth_px = GRAVITY / model.rhow[k] * 0.5*((getTHV(i, k, model) - getTHV(i-1, k, model))/model.thvb[k] + (getTHV(i, k-1, model) - getTHV(i-1, k-1, model))/model.thvb[k-1]) * model.rdx;

            #if defined(AB3)
                // g_rhothbpth_px_m = GRAVITY / model.rhow[k] * 0.5*((model.thm[i][k] - model.thm[i-1][k])/model.thb[k] + (model.thm[i][k-1] - model.thm[i-1][k-1])/model.thb[k-1]) * model.rdx;
                g_rhothbpth_px_m = GRAVITY / model.rhow[k] * 0.5*((getTHVm(i, k, model) - getTHVm(i-1, k, model))/model.thvbm[k] + (getTHVm(i, k-1, model) - getTHVm(i-1, k-1, model))/model.thvbm[k-1]) * model.rdx;
                model.zetap[i][k] += (1.5*DT*g_rhothbpth_px) - (0.5*DT*g_rhothbpth_px_m);
            #else
                model.zetap[i][k] += model.d2t * g_rhothbpth_px;
            #endif

            #if defined(WATER)
                // double g_rhopqv_px = GRAVITY / model.rhow[k] * (0.5*(model.qv[i][k] + model.qv[i][k-1]) - 0.5*(model.qv[i-1][k] + model.qv[i-1][k-1])) * model.rdx;
                double g_rhopqc_px = GRAVITY / model.rhow[k] * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * model.rdx;
                double g_rhopqr_px = GRAVITY / model.rhow[k] * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * model.rdx;
                #if defined(AB3)
                    double g_rhopqc_px_m = GRAVITY / model.rhow[k] * (0.5*(model.qcm[i][k] + model.qcm[i][k-1]) - 0.5*(model.qcm[i-1][k] + model.qcm[i-1][k-1])) * model.rdx;
                    double g_rhopqr_px_m = GRAVITY / model.rhow[k] * (0.5*(model.qrm[i][k] + model.qrm[i][k-1]) - 0.5*(model.qrm[i-1][k] + model.qrm[i-1][k-1])) * model.rdx; 
                    model.zetap[i][k] += 1.5*DT*(- g_rhopqc_px - g_rhopqr_px) - 0.5*DT*(- g_rhopqc_px_m - g_rhopqr_px_m);
                #else
                    model.zetap[i][k] += model.d2t * (- g_rhopqc_px - g_rhopqr_px);
                #endif
            #endif
        }
    }
    return;
}
