#include "Declare.hpp"

void vvm::Bouyancy(vvm &model) {
    double g_rhothbpth_px = 0.;
    double g_rhothbpth_px_m = 0.;
    for (int k = 2; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            g_rhothbpth_px = GRAVITY / model.rhow[k] * 0.5*((model.th[i][k] - model.th[i-1][k])/model.thb[k] + (model.th[i][k-1] - model.th[i-1][k-1])/model.thb[k-1]) * model.rdx;
            g_rhothbpth_px_m = GRAVITY / model.rhow[k] * 0.5*((model.thm[i][k] - model.thm[i-1][k])/model.thb[k] + (model.thm[i][k-1] - model.thm[i-1][k-1])/model.thb[k-1]) * model.rdx;

            #if defined(AB3)
                model.zetap[i][k] += (1.5*DT*g_rhothbpth_px / model.rhow[k]) - (0.5*DT*g_rhothbpth_px_m / model.rhow[k]);
            #else
                model.zetap[i][k] += model.d2t * g_rhothbpth_px / model.rhow[k];
            #endif

            #if defined(WATER)
                double g_rhopqv_px = GRAVITY / model.rhow[k] * (0.5*(model.qv[i][k] + model.qv[i][k-1]) - 0.5*(model.qv[i-1][k] + model.qv[i-1][k-1])) * model.rdx;
                double g_rhopqc_px = GRAVITY / model.rhow[k] * (0.5*(model.qc[i][k] + model.qc[i][k-1]) - 0.5*(model.qc[i-1][k] + model.qc[i-1][k-1])) * model.rdx;
                double g_rhopqr_px = GRAVITY / model.rhow[k] * (0.5*(model.qr[i][k] + model.qr[i][k-1]) - 0.5*(model.qr[i-1][k] + model.qr[i-1][k-1])) * model.rdx;
                #if defined(AB3)
                    double g_rhopqv_px_m = GRAVITY / model.rhow[k] * (0.5*(model.qvm[i][k] + model.qvm[i][k-1]) - 0.5*(model.qvm[i-1][k] + model.qvm[i-1][k-1])) * model.rdx;
                    double g_rhopqc_px_m = GRAVITY / model.rhow[k] * (0.5*(model.qcm[i][k] + model.qcm[i][k-1]) - 0.5*(model.qcm[i-1][k] + model.qcm[i-1][k-1])) * model.rdx;
                    double g_rhopqr_px_m = GRAVITY / model.rhow[k] * (0.5*(model.qrm[i][k] + model.qrm[i][k-1]) - 0.5*(model.qrm[i-1][k] + model.qrm[i-1][k-1])) * model.rdx; 
                    model.zetap[i][k] += 1.5*DT*(0.61 * g_rhopqv_px - g_rhopqc_px - g_rhopqr_px) - 0.5*DT*(0.61*g_rhopqv_px_m - g_rhopqc_px_m - g_rhopqr_px_m);
                #else
                    model.zetap[i][k] += model.d2t * (0.61 * g_rhopqv_px - g_rhopqc_px - g_rhopqr_px);
                #endif
            #endif
        }
    }
    return;
}
