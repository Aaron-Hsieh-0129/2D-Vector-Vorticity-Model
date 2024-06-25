/**
 * @file Advection.cpp
 * @author Aaron Hsieh (b08209006@ntu.edu.tw)
 * @brief 
 * @version 0.1
 * @date 2024-03-13
 * 
 * @copyright Copyright (c) 2024
 * 
*/

#include "Declare.hpp"
#include <cmath>

// TODO: Arakawa Jacobian for streamfunction
// TODO: Fix the blow up in AB3.
double PLUS(double var) {
    return 0.5 * (var + std::fabs(var));
}

double MINU(double var) {
    return 0.5 * (var - std::fabs(var));
}

void vvm::Advection_thermo(double **past, double **now, double **future, double ***dvar, vvm &model) {
    double prhouvar_px_rho = 0., prhowvar_pz_rho = 0.;
    double *flux_ucont, *flux_wcont;
    double **flux_u = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_ucont);
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);

    for (int k = 1; k <= model.nz-1; k++) {
        for (int i = 1; i <= model.nx-1; i++) {
            flux_u[i][k] = model.rhou[k] * model.u[i][k] * (now[i][k] + now[i-1][k]);
            flux_w[i][k] = model.rhow[k] * model.w[i][k] * (now[i][k] + now[i][k-1]);
            #if defined(AB3)
                if (i >= 2 && i <= model.nx-3 && k >= 2 && k <= model.nz-3) {
                    flux_u[i][k] += -ALPHA/3. * 
                                    (model.rhou[k]*PLUS(model.u[i][k]) * (now[i][k] - now[i-1][k])
                                        - std::sqrt(model.rhou[k]*PLUS(model.u[i][k]) * model.rhou[k]*PLUS(model.u[i-1][k]))*(now[i-1][k] - now[i-2][k])
                                   - model.rhou[k]*MINU(model.u[i][k]) * (now[i][k] - now[i-1][k])
                                        - std::sqrt(std::fabs(model.rhou[k]*MINU(model.u[i][k]) * model.rhou[k]*MINU(model.u[i+1][k])))*(now[i+1][k] - now[i][k]));

                    flux_w[i][k] += -ALPHA/3. * 
                                    (model.rhow[k]*PLUS(model.w[i][k]) * (now[i][k] - now[i][k-1])
                                        - std::sqrt(model.rhow[k]*PLUS(model.w[i][k]) * model.rhow[k-1]*PLUS(model.w[i][k-1]))*(now[i][k-1] - now[i][k-2])
                                   - model.rhow[k]*MINU(model.w[i][k]) * (now[i][k] - now[i][k-1])
                                        - std::sqrt(std::fabs(model.rhow[k]*MINU(model.w[i][k]) * model.rhow[k+1]*MINU(model.w[i][k+1])))*(now[i][k+1] - now[i][k]));
                }
            #endif
        }
    }
    model.BoundaryProcess2D_center(flux_u, model.nx, model.nz);
    model.BoundaryProcess2D_center(flux_w, model.nx, model.nz);
    for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][1] = flux_w[i][model.nz-1] = 0.;

    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhouvar_px_rho = (flux_u[i+1][k] - flux_u[i][k]) * model.r2dx / model.rhou[k];
            prhowvar_pz_rho = (flux_w[i][k+1] - flux_w[i][k]) * model.r2dz / model.rhou[k];

            #if defined(AB3)
                dvar[i][k][(model.step+1)%2] = -prhouvar_px_rho - prhowvar_pz_rho;
                if (model.step == 0) dvar[i][k][0] = dvar[i][k][1];
                future[i][k] = now[i][k] + 1.5*model.dt*dvar[i][k][(model.step+1)%2] - 0.5*model.dt*dvar[i][k][model.step%2];
            #else
                future[i][k] = past[i][k] + model.d2t * (-prhouvar_px_rho - prhowvar_pz_rho);
            #endif
        }
    }

    vvm::deallocate2DContinuousArray(flux_u, flux_ucont);
    vvm::deallocate2DContinuousArray(flux_w, flux_wcont);
    return;
}

void vvm::Advection_zeta(vvm &model) {
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            model.U_w[i][k] = 0.25 * (model.rhou[k] * (model.u[i+1][k]+model.u[i][k]) + model.rhou[k-1] * (model.u[i+1][k-1]+model.u[i][k-1]));
            model.W_u[i][k] = 0.25 * (model.rhow[k+1] * (model.w[i][k+1]+model.w[i-1][k+1]) + model.rhow[k] * (model.w[i][k]+model.w[i-1][k]));
        }
    }
    model.BoundaryProcess2D_center(model.U_w, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.W_u, model.nx, model.nz);
    for (int i = 0; i < model.nx; i++) model.W_u[i][0] = model.W_u[i][model.nz-1] = 0.;

    double *flux_ucont, *flux_wcont;
    double **flux_u = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_ucont);
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);

    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 0; i <= model.nx-2; i++) {
            flux_u[i][k] = model.U_w[i][k] * (model.zeta[i+1][k] + model.zeta[i][k]);
            flux_w[i][k] = model.W_u[i][k] * (model.zeta[i][k+1] + model.zeta[i][k]);
            #if defined(AB3)
                if (i >= 2 && i <= model.nx-3 && k >= 2 && k <= model.nz-3) {
                    flux_u[i][k] += -ALPHA/3. * (PLUS(model.U_w[i][k])*(model.zeta[i+1][k]-model.zeta[i][k]) 
                                                    - std::sqrt(PLUS(model.U_w[i][k]) * PLUS(model.U_w[i-1][k]))*(model.zeta[i][k]-model.zeta[i-1][k]) - 
                                                 MINU(model.U_w[i][k])*(model.zeta[i+1][k]-model.zeta[i][k]) 
                                                    - std::sqrt(std::fabs(MINU(model.U_w[i][k]) * MINU(model.U_w[i+1][k])))*(model.zeta[i+2][k]-model.zeta[i+1][k]));
                    flux_w[i][k] += -ALPHA/3. * (PLUS(model.W_u[i][k])*(model.zeta[i][k+1]-model.zeta[i][k]) 
                                                    - std::sqrt(PLUS(model.W_u[i][k]) * PLUS(model.W_u[i][k-1]))*(model.zeta[i][k]-model.zeta[i][k-1]) - 
                                                 MINU(model.W_u[i][k])*(model.zeta[i][k+1]-model.zeta[i][k]) 
                                                    - std::sqrt(std::fabs(MINU(model.W_u[i][k]) * MINU(model.W_u[i][k+1])))*(model.zeta[i][k+2]-model.zeta[i][k+1]));
                }
            #endif
        }
    }
    model.BoundaryProcess2D_center(flux_u, model.nx, model.nz);
    model.BoundaryProcess2D_center(flux_w, model.nx, model.nz);
    for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][model.nz-1] = 0.;

    double prhouzeta_px_rho = 0., prhowzeta_pz_rho = 0.;
    for (int k = 2; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhouzeta_px_rho = (flux_u[i][k] - flux_u[i-1][k]) * model.r2dx / model.rhow[k];
            prhowzeta_pz_rho = (flux_w[i][k] - flux_w[i][k-1]) * model.r2dz / model.rhow[k];

            #if defined(AB3)
                model.dzeta_advect[i][k][(model.step+1)%2] = -prhouzeta_px_rho - prhowzeta_pz_rho;
                if (model.step == 0) model.dzeta_advect[i][k][0] = model.dzeta_advect[i][k][1];
                model.zetap[i][k] = model.zeta[i][k] + 1.5*model.dt*model.dzeta_advect[i][k][(model.step+1)%2] - 0.5*model.dt*model.dzeta_advect[i][k][model.step%2];
            #else
                model.zetap[i][k] = model.zetam[i][k] + model.d2t * (-prhouzeta_px_rho - prhowzeta_pz_rho);
            #endif
        }
    }

    vvm::deallocate2DContinuousArray(flux_u, flux_ucont);
    vvm::deallocate2DContinuousArray(flux_w, flux_wcont);
    return;
}

#if defined(WATER)
void vvm::Advection_qrVT(vvm &model) {
    double prhoVTqr_pz_rho = 0.;
    double *flux_wcont;
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);


    double VT = 0.;
    #if defined(AB3)
        double VT_u = 0., VT_d = 0.;
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            VT = 1E-2 * (3634 * pow(1E-3*model.rhow[k] * 0.5*(model.qr[i][k]+model.qr[i][k-1]), 0.1346) * pow(model.rhow[k]/model.rhow[1], -0.5));

            flux_w[i][k] = model.rhow[k] * VT * (model.qr[i][k] + model.qr[i][k-1]);
            #if defined(AB3)
                if (i >= 2 && i <= model.nx-2 && k >= 2 && k <= model.nz-2) {
                    VT_u = 1E-2 * (3634 * pow(1E-3*model.rhow[k+1] * 0.5*(model.qr[i][k+1]+model.qr[i][k]), 0.1346) * pow(model.rhow[k+1]/model.rhow[1], -0.5));
                    VT_d = 1E-2 * (3634 * pow(1E-3*model.rhow[k-1] * 0.5*(model.qr[i][k-1]+model.qr[i][k-2]), 0.1346) * pow(model.rhow[k-1]/model.rhow[1], -0.5));

                    flux_w[i][k] += -ALPHA/3. * 
                                    (model.rhow[k]*PLUS(VT) * (model.qr[i][k] - model.qr[i][k-1])
                                        - std::sqrt(model.rhow[k]*PLUS(VT) * model.rhow[k-1]*PLUS(VT_d))*(model.qr[i][k-1] - model.qr[i][k-2])
                                   - model.rhow[k]*MINU(VT) * (model.qr[i][k] - model.qr[i][k-1])
                                        - std::sqrt(std::fabs(model.rhow[k]*MINU(VT) * model.rhow[k+1]*MINU(VT_u)))*(model.qr[i][k+1] - model.qr[i][k]));
                }
            #endif
        }
    }
    model.BoundaryProcess2D_center(flux_w, model.nx, model.nz);
    for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][1] = flux_w[i][model.nz-1] = 0.;

    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhoVTqr_pz_rho = (flux_w[i][k+1] - flux_w[i][k]) * model.r2dz / model.rhou[k];

            #if defined(AB3)
                model.dqr_VT[i][k][(model.step+1)%2] = prhoVTqr_pz_rho;
                if (model.step == 0) model.dqr_VT[i][k][0] = model.dqr_VT[i][k][1];
                model.qrp[i][k] += 1.5*model.dt*model.dqr_VT[i][k][(model.step+1)%2] - 0.5*model.dt*model.dqr_VT[i][k][model.step%2];
            #else
                model.qrp[i][k] = model.qrm[i][k] + model.d2t * (prhoVTqr_pz_rho);
            #endif
        }
    }

    vvm::deallocate2DContinuousArray(flux_w, flux_wcont);
    return;
}
#endif
