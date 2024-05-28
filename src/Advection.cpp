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
#include <iostream>

// TODO: Arakawa Jacobian for streamfunction
// TODO: Fix the blow up in AB3.
double PLUS(double var) {
    return 0.5 * (var + std::abs(var));
}

double MINUS(double var) {
    return std::abs(0.5 * (var - std::abs(var)));
}

double U_zeta(int i, int k, vvm &model) {
    return 0.25 * (model.rhou[k] * (model.u[i+1][k]+model.u[i][k]) + model.rhou[k-1] * (model.u[i+1][k-1]+model.u[i][k-1])); 
}

double W_zeta(int i, int k, vvm &model) {
    return 0.25 * (model.rhow[k+1] * (model.w[i][k+1]+model.w[i-1][k+1]) + model.rhow[k] * (model.w[i][k]+model.w[i-1][k]));
}

double Um_zeta(int i, int k, vvm &model) {
    return 0.25 * (model.rhou[k] * (model.um[i+1][k]+model.um[i][k]) + model.rhou[k-1] * (model.um[i+1][k-1]+model.um[i][k-1])); 
}

double Wm_zeta(int i, int k, vvm &model) {
    return 0.25 * (model.rhow[k+1] * (model.wm[i][k+1]+model.wm[i-1][k+1]) + model.rhow[k] * (model.wm[i][k]+model.wm[i-1][k]));
}
void vvm::Advection_zeta(vvm &model) {
    double U_ik_rr = 0., U_ik_r = 0., U_ik_l = 0., U_ik_ll = 0.;
    double W_ik_uu = 0., W_ik_u = 0., W_ik_d = 0., W_ik_dd = 0.;
    double F_ik_r = 0., F_ik_l = 0., G_ik_u = 0., G_ik_d = 0.;
    double prhouzeta_px_rho = 0., prhowzeta_pz_rho = 0.;
    double prhouzeta_px_rho_m = 0., prhowzeta_pz_rho_m = 0.;
    for (int k = 2; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            U_ik_r = U_zeta(i, k, model);
            U_ik_l = U_zeta(i-1, k, model);
            F_ik_r = U_ik_r * (model.zeta[i+1][k] + model.zeta[i][k]);
            F_ik_l = U_ik_l * (model.zeta[i][k] + model.zeta[i-1][k]);

            W_ik_u = W_zeta(i, k, model);
            W_ik_d = W_zeta(i, k-1, model);
            G_ik_u = W_ik_u * (model.zeta[i][k+1] + model.zeta[i][k]);
            G_ik_d = W_ik_d * (model.zeta[i][k] + model.zeta[i][k-1]);
            #if defined(AB3)
            if (i != 1 && i != model.nx-2 && k != 2 && k != model.nz-2 && ALPHA != 0.) {
                U_ik_rr = U_zeta(i+1, k, model);
                U_ik_ll = U_zeta(i-2, k, model);
                F_ik_r += ALPHA/3. * (PLUS(U_ik_r)*(model.zeta[i+1][k]-model.zeta[i][k]) - std::pow(PLUS(U_ik_r), 0.5)*std::pow(PLUS(U_ik_l), 0.5)*(model.zeta[i][k]-model.zeta[i-1][k]) - 
                                      MINUS(U_ik_r)*(model.zeta[i+1][k]-model.zeta[i][k]) - std::pow(MINUS(U_ik_r), 0.5)*std::pow(MINUS(U_ik_rr), 0.5)*(model.zeta[i+2][k]-model.zeta[i+1][k]));
                F_ik_l += ALPHA/3. * (PLUS(U_ik_l)*(model.zeta[i][k]-model.zeta[i-1][k]) - std::pow(PLUS(U_ik_l), 0.5)*std::pow(PLUS(U_ik_ll), 0.5)*(model.zeta[i-1][k]-model.zeta[i-2][k]) - 
                                      MINUS(U_ik_l)*(model.zeta[i][k]-model.zeta[i-1][k]) - std::pow(MINUS(U_ik_l), 0.5)*std::pow(MINUS(U_ik_r), 0.5)*(model.zeta[i+1][k]-model.zeta[i][k]));


                W_ik_uu = W_zeta(i, k+1, model);
                W_ik_dd = W_zeta(i, k-2, model);
                G_ik_u += ALPHA/3. * (PLUS(W_ik_u)*(model.zeta[i][k+1]-model.zeta[i][k]) - std::pow(PLUS(W_ik_u), 0.5)*std::pow(PLUS(W_ik_d), 0.5)*(model.zeta[i][k]-model.zeta[i][k-1]) - 
                                      MINUS(W_ik_u)*(model.zeta[i][k+1]-model.zeta[i][k]) - std::pow(MINUS(W_ik_u), 0.5)*std::pow(MINUS(W_ik_uu), 0.5)*(model.zeta[i][k+2]-model.zeta[i][k+1]));
                G_ik_d += ALPHA/3. * (PLUS(W_ik_d)*(model.zeta[i][k]-model.zeta[i][k-1]) - std::pow(PLUS(W_ik_d), 0.5)*std::pow(PLUS(W_ik_dd), 0.5)*(model.zeta[i][k-1]-model.zeta[i][k-2]) - 
                                      MINUS(W_ik_d)*(model.zeta[i][k]-model.zeta[i][k-1]) - std::pow(MINUS(W_ik_d), 0.5)*std::pow(MINUS(W_ik_u), 0.5)*(model.zeta[i][k+1]-model.zeta[i][k]));
            }
            #endif
            prhouzeta_px_rho = (F_ik_r - F_ik_l) * model.r2dx / model.rhow[k];
            prhowzeta_pz_rho = (G_ik_u - G_ik_d) * model.r2dz / model.rhow[k];

            #if defined(AB3)
            U_ik_r = Um_zeta(i, k, model);
            U_ik_l = Um_zeta(i-1, k, model);
            F_ik_r = U_ik_r * (model.zetam[i+1][k] + model.zetam[i][k]);
            F_ik_l = U_ik_l * (model.zetam[i][k] + model.zetam[i-1][k]);

            W_ik_u = Wm_zeta(i, k, model);
            W_ik_d = Wm_zeta(i, k-1, model);
            G_ik_u = W_ik_u * (model.zetam[i][k+1] + model.zetam[i][k]);
            G_ik_d = W_ik_d * (model.zetam[i][k] + model.zetam[i][k-1]);
            if (i != 1 && i != model.nx-2 && k != 2 && k != model.nz-2 && ALPHA != 0.) {
                U_ik_rr = Um_zeta(i+1, k, model);
                U_ik_ll = Um_zeta(i-2, k, model);
                F_ik_r += ALPHA/3. * (PLUS(U_ik_r)*(model.zetam[i+1][k]-model.zetam[i][k]) - std::pow(PLUS(U_ik_r), 0.5)*std::pow(PLUS(U_ik_l), 0.5)*(model.zetam[i][k]-model.zetam[i-1][k]) - 
                                      MINUS(U_ik_r)*(model.zetam[i+1][k]-model.zetam[i][k]) - std::pow(MINUS(U_ik_r), 0.5)*std::pow(MINUS(U_ik_rr), 0.5)*(model.zetam[i+2][k]-model.zetam[i+1][k]));
                F_ik_l += ALPHA/3. * (PLUS(U_ik_l)*(model.zetam[i][k]-model.zetam[i-1][k]) - std::pow(PLUS(U_ik_l), 0.5)*std::pow(PLUS(U_ik_ll), 0.5)*(model.zetam[i-1][k]-model.zetam[i-2][k]) - 
                                      MINUS(U_ik_l)*(model.zetam[i][k]-model.zetam[i-1][k]) - std::pow(MINUS(U_ik_l), 0.5)*std::pow(MINUS(U_ik_r), 0.5)*(model.zetam[i+1][k]-model.zetam[i][k]));


                W_ik_uu = W_zeta(i, k+1, model);
                W_ik_dd = W_zeta(i, k-2, model);
                G_ik_u += ALPHA/3. * (PLUS(W_ik_u)*(model.zetam[i][k+1]-model.zetam[i][k]) - std::pow(PLUS(W_ik_u), 0.5)*std::pow(PLUS(W_ik_d), 0.5)*(model.zetam[i][k]-model.zetam[i][k-1]) - 
                                      MINUS(W_ik_u)*(model.zetam[i][k+1]-model.zetam[i][k]) - std::pow(MINUS(W_ik_u), 0.5)*std::pow(MINUS(W_ik_uu), 0.5)*(model.zetam[i][k+2]-model.zetam[i][k+1]));
                G_ik_d += ALPHA/3. * (PLUS(W_ik_d)*(model.zetam[i][k]-model.zetam[i][k-1]) - std::pow(PLUS(W_ik_d), 0.5)*std::pow(PLUS(W_ik_dd), 0.5)*(model.zetam[i][k-1]-model.zetam[i][k-2]) - 
                                      MINUS(W_ik_d)*(model.zetam[i][k]-model.zetam[i][k-1]) - std::pow(MINUS(W_ik_d), 0.5)*std::pow(MINUS(W_ik_u), 0.5)*(model.zetam[i][k+1]-model.zetam[i][k]));
            }
            prhouzeta_px_rho_m = (F_ik_r - F_ik_l) * model.r2dx / model.rhow[k];
            prhowzeta_pz_rho_m = (G_ik_u - G_ik_d) * model.r2dz / model.rhow[k];
            #endif

            #if defined(AB3)
                model.zetap[i][k] = model.zeta[i][k] + 1.5*DT*(-prhouzeta_px_rho - prhowzeta_pz_rho) - 0.5*DT*(-prhouzeta_px_rho_m - prhowzeta_pz_rho_m);
            #else
                model.zetap[i][k] = model.zetam[i][k] + model.d2t * (-prhouzeta_px_rho - prhowzeta_pz_rho);
            #endif
        }
    }
    model.BoundaryProcess2D_westdown(model.zetap, model.nx, model.nz);
    return;
}

void vvm::Advection_thermo(double **past, double **now, double **future, vvm &model) {
    double F_ik_r = 0., F_ik_l = 0., G_ik_u = 0., G_ik_d = 0.;
    double pvar_px = 0., prhowvar_pz_rho = 0.;
    double pvar_px_m = 0., prhowvar_pz_rho_m = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            F_ik_r = model.u[i+1][k] * (now[i+1][k] + now[i][k]);
            F_ik_l = model.u[i][k] * (now[i][k] + now[i-1][k]);
            G_ik_u = model.rhow[k+1] * model.w[i][k+1] * (now[i][k+1] + now[i][k]);
            G_ik_d = model.rhow[k] * model.w[i][k] * (now[i][k] + now[i][k-1]);
            #if defined(AB3)
            if (i != 1 && i != model.nx-2 && k != 1 && k != model.nz-2 && ALPHA != 0.) {
                F_ik_r += -ALPHA/3. * (PLUS(model.u[i+1][k]) * (now[i+1][k] - now[i][k]) - std::pow(PLUS(model.u[i+1][k]), 0.5)*std::pow(PLUS(model.u[i][k]), 0.5)*(now[i][k] - now[i-1][k]) - 
                                       MINUS(model.u[i+1][k]) * (now[i+1][k] - now[i][k]) - std::pow(MINUS(model.u[i+1][k]), 0.5)*std::pow(MINUS(model.u[i+2][k]), 0.5)*(now[i+2][k] - now[i+1][k]));
                F_ik_l += -ALPHA/3. * (PLUS(model.u[i][k]) * (now[i][k] - now[i-1][k]) - std::pow(PLUS(model.u[i][k]), 0.5)*std::pow(PLUS(model.u[i-1][k]), 0.5)*(now[i-1][k] - now[i-2][k]) - 
                                       MINUS(model.u[i][k]) * (now[i][k] - now[i-1][k]) - std::pow(MINUS(model.u[i][k]), 0.5)*std::pow(MINUS(model.u[i+1][k]), 0.5)*(now[i+1][k] - now[i][k]));
                G_ik_u += -ALPHA/3. * (model.rhow[k+1]*PLUS(model.w[i][k+1]) * (now[i][k+1] - now[i][k]) - std::pow(model.rhow[k+1]*PLUS(model.w[i][k+1]), 0.5)*std::pow(model.rhow[k]*PLUS(model.w[i][k]), 0.5)*(now[i][k] - now[i][k-1]) - 
                                       model.rhow[k+1]*MINUS(model.w[i][k+1]) * (now[i][k+1] - now[i][k]) - std::pow(model.rhow[k+1]*MINUS(model.w[i][k+1]), 0.5)*std::pow(model.rhow[k+2]*MINUS(model.w[i][k+2]), 0.5)*(now[i][k+2] - now[i][k+1]));
                G_ik_d += -ALPHA/3. * (model.rhow[k]*PLUS(model.w[i][k]) * (now[i][k] - now[i][k-1]) - std::pow(model.rhow[k]*PLUS(model.w[i][k]), 0.5)*std::pow(model.rhow[k]*PLUS(model.w[i][k-1]), 0.5)*(now[i][k-1] - now[i][k-2]) - 
                                       model.rhow[k]*MINUS(model.w[i][k]) * (now[i][k] - now[i][k-1]) - std::pow(model.rhow[k]*MINUS(model.w[i][k]), 0.5)*std::pow(model.rhow[k+1]*MINUS(model.w[i][k+1]), 0.5)*(now[i][k+1] - now[i][k]));
            }
            #endif
            pvar_px = (F_ik_r - F_ik_l) * model.r2dx;
            prhowvar_pz_rho = (G_ik_u - G_ik_d) * model.r2dz / model.rhou[k];

            #if defined(AB3)
            F_ik_r = F_ik_l = G_ik_u = G_ik_d = 0.;
            F_ik_r = model.um[i+1][k] * (past[i+1][k] + past[i][k]);
            F_ik_l = model.um[i][k] * (past[i][k] + past[i-1][k]);
            G_ik_u = model.rhow[k+1] * model.wm[i][k+1] * (past[i][k+1] + past[i][k]);
            G_ik_d = model.rhow[k] * model.wm[i][k] * (past[i][k] + past[i][k-1]);
            if (i != 1 && i != model.nx-2 && k != 1 && k != model.nz-2 && ALPHA != 0.) {
                F_ik_r += -ALPHA/3. * (PLUS(model.um[i+1][k]) * (past[i+1][k] - past[i][k]) - std::pow(PLUS(model.um[i+1][k]), 0.5)*std::pow(PLUS(model.um[i][k]), 0.5)*(past[i][k] - past[i-1][k]) - 
                                       MINUS(model.um[i+1][k]) * (past[i+1][k] - past[i][k]) - std::pow(MINUS(model.um[i+1][k]), 0.5)*std::pow(MINUS(model.um[i+2][k]), 0.5)*(past[i+2][k] - past[i+1][k]));
                F_ik_l += -ALPHA/3. * (PLUS(model.um[i][k]) * (past[i][k] - past[i-1][k]) - std::pow(PLUS(model.um[i][k]), 0.5)*std::pow(PLUS(model.um[i-1][k]), 0.5)*(past[i-1][k] - past[i-2][k]) - 
                                       MINUS(model.um[i][k]) * (past[i][k] - past[i-1][k]) - std::pow(MINUS(model.um[i][k]), 0.5)*std::pow(MINUS(model.um[i+1][k]), 0.5)*(past[i+1][k] - past[i][k]));
                G_ik_u += -ALPHA/3. * (model.rhow[k+1]*PLUS(model.wm[i][k+1]) * (past[i][k+1] - past[i][k]) - std::pow(model.rhow[k+1]*PLUS(model.wm[i][k+1]), 0.5)*std::pow(model.rhow[k]*PLUS(model.wm[i][k]), 0.5)*(past[i][k] - past[i][k-1]) - 
                                       model.rhow[k+1]*MINUS(model.wm[i][k+1]) * (past[i][k+1] - past[i][k]) - std::pow(model.rhow[k+1]*MINUS(model.wm[i][k+1]), 0.5)*std::pow(model.rhow[k+2]*MINUS(model.wm[i][k+2]), 0.5)*(past[i][k+2] - past[i][k+1]));
                G_ik_d += -ALPHA/3. * (model.rhow[k]*PLUS(model.wm[i][k]) * (past[i][k] - past[i][k-1]) - std::pow(model.rhow[k]*PLUS(model.wm[i][k]), 0.5)*std::pow(model.rhow[k]*PLUS(model.wm[i][k-1]), 0.5)*(past[i][k-1] - past[i][k-2]) - 
                                       model.rhow[k]*MINUS(model.wm[i][k]) * (past[i][k] - past[i][k-1]) - std::pow(model.rhow[k]*MINUS(model.wm[i][k]), 0.5)*std::pow(model.rhow[k+1]*MINUS(model.wm[i][k+1]), 0.5)*(past[i][k+1] - past[i][k]));
            }
            pvar_px_m = (F_ik_r - F_ik_l) * model.r2dx;
            prhowvar_pz_rho_m = (G_ik_u - G_ik_d) * model.r2dz / model.rhou[k];
            #endif

            #if defined(AB3)
                future[i][k] = now[i][k] + 1.5*DT*(-pvar_px - prhowvar_pz_rho) - 0.5*DT*(-pvar_px_m - prhowvar_pz_rho_m);
            #else
                future[i][k] = past[i][k] + model.d2t * (-pvar_px - prhowvar_pz_rho);
            #endif
        }
    }
    return;
}

#if defined(WATER)
void vvm::Advection_qrVT(vvm &model) {
    #if defined(VTconst)
        double VT_u = 6., VT_d = 6.;
        double rhoVT_u = 6., rhoVT_d = 6.;
        double rhoVT_uu = 6., rhoVT_dd = 6.;
    #else
        double rhoVT_u = 0., rhoVT_uu = 0., rhoVT_d = 0., rhoVT_dd = 0.;
    #endif
    double G_ik_u = 0., G_ik_d = 0.;
    double prhoVTqr_pz_rho = 0.;
    double prhoVTqr_pz_rho_m = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            #if defined(VTconst)
                G_ik_u = model.rhow[k+1] * VT_u * (model.qr[i][k+1] + model.qr[i][k]);
                G_ik_d = model.rhow[k] * VT_d * (model.qr[i][k] + model.qr[i][k-1]);
            #else
                rhoVT_u = model.rhow[k+1] * 1E-2 * (3634 * pow(1E-3*model.rhow[k+1] * 0.5*(model.qr[i][k+1]+model.qr[i][k]), 0.1346) * pow(model.rhow[k+1]/model.rhow[1], -0.5));
                rhoVT_d = model.rhow[k] * 1E-2 * (3634 * pow(1E-3*model.rhow[k] * 0.5*(model.qr[i][k]+model.qr[i][k-1]), 0.1346) * pow(model.rhow[k]/model.rhow[1], -0.5));
                G_ik_u = rhoVT_u * (model.qr[i][k+1] + model.qr[i][k]);
                G_ik_d = rhoVT_d * (model.qr[i][k] + model.qr[i][k-1]);
            #endif

            #if defined(AB3)
                if (i != 1 && i != model.nx-2 && k != 1 && k != model.nz-2) {
                    #if defined(VTconst)
                        rhoVT_u = model.rhow[k+1] * VT_u;
                        rhoVT_d = model.rhow[k] * VT_d;
                    #else
                        rhoVT_uu = model.rhow[k+2] * 1E-2 * (36.34 * pow(1E-3*model.rhow[k+2] * 0.5*(model.qr[i][k+2]+model.qr[i][k+1]), 0.1346) * pow(model.rhow[k+2]/model.rhow[1], -0.5));
                        rhoVT_dd = model.rhow[k-1] * 1E-2 * (36.34 * pow(1E-3*model.rhow[k-1] * 0.5*(model.qr[i][k-1]+model.qr[i][k-2]), 0.1346) * pow(model.rhow[k-1]/model.rhow[1], -0.5));
                    #endif
                    G_ik_u += -ALPHA/3. * (PLUS(rhoVT_u) * (model.qr[i][k+1] - model.qr[i][k]) - std::pow(PLUS(rhoVT_u), 0.5)*std::pow(PLUS(rhoVT_d), 0.5)*(model.qr[i][k] - model.qr[i][k-1]) - 
                                           MINUS(rhoVT_u) * (model.qr[i][k+1] - model.qr[i][k]) - std::pow(MINUS(rhoVT_u), 0.5)*std::pow(MINUS(rhoVT_uu), 0.5)*(model.qr[i][k+2] - model.qr[i][k+1]));
                    G_ik_d += -ALPHA/3. * (PLUS(rhoVT_d) * (model.qr[i][k] - model.qr[i][k-1]) - std::pow(PLUS(rhoVT_d), 0.5)*std::pow(PLUS(rhoVT_dd), 0.5)*(model.qr[i][k-1] - model.qr[i][k-2]) - 
                                           MINUS(rhoVT_d) * (model.qr[i][k] - model.qr[i][k-1]) - std::pow(MINUS(rhoVT_d), 0.5)*std::pow(MINUS(rhoVT_u), 0.5)*(model.qr[i][k+1] - model.qr[i][k]));
                }
            #endif
            prhoVTqr_pz_rho = (G_ik_u - G_ik_d) * model.r2dz / model.rhou[k];

            #if defined(AB3)
                #if defined(VTconst)
                    G_ik_u = model.rhow[k+1] * VT_u * (model.qrm[i][k+1] + model.qrm[i][k]);
                    G_ik_d = model.rhow[k] * VT_d * (model.qrm[i][k] + model.qrm[i][k-1]);
                #else
                    rhoVT_u = model.rhow[k+1] * 1E-2 * (3634 * pow(1E-3*model.rhow[k+1] * 0.5*(model.qrm[i][k+1]+model.qrm[i][k]), 0.1346) * pow(model.rhow[k+1]/model.rhow[1], -0.5));
                    rhoVT_d = model.rhow[k] * 1E-2 * (3634 * pow(1E-3*model.rhow[k] * 0.5*(model.qrm[i][k]+model.qrm[i][k-1]), 0.1346) * pow(model.rhow[k]/model.rhow[1], -0.5));
                    G_ik_u = rhoVT_u * (model.qrm[i][k+1] + model.qrm[i][k]);
                    G_ik_d = rhoVT_d * (model.qrm[i][k] + model.qrm[i][k-1]);
                #endif

                if (i != 1 && i != model.nx-2 && k != 1 && k != model.nz-2 && ALPHA != 0.) {
                    #if defined(VTconst)
                        rhoVT_u = model.rhow[k+1] * VT_u;
                        rhoVT_d = model.rhow[k] * VT_d;
                    #else
                        rhoVT_uu = model.rhow[k+2] * 1E-2 * (36.34 * pow(1E-3*model.rhow[k+2] * 0.5*(model.qrm[i][k+2]+model.qrm[i][k+1]), 0.1346) * pow(model.rhow[k+2]/model.rhow[1], -0.5));
                        rhoVT_dd = model.rhow[k-1] * 1E-2 * (36.34 * pow(1E-3*model.rhow[k-1] * 0.5*(model.qrm[i][k-1]+model.qrm[i][k-2]), 0.1346) * pow(model.rhow[k-1]/model.rhow[1], -0.5));
                    #endif
                    G_ik_u += -ALPHA/3. * (PLUS(rhoVT_u) * (model.qrm[i][k+1] - model.qrm[i][k]) - std::pow(PLUS(rhoVT_u), 0.5)*std::pow(PLUS(rhoVT_d), 0.5)*(model.qrm[i][k] - model.qrm[i][k-1]) - 
                                           MINUS(rhoVT_u) * (model.qrm[i][k+1] - model.qrm[i][k]) - std::pow(MINUS(rhoVT_u), 0.5)*std::pow(MINUS(rhoVT_uu), 0.5)*(model.qrm[i][k+2] - model.qrm[i][k+1]));
                    G_ik_d += -ALPHA/3. * (PLUS(rhoVT_d) * (model.qrm[i][k] - model.qrm[i][k-1]) - std::pow(PLUS(rhoVT_d), 0.5)*std::pow(PLUS(rhoVT_dd), 0.5)*(model.qrm[i][k-1] - model.qrm[i][k-2]) - 
                                           MINUS(rhoVT_d) * (model.qrm[i][k] - model.qrm[i][k-1]) - std::pow(MINUS(rhoVT_d), 0.5)*std::pow(MINUS(rhoVT_u), 0.5)*(model.qrm[i][k+1] - model.qrm[i][k]));
                }
                prhoVTqr_pz_rho_m = (G_ik_u - G_ik_d) * model.r2dz / model.rhou[k];
            #endif

            #if defined(AB3)
                model.qrp[i][k] += 1.5*DT*(prhoVTqr_pz_rho) - 0.5*DT*(prhoVTqr_pz_rho_m);
            #else
                model.qrp[i][k] = model.qrm[i][k] + model.d2t * (prhoVTqr_pz_rho);
            #endif
        }
    }
}
#endif
