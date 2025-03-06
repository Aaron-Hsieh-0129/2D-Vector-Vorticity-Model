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
// TODO: Fix the blow up in AB2.
double PLUS(double var) {
    return 0.5 * (var + std::fabs(var));
}

double MINU(double var) {
    return 0.5 * (var - std::fabs(var));
}

/*
void vvm::Advection_thermo(double **past, double **now, double **future, double ***dvar, vvm &model) {
    double prhouvar_px_rho = 0., prhowvar_pz_rho = 0.;
    double *flux_ucont, *flux_wcont;
    double **flux_u = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_ucont);
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);

    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-1; k++) {
        for (int i = 1; i <= model.nx-1; i++) {
            flux_u[i][k] = model.rhou[k] * model.u[i][k] * (now[i][k] + now[i-1][k]);
            flux_w[i][k] = model.rhow[k] * model.w[i][k] * (now[i][k] + now[i][k-1]);
            #if defined(AB2)
                if (i >= 2 && i <= model.nx-3 && k >= 2 && k <= model.nz-3) {
                    flux_u[i][k] += -1./3. * 
                                    (model.rhou[k]*PLUS(model.u[i][k]) * (now[i][k] - now[i-1][k])
                                        - std::sqrt(model.rhou[k]*PLUS(model.u[i][k]) * model.rhou[k]*PLUS(model.u[i-1][k]))*(now[i-1][k] - now[i-2][k])
                                   - model.rhou[k]*MINU(model.u[i][k]) * (now[i][k] - now[i-1][k])
                                        - std::sqrt(std::fabs(model.rhou[k]*MINU(model.u[i][k]) * model.rhou[k]*MINU(model.u[i+1][k])))*(now[i+1][k] - now[i][k]));

                    flux_w[i][k] += -1./3. * 
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
    // for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][1] = flux_w[i][model.nz-1] = 0.;
    for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][model.nz-1] = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhouvar_px_rho = (flux_u[i+1][k] - flux_u[i][k]) * model.r2dx / model.rhou[k];
            prhowvar_pz_rho = (flux_w[i][k+1] - flux_w[i][k]) * model.r2dz / model.rhou[k];

            #if defined(AB2)
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
*/

#include <iostream>

void vvm::Advection_thermo(double **past, double **now, double **future, double ***dvar, vvm &model) {
    double prhouvar_px_rho = 0., prhowvar_pz_rho = 0.;
    double *flux_ucont, *flux_wcont;
    double **flux_u = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_ucont);
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);

    // Local variables to avoid struct access in OpenACC
    double* rhou = model.rhou;  // 1D array
    double* rhow = model.rhow;  // 1D array
    double* ucont = model.ucont;  // Assuming contiguous version of u
    double* wcont = model.wcont;  // Assuming contiguous version of w
    double* pastcont = past[0];   // Assuming past is contiguous
    double* nowcont = now[0];     // Assuming now is contiguous
    double* futurecont = future[0];  // Assuming future is contiguous
    double* dvarcont = dvar[0][0];   // Assuming dvar is contiguous across all dimensions
    double r2dx = model.r2dx, r2dz = model.r2dz, dt = model.dt, d2t = model.d2t;
    int nx = model.nx, nz = model.nz, step = model.step;

    #pragma acc data create(flux_ucont[0:nx * nz], flux_wcont[0:nx * nz]) \
                     copyin(rhou[0:nz], rhow[0:nz], ucont[0:nx * nz], wcont[0:nx * nz], \
                            pastcont[0:nx * nz], nowcont[0:nx * nz]) \
                     copy(futurecont[0:nx * nz], dvarcont[0:nx * nz * 2]) \
                     copyin(r2dx, r2dz, dt, d2t, step)
    {

        // First loop: Compute fluxes
        #pragma acc parallel loop collapse(2) private(prhouvar_px_rho, prhowvar_pz_rho)
        for (int k = 1; k <= nz-1; k++) {
            for (int i = 1; i <= nx-1; i++) {
                int idx = i * nz + k;
                int idx_im1 = (i-1) * nz + k;
                int idx_km1 = i * nz + (k-1);

                flux_ucont[idx] = rhou[k] * ucont[idx] * (nowcont[idx] + nowcont[idx_im1]);
                flux_wcont[idx] = rhow[k] * wcont[idx] * (nowcont[idx] + nowcont[idx_km1]);

                #if defined(AB2)
                if (i >= 2 && i <= nx-3 && k >= 2 && k <= nz-3) {
                    int idx_im2 = (i-2) * nz + k;
                    int idx_ip1 = (i+1) * nz + k;
                    int idx_km2 = i * nz + (k-2);
                    int idx_kp1 = i * nz + (k+1);

                    double u_plus = (ucont[idx] > 0 ? ucont[idx] : 0.);
                    double u_minus = (ucont[idx] < 0 ? ucont[idx] : 0.);
                    double w_plus = (wcont[idx] > 0 ? wcont[idx] : 0.);
                    double w_minus = (wcont[idx] < 0 ? wcont[idx] : 0.);

                    flux_ucont[idx] += -1./3. * (
                        rhou[k] * u_plus * (nowcont[idx] - nowcont[idx_im1])
                        - std::sqrt(rhou[k] * u_plus * rhou[k] * (ucont[idx_im1] > 0 ? ucont[idx_im1] : 0.)) * (nowcont[idx_im1] - nowcont[idx_im2])
                        - rhou[k] * u_minus * (nowcont[idx] - nowcont[idx_im1])
                        - std::sqrt(std::fabs(rhou[k] * u_minus * rhou[k] * (ucont[idx_ip1] < 0 ? ucont[idx_ip1] : 0.))) * (nowcont[idx_ip1] - nowcont[idx])
                    );

                    flux_wcont[idx] += -1./3. * (
                        rhow[k] * w_plus * (nowcont[idx] - nowcont[idx_km1])
                        - std::sqrt(rhow[k] * w_plus * rhow[k-1] * (wcont[idx_km1] > 0 ? wcont[idx_km1] : 0.)) * (nowcont[idx_km1] - nowcont[idx_km2])
                        - rhow[k] * w_minus * (nowcont[idx] - nowcont[idx_km1])
                        - std::sqrt(std::fabs(rhow[k] * w_minus * rhow[k+1] * (wcont[idx_kp1] < 0 ? wcont[idx_kp1] : 0.))) * (nowcont[idx_kp1] - nowcont[idx])
                    );
                }
                #endif
            }
        }

        // Boundary conditions (simplified for GPU)
        #pragma acc parallel loop
        for (int i = 0; i < nx; i++) {
            int idx_bot = i * nz + 0;
            int idx_top = i * nz + (nz-1);
            flux_wcont[idx_bot] = 0.;
            flux_wcont[idx_top] = 0.;
        }

        // Second loop: Update future and dvar
        #pragma acc parallel loop collapse(2) private(prhouvar_px_rho, prhowvar_pz_rho)
        for (int k = 1; k <= nz-2; k++) {
            for (int i = 1; i <= nx-2; i++) {
                int idx = i * nz + k;
                int idx_ip1 = (i+1) * nz + k;
                int idx_kp1 = i * nz + (k+1);

                prhouvar_px_rho = (flux_ucont[idx_ip1] - flux_ucont[idx]) * r2dx / rhou[k];
                prhowvar_pz_rho = (flux_wcont[idx_kp1] - flux_wcont[idx]) * r2dz / rhou[k];

                #if defined(AB2)
                    int dvar_idx = (i * nz + k) * 2 + ((step + 1) % 2);
                    dvarcont[dvar_idx] = -prhouvar_px_rho - prhowvar_pz_rho;
                    if (step == 0) dvarcont[(i * nz + k) * 2] = dvarcont[(i * nz + k) * 2 + 1];
                    futurecont[idx] = nowcont[idx] + 1.5 * dt * dvarcont[dvar_idx] - 0.5 * dt * dvarcont[(i * nz + k) * 2 + (step % 2)];
                #else
                    futurecont[idx] = pastcont[idx] + d2t * (-prhouvar_px_rho - prhowvar_pz_rho);
                #endif
            }
        }
    }

    vvm::deallocate2DContinuousArray(flux_u, flux_ucont);
    vvm::deallocate2DContinuousArray(flux_w, flux_wcont);
}


void vvm::Advection_zeta(vvm &model) {
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
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

    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 0; i <= model.nx-2; i++) {
            flux_u[i][k] = model.U_w[i][k] * (model.zeta[i+1][k] + model.zeta[i][k]);
            flux_w[i][k] = model.W_u[i][k] * (model.zeta[i][k+1] + model.zeta[i][k]);
            #if defined(AB2)
                if (i >= 2 && i <= model.nx-3 && k >= 2 && k <= model.nz-3) {
                    flux_u[i][k] += -1./3. * (PLUS(model.U_w[i][k])*(model.zeta[i+1][k]-model.zeta[i][k]) 
                                                    - std::sqrt(PLUS(model.U_w[i][k]) * PLUS(model.U_w[i-1][k]))*(model.zeta[i][k]-model.zeta[i-1][k]) - 
                                                 MINU(model.U_w[i][k])*(model.zeta[i+1][k]-model.zeta[i][k]) 
                                                    - std::sqrt(std::fabs(MINU(model.U_w[i][k]) * MINU(model.U_w[i+1][k])))*(model.zeta[i+2][k]-model.zeta[i+1][k]));
                    flux_w[i][k] += -1./3. * (PLUS(model.W_u[i][k])*(model.zeta[i][k+1]-model.zeta[i][k]) 
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
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 2; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhouzeta_px_rho = (flux_u[i][k] - flux_u[i-1][k]) * model.r2dx / model.rhow[k];
            prhowzeta_pz_rho = (flux_w[i][k] - flux_w[i][k-1]) * model.r2dz / model.rhow[k];

            #if defined(AB2)
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
#if defined(KESSLER_MICROPHY)
void vvm::Advection_qrVT(vvm &model) {
    double prhoVTqr_pz_rho = 0.;
    double *flux_wcont;
    double **flux_w = vvm::allocate2DContinuousArray(model.nx, model.nz, flux_wcont);


    double VT = 0.;
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            VT = 1E-2 * (3634 * pow(1E-3*model.rhow[k] * 0.5*(model.qr[i][k]+model.qr[i][k-1]), 0.1346) * pow(model.rhow[k]/model.rhow[1], -0.5));
            flux_w[i][k] = model.rhow[k] * VT * (model.qr[i][k] + model.qr[i][k-1]);
        }
    }
    model.BoundaryProcess2D_center(flux_w, model.nx, model.nz);
    // for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][1] = flux_w[i][model.nz-1] = 0.;
    for (int i = 0; i < model.nx; i++) flux_w[i][0] = flux_w[i][model.nz-1] = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            prhoVTqr_pz_rho = (flux_w[i][k+1] - flux_w[i][k]) * model.r2dz / model.rhou[k];

            #if defined(AB2)
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
#endif

