#include "Declare.hpp"
#include <algorithm>
#include <iostream>

void vvm::Turbulence::RKM_RKH(vvm &model) {
    double Rzeta = 0.;
    double Rotat = 0.;
    double Ri = 0.;

    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            Rzeta = 0.5*((model.w[i+1][k+1] + model.w[i+1][k]) - (model.w[i-1][k+1] + model.w[i-1][k])) * model.r2dx + 
                    0.5*((model.u[i+1][k+1] + model.u[i][k+1]) - (model.u[i+1][k-1] + model.u[i][k-1])) * model.r2dz;

            Rotat = std::pow((model.u[i+1][k] - model.u[i][k]) * model.rdx, 2) + std::pow((model.w[i][k+1] - model.w[i][k]) * model.rdz, 2);
            
            Ri = (model.GRAVITY / model.thb[k] * (model.thp[i][k+1] - model.thp[i][k-1]) * model.r2dz) / (std::pow(Rzeta, 2) + 2. * Rotat);

            if (Ri < 0) {
                model.RKM[i][k] = model.lambda2[k] * std::sqrt(std::pow(Rzeta, 2) + 2. * Rotat) * std::sqrt(1. - 16. * Ri);
                model.RKH[i][k] = model.lambda2[k] * std::sqrt(std::pow(Rzeta, 2) + 2. * Rotat) * 1.4 * std::sqrt(1. - 40.*Ri);
            }
            else if (0 < Ri && Ri < 0.25) {
                model.RKM[i][k] = model.lambda2[k] * std::sqrt(std::pow(Rzeta, 2) + 2. * Rotat) * std::pow(1. - 4.*Ri, 4);
                model.RKH[i][k] = model.lambda2[k] * std::sqrt(std::pow(Rzeta, 2) + 2. * Rotat) * 1.4 * (1.-1.2*Ri)*std::pow(1.-4.*Ri, 4);
            }
            else {
                model.RKM[i][k] = 0.;
                model.RKH[i][k] = 0.;
            }

            // Diffusion should be larger than 1
            model.RKM[i][k] = std::max(model.RKM[i][k], 1.);
            model.RKH[i][k] = std::max(model.RKH[i][k], 1.);
            
            // Diffusion should be smaller than 0.8 * dx^2 / dt
            model.RKM[i][k] = std::min(model.RKM[i][k], 0.8 * model.dx * model.dz / model.dt);
            model.RKH[i][k] = std::min(model.RKH[i][k], 0.8 * model.dx * model.dz / model.dt);
        }
    }
    model.BoundaryProcess2D_center(model.RKM, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.RKH, model.nx, model.nz);

    vvm::Turbulence::Mparam(model, model.zeta, model.zetap);
    vvm::Turbulence::Hparam(model, model.th, model.thp);
    #if defined(WATER)
        vvm::Turbulence::Hparam(model, model.qv, model.qvp);
        vvm::Turbulence::Hparam(model, model.qc, model.qcp);
        vvm::Turbulence::Hparam(model, model.qr, model.qrp);
    #endif
    // std::cout << "max RKM: " << *std::max_element(model.RKMcont, model.RKMcont + model.nx*model.nz) << " max RKH: " << *std::max_element(model.RKHcont, model.RKHcont + model.nx*model.nz) << std::endl;
    // std::cout << "min RKM: " << *std::min_element(model.RKMcont, model.RKMcont + model.nx*model.nz) << " min RKH: " << *std::min_element(model.RKHcont, model.RKHcont + model.nx*model.nz) << std::endl;
    return;
}

void vvm::Turbulence::Mparam(vvm &model, double **var_now, double **var_future) {
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            // var_future[i][k] += 1. / std::pow(model.rhow[k], 2) * model.rdx2 * model.dt * 
            //                     (model.rhow[k] * 0.5 * (model.RKM[i][k] + model.RKM[i][k-1]) * (model.rhow[k]*var_now[i+1][k] - model.rhow[k]*var_now[i][k]) - 
            //                      model.rhow[k] * 0.5 * (model.RKM[i-1][k] + model.RKM[i-1][k-1]) * (model.rhow[k]*var_now[i][k] - model.rhow[k]*var_now[i-1][k]))
            //                   + 1. / std::pow(model.rhow[k], 2) * model.rdz2 * model.dt * 
            //                     (model.rhou[k] * 0.5 * (model.RKM[i][k] + model.RKM[i-1][k]) * (model.rhow[k+1]*var_now[i][k+1] - model.rhow[k]*var_now[i][k]) - 
            //                      model.rhou[k-1] * 0.5 * (model.RKM[i][k-1] + model.RKM[i-1][k-1]) * (model.rhow[k]*var_now[i][k] - model.rhow[k-1]*var_now[i][k-1]));
            
            var_future[i][k] += model.rdx2 * model.dt * 
                                (0.5 * (model.RKM[i][k] + model.RKM[i][k-1]) * (var_now[i+1][k] - var_now[i][k]) - 
                                 0.5 * (model.RKM[i-1][k] + model.RKM[i-1][k-1]) * (var_now[i][k] - var_now[i-1][k]))
                              + model.rdz2 * model.dt * 
                                (0.5 * (model.RKM[i][k] + model.RKM[i-1][k]) * (var_now[i][k+1] - var_now[i][k]) - 
                                 0.5 * (model.RKM[i][k-1] + model.RKM[i-1][k-1]) * (var_now[i][k] - var_now[i][k-1]));
        }
    }
    return;
}

void vvm::Turbulence::Hparam(vvm &model, double **var_now, double **var_future) {
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            // var_future[i][k] += 1. / std::pow(model.rhou[k], 2) * model.rdx2 * model.dt * 
            //                     (model.rhou[k] * 0.5 * (model.RKH[i+1][k] + model.RKH[i][k]) * (model.rhou[k]*var_now[i+1][k] - model.rhou[k]*var_now[i][k]) - 
            //                      model.rhou[k] * 0.5 * (model.RKH[i][k] + model.RKH[i-1][k]) * (model.rhou[k]*var_now[i][k] - model.rhou[k]*var_now[i-1][k]))
            //                   + 1. / std::pow(model.rhou[k], 2) * model.rdz2 * model.dt * 
            //                     (model.rhow[k+1] * 0.5 * (model.RKH[i][k+1] + model.RKH[i][k]) * (model.rhou[k+1]*var_now[i][k+1] - model.rhou[k]*var_now[i][k]) - 
            //                      model.rhow[k] * 0.5 * (model.RKH[i][k] + model.RKH[i][k-1]) * (model.rhou[k]*var_now[i][k] - model.rhou[k-1]*var_now[i][k-1]));

            var_future[i][k] += model.rdx2 * model.dt * 
                                (0.5 * (model.RKH[i+1][k] + model.RKH[i][k]) * (var_now[i+1][k] - var_now[i][k]) - 
                                 0.5 * (model.RKH[i][k] + model.RKH[i-1][k]) * (var_now[i][k] - var_now[i-1][k]))
                              + model.rdz2 * model.dt * 
                                (0.5 * (model.RKH[i][k+1] + model.RKH[i][k]) * (var_now[i][k+1] - var_now[i][k]) - 
                                 0.5 * (model.RKH[i][k] + model.RKH[i][k-1]) * (var_now[i][k] - var_now[i][k-1]));
        }
    }
    return;
}

