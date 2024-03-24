#include "Declare.hpp"

void vvm::Diffusion(double var_in[][NZ], double var_out[][NZ], vvm &model) {
    for (int k = 1; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            var_out[i][k] += model.d2t * Kx * model.rdx2 * (var_in[i+1][k] - 2. * var_in[i][k] + var_in[i-1][k]) + 
                             model.d2t * Kz * model.rdz2 * (var_in[i][k+1] - 2. * var_in[i][k] + var_in[i][k-1]);
        }
    }
    return;
}

void vvm::TimeFilter(double previous[][NZ], double now[][NZ], double future[][NZ], vvm &model) {
    for (int i = 0; i <= model.nx-1; i++) {
        for (int k = 0; k <= model.nz-1; k++) {
            now[i][k] += TIMETS * (future[i][k] - 2.*now[i][k] + previous[i][k]);
        }
    }
    return;
}
