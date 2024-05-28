#include "Declare.hpp"

void vvm::BoundaryProcess1D_center(double *var, int nz) {
    var[0] = var[1];
    var[nz-1] = var[nz-2];
    return;
}

void vvm::BoundaryProcess2D_center(double **var, int nx, int nz) {
    for (int k = 1; k <= nz-2; k++) {
        var[0][k] = var[nx-2][k];
        var[nx-1][k] = var[1][k];
    }
    for (int i = 0; i <= nx-1; i++) {
        var[i][0] = var[i][1];
        var[i][nz-1] = var[i][nz-2];
    }
    return;
}

void vvm::BoundaryProcess2D_westdown(double **var, int nx, int nz) {
    for (int k = 1; k <= nz-2; k++) {
        var[0][k] = var[nx-2][k];
        var[nx-1][k] = var[1][k];
    }
    for (int i = 0; i <= nx-1; i++) {
        var[i][0] = var[i][1] = var[i][nz-1] = 0.;
    }
    return;
}

void vvm::BoundaryProcess2D_all(vvm &model) {
    model.BoundaryProcess2D_westdown(model.zetap, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.thp, model.nx, model.nz);
    model.BoundaryProcess2D_westdown(model.w, model.nx, model.nz);
    model.BoundaryProcess2D_center(model.u, model.nx, model.nz);
    #if defined(WATER)
        model.BoundaryProcess2D_center(model.qvp, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.qcp, model.nx, model.nz);
        model.BoundaryProcess2D_center(model.qrp, model.nx, model.nz);
    #endif
}
