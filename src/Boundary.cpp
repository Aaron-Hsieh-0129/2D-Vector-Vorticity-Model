#include "Declare.hpp"

void vvm::BoundaryProcess1D_center(double var[]) {
    var[0] = var[1];
    var[NZ-1] = var[NZ-2];
    return;
}

void vvm::BoundaryProcess2D_center(double var[][NZ]) {
    for (int k = 1; k <= NZ-2; k++) {
        var[0][k] = var[NX-2][k];
        var[NX-1][k] = var[1][k];
    }
    for (int i = 0; i <= NX-1; i++) {
        var[i][0] = var[i][1];
        var[i][NZ-1] = var[i][NZ-2];
    }
    return;
}

void vvm::BoundaryProcess2D_westdown(double var[][NZ]) {
    for (int k = 1; k <= NZ-2; k++) {
        var[0][k] = var[NX-2][k];
        var[NX-1][k] = var[1][k];
    }
    for (int i = 0; i <= NX-1; i++) {
        var[i][0] = var[i][1] = var[i][NZ-1] = 0.;
    }
    return;
}
