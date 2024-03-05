#include "Declare.hpp"

vvmArray::vvmArray() {
    for (int k = 0; k <= nz - 1; k++) {
        for (int i = 0; i <= nx - 1; i++) {
            zetap[i][k] = zeta[i][k] = zetam[i][k] = 0.;
            thp[i][k] = th[i][k] = thm[i][k] = 0.;
            qvp[i][k] = qv[i][k] = qvm[i][k] = 0.;
            qcp[i][k] = qc[i][k] = qcm[i][k] = 0.;
            qrp[i][k] = qr[i][k] = qrm[i][k] = 0.;
            w[i][k] = 0.;
            u[i][k] = 0.;
            #if defined(TROPICALFORCING)
                init_th_forcing[i][k] = 0.;
            #endif
            #if defined(WATER)
                evaporation[i][k] = 0.;
                accretion[i][k] = 0.;
                autoconversion[i][k] = 0.;
            #endif
        }
        tb[k] = 0.;
        tb_zeta[k] = 0.;
        rhou[k] = 0.;
        rhow[k] = 0.;
        pib[k] = 0.;
        qvb[k] = 0.;
        qvsb[k] = 0.;
        tvb[k] = 0.;
        pb[k] = 0.;
    }

    for (int i = 0; i <= nx - 1; i++) {
        xi[i] = 0.;
        uxi[i] = 0.;
        addflx[i] = 1.;
        qrAcc[i] = 0.;
    }

    A = Eigen::SparseMatrix<double>((nx-2)*(nz-3), (nx-2)*(nz-3));
    G = Eigen::SparseMatrix<double>(nx-2, nx-2);
    return;
}

void vvmArray::BoundaryProcess(double tmp[][nz]) {
  for (int k = 1; k <= nz - 2; k++) {
    tmp[0][k] = tmp[nx - 2][k];
    tmp[nx - 1][k] = tmp[1][k];
  }
  for (int i = 0; i <= nx - 1; i++) {
    tmp[i][0] = tmp[i][1];
    tmp[i][nz - 1] = tmp[i][nz - 2];
  }
}

void vvmArray::BoundaryProcessZETA(double tmp[][nz]) {
    for (int k = 1; k <= nz - 2; k++) {
        tmp[0][k] = tmp[nx - 2][k];
        tmp[nx - 1][k] = tmp[1][k];
    }
    for (int i = 0; i <= nx - 1; i++) {
        tmp[i][0] = tmp[i][1];
        tmp[i][nz - 1] = 0.;
    }
}

void vvmArray::BoundaryProcessDouble(double tmp[][nz]) {
    for (int k = 1; k <= nz - 2; k++) {
        tmp[0][k] = tmp[nx - 2][k];
        tmp[nx - 1][k] = tmp[1][k];
    }
    for (int i = 0; i <= nx - 1; i++) {
        tmp[i][0] = tmp[i][nz - 2];
        tmp[i][nz - 1] = tmp[i][1];
    }
}
