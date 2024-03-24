#include "Declare.hpp"


vvm::vvm() {
    // 0D variables
    ubarTopp = ubarTop = ubarTopm = 0.;

    // 1D variables
    for (int k = 0; k < nz; k++) {
        thb[k] = thb_zeta[k] = 0.;
        rhou[k] = rhow[k] = 0.;
        pib[k] = 0.;
        tvb[k] = 0.;
        pb[k] = pib[k] = 0.;
        qvb[k] = qvsb[k] = 0.;

        #if defined(TROPICALFORCING)
            Q1LS[k] = Q2LS[k] = 0.;
        #endif
    }
    for (int i = 0; i < NX; i++) {
        xi[i] = 0.;
        uxi[i] = 0.;
        #if defined(WATER)
            qrAcc[i] = 0.;
        #endif
    }

    // 2D variables
    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            zetap[i][k] = zeta[i][k] = zetam[i][k] = 0.;
            thp[i][k] = th[i][k] = thm[i][k] = 0.;
            u[i][k] = um[i][k] = 0.;
            w[i][k] = wm[i][k] = 0.;
            #if defined(STREAMFUNCTION)
                psi[i][k] = 0.;
            #endif

            #if defined(WATER)
                qvp[i][k] = qv[i][k] = qvm[i][k] = 0.;
                qcp[i][k] = qc[i][k] = qcm[i][k] = 0.;
                qrp[i][k] = qr[i][k] = qrm[i][k] = 0.;
                evaporation[i][k] = 0.;
                accretion[i][k] = 0.;
                autoconversion[i][k] = 0.;
            #endif

            #if defined(TROPICALFORCING)
                init_th_forcing[i][k] = 0.;
            #endif
        }
    }

    #ifndef PETSC
        // Metrices to solve the Poisson equations
        A = Eigen::SparseMatrix<double>((nx-2)*(nz-3), (nx-2)*(nz-3));
        G = Eigen::SparseMatrix<double>(nx-2, nx-2);
    #endif
    return;
}

