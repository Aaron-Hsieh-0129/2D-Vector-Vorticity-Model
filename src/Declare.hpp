#include <iostream>
#include <cmath>
#include "Eigen/Sparse"
#include "Define.hpp"

typedef Eigen::Triplet<long double> T;
class vvmArray {
	public:
		// constructor
		vvmArray();

		long double tb[nz], tb_zeta[nz], rhou[nz], rhow[nz], pib[nz], qvb[nz], qvsb[nz], tvb[nz], pb[nz];
		#if defined(TROPICALFORCING)
			long double Q1LS[nz], Q2LS[nz];
			long double init_th_forcing[nx][nz];
			bool status_for_adding_forcing = false;
		#endif
		long double zetap[nx][nz], zeta[nx][nz], zetam[nx][nz];
		long double thp[nx][nz], th[nx][nz], thm[nx][nz];
		long double qvp[nx][nz], qv[nx][nz], qvm[nx][nz];
		long double qcp[nx][nz], qc[nx][nz], qcm[nx][nz];
		long double qrp[nx][nz], qr[nx][nz], qrm[nx][nz];
		long double u[nx][nz];

		long double ubarTopp = 0., ubarTop = 0., ubarTopm = 0.;
		long double w[nx][nz];
		long double xi[nx], uxi[nx];
		long double addflx[nx];
		long double qrAcc[nx];

		// Metrices to solve the Poisson equation
		Eigen::SparseMatrix<long double> A;
		Eigen::SparseMatrix<long double> G;

		static void BoundaryProcess(long double tmp[][nz]);
		static void BoundaryProcessZETA(long double tmp[][nz]);
		static void BoundaryProcessDouble(long double tmp[][nz]);
};
