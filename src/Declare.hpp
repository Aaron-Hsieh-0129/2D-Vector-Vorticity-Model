#include <iostream>
#include <cmath>
#include "../include/Eigen/Sparse"
#include "Define.hpp"

typedef Eigen::Triplet<double> T;
class vvmArray {
	public:
		// constructor
		vvmArray();

		double tb[nz], tb_zeta[nz], rhou[nz], rhow[nz], pib[nz], qvb[nz], qvsb[nz], tvb[nz], pb[nz];
		#if defined(TROPICALFORCING)
			double Q1LS[nz], Q2LS[nz];
			double init_th_forcing[nx][nz];
			bool status_for_adding_forcing = false;
		#endif
		double zetap[nx][nz], zeta[nx][nz], zetam[nx][nz];
		double thp[nx][nz], th[nx][nz], thm[nx][nz];
		double qvp[nx][nz], qv[nx][nz], qvm[nx][nz];
		double qcp[nx][nz], qc[nx][nz], qcm[nx][nz];
		double qrp[nx][nz], qr[nx][nz], qrm[nx][nz];
		double u[nx][nz];
		#if defined(WATER)
			double evaporation[nx][nz], accretion[nx][nz], autoconversion[nx][nz];
		#endif

		double ubarTopp = 0., ubarTop = 0., ubarTopm = 0.;
		double w[nx][nz];
		double xi[nx], uxi[nx];
		double addflx[nx];
		double qrAcc[nx];

		// Metrices to solve the Poisson equation
		Eigen::SparseMatrix<double> A;
		Eigen::SparseMatrix<double> G;

		static void BoundaryProcess(double tmp[][nz]);
		static void BoundaryProcessZETA(double tmp[][nz]);
		static void BoundaryProcessDouble(double tmp[][nz]);
};
