#include <iostream>
#include <cmath>
#include "../include/Eigen/Sparse"
#include "Define.hpp"

typedef Eigen::Triplet<double> T;
class vvmArray {
	public:
		// constructor
		vvmArray();

		double rdx = 1. / (double) dx, r2dx = 1. / (2. * (double) dx);
		double rdz = 1. / (double) dz, r2dz = 1. / (2. * (double) dz);
		double rdx2 = rdx * rdx, rdz2 = rdz * rdz;
		int nx = NX, nz = NZ;


		double tb[NZ], tb_zeta[NZ], rhou[NZ], rhow[NZ], pib[NZ], qvb[NZ], qvsb[NZ], tvb[NZ], pb[NZ];
		#if defined(TROPICALFORCING)
			double Q1LS[NZ], Q2LS[NZ];
			double init_th_forcing[NX][NZ];
			bool status_for_adding_forcing = false;
		#endif
		double zetap[NX][NZ], zeta[NX][NZ], zetam[NX][NZ];
		double thp[NX][NZ], th[NX][NZ], thm[NX][NZ];
		double qvp[NX][NZ], qv[NX][NZ], qvm[NX][NZ];
		double qcp[NX][NZ], qc[NX][NZ], qcm[NX][NZ];
		double qrp[NX][NZ], qr[NX][NZ], qrm[NX][NZ];
		double u[NX][NZ], w[NX][NZ];
		#if defined(STREAMFUNCTION)
			double psi[NX][NZ];
		#endif

		#if defined(WATER)
			double evaporation[NX][NZ], accretion[NX][NZ], autoconversion[NX][NZ];
		#endif

		double ubarTopp = 0., ubarTop = 0., ubarTopm = 0.;
		double xi[NX], uxi[NX];
		double addflx[NX];
		double qrAcc[NX];

		// Metrices to solve the Poisson equation
		Eigen::SparseMatrix<double> A;
		Eigen::SparseMatrix<double> G;

		static void BoundaryProcess(double tmp[][NZ]);
		static void BoundaryProcessZETA(double tmp[][NZ]);
		static void BoundaryProcessDouble(double tmp[][NZ]);
};
