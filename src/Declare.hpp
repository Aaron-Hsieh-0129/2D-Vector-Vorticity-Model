#include <iostream>
#include <cmath>
#include "Const.hpp"

class vvmArray {
	public:
			// constructor
			vvmArray();

			double tb[nz], tb_zeta[nz], rhou[nz], rhow[nz], pib[nz], qvb[nz], qvsb[nz], tvb[nz], pb[nz];
			double zetap[nx][nz], zeta[nx][nz], zetam[nx][nz];
			double thp[nx][nz], th[nx][nz], thm[nx][nz];
			double qvp[nx][nz], qv[nx][nz], qvm[nx][nz];
			double qcp[nx][nz], qc[nx][nz], qcm[nx][nz];
			double qrp[nx][nz], qr[nx][nz], qrm[nx][nz];
			double u[nx][nz];
			double w[nx][nz];
			double xi[nx], uxi[nx];
			double addflx[nx];

			static void BoundaryProcess(double tmp[][nz]);
			static void BoundaryProcessZETA(double tmp[][nz]);
			static void BoundaryProcessDouble(double tmp[][nz]);
};