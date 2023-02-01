#include "Outputfile.hpp"

using namespace netCDF;
using std::vector;

void Output::printInit(vvmArray & myArray) {
	double z;
	std::cout << "z          tb        rhou       rhow       qvb	 RH      pib" << std::endl;
	for (int k = 0; k <= nz-1;k++){
		z = (double) (k - 0.5) * dz ;
		std::cout << std::fixed << std::setprecision(2) << z << "    " << myArray.tb[k] << "    " << myArray.rhou[k] << "     " 
		<< myArray.rhow[k] << "   	 " << myArray.qvb[k] * 1000 << "    " << myArray.qvb[k] / myArray.qvsb[k] << "    "
		<< myArray.pib[k] << std::endl;
	}
	std::fstream initout;
	std::string initName = "../outputs/init.txt";
	initout.open(initName, std::ios::out);
	for (int k = 0; k <= nz-1; k++) {
		z = (double) (k - 0.5) * dz ;
		initout << z << "    " << myArray.tb[k] << "    " << myArray.rhou[k] << "     " 
		<< myArray.rhow[k] << "   	 " << myArray.qvb[k] << "    " << myArray.qvsb[k] << "    " << myArray.qvb[k] / myArray.qvsb[k] << "    "
		<< myArray.pib[k] << std::endl;
	}
	return;
}


void Output::output_zeta(int n, vvmArray & myArray) {
	std::fstream foutzeta;
	std::string zetaName = "../outputs/zeta/zeta_" + std::to_string(n) + ".txt";
	foutzeta.open(zetaName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutzeta << myArray.zeta[i][k] << " ";
		}
	}
}

void Output::output_th(int n, vvmArray & myArray) {
	std::fstream foutth;
	std::string thName = "../outputs/th/th_" + std::to_string(n) + ".txt";
	foutth.open(thName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutth << myArray.th[i][k] << " ";
		}
	}
}

void Output::output_u(int n, vvmArray & myArray) {
	std::fstream foutu;
	std::string uName = "../outputs/u/u_" + std::to_string(n) + ".txt";
	foutu.open(uName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutu << myArray.u[i][k] << " ";
		}
	}
}

void Output::output_w(int n, vvmArray & myArray) {
	std::fstream foutw;
	std::string wName = "../outputs/w/w_" + std::to_string(n) + ".txt";
	foutw.open(wName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutw << myArray.w[i][k] << " ";
		}
	}
}

void Output::output_qv(int n, vvmArray & myArray) {
	std::fstream foutqv;
	std::string qvName = "../outputs/qv/qv_" + std::to_string(n) + ".txt";
	foutqv.open(qvName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqv << myArray.qv[i][k] + myArray.qvb[k] << " ";
		}
	}
}

void Output::output_qc(int n, vvmArray & myArray) {
	std::fstream foutqc;
	std::string qcName = "../outputs/qc/qc_" + std::to_string(n) + ".txt";
	foutqc.open(qcName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqc << myArray.qc[i][k] << " ";
		}
	}
}

void Output::output_qr(int n, vvmArray & myArray) {
	std::fstream foutqr;
	std::string qrName = "../outputs/qr/qr_" + std::to_string(n) + ".txt";
	foutqr.open(qrName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqr << myArray.qr[i][k] << " ";
		}
	}
}

void Output::output_nc(int n, vvmArray &myArray) {
    std::string ncName = "../outputs/nc/" + std::to_string(n) + ".nc";

    NcFile dataFile(ncName, NcFile::replace);
    // Create netCDF dimensions
    NcDim p = dataFile.addDim("p", 6);
    NcDim xDim = dataFile.addDim("x", nx);
    NcDim zDim = dataFile.addDim("z", nz);

    vector<NcDim> xzDim;
    xzDim.push_back(xDim);
    xzDim.push_back(zDim);


    NcVar thData = dataFile.addVar("th", ncDouble, xzDim);
    NcVar zetaData = dataFile.addVar("zeta", ncDouble, xzDim);
    NcVar uData = dataFile.addVar("u", ncDouble, xzDim);
    NcVar wData = dataFile.addVar("w", ncDouble, xzDim);
    NcVar qcData = dataFile.addVar("qc", ncDouble, xzDim);
    NcVar qrData = dataFile.addVar("qr", ncDouble, xzDim);
    NcVar qvData = dataFile.addVar("qv", ncDouble, xzDim);

	thData.putVar(myArray.th);
	zetaData.putVar(myArray.zeta);
	uData.putVar(myArray.u);
	wData.putVar(myArray.w);
	qcData.putVar(myArray.qc);
	qrData.putVar(myArray.qr);
	qvData.putVar(myArray.qr);
}