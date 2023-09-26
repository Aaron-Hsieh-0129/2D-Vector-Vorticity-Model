#include "Outputfile.hpp"

using namespace netCDF;
using std::vector;
using std::string;

void Output::printInit(vvmArray &model) {
	double z;
	std::cout << "z          tb        rhou       rhow       qvb	 RH      pib" << std::endl;
	for (int k = 0; k <= nz-1;k++){
		z = (double) (k - 0.5) * dz ;
		std::cout << std::fixed << std::setprecision(2) << z << "    " << model.tb[k] << "    " << model.rhou[k] << "     " 
		<< model.rhow[k] << "   	 " << model.qvb[k] * 1000 << "    " << model.qvb[k] / model.qvsb[k] << "    "
		<< model.pib[k] << std::endl;
	}
	std::fstream initout;
	string initName = OUTPUTPATH + (string) "init.txt";
	initout.open(initName, std::ios::out);
	for (int k = 0; k <= nz-1; k++) {
		z = (double) (k - 0.5) * dz ;
		initout << z << "    " << model.tb[k] << "    " << model.rhou[k] << "     " 
		<< model.rhow[k] << "   	 " << model.qvb[k] << "    " << model.qvsb[k] << "    " << model.qvb[k] / model.qvsb[k] << "    "
		<< model.pib[k] << std::endl;
	}
	return;
}


void Output::output_zeta(int n, vvmArray &model) {
	std::fstream foutzeta;
	string zetaName = OUTPUTPATH + (string) "/zeta/zeta_" + std::to_string(n) + (string) ".txt";
	foutzeta.open(zetaName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutzeta << model.zeta[i][k] << " ";
		}
	}
}

void Output::output_th(int n, vvmArray &model) {
	std::fstream foutth;
	string thName = OUTPUTPATH + (string) "th/th_" + std::to_string(n) + (string) ".txt";
	foutth.open(thName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutth << model.th[i][k] << " ";
		}
	}
}

void Output::output_u(int n, vvmArray &model) {
	std::fstream foutu;
	string uName = OUTPUTPATH + (string) "u/u_" + std::to_string(n) + (string) ".txt";
	foutu.open(uName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutu << model.u[i][k] << " ";
		}
	}
}

void Output::output_w(int n, vvmArray &model) {
	std::fstream foutw;
	string wName = OUTPUTPATH + (string) "w/w_" + std::to_string(n) + (string) ".txt";
	foutw.open(wName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutw << model.w[i][k] << " ";
		}
	}
}

void Output::output_qv(int n, vvmArray &model) {
	std::fstream foutqv;
	string qvName = OUTPUTPATH + (string) "qv/qv_" + std::to_string(n) + (string) ".txt";
	foutqv.open(qvName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqv << model.qv[i][k] + model.qvb[k] << " ";
		}
	}
}

void Output::output_qc(int n, vvmArray &model) {
	std::fstream foutqc;
	string qcName = OUTPUTPATH + (string) "qc/qc_" + std::to_string(n) + (string) ".txt";
	foutqc.open(qcName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqc << model.qc[i][k] << " ";
		}
	}
}

void Output::output_qr(int n, vvmArray &model) {
	std::fstream foutqr;
	string qrName = OUTPUTPATH + (string) "qr/qr_" + std::to_string(n) + (string) ".txt";
	foutqr.open(qrName, std::ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqr << model.qr[i][k] << " ";
		}
	}
}

void Output::output_nc(int n, vvmArray &model) {
	string ncName = OUTPUTPATH + (string) "nc/" + std::to_string(n) + (string) ".nc";

	NcFile dataFile(ncName, NcFile::replace);

	// Create netCDF dimensions
	NcDim xDim = dataFile.addDim("x", nx);
	NcDim zDim = dataFile.addDim("z", nz);
	vector<NcDim> xzDim, zNcDim;
	xzDim.push_back(xDim);
	xzDim.push_back(zDim);
	zNcDim.push_back(zDim);

	NcVar thData = dataFile.addVar("th", ncDouble, xzDim);
	NcVar zetaData = dataFile.addVar("zeta", ncDouble, xzDim);
	NcVar uData = dataFile.addVar("u", ncDouble, xzDim);
	NcVar wData = dataFile.addVar("w", ncDouble, xzDim);
	NcVar qcData = dataFile.addVar("qc", ncDouble, xzDim);
	NcVar qrData = dataFile.addVar("qr", ncDouble, xzDim);
	NcVar qvData = dataFile.addVar("qv", ncDouble, xzDim);

	thData.putVar(model.th);
	zetaData.putVar(model.zeta);
	uData.putVar(model.u);
	wData.putVar(model.w);
	qcData.putVar(model.qc);
	qrData.putVar(model.qr);
	qvData.putVar(model.qr);

	if (n == 0) {
		NcVar tbData = dataFile.addVar("tb", ncDouble, zNcDim);
		NcVar rhoData = dataFile.addVar("rho", ncDouble, zNcDim);
		NcVar qvbData = dataFile.addVar("qvb", ncDouble, zNcDim);
		NcVar qvsbData = dataFile.addVar("qvsb", ncDouble, zNcDim);

		tbData.putVar(model.tb);
		rhoData.putVar(model.rhou);
		qvbData.putVar(model.qvb);
		qvsbData.putVar(model.qvsb);
	}
}

void Output::create_directory(string directory_name) {
    string str = "mkdir -p " + directory_name;
    const char *command = str.c_str();
    const int dir_err = system(command);
    if (-1 == dir_err) {
        std::cout << "Error on creating directory!\n" << std::endl;
        return;
    }
    return;
}

void Output::create_all_directory() {
    // data directory
    #ifdef TXTOUTPUT
        create_directory(OUTPUTPATH + (string) "u");
        create_directory(OUTPUTPATH + (string) "w");
        create_directory(OUTPUTPATH + (string) "zeta");
        create_directory(OUTPUTPATH + (string) "th");
        create_directory(OUTPUTPATH + (string) "qc");
        create_directory(OUTPUTPATH + (string) "qr");
        create_directory(OUTPUTPATH + (string) "qv");
    #endif
    #ifdef NCOUTPUT
        create_directory(OUTPUTPATH + (string) "nc");
    #endif

    // plot directory
    create_directory(OUTPUTPATH + (string) "graphs/u");
    create_directory(OUTPUTPATH + (string) "graphs/w");
    create_directory(OUTPUTPATH + (string) "graphs/zeta");
    create_directory(OUTPUTPATH + (string) "graphs/th");
    create_directory(OUTPUTPATH + (string) "graphs/qc");
    create_directory(OUTPUTPATH + (string) "graphs/qr");
    create_directory(OUTPUTPATH + (string) "graphs/qv");
    create_directory(OUTPUTPATH + (string) "graphs/qc+qr+th+u+w");
    create_directory(OUTPUTPATH + (string) "graphs/qv+qc");
    create_directory(OUTPUTPATH + (string) "graphs/qc+qr");
    create_directory(OUTPUTPATH + (string) "graphs/qr+th+u+w");
    create_directory(OUTPUTPATH + (string) "graphs/qc+qr+th+u+w");
}