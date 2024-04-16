#include "Declare.hpp"
#include <iostream>

#if defined(WATER)
void vvm::MicroPhysics::condensation(vvm &model) {
    double qvs = 0., phi = 0., C = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        double pc = 380. / (pow(model.pib[k], C_p / Rd) * P0); 	// coefficient
        for (int i = 1; i <= model.nx-2; i++) {
            qvs = pc * exp(17.27 * (model.pib[k] * model.thp[i][k] - 273.) / (model.pib[k] * model.thp[i][k] - 36.));
            phi = qvs * (17.27 * 237. * Lv) / (C_p * pow(model.thp[i][k] * model.pib[k] - 36., 2));
            C = (model.qvp[i][k] - qvs) / (1 + phi); 
    
            // C should less than qc (C can be sink for qc and source for qv, so it should not excess qc)
            if (fabs(C) > model.qcp[i][k] && C < 0) C = -model.qcp[i][k];

            model.qvp[i][k] = model.qvp[i][k] - C;
            model.qcp[i][k] = model.qcp[i][k] + C;
            model.thp[i][k] = model.thp[i][k] + Lv / (C_p * model.pib[k]) * C;
        }
    }
    return;
}

// autoconversion of qc to qr
void vvm::MicroPhysics::autoconversion(vvm & model) {
    double autort = 0.001, autotr = 0.001; // autocon rate [1/sec], autocon threshold [kg/kg]
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {

            double qcplus = std::max(0., model.qcp[i][k]);
            double ar = autort * (qcplus - autotr);
            ar = std::max(0., ar);
            #if defined(AB3)
                double arcrdt = std::min(ar * DT, qcplus);
            #else
                double arcrdt = std::min(ar * model.d2t, qcplus);
            #endif
            model.qcp[i][k] = model.qcp[i][k] - arcrdt;
            model.qrp[i][k] = model.qrp[i][k] + arcrdt;
            model.autoconversion[i][k] = arcrdt;
        }
    }
    return;
}


void vvm::MicroPhysics::accretion(vvm &model) {
    double autort = 0.001, autotr = 0.001; // autocon rate [1/sec], autocon threshold [kg/kg]
    double qcplus = 0., ar = 0., arcrdt = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            qcplus = std::max(0., model.qcp[i][k]);
            ar = autort * (qcplus - autotr);
            ar = std::max(0., ar);
            #if defined(AB3)
                arcrdt = std::min(ar * DT, qcplus);
            #else
                arcrdt = std::min(ar * model.d2t, qcplus);
            #endif
            model.qcp[i][k] = model.qcp[i][k] - arcrdt;
            model.qrp[i][k] = model.qrp[i][k] + arcrdt;
            model.accretion[i][k] = arcrdt;
        }
    }
    return;
}

// evaporation of rain water
void vvm::MicroPhysics::evaporation(vvm &model) {
    double pc = 0.;
    double qrplus = 0., qvplus = 0., qvs = 0., coef = 0., deficit = 0.;
    double er = 0., erdt = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        pc = 380. / (pow(model.pib[k], C_p / Rd) * P0);	 // coefficient
        coef = 1.6 + 30.39 * pow((model.rhou[k] * qrplus), 0.2046);	// ventilation coef.
        for (int i = 1; i <= model.nx-2; i++) {
            qrplus = std::max(0., model.qrp[i][k]);
            qvplus = std::max(0., model.qvp[i][k]);
            qvs = pc * exp(17.27 * (model.pib[k] * model.thp[i][k] - 273.) / (model.pib[k] * model.thp[i][k] - 36.));	// Tetens equation

            deficit = std::max((1. - qvplus / qvs), 0.);					// saturation dificit (RH < 100%)
            er = coef * deficit * (pow(model.rhou[k] * qrplus, 0.525)) / ((2.03e4 + 9.584e6 / (model.pb[k] * qvs)) * model.rhou[k]);
            #if defined(AB3)
                double erdt = std::min(qrplus, std::max(0., er * DT));
            #else
                double erdt = std::min(qrplus, std::max(0., er * model.d2t));
            #endif
            
            if (erdt < 0.) {
                std::cout << "Evaporation Wrong" << std::endl;
            }

            model.qrp[i][k] = model.qrp[i][k] - erdt;
            model.qvp[i][k] = model.qvp[i][k] + erdt;
            model.thp[i][k] = model.thp[i][k] - Lv * erdt / (C_p * model.pib[k]);
            model.evaporation[i][k] = erdt;
        }
    }
    return;
}

void vvm::MicroPhysics::NegativeValueProcess(double var[][NZ]) {
    double positive = 0.;
    double negative = 0.;
    for (int k = 1; k <= NZ-2; k++) {
        for (int i = 1; i <= NX-2; i++) {
            if (var[i][k] >= 0.) positive += var[i][k];
            else {
                negative += var[i][k];
                var[i][k] = 0.;
            }
        }
    }

    if (positive == 0. || std::abs(negative) > positive) return;

    double correctionRatio = 1. - std::abs(negative/positive);
    for (int k = 1; k <= NZ-2; k++) {
        for (int i = 1; i <= NX-2; i++) {
            if (var[i][k] > 0) var[i][k] = var[i][k] * correctionRatio;
        }
    }
    return;
}
#endif
