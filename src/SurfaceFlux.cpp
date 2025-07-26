#include "Declare.hpp"

#if defined(SFCFLX)

void vvm::FluxCoefficient(double thvm, double thvsm, double speed1, double zr, double zrough, double &ustar, double &ventfc, double &vent2, double &molen) {
    double vk = 0.4; // von Karman constant
    double grav = 9.80665;
    double bus = 0.74;
    double crit = 0.003;
    int maxit = 5;

    double tem1 = std::log(zr/zrough);
    double cuni = tem1 / vk;
    double ctni = cuni * bus;
    double cu = 1. / cuni;
    double ct = 1. / ctni;

    bool stable = thvsm < 0;
    bool stopit = false;
    double speedm = std::max(speed1, 1e-3);

    int it = 0;

    double zeta = 0.;
    double tem2 = 0., tem3 = 0.;
    double cui = 0., cti = 0.;
    double x = 0., y = 0.;
    if (!stable) speedm = std::max(speedm, 0.);
    while (!stopit) {
        it++;
        zeta = -zr * ct * vk * grav * thvsm / (thvm * cu * cu * speedm * speedm);

        if (stable) {
            if (zeta >= 2.45) {
                stopit = true;
                zeta = 2.45;
            }
            tem2 = tem1 + 4.7 * zeta;
            tem3 = tem1 + 4.7 / bus * zeta;

            cui = tem2 / vk;
            cti = bus * tem3 / vk;
        }
        else {
            x = std::pow(1. - 15.*zeta, 0.25);
            y = std::pow(1. - 9.*zeta, 0.25);
            tem2 = tem1 - (std::log(0.5*(1+x*x)) + 2.*std::log(0.5*(1+x)) - 2.*std::atan(x) + M_PI / 2.);
            tem3 = tem1 - 2. * std::log(0.5*(1+y*y));
            cui = tem2 / vk;
            cui = std::max(cui, 0.5*cuni);
            cti = bus * tem3 / vk;
            cti = std::max(cti, 0.3*ctni);
        }

        stopit = stopit || it == maxit;
        if (stopit) {
            cu = 1. / cui;
            ct = 1. / cti;
        }

        // Check for convergence
        double custar = cu;
        double ctstar = ct;
        cu = 1. / cui;
        ct = 1. / cti;
        stopit = std::abs(cu/custar - 1.) <= crit && std::abs(ct/ctstar-1.) <= crit;
    }

    if (stable) {
        ustar = cu * speedm;
        ventfc = cu * ustar;
        vent2 = ct * ustar;
    }
    else {
        ustar = cu * speedm;
        ventfc = cu * ustar;
        vent2 = ct * ustar;

        if (cti < 0.3*ctni) vent2 = std::max(vent2, 0.0019*std::pow(thvsm, 1./3.));
    }

    zeta = -zr * ct * vk * grav * thvsm / (thvm * cu * cu * speedm * speedm);
    zeta = std::max(std::abs(zeta), 1e-6) * std::copysign(1., zeta);
    molen = zr / std::min(zeta, 2.45);
}


void vvm::SurfaceFlux(vvm &model) { 
    double ztmp = 0.5 * model.dz / model.flex_height_coef_th[1];
    std::vector<double> tg(model.nx), ts(model.nx), Q(model.nx), T(model.nx), z_rough(model.nx), speedm(model.nx), thvsm(model.nx);
    std::vector<double> ustar(model.nx), vent(model.nx), vent2(model.nx), molen(model.nx);
    std::vector<double> WT(model.nx), WTH(model.nx), WQ(model.nx), UW(model.nx);
    double Lvm = 3.36E5;
    double GWET = 1.;
    double DELTA = 0.608;

    for (int i = 1; i < model.nx-1; i++) {
        tg[i] = 303. * model.pib[1];
        double es1 = 611.2 * std::exp(17.67 * (tg[i]-273.15) / (tg[i]-273.15+243.5));
        double qvs_sfc = es1 * 0.622 / (model.pb_lev[1] - 0.378 * es1);
        ts[i] = model.Cp * tg[i] + model.GRAVITY * model.z_zeta[1];
        Q[i] = model.qv[i][1] + model.qc[i][1] + model.qitot[i][1];
        T[i] = model.Cp * model.th[i][1] * model.pib[1] - model.Lv * model.qc[i][1] + model.GRAVITY * model.z[1] - (model.Lv+Lvm) * model.qitot[i][1];

        z_rough[i] = 2e-4;
        speedm[i] = std::max(0.5*std::abs(model.u[i+1][1]+model.u[i][1]), 1e-3);
        thvsm[i] = tg[i] / model.pib[1] - model.th[i][1] + std::abs(GWET)*model.thb[1]*(DELTA*(qvs_sfc-model.qv[i][1]));

        FluxCoefficient(model.thb[1], thvsm[i], speedm[i], ztmp, z_rough[i], ustar[i], vent[i], vent2[i], molen[i]);

        WT[i] = vent2[i] * (ts[i] - T[i]);
        WQ[i] = vent2[i] * std::abs(GWET) * (qvs_sfc - Q[i]);

        model.heatflux[i] = WT[i] * model.rhow[1] / (model.Cp * model.pib[1]);
        model.waterflux[i] = WQ[i] * model.rhow[1];
        model.momentumflux[i] = -0.5*(vent[i]+vent[i-1]) * model.u[i][1] * model.rhow[1];
        
        model.thp[i][1] += model.heatflux[i] / model.rhou[1] * model.flex_height_coef_th[1] * model.rdz * model.dt;
        model.qvp[i][1] +=  model.waterflux[i] / model.rhou[1] * model.flex_height_coef_th[1] * model.rdz * model.dt;
        model.zeta[i][1] +=  model.momentumflux[i] / model.rhou[1] * model.flex_height_coef_zeta[1] * model.flex_height_coef_th[1] * model.rdz2 * model.dt;
    }
    return;
}
#endif
