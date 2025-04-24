#include "Declare.hpp"

void vvm::SurfaceFlux(vvm &model) { 
    double cdh = 7E-3;
    for (int k = 1; k <=2; k++) {
        for (int i = 1; i <= model.nx/2; i++) {
            double thdif = std::max(model.th_ground[i]-model.th[i][k], 0.);
            double avgu = std::max(0.5*std::abs(model.u[i+1][k]+model.u[i][k]), 2.); // Note: which time step?
            double wnetc = 2. * std::sqrt(thdif); // convective velocity adjustment (Forvell)
            double vel = std::sqrt(avgu*avgu+wnetc*wnetc);

            double qvdif = std::max(model.qvs_ground[i]-model.qv[i][k], 0.);
            double wnetqv = 2. * std::sqrt(qvdif);
            double velqv = std::sqrt(avgu*avgu+wnetqv*wnetqv);


            model.heatflux[i] = cdh * vel * model.addflux[i] * thdif / model.dz;
            model.waterflux[i] = cdh * velqv * model.addflux[i] * qvdif / model.dz;
            
            model.thp[i][k] += model.heatflux[i] * model.dt;
            model.qvp[i][k] += model.waterflux[i] * model.dt;
        }
    }
    return;
}
