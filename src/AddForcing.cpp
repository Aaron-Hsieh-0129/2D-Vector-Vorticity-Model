#include "Declare.hpp"

#if defined(WATER)


void vvm::AddForcing(vvm &model) {
    double dt = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            #ifndef AB2
                dt = model.d2t;
            #else
                dt = model.dt;
            #endif
            #if defined(TROPICALFORCING)
                model.thp[i][k] += dt * model.Q1LS[k];
                model.qvp[i][k] += dt * model.Q2LS[k];
            #endif

            #if defined(RTERRTMGP)
                model.thp[i][k] += dt * model.radiation_heating_rate[i][k] / 86400.; 
            #endif
            #if defined(RADIATIONCOOLING) 
                model.thp[i][k] += dt * (-2. / 86400.);
            #endif

            if (model.status_for_adding_forcing == true) model.thp[i][k] += dt * model.init_th_forcing[i][k];
        }
    }
}
#endif
