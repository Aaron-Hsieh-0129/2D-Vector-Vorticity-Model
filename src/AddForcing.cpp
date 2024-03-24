#include "Declare.hpp"

#if defined(WATER)
#if defined(TROPICALFORCING)


void vvm::AddForcing(vvm &model) {
    double dt = 0.;
    for (int k = 1; k <= model.nz-2; k++) {
        for (int i = 1; i <= model.nx-2; i++) {
            #ifndef AB3
                dt = model.d2t;
            #else
                dt = DT;
            #endif
            model.thp[i][k] += dt * model.Q1LS[k];
            if (model.status_for_adding_forcing == true) model.thp[i][k] += dt * model.init_th_forcing[i][k];
            #if defined(RADIATIONCOOLING) 
                model.thp[i][k] += dt * (-2. / 86400.);
            #endif

            model.qvp[i][k] += dt * model.Q2LS[k];
        }
    }
    model.BoundaryProcess2D_center(model.thp);
    model.BoundaryProcess2D_center(model.qvp);
}
#endif
#endif
