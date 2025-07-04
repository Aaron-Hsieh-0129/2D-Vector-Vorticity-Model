#include "Declare.hpp"

#if defined(TIMEFILTER)
void vvm::NumericalProcess::TimeFilter(double **previous, double **now, double **future, vvm &model) {
    for (int i = 0; i <= model.nx-1; i++) {
        for (int k = 0; k <= model.nz-1; k++) {
            now[i][k] += model.TIMETS * (future[i][k] - 2.*now[i][k] + previous[i][k]);
        }
    }
    return;
}

void vvm::NumericalProcess::timeFilterAll(vvm &model) {
    #if defined(TIMEFILTER)
        TimeFilter(model.zetam, model.zeta, model.zetap, model);
        TimeFilter(model.thm, model.th, model.thp, model);
        #if defined(WATER)
            TimeFilter(model.qvm, model.qv, model.qvp, model);
            TimeFilter(model.qcm, model.qc, model.qcp, model);
            TimeFilter(model.qrm, model.qr, model.qrp, model);
        #endif
    #endif
}
#endif

#if defined(DIFFUSION_VVM)
void vvm::NumericalProcess::Diffusion(double **var_in, double **var_out, vvm &model) {
    #ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int k = 1; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            var_out[i][k] += model.d2t * model.Kx * model.rdx2 * (var_in[i+1][k] - 2. * var_in[i][k] + var_in[i-1][k]) + 
                             model.d2t * model.Kz * model.rdz2 * (var_in[i][k+1] - 2. * var_in[i][k] + var_in[i][k-1]);
        }
    }
    return;
}


void vvm::NumericalProcess::DiffusionAll(vvm &model) {
    Diffusion(model.zetam, model.zetap, model);
    Diffusion(model.thm, model.thp, model);
    #if defined(WATER)
        Diffusion(model.qvm, model.qvp, model);
        Diffusion(model.qcm, model.qcp, model);
        Diffusion(model.qrm, model.qrp, model);
        #if defined(P3_MICROPHY)
            Diffusion(model.ncm, model.ncp, model);
            Diffusion(model.nrm, model.nrp, model);
            Diffusion(model.nim, model.nip, model);
            Diffusion(model.qitotm, model.qitotp, model);
            Diffusion(model.qirimm, model.qirimp, model);
            Diffusion(model.birimm, model.birimp, model);
        #endif
    #endif
}
#endif


void vvm::NumericalProcess::GravityWaveDampingExponential(vvm &model) {
    for (int k = model.k_diff_start; k < model.nz-1; k++) {
        for (int i = 1; i < model.nx-1; i++) {
            model.thp[i][k] -= model.dt * model.nudge_tau[k] * (model.thp[i][k] - model.thb_init[k]);
            model.zetap[i][k] -= model.dt * model.nudge_tau[k] * (model.zetap[i][k]);
            #if defined(WATER)
                model.qvp[i][k] -= model.dt * model.nudge_tau[k] * (model.qvp[i][k] - model.qvb0[k]);
                model.qcp[i][k] -= model.dt * model.nudge_tau[k] * model.qcp[i][k];
                model.qrp[i][k] -= model.dt * model.nudge_tau[k] * model.qrp[i][k];
                #if defined(P3_MICROPHY)
                    model.ncp[i][k] -= model.dt * model.nudge_tau[k] * model.ncp[i][k];
                    model.nrp[i][k] -= model.dt * model.nudge_tau[k] * model.nrp[i][k];
                    model.nip[i][k] -= model.dt * model.nudge_tau[k] * model.nip[i][k];
                    model.qitotp[i][k] -= model.dt * model.nudge_tau[k] * model.qitotp[i][k];
                    model.qirimp[i][k] -= model.dt * model.nudge_tau[k] * model.qirimp[i][k];
                    model.birimp[i][k] -= model.dt * model.nudge_tau[k] * model.birimp[i][k];
                #endif
            #endif
        }
    }
}


void vvm::NumericalProcess::Nudge_qv(vvm &model) {
    if (model.moisture_nudge_time == 0) return;

    for (int i = 0; i <= model.nx-1; i++) {
        for (int k = 0; k <= model.nz-1; k++) {
            model.qvp[i][k] = model.qvp[i][k] + model.dt * (model.qvb0[k] - model.qv[i][k]) / model.moisture_nudge_time;
        }
    }
    return;
}

void vvm::NumericalProcess::NegativeValueProcess(double **var, int nx, int nz) {
    double positive = 0.;
    double negative = 0.;
    for (int k = 1; k <= nz-2; k++) {
        for (int i = 1; i <= nx-2; i++) {
            if (var[i][k] >= 0.) positive += var[i][k];
            else {
                negative += var[i][k];
                var[i][k] = 0.;
            }
        }
    }

    if (positive == 0. || std::abs(negative) > positive) return;

    double correctionRatio = 1. - std::abs(negative/positive);
    for (int k = 1; k <= nz-2; k++) {
        for (int i = 1; i <= nx-2; i++) {
            if (var[i][k] > 0) var[i][k] = var[i][k] * correctionRatio;
        }
    }
    return;
}


// Function to perform interpolation/extrapolation
void vvm::NumericalProcess::interpolate(const std::vector<double>& known_heights,
                 const std::vector<std::vector<double>>& known_data_fields,
                 const std::vector<double>& new_heights,
                 std::vector<std::vector<double>>& interpolated_data_fields) {
    // Ensure input arrays have consistent sizes
    assert(!known_heights.empty());
    assert(known_heights.size() >= 2); // Need at least 2 points for interpolation/extrapolation
    for (const auto& field : known_data_fields) {
        assert(field.size() == known_heights.size());
    }

    // Resize output vectors to match new_heights size and number of data fields
    interpolated_data_fields.resize(known_data_fields.size());
    for (auto& field : interpolated_data_fields) {
        field.resize(new_heights.size());
    }

    for (size_t k = 0; k < new_heights.size(); ++k) {
        double h = new_heights[k];

        // Handle extrapolation for heights below the minimum known height
        if (h < known_heights[0]) {
            // Use the first two points to compute the slope for extrapolation
            double h0 = known_heights[0];
            double h1 = known_heights[1];

            for (size_t j = 0; j < known_data_fields.size(); ++j) {
                double v0 = known_data_fields[j][0];
                double v1 = known_data_fields[j][1];
                double slope = (v1 - v0) / (h1 - h0);
                interpolated_data_fields[j][k] = v0 + slope * (h - h0);
            }
            continue;
        }

        // Handle extrapolation for heights above the maximum known height
        if (h > known_heights.back()) {
            // Use the last two points to compute the slope for extrapolation
            size_t n = known_heights.size() - 1;
            double h0 = known_heights[n - 1];
            double h1 = known_heights[n];

            for (size_t j = 0; j < known_data_fields.size(); ++j) {
                double v0 = known_data_fields[j][n - 1];
                double v1 = known_data_fields[j][n];
                double slope = (v1 - v0) / (h1 - h0);
                interpolated_data_fields[j][k] = v1 + slope * (h - h1);
            }
            continue;
        }

        // Interpolation for heights within the range
        // Find the first height greater than h using binary search
        auto it = std::upper_bound(known_heights.begin(), known_heights.end(), h);

        // If h equals the last known height
        if (it == known_heights.end()) {
            for (size_t j = 0; j < known_data_fields.size(); ++j) {
                interpolated_data_fields[j][k] = known_data_fields[j].back();
            }
            continue;
        }

        // Find the bracketing indices
        size_t i = it - known_heights.begin() - 1;
        double h0 = known_heights[i];
        double h1 = known_heights[i + 1];
        double factor = (h - h0) / (h1 - h0);

        // Linear interpolation for each data field
        for (size_t j = 0; j < known_data_fields.size(); ++j) {
            interpolated_data_fields[j][k] = known_data_fields[j][i] + 
                factor * (known_data_fields[j][i + 1] - known_data_fields[j][i]);
        }
    }
}
