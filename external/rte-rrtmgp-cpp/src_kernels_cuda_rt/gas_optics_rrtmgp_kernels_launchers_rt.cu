#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "gas_optics_rrtmgp_kernels_cuda_rt.h"
#include "tools_gpu.h"
#include "tuner.h"


namespace
{
    #include "gas_optics_rrtmgp_kernels_rt.cu"
}


namespace Gas_optics_rrtmgp_kernels_cuda_rt
{
    void reorder123x321(
            const int ni, const int nj, const int nk,
            const Float* arr_in, Float* arr_out)
    {
        Tuner_map& tunings = Tuner::get_map();

        dim3 grid{ni, nj, nk};
        dim3 block;

        if (tunings.count("reorder123x321_kernel_rt") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                "reorder123x321_kernel_rt",
                dim3{ni, nj, nk},
                {1, 2, 4, 8, 16, 24, 32, 48, 64, 96},
                {1, 2, 4, 8, 16, 24, 32, 48, 64, 96},
                {1, 2, 4, 8, 16, 24, 32, 48, 64, 96},
                reorder123x321_kernel,
                ni, nj, nk, arr_in, arr_out);

            tunings["reorder123x321_kernel_rt"].first = grid;
            tunings["reorder123x321_kernel_rt"].second = block;
        }
        else
        {
            grid = tunings["reorder123x321_kernel_rt"].first;
            block = tunings["reorder123x321_kernel_rt"].second;
        }

        reorder123x321_kernel<<<grid, block>>>(
                ni, nj, nk, arr_in, arr_out);
    }


    void reorder12x21(const int ni, const int nj,
                      const Float* arr_in, Float* arr_out)
    {
        const int block_i = 32;
        const int block_j = 16;

        const int grid_i = ni/block_i + (ni%block_i > 0);
        const int grid_j = nj/block_j + (nj%block_j > 0);

        dim3 grid_gpu(grid_i, grid_j);
        dim3 block_gpu(block_i, block_j);

        reorder12x21_kernel<<<grid_gpu, block_gpu>>>(
                ni, nj, arr_in, arr_out);
    }


    void zero_array(const int ni, const int nj, const int nk, const int nn, Float* arr)
    {
        const int block_i = 32;
        const int block_j = 16;
        const int block_k = 1;

        const int grid_i = ni/block_i + (ni%block_i > 0);
        const int grid_j = nj/block_j + (nj%block_j > 0);
        const int grid_k = nk/block_k + (nk%block_k > 0);

        dim3 grid_gpu(grid_i, grid_j, grid_k);
        dim3 block_gpu(block_i, block_j, block_k);

        zero_array_kernel<<<grid_gpu, block_gpu>>>(
                ni, nj, nk, nn, arr);
    }


    void zero_array(const int ni, const int nj, const int nk, Float* arr)
    {
        const int block_i = 32;
        const int block_j = 16;
        const int block_k = 1;

        const int grid_i = ni/block_i + (ni%block_i > 0);
        const int grid_j = nj/block_j + (nj%block_j > 0);
        const int grid_k = nk/block_k + (nk%block_k > 0);

        dim3 grid_gpu(grid_i, grid_j, grid_k);
        dim3 block_gpu(block_i, block_j, block_k);

        zero_array_kernel<<<grid_gpu, block_gpu>>>(
                ni, nj, nk, arr);

    }


    void zero_array(const int ni, const int nj, Float* arr)
    {
        const int block_i = 32;
        const int block_j = 16;

        const int grid_i = ni/block_i + (ni%block_i > 0);
        const int grid_j = nj/block_j + (nj%block_j > 0);

        dim3 grid_gpu(grid_i, grid_j, 1);
        dim3 block_gpu(block_i, block_j, 1);

        zero_array_kernel<<<grid_gpu, block_gpu>>>(
                ni, nj, arr);

    }

    void zero_array(const int ni, int* arr)
    {
        const int block_i = 32;

        const int grid_i = ni/block_i + (ni%block_i > 0);

        dim3 grid_gpu(grid_i);
        dim3 block_gpu(block_i);

        zero_array_kernel<<<grid_gpu, block_gpu>>>(
                ni, arr);
    }

    void interpolation(
            const int ncol, const int nlay,
            const int ngas, const int nflav, const int neta, const int npres, const int ntemp,
            const int* flavor,
            const Float* press_ref_log,
            const Float* temp_ref,
            Float press_ref_log_delta,
            Float temp_ref_min,
            Float temp_ref_delta,
            Float press_ref_trop_log,
            const Float* vmr_ref,
            const Float* play,
            const Float* tlay,
            Float* col_gas,
            int* jtemp,
            Float* fmajor, Float* fminor,
            Float* col_mix,
            Bool* tropo,
            int* jeta,
            int* jpress)
    {
        Tuner_map& tunings = Tuner::get_map();
        Float tmin = std::numeric_limits<Float>::min();

        dim3 grid(ncol, nlay, nflav), block;
        if (tunings.count("interpolation_kernel_rt") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                    "interpolation_kernel_rt",
                    dim3{ncol, nlay, nflav},
                    {1, 2, 4, 8, 16, 32, 64, 128, 256}, {1}, {1, 2, 4, 8, 16, 32, 64, 128, 256},
                    interpolation_kernel,
                    ncol, nlay, ngas, nflav, neta, npres, ntemp, tmin,
                    flavor, press_ref_log, temp_ref,
                    press_ref_log_delta, temp_ref_min,
                    temp_ref_delta, press_ref_trop_log,
                    vmr_ref, play, tlay,
                    col_gas, jtemp, fmajor,
                    fminor, col_mix, tropo,
                    jeta, jpress);
            tunings["interpolation_kernel_rt"].first = grid;
            tunings["interpolation_kernel_rt"].second = block;
        }
        else
        {
            grid = tunings["interpolation_kernel_rt"].first;
            block = tunings["interpolation_kernel_rt"].second;
        }

        interpolation_kernel<<<grid, block>>>(
                ncol, nlay, ngas, nflav, neta, npres, ntemp, tmin,
                flavor, press_ref_log, temp_ref,
                press_ref_log_delta, temp_ref_min,
                temp_ref_delta, press_ref_trop_log,
                vmr_ref, play, tlay,
                col_gas, jtemp, fmajor,
                fminor, col_mix, tropo,
                jeta, jpress);

    }


    void minor_scalings(
            const int ncol, const int nlay, const int nflav, const int ngpt,
            const int nminorlower, const int nminorupper,
            const int idx_h2o,
            const int* gpoint_flavor,
            const int* minor_limits_gpt_lower,
            const int* minor_limits_gpt_upper,
            const Bool* minor_scales_with_density_lower,
            const Bool* minor_scales_with_density_upper,
            const Bool* scale_by_complement_lower,
            const Bool* scale_by_complement_upper,
            const int* idx_minor_lower,
            const int* idx_minor_upper,
            const int* idx_minor_scaling_lower,
            const int* idx_minor_scaling_upper,
            const Float* play,
            const Float* tlay,
            const Float* col_gas,
            const Bool* tropo,
            Float* scalings_lower,
            Float* scalings_upper)
    {
        Tuner_map& tunings = Tuner::get_map();

        // lower atmosphere
        int idx_tropo=1;

        dim3 grid_1(ncol, nlay, nminorlower), block_1;
        if (tunings.count("minor_scalings_lower_kernel_rt") == 0)
        {
            std::tie(grid_1, block_1) = tune_kernel(
                    "minor_scalings_lower_kernel_rt",
                    dim3{ncol, nlay, nminorlower},
                    {8, 16, 24, 32, 48, 64, 128}, {1}, {1, 2, 4, 8, 16, 32},
                    scaling_kernel,
                    ncol, nlay, nflav, nminorlower,
                    idx_h2o, idx_tropo,
                    gpoint_flavor,
                    minor_limits_gpt_lower,
                    minor_scales_with_density_lower,
                    scale_by_complement_lower,
                    idx_minor_lower,
                    idx_minor_scaling_lower,
                    play, tlay, col_gas,
                    tropo, scalings_lower);

            tunings["minor_scalings_lower_kernel_rt"].first = grid_1;
            tunings["minor_scalings_lower_kernel_rt"].second = block_1;
        }
        else
        {
            grid_1 =  tunings["minor_scalings_lower_kernel_rt"].first;
            block_1 = tunings["minor_scalings_lower_kernel_rt"].second;
        }

        scaling_kernel<<<grid_1, block_1>>>(
                ncol, nlay, nflav, nminorlower,
                idx_h2o, idx_tropo,
                gpoint_flavor,
                minor_limits_gpt_lower,
                minor_scales_with_density_lower,
                scale_by_complement_lower,
                idx_minor_lower,
                idx_minor_scaling_lower,
                play, tlay, col_gas,
                tropo, scalings_lower);

        // upper atmosphere
        idx_tropo=0;

        dim3 grid_2(ncol, nlay, nminorupper), block_2;
        if (tunings.count("minor_scalings_upper_kernel_rt") == 0)
        {
            std::tie(grid_2, block_2) = tune_kernel(
                    "minor_scalings_upper_kernel_rt",
                    dim3{ncol, nlay, nminorupper},
                    {8, 16, 24, 32, 48, 64, 128}, {1}, {1, 2, 4, 8, 16, 32},
                    scaling_kernel,
                    ncol, nlay, nflav, nminorupper,
                    idx_h2o, idx_tropo,
                    gpoint_flavor,
                    minor_limits_gpt_upper,
                    minor_scales_with_density_upper,
                    scale_by_complement_upper,
                    idx_minor_upper,
                    idx_minor_scaling_upper,
                    play, tlay, col_gas,
                    tropo, scalings_upper);

            tunings["minor_scalings_upper_kernel_rt"].first = grid_2;
            tunings["minor_scalings_upper_kernel_rt"].second = block_2;
        }
        else
        {
            grid_2 =  tunings["minor_scalings_upper_kernel_rt"].first;
            block_2 = tunings["minor_scalings_upper_kernel_rt"].second;
        }

        scaling_kernel<<<grid_2, block_2>>>(
                ncol, nlay, nflav, nminorupper,
                idx_h2o, idx_tropo,
                gpoint_flavor,
                minor_limits_gpt_upper,
                minor_scales_with_density_upper,
                scale_by_complement_upper,
                idx_minor_upper,
                idx_minor_scaling_upper,
                play, tlay, col_gas,
                tropo, scalings_upper);
    }


    void combine_abs_and_rayleigh(
            const int ncol, const int nlay,
            const Float* tau_abs, const Float* tau_rayleigh,
            Float* tau, Float* ssa, Float* g)
    {
        Tuner_map& tunings = Tuner::get_map();

        Float tmin = std::numeric_limits<Float>::epsilon();

        dim3 grid{ncol, nlay, 1}, block;
        if (tunings.count("combine_abs_and_rayleigh_kernel_rt") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                "combine_abs_and_rayleigh_kernel_rt",
                dim3{ncol, nlay, 1},
                {24, 32, 48, 64, 96, 128, 256, 512}, {1, 2, 4}, {1},
                combine_abs_and_rayleigh_kernel,
                ncol, nlay, tmin,
                tau_abs, tau_rayleigh,
                tau, ssa, g);

            tunings["combine_abs_and_rayleigh_kernel_rt"].first = grid;
            tunings["combine_abs_and_rayleigh_kernel_rt"].second = block;
        }
        else
        {
            grid = tunings["combine_abs_and_rayleigh_kernel_rt"].first;
            block = tunings["combine_abs_and_rayleigh_kernel_rt"].second;
        }

        combine_abs_and_rayleigh_kernel<<<grid, block>>>(
                ncol, nlay, tmin,
                tau_abs, tau_rayleigh,
                tau, ssa, g);
    }


    void compute_tau_rayleigh(
            const int ncol, const int nlay, const int nbnd, const int ngpt, const int igpt,
            const int ngas, const int nflav, const int neta, const int npres, const int ntemp,
            const int* gpoint_flavor,
            const int* gpoint_bands,
            const int* band_lims_gpt,
            const Float* krayl,
            int idx_h2o, const Float* col_dry, const Float* col_gas,
            const Float* fminor, const int* jeta,
            const Bool* tropo, const int* jtemp,
            Float* tau_rayleigh)
    {
        Tuner_map& tunings = Tuner::get_map();

        dim3 grid{ncol, nlay, 1}, block;
        if (tunings.count("compute_tau_rayleigh_kernel_rt") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                "compute_tau_rayleigh_kernel_rt",
                dim3{ncol, nlay, 1},
                {24, 32, 64, 128, 256, 512}, {1, 2}, {1},
                compute_tau_rayleigh_kernel,
                ncol, nlay, nbnd, ngpt,
                ngas, nflav, neta, npres, ntemp,
                igpt,
                gpoint_flavor,
                gpoint_bands,
                band_lims_gpt,
                krayl,
                idx_h2o, col_dry, col_gas,
                fminor, jeta,
                tropo, jtemp,
                tau_rayleigh);

            tunings["compute_tau_rayleigh_kernel_rt"].first = grid;
            tunings["compute_tau_rayleigh_kernel_rt"].second = block;
        }
        else
        {
            grid = tunings["compute_tau_rayleigh_kernel_rt"].first;
            block = tunings["compute_tau_rayleigh_kernel_rt"].second;
        }

        compute_tau_rayleigh_kernel<<<grid, block>>>(
                ncol, nlay, nbnd, ngpt,
                ngas, nflav, neta, npres, ntemp,
                igpt,
                gpoint_flavor,
                gpoint_bands,
                band_lims_gpt,
                krayl,
                idx_h2o, col_dry, col_gas,
                fminor, jeta,
                tropo, jtemp,
                tau_rayleigh);
    }


    void compute_tau_absorption(
            const int ncol, const int nlay, const int nband, const int ngpt, const int igpt,
            const int ngas, const int nflav, const int neta, const int npres, const int ntemp,
            const int nminorlower, const int nminorklower,
            const int nminorupper, const int nminorkupper,
            const int idx_h2o,
            const int* gpoint_flavor,
            const int* band_lims_gpt,
            const Float* kmajor,
            const Float* kminor_lower,
            const Float* kminor_upper,
            const int* minor_limits_gpt_lower,
            const int* minor_limits_gpt_upper,
            const int* first_last_minor_lower,
            const int* first_last_minor_upper,
            const Bool* minor_scales_with_density_lower,
            const Bool* minor_scales_with_density_upper,
            const Bool* scale_by_complement_lower,
            const Bool* scale_by_complement_upper,
            const int* idx_minor_lower,
            const int* idx_minor_upper,
            const int* idx_minor_scaling_lower,
            const int* idx_minor_scaling_upper,
            const int* kminor_start_lower,
            const int* kminor_start_upper,
            const Bool* tropo,
            const Float* col_mix, const Float* fmajor,
            const Float* fminor, const Float* play,
            const Float* tlay, const Float* col_gas,
            const int* jeta, const int* jtemp,
            const int* jpress,
            const Float* scalings_lower,
            const Float* scalings_upper,
            Float* tau)
    {
        Tuner_map& tunings = Tuner::get_map();

        dim3 grid_maj{nlay, ncol, 1};
        dim3 block_maj;

        if (tunings.count("gas_optical_depths_major_kernel_rt") == 0)
        {
            Float* tau_tmp = Tools_gpu::allocate_gpu<Float>(nlay*ncol);

            std::tie(grid_maj, block_maj) = tune_kernel(
                    "gas_optical_depths_major_kernel_rt",
                    dim3{nlay, ncol, 1},
                    {1, 2}, {64, 96, 128, 256, 512, 768, 1024}, {1},
                    gas_optical_depths_major_kernel,
                    ncol, nlay, nband, ngpt,
                    nflav, neta, npres, ntemp,
                    igpt,
                    gpoint_flavor, band_lims_gpt,
                    kmajor, col_mix, fmajor, jeta,
                    tropo, jtemp, jpress,
                    tau_tmp);

            Tools_gpu::free_gpu<Float>(tau_tmp);

            tunings["gas_optical_depths_major_kernel_rt"].first = grid_maj;
            tunings["gas_optical_depths_major_kernel_rt"].second = block_maj;
        }
        else
        {
            grid_maj = tunings["gas_optical_depths_major_kernel_rt"].first;
            block_maj = tunings["gas_optical_depths_major_kernel_rt"].second;
        }

        gas_optical_depths_major_kernel<<<grid_maj, block_maj>>>(
            ncol, nlay, nband, ngpt,
            nflav, neta, npres, ntemp,
            igpt,
            gpoint_flavor, band_lims_gpt,
            kmajor, col_mix, fmajor, jeta,
            tropo, jtemp, jpress,
            tau);


        const int nscale_lower = nminorlower;
        const int nscale_upper = nminorupper;

        // Lower
        int idx_tropo = 1;

        dim3 grid_min_1(nlay, ncol, 1), block_min_1;
        if (tunings.count("gas_optical_depths_minor_kernel_lower_rt") == 0)
        {
            Float* tau_tmp = Tools_gpu::allocate_gpu<Float>(nlay*ncol);

            std::tie(grid_min_1, block_min_1) = tune_kernel(
                        "gas_optical_depths_minor_kernel_lower_rt",
                        dim3{nlay, ncol, 1},
                        {1}, {32, 48, 64, 96, 128, 256, 384, 512}, {1},
                        gas_optical_depths_minor_kernel,
                        ncol, nlay, ngpt, igpt,
                        ngas, nflav, ntemp, neta,
                        nscale_lower,
                        nminorlower,
                        nminorklower,
                        idx_h2o, idx_tropo,
                        gpoint_flavor,
                        kminor_lower,
                        minor_limits_gpt_lower,
                        first_last_minor_lower,
                        minor_scales_with_density_lower,
                        scale_by_complement_lower,
                        idx_minor_lower,
                        idx_minor_scaling_lower,
                        kminor_start_lower,
                        play, tlay, col_gas,
                        fminor, jeta, jtemp,
                        tropo, scalings_lower,
                        tau_tmp);
            Tools_gpu::free_gpu<Float>(tau_tmp);

            tunings["gas_optical_depths_minor_kernel_lower_rt"].first = grid_min_1;
            tunings["gas_optical_depths_minor_kernel_lower_rt"].second = block_min_1;
        }
        else
        {
            grid_min_1 = tunings["gas_optical_depths_minor_kernel_lower_rt"].first;
            block_min_1 = tunings["gas_optical_depths_minor_kernel_lower_rt"].second;
        }

        gas_optical_depths_minor_kernel<<<grid_min_1, block_min_1>>>(
                ncol, nlay, ngpt, igpt,
                ngas, nflav, ntemp, neta,
                nscale_lower,
                nminorlower,
                nminorklower,
                idx_h2o, idx_tropo,
                gpoint_flavor,
                kminor_lower,
                minor_limits_gpt_lower,
                first_last_minor_lower,
                minor_scales_with_density_lower,
                scale_by_complement_lower,
                idx_minor_lower,
                idx_minor_scaling_lower,
                kminor_start_lower,
                play, tlay, col_gas,
                fminor, jeta, jtemp,
                tropo, scalings_lower, tau);

        // Upper
        idx_tropo = 0;

        dim3 grid_min_2(nlay, ncol, 1), block_min_2;
        if (tunings.count("gas_optical_depths_minor_kernel_upper_rt") == 0)
        {
            Float* tau_tmp = Tools_gpu::allocate_gpu<Float>(nlay*ncol);
            std::tie(grid_min_2, block_min_2) = tune_kernel(
                   "gas_optical_depths_minor_kernel_upper_rt",
                   dim3{nlay, ncol, 1},
                   {1}, {32, 48, 64, 96, 128, 256, 384, 512}, {1},
                   gas_optical_depths_minor_kernel,
                   ncol, nlay, ngpt, igpt,
                   ngas, nflav, ntemp, neta,
                   nscale_upper,
                   nminorupper,
                   nminorkupper,
                   idx_h2o, idx_tropo,
                   gpoint_flavor,
                   kminor_upper,
                   minor_limits_gpt_upper,
                   first_last_minor_upper,
                   minor_scales_with_density_upper,
                   scale_by_complement_upper,
                   idx_minor_upper,
                   idx_minor_scaling_upper,
                   kminor_start_upper,
                   play, tlay, col_gas,
                   fminor, jeta, jtemp,
                   tropo, scalings_upper,
                   tau_tmp);
            Tools_gpu::free_gpu<Float>(tau_tmp);

            tunings["gas_optical_depths_minor_kernel_upper_rt"].first = grid_min_2;
            tunings["gas_optical_depths_minor_kernel_upper_rt"].second = block_min_2;
        }
        else
        {
            grid_min_2 = tunings["gas_optical_depths_minor_kernel_upper_rt"].first;
            block_min_2 = tunings["gas_optical_depths_minor_kernel_upper_rt"].second;
        }

        gas_optical_depths_minor_kernel<<<grid_min_2, block_min_2>>>(
                ncol, nlay, ngpt, igpt,
                ngas, nflav, ntemp, neta,
                nscale_upper,
                nminorupper,
                nminorkupper,
                idx_h2o, idx_tropo,
                gpoint_flavor,
                kminor_upper,
                minor_limits_gpt_upper,
                first_last_minor_upper,
                minor_scales_with_density_upper,
                scale_by_complement_upper,
                idx_minor_upper,
                idx_minor_scaling_upper,
                kminor_start_upper,
                play, tlay, col_gas,
                fminor, jeta, jtemp,
                tropo, scalings_upper, tau);

    }



    void Planck_source(
            const int ncol, const int nlay, const int nbnd, const int ngpt, const int igpt,
            const int nflav, const int neta, const int npres, const int ntemp,
            const int nPlanckTemp,
            const Float* tlay,
            const Float* tlev,
            const Float* tsfc,
            const int sfc_lay,
            const Float* fmajor,
            const int* jeta,
            const Bool* tropo,
            const int* jtemp,
            const int* jpress,
            const int* gpoint_bands,
            const int* band_lims_gpt,
            const Float* pfracin,
            const Float temp_ref_min, const Float totplnk_delta,
            const Float* totplnk,
            const int* gpoint_flavor,
            Float* sfc_src,
            Float* lay_src,
            Float* lev_src_inc,
            Float* lev_src_dec,
            Float* sfc_src_jac)
    {
        Tuner_map& tunings = Tuner::get_map();

        const Float delta_Tsurf = Float(1.);

        dim3 grid(ncol, nlay, 1), block;
        if (tunings.count("Planck_source_kernel_rt") == 0)
        {
            std::tie(grid, block) = tune_kernel(
                    "Planck_source_kernel_rt",
                    dim3{ncol, nlay, 1},
                    {16, 32, 48, 64, 96, 128, 256, 512}, {1, 2, 4, 8}, {1},
                    Planck_source_kernel,
                    ncol, nlay, nbnd, ngpt,
                    nflav, neta, npres, ntemp, nPlanckTemp, igpt,
                    tlay, tlev, tsfc, sfc_lay,
                    fmajor, jeta, tropo, jtemp,
                    jpress, gpoint_bands, band_lims_gpt,
                    pfracin, temp_ref_min, totplnk_delta,
                    totplnk, gpoint_flavor,
                    delta_Tsurf, sfc_src, lay_src,
                    lev_src_inc, lev_src_dec,
                    sfc_src_jac);

            tunings["Planck_source_kernel_rt"].first = grid;
            tunings["Planck_source_kernel_rt"].second = block;
        }
        else
        {
            grid = tunings["Planck_source_kernel_rt"].first;
            block = tunings["Planck_source_kernel_rt"].second;
        }

        Planck_source_kernel<<<grid, block>>>(
                ncol, nlay, nbnd, ngpt,
                nflav, neta, npres, ntemp, nPlanckTemp, igpt,
                tlay, tlev, tsfc, sfc_lay,
                fmajor, jeta, tropo, jtemp,
                jpress, gpoint_bands, band_lims_gpt,
                pfracin, temp_ref_min, totplnk_delta,
                totplnk, gpoint_flavor,
                delta_Tsurf,
                sfc_src, lay_src,
                lev_src_inc, lev_src_dec,
                sfc_src_jac);
    }
}

