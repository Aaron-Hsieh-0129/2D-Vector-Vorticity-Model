#include "Declare.hpp"
#include <iostream>
#if defined(GPU_POISSON)
#include <amgx_c.h>
#include <cuda_runtime.h>
#include <mpi.h>

void vvm::PoissonSolver::InitAMGX(vvm &model) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set and verify device before AMGX initialization
    if (cudaSetDevice(gpu_id) != cudaSuccess) {
        std::cerr << "Rank " << rank << ": Failed to set GPU " << gpu_id << " before AMGX init" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int device;
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": Expected GPU " << gpu_id << " but got " << device << " before AMGX init" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::cout << "Rank " << rank << ": Device set to GPU " << device << " before AMGX_initialize" << std::endl;

    // Initialize AMGX
    AMGX_initialize();

    // Immediately check and reset device after AMGX_initialize
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": GPU changed to " << device << " after AMGX_initialize, resetting to " << gpu_id << std::endl;
        if (cudaSetDevice(gpu_id) != cudaSuccess) {
            std::cerr << "Rank " << rank << ": Failed to reset GPU to " << gpu_id << " after AMGX_initialize" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cudaGetDevice(&device);
        std::cout << "Rank " << rank << ": Reset device to GPU " << device << " after AMGX_initialize" << std::endl;
    }

    // Proceed with AMGX setup
    std::string config_w = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"ilu\", \"solver\": \"ILU0\"}, \"scope\": \"main\", \"solver\": \"BICGSTAB\", \"tolerance\": 1e-12, \"max_iters\": 10000, \"monitor_residual\": 1, \"print_solve_stats\": 0}}";
    AMGX_config_create(&model.cfg_w, config_w.c_str());
    AMGX_resources_create_simple(&model.rsc_w, model.cfg_w);

    // Check device after resource creation
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": GPU changed to " << device << " after resource creation for w" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    AMGX_matrix_create(&model.A, model.rsc_w, AMGX_mode_dDDI);
    AMGX_vector_create(&model.b_vec_w, model.rsc_w, AMGX_mode_dDDI);
    AMGX_vector_create(&model.x_vec_w, model.rsc_w, AMGX_mode_dDDI);
    AMGX_solver_create(&model.solver_w, model.rsc_w, AMGX_mode_dDDI, model.cfg_w);

    std::string config_u = "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"scope\": \"main\", \"solver\": \"CG\", \"tolerance\": 1e-12, \"max_iters\": 10000, \"monitor_residual\": 1, \"print_solve_stats\": 0}}";
    AMGX_config_create(&model.cfg_u, config_u.c_str());
    AMGX_resources_create_simple(&model.rsc_u, model.cfg_u);

    // Check device after resource creation for u
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": GPU changed to " << device << " after resource creation for u" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    AMGX_matrix_create(&model.G, model.rsc_u, AMGX_mode_dDDI);
    AMGX_vector_create(&model.h_vec_u, model.rsc_u, AMGX_mode_dDDI);
    AMGX_vector_create(&model.y_vec_u, model.rsc_u, AMGX_mode_dDDI);
    AMGX_solver_create(&model.solver_u, model.rsc_u, AMGX_mode_dDDI, model.cfg_u);

    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": GPU changed to " << device << " after solver creation" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::cout << "Rank " << rank << ": Initialized AMGX on GPU " << gpu_id << std::endl;
}

void vvm::PoissonSolver::CleanupAMGX(vvm &model) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Ensure all solver operations are complete by destroying solvers first
    if (model.solver_w) {
        AMGX_solver_destroy(model.solver_w);
        model.solver_w = nullptr;
    }
    if (model.solver_u) {
        AMGX_solver_destroy(model.solver_u);
        model.solver_u = nullptr;
    }

    // Destroy vectors
    if (model.x_vec_w) {
        AMGX_vector_destroy(model.x_vec_w);
        model.x_vec_w = nullptr;
    }
    if (model.b_vec_w) {
        AMGX_vector_destroy(model.b_vec_w);
        model.b_vec_w = nullptr;
    }
    if (model.y_vec_u) {
        AMGX_vector_destroy(model.y_vec_u);
        model.y_vec_u = nullptr;
    }
    if (model.h_vec_u) {
        AMGX_vector_destroy(model.h_vec_u);
        model.h_vec_u = nullptr;
    }

    // Destroy matrices
    if (model.A) {
        AMGX_matrix_destroy(model.A);
        model.A = nullptr;
    }
    if (model.G) {
        AMGX_matrix_destroy(model.G);
        model.G = nullptr;
    }

    // Destroy resources
    if (model.rsc_w) {
        AMGX_resources_destroy(model.rsc_w);
        model.rsc_w = nullptr;
    }
    if (model.rsc_u) {
        AMGX_resources_destroy(model.rsc_u);
        model.rsc_u = nullptr;
    }

    // Destroy configurations
    if (model.cfg_w) {
        AMGX_config_destroy(model.cfg_w);
        model.cfg_w = nullptr;
    }
    if (model.cfg_u) {
        AMGX_config_destroy(model.cfg_u);
        model.cfg_u = nullptr;
    }

    // Finalize AMGX last
    AMGX_finalize();

    std::cout << "Rank " << rank << ": Cleaned up AMGX on GPU " << vvm::PoissonSolver::gpu_id << std::endl;
}

void vvm::PoissonSolver::InitPoissonMatrix(vvm &model) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int device;
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": Expected GPU " << gpu_id << " but got " << device << " in InitPoissonMatrix" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // For A (w)
    int size_w = (model.nx - 2) * (model.nz - 3);
    int max_nnz_w = 5 * size_w;
    model.row_ptr_w = new int[size_w + 1];
    model.col_idx_w = new int[max_nnz_w];
    model.values_w = new double[max_nnz_w];
    int nnz_w = 0, k = 1;

    for (int idx = 0; idx < size_w; idx++) {
        if (idx % (model.nx - 2) == 0) k++;
        model.row_ptr_w[idx] = nnz_w;

        // Down term (idx - (nx-2))
        if (idx >= (model.nx - 2)) {
            model.col_idx_w[nnz_w] = idx - (model.nx - 2);
            model.values_w[nnz_w] = model.rhow[k] / model.rhou[k - 1];
            nnz_w++;
        }
        // Left
        if (idx % (model.nx - 2) != 0) {
            model.col_idx_w[nnz_w] = idx - 1;
            model.values_w[nnz_w] = 1.0;
            nnz_w++;
        }
        // Boundary down (idx - (nx-3))
        if (idx >= (model.nx - 3) && (idx - (model.nx - 3)) % (model.nx - 2) == 0) {
            model.col_idx_w[nnz_w] = idx - (model.nx - 3);
            model.values_w[nnz_w] = 1.0;
            nnz_w++;
        }
        // Diagonal
        model.col_idx_w[nnz_w] = idx;
        model.values_w[nnz_w] = -(2. + (model.rhow[k] / model.rhou[k]) + (model.rhow[k] / model.rhou[k - 1]));
        nnz_w++;
        // Right
        if ((idx + 1) % (model.nx - 2) != 0) {
            model.col_idx_w[nnz_w] = idx + 1;
            model.values_w[nnz_w] = 1.0;
            nnz_w++;
        }
        // Boundary up (idx + (nx-3))
        if (idx % (model.nx - 2) == 0 && idx + (model.nx - 3) < size_w) {
            model.col_idx_w[nnz_w] = idx + (model.nx - 3);
            model.values_w[nnz_w] = 1.0;
            nnz_w++;
        }
        // Up term (idx + (nx-2))
        if (idx < (model.nx - 2) * (model.nz - 4)) {
            model.col_idx_w[nnz_w] = idx + (model.nx - 2);
            model.values_w[nnz_w] = model.rhow[k] / model.rhou[k];
            nnz_w++;
        }
    }
    model.row_ptr_w[size_w] = nnz_w;
    model.nnz_w = nnz_w;

    // For G (u) - Unchanged from original
    int size_u = model.nx - 2;
    int max_nnz_u = 3 * size_u + 2;
    model.row_ptr_u = new int[size_u + 1];
    model.col_idx_u = new int[max_nnz_u];
    model.values_u = new double[max_nnz_u];
    int nnz_u = 0;

    for (int i = 0; i < size_u; i++) {
        model.row_ptr_u[i] = nnz_u;
        model.col_idx_u[nnz_u] = i;
        model.values_u[nnz_u] = -2.0;
        nnz_u++;
        if (i != size_u - 1) {
            model.col_idx_u[nnz_u] = i + 1;
            model.values_u[nnz_u] = 1.0;
            nnz_u++;
        }
        if (i != 0) {
            model.col_idx_u[nnz_u] = i - 1;
            model.values_u[nnz_u] = 1.0;
            nnz_u++;
        }
        if (i == 0 && size_u > 1) {
            model.col_idx_u[nnz_u] = size_u - 1;
            model.values_u[nnz_u] = 1.0;
            nnz_u++;
        }
        else if (i == size_u - 1 && size_u > 1) {
            model.col_idx_u[nnz_u] = 0;
            model.values_u[nnz_u] = 1.0;
            nnz_u++;
        }
    }
    model.row_ptr_u[size_u] = nnz_u;
    model.nnz_u = nnz_u;

    AMGX_matrix_upload_all(model.A, size_w, model.nnz_w, 1, 1, model.row_ptr_w, model.col_idx_w, model.values_w, nullptr);
    AMGX_solver_setup(model.solver_w, model.A);

    AMGX_matrix_upload_all(model.G, size_u, model.nnz_u, 1, 1, model.row_ptr_u, model.col_idx_u, model.values_u, nullptr);
    AMGX_solver_setup(model.solver_u, model.G);
}

void vvm::PoissonSolver::cal_w(vvm &model, int p, int i, int j) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int device;
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": Expected GPU " << gpu_id << " but got " << device << " in cal_w" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size = (model.nx - 2) * (model.nz - 3);
    double* b = new double[size]();
    double* x = new double[size]();

    int count = 0;
    for (int k = 2; k <= model.nz - 2; k++) {
        for (int i = 1; i <= model.nx - 2; i++) {
            b[count] = model.rhow[k] * model.rhow[k] * (model.zetap[i + 1][k] - model.zetap[i][k]) * model.dx;
            count++;
        }
    }

    AMGX_vector_upload(model.b_vec_w, size, 1, b);
    double* initial_guess = new double[size]();
    AMGX_vector_upload(model.x_vec_w, size, 1, initial_guess);
    AMGX_solver_solve(model.solver_w, model.b_vec_w, model.x_vec_w);
    AMGX_vector_download(model.x_vec_w, x);

    AMGX_SOLVE_STATUS status;
    AMGX_solver_get_status(model.solver_w, &status);
    if (status != AMGX_SOLVE_SUCCESS) {
        std::cout << "Rank " << rank << ": W Solve Warning!!!!!!!!!!!! p, i, j = " << p << ", " << i << ", " << j << "\n";
        int iterations;
        AMGX_solver_get_iterations_number(model.solver_w, &iterations);
        std::cout << "Rank " << rank << ": W: #iterations: " << iterations << ", residual printed in AMGX stats\n";
    }

    int cnt = 0;
    for (int k = 2; k <= model.nz - 2; k++) {
        for (int i = 1; i <= model.nx - 2; i++) {
            model.w[i][k] = x[cnt] / model.rhow[k];
            cnt++;
        }
    }
    model.BoundaryProcess2D_westdown(model.w, model.nx, model.nz);

    delete[] b;
    delete[] x;
    delete[] initial_guess;

    // std::cout << "Rank " << rank << ": Completed cal_w on GPU " << gpu_id << std::endl;
}

void vvm::PoissonSolver::cal_u(vvm &model) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int device;
    cudaGetDevice(&device);
    if (device != gpu_id) {
        std::cerr << "Rank " << rank << ": Expected GPU " << gpu_id << " but got " << device << " in cal_u" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size = model.nx - 2;
    double* h = new double[size]();
    double* y = new double[size]();

    double tmp = 0.0;
    for (int i = 1; i <= model.nx - 2; i++) {
        h[i - 1] = -(0. - model.rhow[model.nz - 2] * model.w[i][model.nz - 2]) / model.rhou[model.nz - 2] * model.dx;
        tmp += h[i - 1];
    }
    tmp /= (model.nx - 2);
    for (int i = 0; i < size; i++) {
        h[i] -= tmp;
    }

    AMGX_vector_upload(model.h_vec_u, size, 1, h);
    double* initial_guess = new double[size]();
    AMGX_vector_upload(model.y_vec_u, size, 1, initial_guess);
    AMGX_solver_solve(model.solver_u, model.h_vec_u, model.y_vec_u);
    AMGX_vector_download(model.y_vec_u, y);

    AMGX_SOLVE_STATUS status;
    AMGX_solver_get_status(model.solver_u, &status);
    if (status != AMGX_SOLVE_SUCCESS) {
        std::cout << "Rank " << rank << ": U Solve Warning!!!!!!!!!!!!\n";
        int iterations;
        AMGX_solver_get_iterations_number(model.solver_u, &iterations);
        std::cout << "Rank " << rank << ": U: #iterations: " << iterations << ", residual printed in AMGX stats\n";
    }

    for (int i = 1; i <= model.nx - 2; i++) {
        model.xi[i] = y[i - 1];
    }
    model.xi[0] = model.xi[model.nx - 2];
    model.xi[model.nx - 1] = model.xi[1];

    for (int i = 1; i <= model.nx - 2; i++) {
        model.uxi[i] = (model.xi[i] - model.xi[i - 1]) * model.rdx;
    }
    // mean uxi
    tmp = 0.;
    for (int i = 1; i < model.nx-1; i++) {
        tmp += model.uxi[i];
    }
    tmp /= (model.nx - 2);
    for (int i = 1; i <= model.nx - 2; i++) {
        model.uxi[i] -= tmp;
    }
    model.uxi[0] = model.uxi[model.nx - 2];
    model.uxi[model.nx - 1] = model.uxi[1];

    for (int i = 1; i <= model.nx - 2; i++) {
        model.u[i][model.nz - 2] = model.uxi[i] + model.ubarTopp;
    }
    model.u[0][model.nz - 2] = model.u[model.nx - 2][model.nz - 2];
    model.u[model.nx - 1][model.nz - 2] = model.u[1][model.nz - 2];

    for (int i = 1; i <= model.nx - 2; i++) {
        double area = 0.0;
        for (int k = model.nz - 3; k >= 1; k--) {
            area += ((model.w[i][k + 1] - model.w[i - 1][k + 1]) * model.rdx - model.rhow[k + 1] * model.zetap[i][k + 1]) * -model.dz;
            model.u[i][k] = area + model.u[i][model.nz - 2];
        }
    }
    model.BoundaryProcess2D_center(model.u, model.nx, model.nz);

    delete[] h;
    delete[] y;
    delete[] initial_guess;
}

void vvm::PoissonSolver::pubarTop_pt(vvm &model) {
    double rhouwDown = 0.;
    double prhouwb_pz_rhob = 0.;
    for (int i = 1; i < model.nx-1; i++) {
        rhouwDown += 0.25 * (model.rhou[model.nz-2]*model.u[i][model.nz-2] + model.rhou[model.nz-3]*model.u[i][model.nz-3]) * (model.w[i][model.nz-2] + model.w[i-1][model.nz-2]);
    }
    rhouwDown /= ((double) (model.nx - 2.));

    prhouwb_pz_rhob = (- rhouwDown) * model.rdz / model.rhou[model.nz-2];
    model.dubarTop_advect[(model.step+1)%2] = -prhouwb_pz_rhob;
    if (model.step == 0) model.dubarTop_advect[0] = model.dubarTop_advect[1];
    model.ubarTopp = model.ubarTop + 1.5*model.dt*model.dubarTop_advect[(model.step+1)%2] - 0.5*model.dt*model.dubarTop_advect[model.step%2];
    return;
}
#endif
