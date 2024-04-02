#include "Declare.hpp"
#include <ios>
#include <petsc.h>
#include <petscksp.h>
#include <petscvec.h>
#include <iostream>
#include <petscviewer.h>


#if defined(STREAMFUNCTION)
void vvm::PoissonSolver::calpsiuw(vvm &model) {
    #if defined(PETSC)
        Vec x, b;
        Mat A;
        KSP ksp;
       
        int j[5];
        double v[5];
        
        PetscInt m = (model.nx-2)*(model.nz-3);
    
        // Create vectors
        VecCreate(PETSC_COMM_WORLD, &x);
        VecSetSizes(x, PETSC_DECIDE, m);
        VecSetFromOptions(x);
        VecDuplicate(x, &b);

        // Create the matrix
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, m);
        MatSetFromOptions(A);
        MatSetUp(A);

        int idx_i = 0, k = 1;
        for (int i = 0; i < (model.nx-2)*(model.nz-3); i++) {    
            // Height
            if (i % (model.nx-2) == 0) k++;

            // Diagonal
            v[0] = -(2./model.rhow[k] + 1./model.rhou[k] + 1./model.rhou[k-1]);
            j[0] = i;
            MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);

            
            // D
            // Fill left right except the point near the boundaries
            if (i % (model.nx-2) != 0 && i % (model.nx-2) != (model.nx-2)-1) {
                v[0] = 1./model.rhow[k]; v[1] = 1./model.rhow[k]; 
                j[0] = i-1; j[1] = i+1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
            // Fill the leftest point
            if (i % (model.nx-2) == 0) {
                v[0] = 1./model.rhow[k]; v[1] = 1./model.rhow[k];
                j[0] = i+1, j[1] = i+(model.nx-2)-1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
        
            // Fill the rightest point
            if (i % (model.nx-2) == (model.nx-2)-1) {
                v[0] = 1./model.rhow[k]; v[1] = 1./model.rhow[k];
                j[0] = i-(model.nx-2)+1, j[1] = i-1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
            
            // E
            if (i < (model.nx-2)*(model.nz-3-1)) {
                v[0] = 1./model.rhou[k];
                j[0] = i+(model.nx-2);
                MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);
            }

            // F
            if (i >= model.nx-2) {
                v[0] = 1./model.rhou[k-1];
                j[0] = i-(model.nx-2);
                MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);
            }
            
            
            idx_i = (i % (model.nx-2)) + 1;
            
            double bval = model.zetap[idx_i][k] * dx * dx;
            VecSetValues(b, 1, &i, &bval, INSERT_VALUES);
        }

        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b); VecAssemblyEnd(b);
        // VecView(b, PETSC_VIEWER_STDOUT_WORLD);

        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSetTolerances(ksp, 1E-35, 1E-15, 1E2, 10000);
        KSPSolve(ksp, b, x);
        // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
        int iterNum = 0.;
        double normError = 0.;
        KSPGetIterationNumber(ksp, &iterNum);
        KSPGetResidualNorm(ksp, &normError);
        std::cout << "Solving psi: " << ", Norm of error = " << std::scientific << normError << ", Iterations = " << iterNum << std::endl;

        PetscScalar *x_arr;
        VecGetArray(x, &x_arr);

        
        int K = 1;
        int I = 0;
        for (PetscInt i = 0; i < m; i++) {
            I = (i % (model.nx-2)) + 1;
            if (i % (model.nx-2) == 0) K++;
            model.psi[I][K] = x_arr[i];
        }

        MatDestroy(&A);
        VecDestroy(&b);
        VecDestroy(&x);
        KSPDestroy(&ksp);

    #else
        // eigen solver
        Eigen::VectorXd x((model.nx-2)*(model.nz-3)), b((model.nx-2)*(model.nz-3));

        // b
        int count = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                // b(count) = model.rhow[k] * model.zetap[i][k] * dx * dx;
                b(count) = model.zetap[i][k] * dx * dx;
                count++;
            }
        }

        auto A = model.A;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
        solver.setTolerance(1e-31);
        solver.setMaxIterations(10000);
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cout << "Decompose Warning!!!!!!!!!!!" << std::endl;
        }
        x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            std::cout << "Solve Warning!!!!!!!!!!!" << std::endl;
        }
        std::cout << "psi:  #iterations:     " << solver.iterations() << ", estimated error: " << std::scientific <<  solver.error() << std::endl;

        int cnt = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                model.psi[i][k] = x[cnt];
                cnt++;
            }
        }
    #endif
    for (int k = 2; k <= model.nz-2; k++) {
        model.psi[0][k] = model.psi[model.nx-2][k];
        model.psi[model.nx-1][k] = model.psi[1][k];
    }
    for (int i = 0; i <= model.nx-1; i++) {
        model.psi[i][0] = model.psi[i][1] = model.psi[i][model.nz-1] = 0.; 
    }
    model.BoundaryProcess2D_westdown(model.psi);

    for (int i = 1; i <= model.nx-2; i++) {
        for (int k = 1; k <= model.nz-2; k++) {
            model.u[i][k] = -(model.rhow[k+1]*model.psi[i][k+1] - model.rhow[k]*model.psi[i][k]) * model.rdz / model.rhou[k];
            model.w[i][k] = (model.psi[i+1][k]-model.psi[i][k]) * model.rdx;
        }
    }
	for (int i = 0; i < model.nz; i++) {
        model.w[i][model.nz-1] = model.w[i][0] = model.w[i][1] = 0.;
    }
    model.BoundaryProcess2D_center(model.u);
    model.BoundaryProcess2D_westdown(model.w);
    
    return;
}
#else


void vvm::PoissonSolver::cal_w(vvm &model) {
    #if defined(PETSC) 
        Vec x, b;
        Mat A;
        KSP ksp;
       
        int j[5];
        double v[5];
        
        PetscInt m = (model.nx-2)*(model.nz-3);
    
        // Create vectors
        VecCreate(PETSC_COMM_WORLD, &x);
        VecSetSizes(x, PETSC_DECIDE, m);
        VecSetFromOptions(x);
        VecDuplicate(x, &b);
        #if defined(POISSONTEST)
            Vec x_ans;
            VecDuplicate(x, &x_ans);
        #endif

        // Create the matrix
        MatCreate(PETSC_COMM_WORLD, &A);
        MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, m);
        MatSetFromOptions(A);
        MatSetUp(A);

        int idx_i = 0, k = 1;
        for (int i = 0; i < (model.nx-2)*(model.nz-3); i++) {    
            // Height
            if (i % (model.nx-2) == 0) k++;

            // Diagonal
            // v[0] = -(2. + (model.rhow[k]/model.rhou[k]) + (model.rhow[k]/model.rhou[k-1])) * model.rdx2 + POISSONPARAM;
            v[0] = -(2. + (model.rhow[k]/model.rhou[k]) + (model.rhow[k]/model.rhou[k-1])) + POISSONPARAM;
            j[0] = i;
            MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);

            
            // D
            // Fill left right except the point near the boundaries
            if (i % (model.nx-2) != 0 && i % (model.nx-2) != (model.nx-2)-1) {
                // v[0] = 1.*model.rdx2; v[1] = 1.*model.rdx2; 
                v[0] = 1.; v[1] = 1.; 
                j[0] = i-1; j[1] = i+1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
            // Fill the leftest point
            if (i % (model.nx-2) == 0) {
                // v[0] = 1.*model.rdx2; v[1] = 1.*model.rdx2;
                v[0] = 1.; v[1] = 1.;
                j[0] = i+1, j[1] = i+(model.nx-2)-1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
        
            // Fill the rightest point
            if (i % (model.nx-2) == (model.nx-2)-1) {
                // v[0] = 1.*model.rdx2; v[1] = 1.*model.rdx2;
                v[0] = 1.; v[1] = 1.;
                j[0] = i-(model.nx-2)+1, j[1] = i-1;
                MatSetValues(A, 1, &i, 2, j, v, INSERT_VALUES);
            }
            
            // E
            if (i < (model.nx-2)*(model.nz-3-1)) {
                // v[0] = model.rhow[k]/model.rhou[k] * model.rdx2;
                v[0] = model.rhow[k]/model.rhou[k];
                j[0] = i+(model.nx-2);
                MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);
            }

            // F
            if (i >= model.nx-2) {
                // v[0] = model.rhow[k]/model.rhou[k-1] * model.rdx2;
                v[0] = model.rhow[k]/model.rhou[k-1];
                j[0] = i-(model.nx-2);
                MatSetValues(A, 1, &i, 1, j, v, INSERT_VALUES);
            }
            
            #if defined(POISSONTEST)
                double x_ansval = exp(sin(i));
                VecSetValues(x_ans, 1, &i, &x_ansval, INSERT_VALUES);
            #else
                idx_i = (i % (model.nx-2)) + 1;
                // double bval = model.rhow[k]*model.rhow[k] * (model.zetap[idx_i+1][k] - model.zetap[idx_i][k]) * model.rdx;
                double bval = model.rhow[k]*model.rhow[k] * (model.zetap[idx_i+1][k] - model.zetap[idx_i][k]) * dx;
                VecSetValues(b, 1, &i, &bval, INSERT_VALUES);
            #endif
        }

        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(b); VecAssemblyEnd(b);
        // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
        #if defined(POISSONTEST)
            VecAssemblyBegin(x_ans); VecAssemblyEnd(x_ans);
            MatMult(A, x_ans, b);
        #endif

        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, A, A);
        KSPSetFromOptions(ksp);
        KSPSetTolerances(ksp, 1E-30, 1E-12, 1E2, 10000);
        KSPSolve(ksp, b, x);

        #if defined(POISSONTEST)
            double error_norm;
            VecAXPY(x, -1, x_ans);
            VecNorm(x, NORM_2, &error_norm);
            printf("Grid %d: Error = %1e\n", m, error_norm);
        #else
            PetscScalar *x_arr;
            VecGetArray(x, &x_arr);
            
            int K = 1;
            int I = 0;
            for (PetscInt i = 0; i < m; i++) {
                I = (i % (model.nx-2)) + 1;
                if (i % (model.nx-2) == 0) K++;
                model.w[I][K] = x_arr[i] / model.rhow[K];
            }

            
            for (int k = 1; k <= NZ-2; k++) {
                model.w[0][k] = model.w[NX-2][k];
                model.w[NX-1][k] = model.w[1][k];
            }
            for (int i = 0; i <= NX-1; i++) {
                model.w[i][0] = model.w[i][1] = model.w[i][NZ-1] = 0.;
            }
        #endif
        int iterNum = 0.;
        double normError = 0.;
        KSPGetIterationNumber(ksp, &iterNum);
        KSPGetResidualNorm(ksp, &normError);
        std::cout << "Solving w: " << ", Norm of error = " << std::scientific << normError << ", Iterations = " << iterNum << std::endl;

        MatDestroy(&A);
        VecDestroy(&b);
        VecDestroy(&x);
        KSPDestroy(&ksp);
    #else
        // eigen solver
        Eigen::VectorXd x((model.nx-2)*(model.nz-3)), b((model.nx-2)*(model.nz-3));
    
        // b
        int count = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                b(count) = model.rhow[k]*model.rhow[k] * (model.zetap[i+1][k] - model.zetap[i][k]) * dx;
                count++;
            }
        }

        // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
        solver.setTolerance(1e-30);
        solver.setMaxIterations(10000);
        Eigen::SparseMatrix<double> A = model.A;
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cout << "WWWWWWWWW Warning!!!!!!!!!!!!" << std::endl;
        }
        x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            std::cout << "W Solve Warning!!!!!!!!!!!!" << std::endl;
        }
        std::cout << "W:  #iterations:     " << solver.iterations() << ", estimated error: " << std::scientific <<  solver.error() << std::endl;

    
        int cnt = 0;
        for (int k = 2; k <= model.nz-2; k++) {
            for (int i = 1; i <= model.nx-2; i++) {
                model.w[i][k] = x[cnt] / model.rhow[k];
                cnt++;
            }
        }
        model.BoundaryProcess2D_westdown(model.w);
    #endif
    return;
}

void vvm::PoissonSolver::cal_u(vvm &model) {
    #if defined(PETSC)
        Vec y, h;
        Mat G;
        KSP ksp;
       
        int j[3];
        double v[3];
        
        PetscInt m = model.nx-2;
    
        // Create vectors
        VecCreate(PETSC_COMM_WORLD, &y);
        VecSetSizes(y, PETSC_DECIDE, m);
        VecSetFromOptions(y);
        VecDuplicate(y, &h);
        #if defined(POISSONTEST)
            Vec y_ans;
            VecDuplicate(y, &y_ans);
        #endif

        // Create the matrix
        MatCreate(PETSC_COMM_WORLD, &G);
        MatSetSizes(G, PETSC_DECIDE, PETSC_DECIDE, m, m);
        MatSetFromOptions(G);
        MatSetUp(G);

        int startind = 0, endind = 0;
        for (int i = 0; i < model.nx-2; i++) {
            if (i == 0) {
                // v[0] = -2.*model.rdx2 + POISSONPARAMU; v[1] = 1.*model.rdx2; v[2] = 1.*model.rdx2;
                v[0] = -2. + POISSONPARAMU; v[1] = 1.; v[2] = 1.;
                j[0] = 0; j[1] = 1; j[2] = (model.nx-2)-1;
                MatSetValues(G, 1, &i, 3, j, v, INSERT_VALUES);
            }
            else if (i == (model.nx-2)-1) {
                // v[0] = 1.*model.rdx2; v[1] = 1.*model.rdx2; v[2] = -2.*model.rdx2 + POISSONPARAMU;
                v[0] = 1.; v[1] = 1.; v[2] = -2. + POISSONPARAMU;
                j[0] = 0; j[1] = (model.nx-2)-2; j[2] = (model.nx-2)-1;
                MatSetValues(G, 1, &i, 3, j, v, INSERT_VALUES);
            }
            else {
                // v[0] = 1.*model.rdx2; v[1] = -2.*model.rdx2 + POISSONPARAMU; v[2] = 1.*model.rdx2;
                v[0] = 1.; v[1] = -2. + POISSONPARAMU; v[2] = 1.;
                j[0] = i-1; j[1] = i; j[2] = i+1;
                MatSetValues(G, 1, &i, 3, j, v, INSERT_VALUES);
            }

            #if defined(POISSONTEST)
                double y_ansval = exp(sin(i*i));
                VecSetValues(y_ans, 1, &i, &y_ansval, INSERT_VALUES);
            #else
                // double hval = -(0. - model.rhow[model.nz-2] * model.w[i+1][model.nz-2]) / model.rhou[model.nz-2] * model.rdz;
                double hval = -(0. - model.rhow[model.nz-2] * model.w[i+1][model.nz-2]) / model.rhou[model.nz-2] * dz;
                VecSetValues(h, 1, &i, &hval, INSERT_VALUES);
            #endif
        }
        MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(h); VecAssemblyEnd(h);
        // MatView(G, PETSC_VIEWER_STDOUT_WORLD);
        // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
        #if defined(POISSONTEST)
            VecAssemblyBegin(y_ans); VecAssemblyEnd(y_ans);
            MatMult(G, y_ans, h);
        #endif

        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, G, G);
        KSPSetFromOptions(ksp);
        KSPSetTolerances(ksp, 1E-15, 1E-9, 1E2, 10000);
        KSPSolve(ksp, h, y);
        // VecView(x, PETSC_VIEWER_STDOUT_WORLD);

        #if defined(POISSONTEST)
            double error_norm;
            VecAXPY(y, -1, y_ans);
            VecNorm(y, NORM_2, &error_norm);
            printf("Grid %d: Error = %1e\n", m, error_norm);
        #else
            PetscScalar *xi;
            VecGetArray(y, &xi);

            for (PetscInt i = 0; i < model.nx-2; i++) {
                model.xi[i+1] = xi[i];
            }
            model.xi[0] = model.xi[model.nx-2];
            model.xi[model.nx-1] = model.xi[1];
        #endif

        int iterNum = 0.;
        double normError = 0.;
        KSPGetIterationNumber(ksp, &iterNum);
        KSPGetResidualNorm(ksp, &normError);
        std::cout << "Solving u: " << ", Norm of error = " << std::scientific << normError << ", Iterations = " << iterNum << std::endl;

        MatDestroy(&G);
        VecDestroy(&h);
        VecDestroy(&y);
        KSPDestroy(&ksp);
    #else
        Eigen::VectorXd y(model.nx-2), h(model.nx-2);
        // h
        for (int i = 1; i <= model.nx-2; i++) {
            h(i-1) = -(0. - model.rhow[model.nz-2] * model.w[i][model.nz-2]) / model.rhou[model.nz-2] * dx;
        }

        Eigen::SparseMatrix<double> G = model.G;
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver_xi;
        // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_xi;
        // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver_xi;
        solver_xi.setTolerance(1e-30);
        solver_xi.setMaxIterations(10000);
        solver_xi.compute(G);
        if (solver_xi.info() != Eigen::Success) {
            std::cout << "UUUUUUUU Warning!!!!!!!!!!!!" << std::endl;
        }
        y = solver_xi.solve(h);
        if (solver_xi.info() != Eigen::Success) {
            std::cout << "U Solve Warning!!!!!!!!!!!!" << std::endl;
        }
        std::cout << "U:  #iterations:     " << solver_xi.iterations() << ", estimated error: " << std::scientific <<  solver_xi.error() << std::endl;

        // xi
        for (int i = 1; i <= model.nx-2; i++) {
            model.xi[i] = y[i-1];
        }
        model.xi[0] = model.xi[model.nx-2];
        model.xi[model.nx-1] = model.xi[1];
    #endif

    // uxi
    for (int i = 1; i <= model.nx-2; i++) {
        model.uxi[i] = (model.xi[i] - model.xi[i-1]) * model.rdx;
    }
    model.uxi[0] = model.uxi[model.nx-2];
    model.uxi[model.nx-1] = model.uxi[1];

    // u_top
    for (int i = 1; i <= model.nx-2; i++) {
        model.u[i][model.nz-2] = model.uxi[i] + model.ubarTopp;
    }
    model.u[0][model.nz-2] = model.u[model.nx-2][model.nz-2];
    model.u[model.nx-1][model.nz-2] = model.u[1][model.nz-2];

    // u
    for (int i = 1; i <= model.nx-2; i++) {
        double area = 0.;
        for (int k = model.nz-3; k >= 1; k--) {
            // area += ((0.5*(model.rhow[k+2]*model.w[i][k+2] + model.rhow[k+1]*model.w[i][k+1]) - 0.5*(model.rhow[k+2]*model.w[i-1][k+2] + model.rhow[k+1]*model.w[i-1][k+1])) / model.rhou[k+1] * model.rdx - 0.5*(model.rhow[k+2]*model.zetap[i][k+2] + model.rhow[k+1]*model.zetap[i][k+1])) * -dz;
            
            area += ((model.w[i][k+1]-model.w[i-1][k+1]) * model.rdx - model.rhow[k+1]*model.zetap[i][k+1]) * -dz;
            // area += ((model.w[i][k+1]-model.w[i-1][k+1]) * model.rdx - model.zetap[i][k+1]) * -dz;
            
            /*
            if (i == model.nz-3 || i == 1) area += ((model.w[i][k+1]-model.w[i-1][k+1]) * model.rdx - model.rhow[k+1]*model.zetap[i][k+1]) * -dz;
            else {
                double fa = ((0.5*(model.w[i][k+2] + model.w[i][k+1]) - 0.5*(model.w[i-1][k+2] + model.w[i-1][k+1])) * model.rdx - 0.5*(model.rhow[k+2]*model.zetap[i][k+2] + model.rhow[k+1]*model.zetap[i][k+1]));
                double fb = ((0.5*(model.w[i][k+1] + model.w[i][k]) - 0.5*(model.w[i-1][k+1] + model.w[i-1][k])) * model.rdx - 0.5*(model.rhow[k+1]*model.zetap[i][k+1] + model.rhow[k]*model.zetap[i][k]));
                double f_mid = ((model.w[i][k+1]-model.w[i-1][k+1]) * model.rdx - model.rhow[k+1]*model.zetap[i][k+1]);
                area += (-dz) / 6. * (fa + 4*f_mid + fb);
            }
            */
            model.u[i][k] = area + model.u[i][model.nz-2];
        }
    }
    model.BoundaryProcess2D_center(model.u);
    return;
}

void vvm::PoissonSolver::pubarTop_pt(vvm &model) {
    double rhouwUp = 0., rhouwDown = 0.;
    double prhouwb_pz_rhob = 0.;
    for (int i = 1; i < model.nx-1; i++) {
        rhouwUp += 0.25 * (model.rhou[model.nz-1]*model.u[i][model.nz-1] + model.rhou[model.nz-2]*model.u[i][model.nz-2]) * (model.w[i][model.nz-1] + model.w[i-1][model.nz-1]);
        rhouwDown += 0.25 * (model.rhou[model.nz-2]*model.u[i][model.nz-2] + model.rhou[model.nz-3]*model.u[i][model.nz-3]) * (model.w[i][model.nz-2] + model.w[i-1][model.nz-2]);
    }
    rhouwUp /= ((double) (model.nx - 2.));
    rhouwDown /= ((double) (model.nx - 2.));
    prhouwb_pz_rhob = (rhouwUp - rhouwDown) * model.rdz / model.rhou[model.nz-2];
    model.ubarTopp = model.ubarTopm + model.d2t * (-prhouwb_pz_rhob);
    return;
}
#endif

#ifndef PETSC
void vvm::PoissonSolver::InitPoissonMatrix(vvm &model) {
    #if defined(STREAMFUNCTION)
        typedef Eigen::Triplet<double> T;
        // ###########################################
        // For solving w
        // A: i = 1~model.nx-2, k = 2~model.nz-2
        // ###########################################
        /*
        int k = 1;
        std::vector<T> coeff;
        for (int idx = 1; idx < (model.nx-2)*(model.nz-3); idx++) {
            // Height
            if (idx % (model.nx-2) == 1) k++;

            // D
            coeff.push_back(T(idx-1, idx-1, -4.));

            // left/right: -1
            if ((idx-1) % (model.nx-2) != 0) coeff.push_back(T(idx-1, idx-2, 1.));
            if (idx % (model.nx-2) != 0) coeff.push_back(T(idx-1, idx, 1.));

            // Boundary
            if ((idx-1) % (model.nx-2) == 0) {
                coeff.push_back(T(idx-1, idx-1+(model.nx-3), 1.));
                coeff.push_back(T(idx-1+(model.nx-3), idx-1, 1.));
            }
        }
        // Last row of A
        coeff.push_back(T((model.nx-2)*(model.nz-3)-1, (model.nx-2)*(model.nz-3)-1, -4.));
        coeff.push_back(T((model.nx-2)*(model.nz-3)-1, (model.nx-2)*(model.nz-3)-1-1, 1.)); // left

        k = 1;
        for (int idx = 1; idx <= (model.nx-2)*(model.nz-4); idx++) {
            // Height
            if (idx % (model.nx-2) == 1) k++;

            // E
            coeff.push_back(T(idx-1, idx+(model.nx-2)-1, 1. + 0.5*(model.rhou[k] - model.rhou[k-1]) / model.rhow[k]));
            
            // F (the k of row should be minus by 1 because it starts from k-1)
            coeff.push_back(T(idx+(model.nx-2)-1, idx-1, 1. - 0.5*(model.rhou[k] - model.rhou[k-1]) / model.rhow[k]));
        }
        */


        std::vector<T> coeff;
        int k = 1; 
        for (int idx = 0; idx < (model.nx-2)*(model.nz-3); idx++) {
            // Height
            if (idx % (model.nx-2) == 0) k++;

            coeff.push_back(T(idx, idx, -2./model.rhow[k] - (1./model.rhou[k] + 1./model.rhou[k-1]) + POISSONPARAM ));

            // left/right = 1
            if (idx % (model.nx-2) != 0) coeff.push_back(T(idx, idx-1, 1./model.rhow[k])); // left
            if ((idx+1) % (model.nx-2) != 0) coeff.push_back(T(idx, idx+1, 1./model.rhow[k])); // right

            // Boundary
            if (idx % (model.nx-2) == 0) {
                coeff.push_back(T(idx, idx+(model.nx-3), 1./model.rhow[k]));
                coeff.push_back(T(idx+(model.nx-3), idx, 1./model.rhow[k]));
            }
        }
        
        k = 1;
        for (int idx = 0; idx < (model.nx-2)*(model.nz-3); idx++) {
            // Height
            if (idx % (model.nx-2) == 0) k++;

            // E
            if (idx < (model.nx-2)*(model.nz-4)) coeff.push_back(T(idx, idx+(model.nx-2), 1./model.rhou[k]));
            
            // F
            if (idx >= (model.nx-2)) coeff.push_back(T(idx, idx-(model.nx-2), 1./model.rhou[k-1]));
        }
        model.A.setFromTriplets(coeff.begin(), coeff.end());
    #else


    // ###########################################
    // For solving w
    // A: i = 1~model.nx-2, k = 1~model.nz-3
    // ###########################################
    typedef Eigen::Triplet<double> T;
    std::vector<T> coeff;
    int k = 1; 
    for (int idx = 0; idx < (model.nx-2)*(model.nz-3); idx++) {
        // Height
        if (idx % (model.nx-2) == 0) k++;

        // coeff.push_back(T(idx, idx, -(2. + (model.rhow[k]/model.rhou[k]) + (model.rhow[k]/model.rhou[k-1]))*model.rdx2 + POISSONPARAM ));
        coeff.push_back(T(idx, idx, -(2. + (model.rhow[k]/model.rhou[k]) + (model.rhow[k]/model.rhou[k-1])) + POISSONPARAM ));
        // coeff.push_back(T(idx, idx, -4.)); // test

        // left/right = 1
        // if (idx % (model.nx-2) != 0) coeff.push_back(T(idx, idx-1, 1.*model.rdx2)); // left
        if (idx % (model.nx-2) != 0) coeff.push_back(T(idx, idx-1, 1.)); // left
        // if ((idx+1) % (model.nx-2) != 0) coeff.push_back(T(idx, idx+1, 1.*model.rdx2)); // right
        if ((idx+1) % (model.nx-2) != 0) coeff.push_back(T(idx, idx+1, 1.)); // right
        // if (idx % (model.nx-2) != 0) coeff.push_back(T(idx, idx-1, 1.)); // left_test
        // if ((idx+1) % (model.nx-2) != 0) coeff.push_back(T(idx, idx+1, 1.)); // right right

        // Boundary
        if (idx % (model.nx-2) == 0) {
            // coeff.push_back(T(idx, idx+(model.nx-3), 1.*model.rdx2));
            coeff.push_back(T(idx, idx+(model.nx-3), 1.));
            // coeff.push_back(T(idx+(model.nx-3), idx, 1.*model.rdx2));
            coeff.push_back(T(idx+(model.nx-3), idx, 1.));
            // coeff.push_back(T(idx, idx+(model.nx-3), 1.)); // test
            // coeff.push_back(T(idx+(model.nx-3), idx, 1.)); // test
        }
    }
    
    k = 1;
    for (int idx = 0; idx < (model.nx-2)*(model.nz-3); idx++) {
        // Height
        if (idx % (model.nx-2) == 0) k++;

        // E
        // if (idx < (model.nx-2)*(model.nz-4)) coeff.push_back(T(idx, idx+(model.nx-2), model.rhow[k]/model.rhou[k] * model.rdx2));
        if (idx < (model.nx-2)*(model.nz-4)) coeff.push_back(T(idx, idx+(model.nx-2), model.rhow[k]/model.rhou[k]));
        // if (idx < (model.nx-2)*(model.nz-3))coeff.push_back(T(idx, idx+(model.nx-2), 1.)); // test
        
        // F
        // if (idx >= (model.nx-2)) coeff.push_back(T(idx, idx-(model.nx-2), model.rhow[k]/model.rhou[k-1] * model.rdx2));
        if (idx >= (model.nx-2)) coeff.push_back(T(idx, idx-(model.nx-2), model.rhow[k]/model.rhou[k-1]));
        // if (idx >= (model.nx-2)*2) coeff.push_back(T(idx, idx-(model.nx-2), 1.)); // test
    }

    model.A.setFromTriplets(coeff.begin(), coeff.end());
    
    /*
    // eigen solver
    Eigen::VectorXd x((model.nx-2)*(model.nz-3)), b((model.nx-2)*(model.nz-3)), x_ans((model.nx-2)*(model.nz-3));
    
    // b
    int count = 0;
    for (int k = 0; k <= model.nz-4; k++) {
        for (int i = 0; i <= model.nx-3; i++) {
            x_ans(count) = exp(sin(i+k));
            count++;
        }
    }
    b = model.A * x_ans;

    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.setTolerance(1e-20);
    solver.setMaxIterations(10000);
    Eigen::SparseMatrix<double> A = model.A;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cout << "WWWWWWWWW Warning!!!!!!!!!!!!" << std::endl;
    }
    x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cout << "W Solve Warning!!!!!!!!!!!!" << std::endl;
    }
    */
    // std::cout << std::setprecision(5) << model.A << std::endl;

    // ###########################################
    // For solving u
    // G: i = 1~model.nx-2
    // ###########################################
    std::vector<T> coeff_xi;

    for (int i = 0; i <= model.nx-3; i++) {
        // D
        // coeff_xi.push_back(T(i, i, -2.*model.rdx2+POISSONPARAM));
        coeff_xi.push_back(T(i, i, -2. + POISSONPARAMU));
        // if (i != model.nx-3) coeff_xi.push_back(T(i, i+1, 1.*model.rdx2));
        if (i != model.nx-3) coeff_xi.push_back(T(i, i+1, 1.));
        // if (i != 0) coeff_xi.push_back(T(i, i-1, 1.*model.rdx2));
        if (i != 0) coeff_xi.push_back(T(i, i-1, 1.));
    }
    // Boundary
    // coeff_xi.push_back(T(0, model.nx-3, 1.*model.rdx2));
    coeff_xi.push_back(T(0, model.nx-3, 1.));
    // coeff_xi.push_back(T(model.nx-3, 0, 1.*model.rdx2));
    coeff_xi.push_back(T(model.nx-3, 0, 1.));
    model.G.setFromTriplets(coeff_xi.begin(), coeff_xi.end());
    
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error()      << std::endl;
    // std::cout << ((x_ans-x).cwiseAbs()).sum() << std::endl;
    /*
    Eigen::VectorXd y(model.nx-2), h(model.nx-2), y_ans(model.nx-2);
    // h
    for (int i = 0; i <= model.nx-3; i++) {
        y_ans(i) = exp(sin(i));
    }
    h = model.G * y_ans;

    Eigen::SparseMatrix<double> G = model.G;
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver_xi;
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_xi;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver_xi;
    // solver_xi.setTolerance(1e-20);
    solver_xi.compute(G);
    if (solver_xi.info() != Eigen::Success) {
        std::cout << "Decompose Warning!!!!!!!!!!!!" << std::endl;
    }
    y = solver_xi.solve(h);
    if (solver_xi.info() != Eigen::Success) {
        std::cout << "U Solve Warning!!!!!!!!!!!!" << std::endl;
    }

    std::cout << ((y_ans-y).cwiseAbs()).sum() << std::endl;
    */	
    #endif
    return;
}
#endif
