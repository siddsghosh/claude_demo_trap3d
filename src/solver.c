#include "solver.h"
#include <stdio.h>
#include <string.h>

/**
 * Solve linear system using PETSc KSP
 */
int solver_solve(Mat K, Vec F, Vec *x) {
    PetscErrorCode ierr;
    KSP ksp;
    PC pc;
    PetscInt its;
    PetscReal rnorm;

    // Create solution vector
    ierr = VecDuplicate(F, x); CHKERRQ(ierr);

    // Create KSP context
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, K, K); CHKERRQ(ierr);

    // Set solver type (CG for symmetric positive definite)
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);

    // Set preconditioner (GAMG for 3D problems)
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCGAMG); CHKERRQ(ierr);

    // Set tolerances
    ierr = KSPSetTolerances(ksp, 1.0e-10, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
    CHKERRQ(ierr);

    // Allow options from command line
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    // Solve
    printf("Solving linear system...\n");
    ierr = KSPSolve(ksp, F, *x); CHKERRQ(ierr);

    // Get convergence info
    ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);

    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);

    if (reason < 0) {
        printf("WARNING: KSP did not converge (reason: %d)\n", reason);
    } else {
        printf("KSP converged in %d iterations, residual norm = %.3e\n", its, rnorm);
    }

    // Cleanup
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

    return 0;
}

/**
 * Write solution to CSV file
 */
int solver_write_csv(const char *filename, const Mesh *mesh, Vec x) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", filename);
        return 1;
    }

    // Write header
    fprintf(fp, "x,y,z,phi\n");

    // Get array access
    PetscErrorCode ierr;
    const PetscScalar *sol_array;
    ierr = VecGetArrayRead(x, &sol_array); CHKERRQ(ierr);

    // Write node data
    for (size_t i = 0; i < mesh->num_nodes; i++) {
        double xc = mesh->nodes[i * 3 + 0];
        double yc = mesh->nodes[i * 3 + 1];
        double zc = mesh->nodes[i * 3 + 2];
        double phi = PetscRealPart(sol_array[i]);

        fprintf(fp, "%.12e,%.12e,%.12e,%.12e\n", xc, yc, zc, phi);
    }

    ierr = VecRestoreArrayRead(x, &sol_array); CHKERRQ(ierr);

    fclose(fp);
    printf("Wrote CSV output to %s\n", filename);

    return 0;
}

/**
 * Write solution to VTK file
 * Legacy VTK format for tetrahedral mesh
 */
int solver_write_vtk(const char *filename, const Mesh *mesh, Vec x) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", filename);
        return 1;
    }

    // VTK header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Poisson FEM Solution\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Points
    fprintf(fp, "POINTS %zu double\n", mesh->num_nodes);
    for (size_t i = 0; i < mesh->num_nodes; i++) {
        fprintf(fp, "%.12e %.12e %.12e\n",
                mesh->nodes[i * 3 + 0],
                mesh->nodes[i * 3 + 1],
                mesh->nodes[i * 3 + 2]);
    }

    // Cells (tetrahedra)
    fprintf(fp, "\nCELLS %zu %zu\n", mesh->num_elements, mesh->num_elements * 5);
    for (size_t i = 0; i < mesh->num_elements; i++) {
        fprintf(fp, "4 %d %d %d %d\n",
                mesh->elements[i * 4 + 0],
                mesh->elements[i * 4 + 1],
                mesh->elements[i * 4 + 2],
                mesh->elements[i * 4 + 3]);
    }

    // Cell types (10 = tetrahedron)
    fprintf(fp, "\nCELL_TYPES %zu\n", mesh->num_elements);
    for (size_t i = 0; i < mesh->num_elements; i++) {
        fprintf(fp, "10\n");
    }

    // Point data (potential)
    fprintf(fp, "\nPOINT_DATA %zu\n", mesh->num_nodes);
    fprintf(fp, "SCALARS phi double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    PetscErrorCode ierr;
    const PetscScalar *sol_array;
    ierr = VecGetArrayRead(x, &sol_array); CHKERRQ(ierr);

    for (size_t i = 0; i < mesh->num_nodes; i++) {
        fprintf(fp, "%.12e\n", PetscRealPart(sol_array[i]));
    }

    ierr = VecRestoreArrayRead(x, &sol_array); CHKERRQ(ierr);

    fclose(fp);
    printf("Wrote VTK output to %s\n", filename);

    return 0;
}
