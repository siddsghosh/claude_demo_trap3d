#include "../src/solver.h"
#include <petsc.h>
#include <stdio.h>
#include <math.h>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) ASSERT(fabs((a)-(b)) < (tol), msg)

/**
 * Test PETSc initialization
 */
int test_petsc_init() {
    printf("Running test_petsc_init...\n");
    // If we got here, PETSc is initialized
    printf("  PASS\n");
    return 0;
}

/**
 * Test solving identity system: Ix = b
 */
int test_solve_identity() {
    printf("Running test_solve_identity...\n");

    int n = 5;
    Mat K;
    Vec F, x;

    // Create identity matrix
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(K);
    MatSetUp(K);

    for (int i = 0; i < n; i++) {
        MatSetValue(K, i, i, 1.0, INSERT_VALUES);
    }

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

    // Create RHS
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, PETSC_DECIDE, n);
    VecSetFromOptions(F);

    for (int i = 0; i < n; i++) {
        VecSetValue(F, i, (double)(i + 1), INSERT_VALUES);
    }

    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    // Solve
    solver_solve(K, F, &x);

    // Check solution: x should equal F
    for (int i = 0; i < n; i++) {
        PetscScalar val;
        VecGetValues(x, 1, &i, &val);
        ASSERT_NEAR(PetscRealPart(val), (double)(i + 1), 1e-10, "Solution of Ix=b should be x=b");
    }

    MatDestroy(&K);
    VecDestroy(&F);
    VecDestroy(&x);

    printf("  PASS\n");
    return 0;
}

/**
 * Test simple 3-node 1D Laplacian problem
 * -d²u/dx² = 0, u(0) = 0, u(1) = 1
 * Solution: u = x
 */
int test_solve_laplacian_1d() {
    printf("Running test_solve_laplacian_1d...\n");

    int n = 3;
    Mat K;
    Vec F, x;

    // Create 1D Laplacian matrix (tridiagonal: 2, -1, -1)
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(K);
    MatSetUp(K);

    double h = 0.5; // dx
    double k = 1.0 / h; // stiffness

    MatSetValue(K, 0, 0, k, INSERT_VALUES);
    MatSetValue(K, 0, 1, -k, INSERT_VALUES);
    MatSetValue(K, 1, 0, -k, INSERT_VALUES);
    MatSetValue(K, 1, 1, 2*k, INSERT_VALUES);
    MatSetValue(K, 1, 2, -k, INSERT_VALUES);
    MatSetValue(K, 2, 1, -k, INSERT_VALUES);
    MatSetValue(K, 2, 2, k, INSERT_VALUES);

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

    // Create RHS (all zeros for Laplace)
    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, PETSC_DECIDE, n);
    VecSetFromOptions(F);
    VecSet(F, 0.0);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    // Apply BCs: u(0) = 0, u(2) = 1
    PetscInt rows[2] = {0, 2};
    MatZeroRows(K, 2, rows, 1.0, NULL, NULL);
    VecSetValue(F, 0, 0.0, INSERT_VALUES);
    VecSetValue(F, 2, 1.0, INSERT_VALUES);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    // Solve
    solver_solve(K, F, &x);

    // Check solution: should be linear u = x/1.0 at x=0, 0.5, 1.0
    double expected[3] = {0.0, 0.5, 1.0};
    for (int i = 0; i < n; i++) {
        PetscScalar val;
        VecGetValues(x, 1, &i, &val);
        ASSERT_NEAR(PetscRealPart(val), expected[i], 1e-8, "1D Laplacian solution incorrect");
    }

    MatDestroy(&K);
    VecDestroy(&F);
    VecDestroy(&x);

    printf("  PASS\n");
    return 0;
}

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    printf("=== Solver Tests ===\n\n");

    int failures = 0;

    failures += test_petsc_init();
    failures += test_solve_identity();
    failures += test_solve_laplacian_1d();

    printf("\n");
    if (failures == 0) {
        printf("All solver tests passed!\n");
    } else {
        printf("%d test(s) failed\n", failures);
    }

    PetscFinalize();
    return failures;
}
