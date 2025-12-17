#include "../src/mesh_io.h"
#include "../src/assembly.h"
#include "../src/boundary.h"
#include "../src/solver.h"
#include <petsc.h>
#include <stdio.h>
#include <math.h>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) ASSERT(fabs((a)-(b)) < (tol), msg)

/**
 * Test: Unit cube with φ=0 at z=0, φ=1 at z=1
 * Solution should be approximately φ=z
 */
int test_unit_cube_linear() {
    printf("Running test_unit_cube_linear...\n");

    // Load mesh
    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/unit_cube.msh", &mesh);
    ASSERT(ret == 0, "Failed to load mesh");

    // Assemble system
    Mat K;
    Vec F;
    ret = assembly_global(&mesh, 0.0, &K, &F);
    ASSERT(ret == 0, "Assembly failed");

    // Apply BCs: z=0 → φ=0, z=1 → φ=1
    BoundaryCondition bcs[2];

    // BC 1: Z_MIN (physical ID 1) → φ=0
    bcs[0].surface_name = "Z_MIN";
    bcs[0].physical_id = 1;
    bcs[0].value = 0.0;
    bcs[0].node_indices = mesh_get_surface_nodes(&mesh, 1, &bcs[0].num_nodes);
    ASSERT(bcs[0].node_indices != NULL, "Failed to get Z_MIN nodes");

    // BC 2: Z_MAX (physical ID 2) → φ=1
    bcs[1].surface_name = "Z_MAX";
    bcs[1].physical_id = 2;
    bcs[1].value = 1.0;
    bcs[1].node_indices = mesh_get_surface_nodes(&mesh, 2, &bcs[1].num_nodes);
    ASSERT(bcs[1].node_indices != NULL, "Failed to get Z_MAX nodes");

    ret = boundary_apply_dirichlet_all(K, F, bcs, 2);
    ASSERT(ret == 0, "Boundary condition application failed");

    // Solve
    Vec x;
    ret = solver_solve(K, F, &x);
    ASSERT(ret == 0, "Solver failed");

    // Check solution: φ should be approximately equal to z coordinate
    const PetscScalar *sol_array;
    VecGetArrayRead(x, &sol_array);

    double max_error = 0.0;
    for (size_t i = 0; i < mesh.num_nodes; i++) {
        double z = mesh.nodes[i * 3 + 2];
        double phi = PetscRealPart(sol_array[i]);
        double error = fabs(phi - z);
        if (error > max_error) max_error = error;
    }

    VecRestoreArrayRead(x, &sol_array);

    printf("  Max error: %.3e\n", max_error);
    ASSERT(max_error < 1e-6, "Solution error too large");

    // Cleanup
    free(bcs[0].node_indices);
    free(bcs[1].node_indices);
    MatDestroy(&K);
    VecDestroy(&F);
    VecDestroy(&x);
    mesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

/**
 * Test output file generation
 */
int test_output_csv() {
    printf("Running test_output_csv...\n");

    // Load mesh
    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/single_tet.msh", &mesh);
    ASSERT(ret == 0, "Failed to load mesh");

    // Create a dummy solution vector
    Vec x;
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, mesh.num_nodes);
    VecSetFromOptions(x);
    VecSet(x, 1.0);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    // Write CSV
    ret = solver_write_csv("test_output.csv", &mesh, x);
    ASSERT(ret == 0, "CSV write failed");

    // Check that file exists and has content
    FILE *fp = fopen("test_output.csv", "r");
    ASSERT(fp != NULL, "Output CSV file not created");

    char line[256];
    int line_count = 0;
    while (fgets(line, sizeof(line), fp)) {
        line_count++;
    }
    fclose(fp);

    ASSERT(line_count == (int)mesh.num_nodes + 1, "CSV should have header + data rows");

    // Cleanup
    VecDestroy(&x);
    mesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

/**
 * Test VTK output
 */
int test_output_vtk() {
    printf("Running test_output_vtk...\n");

    // Load mesh
    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/single_tet.msh", &mesh);
    ASSERT(ret == 0, "Failed to load mesh");

    // Create a dummy solution vector
    Vec x;
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, mesh.num_nodes);
    VecSetFromOptions(x);
    VecSet(x, 1.0);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    // Write VTK
    ret = solver_write_vtk("test_output.vtk", &mesh, x);
    ASSERT(ret == 0, "VTK write failed");

    // Check that file exists
    FILE *fp = fopen("test_output.vtk", "r");
    ASSERT(fp != NULL, "Output VTK file not created");
    fclose(fp);

    // Cleanup
    VecDestroy(&x);
    mesh_free(&mesh);

    printf("  PASS\n");
    return 0;
}

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    printf("=== Integration Tests ===\n\n");

    int failures = 0;

    failures += test_unit_cube_linear();
    failures += test_output_csv();
    failures += test_output_vtk();

    printf("\n");
    if (failures == 0) {
        printf("All integration tests passed!\n");
    } else {
        printf("%d test(s) failed\n", failures);
    }

    PetscFinalize();
    return failures;
}
