#include "../src/mesh_io.h"
#include "../src/boundary.h"
#include <petsc.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)

/**
 * Test surface node identification
 */
int test_identify_boundary_nodes() {
    printf("Running test_identify_boundary_nodes...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/unit_cube.msh", &mesh);
    ASSERT(ret == 0, "Failed to read mesh");

    // Get nodes on surface 1 (Z_MIN)
    size_t num_nodes;
    int *nodes = mesh_get_surface_nodes(&mesh, 1, &num_nodes);

    ASSERT(nodes != NULL, "Should find nodes on surface 1");
    ASSERT(num_nodes > 0, "Surface should have nodes");

    // Check that all nodes have z ≈ 0
    for (size_t i = 0; i < num_nodes; i++) {
        int node_idx = nodes[i];
        double z = mesh.nodes[node_idx * 3 + 2];
        ASSERT(fabs(z) < 1e-6, "Node on z=0 surface should have z≈0");
    }

    free(nodes);
    mesh_free(&mesh);
    printf("  PASS\n");
    return 0;
}

/**
 * Test basic boundary condition application
 */
int test_dirichlet_single() {
    printf("Running test_dirichlet_single...\n");

    // Create a simple 3x3 matrix
    Mat K;
    Vec F;
    MatCreate(PETSC_COMM_WORLD, &K);
    MatSetSizes(K, PETSC_DECIDE, PETSC_DECIDE, 3, 3);
    MatSetFromOptions(K);
    MatSetUp(K);

    VecCreate(PETSC_COMM_WORLD, &F);
    VecSetSizes(F, PETSC_DECIDE, 3);
    VecSetFromOptions(F);

    // Set some values
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double val = (i == j) ? 2.0 : -1.0;
            MatSetValue(K, i, j, val, INSERT_VALUES);
        }
        VecSetValue(F, i, 1.0, INSERT_VALUES);
    }

    MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(F);
    VecAssemblyEnd(F);

    // Apply BC on node 0
    BoundaryCondition bc;
    bc.surface_name = "test";
    bc.physical_id = 1;
    bc.value = 5.0;
    bc.num_nodes = 1;
    int node = 0;
    bc.node_indices = &node;

    boundary_apply_dirichlet(K, F, &bc);

    // Check that row 0 has diagonal = 1 and F[0] = 5.0
    double diag;
    MatGetValue(K, 0, 0, &diag);
    ASSERT(fabs(diag - 1.0) < 1e-10, "Diagonal should be 1.0 after BC");

    PetscScalar val;
    VecGetValues(F, 1, &node, &val);
    ASSERT(fabs(PetscRealPart(val) - 5.0) < 1e-10, "RHS should be 5.0 after BC");

    MatDestroy(&K);
    VecDestroy(&F);

    printf("  PASS\n");
    return 0;
}

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    printf("=== Boundary Condition Tests ===\n\n");

    int failures = 0;

    failures += test_identify_boundary_nodes();
    failures += test_dirichlet_single();

    printf("\n");
    if (failures == 0) {
        printf("All boundary tests passed!\n");
    } else {
        printf("%d test(s) failed\n", failures);
    }

    PetscFinalize();
    return failures;
}
