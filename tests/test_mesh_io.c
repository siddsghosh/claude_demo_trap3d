#include "../src/mesh_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) ASSERT(fabs((a)-(b)) < (tol), msg)

/**
 * Test reading single tetrahedron mesh
 */
int test_read_single_tet() {
    printf("Running test_read_single_tet...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/single_tet.msh", &mesh);
    ASSERT(ret == 0, "Failed to read mesh");

    ASSERT(mesh.num_nodes == 4, "Expected 4 nodes");
    ASSERT(mesh.num_elements == 1, "Expected 1 element");

    // Check node coordinates
    ASSERT_NEAR(mesh.nodes[0], 0.0, 1e-10, "Node 0 x");
    ASSERT_NEAR(mesh.nodes[1], 0.0, 1e-10, "Node 0 y");
    ASSERT_NEAR(mesh.nodes[2], 0.0, 1e-10, "Node 0 z");

    mesh_free(&mesh);
    printf("  PASS\n");
    return 0;
}

/**
 * Test reading unit cube mesh
 */
int test_read_unit_cube() {
    printf("Running test_read_unit_cube...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/unit_cube.msh", &mesh);
    ASSERT(ret == 0, "Failed to read mesh");

    ASSERT(mesh.num_nodes > 0, "Expected nodes");
    ASSERT(mesh.num_elements > 0, "Expected elements");

    // Check that all nodes are in [0,1]^3
    for (size_t i = 0; i < mesh.num_nodes; i++) {
        double x = mesh.nodes[i * 3 + 0];
        double y = mesh.nodes[i * 3 + 1];
        double z = mesh.nodes[i * 3 + 2];

        ASSERT(x >= -1e-10 && x <= 1.0 + 1e-10, "Node x out of bounds");
        ASSERT(y >= -1e-10 && y <= 1.0 + 1e-10, "Node y out of bounds");
        ASSERT(z >= -1e-10 && z <= 1.0 + 1e-10, "Node z out of bounds");
    }

    mesh_free(&mesh);
    printf("  PASS\n");
    return 0;
}

/**
 * Test physical groups
 */
int test_physical_groups() {
    printf("Running test_physical_groups...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/unit_cube.msh", &mesh);
    ASSERT(ret == 0, "Failed to read mesh");

    ASSERT(mesh.num_physical_surfaces >= 2, "Expected at least 2 physical surfaces");

    // Check that physical groups have names
    for (size_t i = 0; i < mesh.num_physical_surfaces; i++) {
        ASSERT(mesh.physical_surface_names[i] != NULL, "Physical surface name is NULL");
        ASSERT(mesh.physical_surface_ids[i] > 0, "Physical surface ID should be positive");
    }

    mesh_free(&mesh);
    printf("  PASS\n");
    return 0;
}

/**
 * Test node ordering in elements
 */
int test_node_ordering() {
    printf("Running test_node_ordering...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("meshes/single_tet.msh", &mesh);
    ASSERT(ret == 0, "Failed to read mesh");

    // Check that all element node indices are valid
    for (size_t i = 0; i < mesh.num_elements; i++) {
        for (int j = 0; j < 4; j++) {
            int node = mesh.elements[i * 4 + j];
            ASSERT(node >= 0 && node < (int)mesh.num_nodes, "Invalid node index");
        }
    }

    mesh_free(&mesh);
    printf("  PASS\n");
    return 0;
}

/**
 * Test invalid file handling
 */
int test_invalid_file() {
    printf("Running test_invalid_file...\n");

    Mesh mesh;
    int ret = mesh_read_gmsh("nonexistent.msh", &mesh);
    ASSERT(ret != 0, "Should fail on nonexistent file");

    printf("  PASS\n");
    return 0;
}

int main() {
    printf("=== Mesh I/O Tests ===\n\n");

    int failures = 0;

    failures += test_read_single_tet();
    failures += test_read_unit_cube();
    failures += test_physical_groups();
    failures += test_node_ordering();
    failures += test_invalid_file();

    printf("\n");
    if (failures == 0) {
        printf("All mesh_io tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
