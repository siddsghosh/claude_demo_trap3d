#include "../src/assembly.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) ASSERT(fabs((a)-(b)) < (tol), msg)

/**
 * Test that shape function gradients sum to zero
 * ∇N₁ + ∇N₂ + ∇N₃ + ∇N₄ = 0
 */
int test_shape_gradients() {
    printf("Running test_shape_gradients...\n");

    // Reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    double nodes[12] = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };

    double grads[4][3];
    double vol = assembly_compute_gradients(nodes, grads);

    // Sum gradients
    for (int j = 0; j < 3; j++) {
        double sum = 0.0;
        for (int i = 0; i < 4; i++) {
            sum += grads[i][j];
        }
        ASSERT_NEAR(sum, 0.0, 1e-12, "Sum of shape gradients should be zero");
    }

    printf("  PASS\n");
    return 0;
}

/**
 * Test element volume calculation
 */
int test_element_volume() {
    printf("Running test_element_volume...\n");

    // Reference tetrahedron: volume = 1/6
    double nodes[12] = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };

    double grads[4][3];
    double vol = assembly_compute_gradients(nodes, grads);

    double expected_vol = 1.0 / 6.0;
    ASSERT_NEAR(vol, expected_vol, 1e-12, "Volume of reference tet should be 1/6");

    printf("  PASS\n");
    return 0;
}

/**
 * Test stiffness matrix symmetry
 * K_ij should equal K_ji
 */
int test_stiffness_symmetry() {
    printf("Running test_stiffness_symmetry...\n");

    double nodes[12] = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };

    double K[4][4];
    assembly_element_stiffness(nodes, K);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ASSERT_NEAR(K[i][j], K[j][i], 1e-12, "Stiffness matrix should be symmetric");
        }
    }

    printf("  PASS\n");
    return 0;
}

/**
 * Test that stiffness matrix row sums are zero
 * This is a property of the Laplacian
 */
int test_stiffness_rowsum() {
    printf("Running test_stiffness_rowsum...\n");

    double nodes[12] = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };

    double K[4][4];
    assembly_element_stiffness(nodes, K);

    for (int i = 0; i < 4; i++) {
        double sum = 0.0;
        for (int j = 0; j < 4; j++) {
            sum += K[i][j];
        }
        ASSERT_NEAR(sum, 0.0, 1e-10, "Row sum of stiffness matrix should be zero");
    }

    printf("  PASS\n");
    return 0;
}

/**
 * Test reference tetrahedron against hand-calculated values
 */
int test_reference_tet() {
    printf("Running test_reference_tet...\n");

    // Reference tetrahedron
    double nodes[12] = {
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    };

    double grads[4][3];
    double vol = assembly_compute_gradients(nodes, grads);

    // For reference tet, gradients are:
    // ∇N₁ = (-1, -1, -1)
    // ∇N₂ = ( 1,  0,  0)
    // ∇N₃ = ( 0,  1,  0)
    // ∇N₄ = ( 0,  0,  1)

    ASSERT_NEAR(grads[0][0], -1.0, 1e-10, "grad N1 x");
    ASSERT_NEAR(grads[0][1], -1.0, 1e-10, "grad N1 y");
    ASSERT_NEAR(grads[0][2], -1.0, 1e-10, "grad N1 z");

    ASSERT_NEAR(grads[1][0],  1.0, 1e-10, "grad N2 x");
    ASSERT_NEAR(grads[1][1],  0.0, 1e-10, "grad N2 y");
    ASSERT_NEAR(grads[1][2],  0.0, 1e-10, "grad N2 z");

    ASSERT_NEAR(grads[2][0],  0.0, 1e-10, "grad N3 x");
    ASSERT_NEAR(grads[2][1],  1.0, 1e-10, "grad N3 y");
    ASSERT_NEAR(grads[2][2],  0.0, 1e-10, "grad N3 z");

    ASSERT_NEAR(grads[3][0],  0.0, 1e-10, "grad N4 x");
    ASSERT_NEAR(grads[3][1],  0.0, 1e-10, "grad N4 y");
    ASSERT_NEAR(grads[3][2],  1.0, 1e-10, "grad N4 z");

    // Compute stiffness
    double K[4][4];
    assembly_element_stiffness(nodes, K);

    // K11 = (∇N₁ · ∇N₁) * V = (1+1+1) * 1/6 = 0.5
    ASSERT_NEAR(K[0][0], 0.5, 1e-10, "K[0][0]");

    // K12 = (∇N₁ · ∇N₂) * V = (-1) * 1/6 = -1/6
    ASSERT_NEAR(K[0][1], -1.0/6.0, 1e-10, "K[0][1]");

    printf("  PASS\n");
    return 0;
}

int main() {
    printf("=== Assembly Tests ===\n\n");

    int failures = 0;

    failures += test_shape_gradients();
    failures += test_element_volume();
    failures += test_stiffness_symmetry();
    failures += test_stiffness_rowsum();
    failures += test_reference_tet();

    printf("\n");
    if (failures == 0) {
        printf("All assembly tests passed!\n");
        return 0;
    } else {
        printf("%d test(s) failed\n", failures);
        return 1;
    }
}
