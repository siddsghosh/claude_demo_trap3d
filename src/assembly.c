#include "assembly.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * Compute shape function gradients for a tetrahedron
 * P1 tetrahedral elements:
 * N1 = 1 - xi - eta - zeta
 * N2 = xi
 * N3 = eta
 * N4 = zeta
 *
 * Returns element volume
 */
double assembly_compute_gradients(const double *nodes, double grads[4][3]) {
    // nodes: [x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3]

    // Compute Jacobian matrix
    // J = [x1-x0  x2-x0  x3-x0]
    //     [y1-y0  y2-y0  y3-y0]
    //     [z1-z0  z2-z0  z3-z0]

    double J[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            J[i][j] = nodes[(j+1)*3 + i] - nodes[i];
        }
    }

    // Compute determinant
    double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
               - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
               + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

    double volume = fabs(det) / 6.0;

    if (fabs(det) < 1e-15) {
        fprintf(stderr, "Warning: Degenerate element detected (det = %.3e)\n", det);
    }

    // Compute inverse Jacobian
    double inv_det = 1.0 / det;
    double J_inv[3][3];

    J_inv[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv_det;
    J_inv[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * inv_det;
    J_inv[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv_det;

    J_inv[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * inv_det;
    J_inv[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv_det;
    J_inv[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * inv_det;

    J_inv[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv_det;
    J_inv[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * inv_det;
    J_inv[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv_det;

    // Shape function gradients in reference element:
    // dN1/dxi = -1, dN1/deta = -1, dN1/dzeta = -1
    // dN2/dxi =  1, dN2/deta =  0, dN2/dzeta =  0
    // dN3/dxi =  0, dN3/deta =  1, dN3/dzeta =  0
    // dN4/dxi =  0, dN4/deta =  0, dN4/dzeta =  1

    double grad_ref[4][3] = {
        {-1.0, -1.0, -1.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0}
    };

    // Transform to physical element: grad_N = J^{-T} * grad_ref
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            grads[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                grads[i][j] += J_inv[k][j] * grad_ref[i][k];
            }
        }
    }

    return volume;
}

/**
 * Compute element stiffness matrix
 * K_ij = integral (grad N_i · grad N_j) dV = (grad N_i · grad N_j) * V
 */
void assembly_element_stiffness(const double *nodes, double K_local[4][4]) {
    double grads[4][3];
    double volume = assembly_compute_gradients(nodes, grads);

    // K_ij = (grad N_i · grad N_j) * volume
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double dot = 0.0;
            for (int k = 0; k < 3; k++) {
                dot += grads[i][k] * grads[j][k];
            }
            K_local[i][j] = dot * volume;
        }
    }
}

/**
 * Compute element load vector
 * For constant source f: F_i = f * V / 4
 */
void assembly_element_load(const double *nodes, double f, double F_local[4]) {
    double grads[4][3];
    double volume = assembly_compute_gradients(nodes, grads);

    // For constant f, integral (f * N_i) dV = f * V / 4
    double val = f * volume / 4.0;
    for (int i = 0; i < 4; i++) {
        F_local[i] = val;
    }
}

/**
 * Assemble global stiffness matrix and load vector
 */
int assembly_global(const Mesh *mesh, double source_term, Mat *K, Vec *F) {
    PetscErrorCode ierr;
    PetscInt n = (PetscInt)mesh->num_nodes;

    // Create matrix and vector
    ierr = MatCreate(PETSC_COMM_WORLD, K); CHKERRQ(ierr);
    ierr = MatSetSizes(*K, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
    ierr = MatSetFromOptions(*K); CHKERRQ(ierr);

    // Preallocate: tetrahedral mesh, typical connectivity ~15-30 per row
    ierr = MatSeqAIJSetPreallocation(*K, 30, NULL); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(*K, 30, NULL, 30, NULL); CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, F); CHKERRQ(ierr);
    ierr = VecSetSizes(*F, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetFromOptions(*F); CHKERRQ(ierr);

    // Assemble element by element
    for (size_t e = 0; e < mesh->num_elements; e++) {
        // Get element nodes
        int elem_nodes[4];
        for (int i = 0; i < 4; i++) {
            elem_nodes[i] = mesh->elements[e * 4 + i];
        }

        // Get node coordinates
        double coords[12];
        for (int i = 0; i < 4; i++) {
            int node_idx = elem_nodes[i];
            coords[i * 3 + 0] = mesh->nodes[node_idx * 3 + 0];
            coords[i * 3 + 1] = mesh->nodes[node_idx * 3 + 1];
            coords[i * 3 + 2] = mesh->nodes[node_idx * 3 + 2];
        }

        // Compute element matrices
        double K_local[4][4];
        double F_local[4];

        assembly_element_stiffness(coords, K_local);
        assembly_element_load(coords, source_term, F_local);

        // Add to global system
        PetscInt indices[4];
        for (int i = 0; i < 4; i++) {
            indices[i] = (PetscInt)elem_nodes[i];
        }

        ierr = MatSetValues(*K, 4, indices, 4, indices, (PetscScalar*)K_local,
                           ADD_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(*F, 4, indices, F_local, ADD_VALUES); CHKERRQ(ierr);
    }

    // Assemble
    ierr = MatAssemblyBegin(*K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*K, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = VecAssemblyBegin(*F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*F); CHKERRQ(ierr);

    printf("Global assembly complete: %d DOFs\n", n);

    return 0;
}
