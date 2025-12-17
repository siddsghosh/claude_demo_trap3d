#include "boundary.h"
#include <stdio.h>

/**
 * Apply Dirichlet boundary condition using row elimination
 * For each constrained node i:
 * - Set row i to [0 ... 0 1 0 ... 0] (identity)
 * - Set F[i] = boundary_value
 */
int boundary_apply_dirichlet(Mat K, Vec F, const BoundaryCondition *bc) {
    PetscErrorCode ierr;

    printf("Applying BC on %s: %zu nodes, value = %.3f V\n",
           bc->surface_name, bc->num_nodes, bc->value);

    for (size_t i = 0; i < bc->num_nodes; i++) {
        PetscInt row = (PetscInt)bc->node_indices[i];

        // Zero out the row
        ierr = MatZeroRows(K, 1, &row, 1.0, NULL, NULL); CHKERRQ(ierr);

        // Set RHS value
        ierr = VecSetValue(F, row, bc->value, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

    return 0;
}

/**
 * Apply all Dirichlet boundary conditions
 */
int boundary_apply_dirichlet_all(Mat K, Vec F, const BoundaryCondition *bcs, size_t num_bcs) {
    // Collect all constrained rows
    PetscInt total_rows = 0;
    for (size_t i = 0; i < num_bcs; i++) {
        total_rows += (PetscInt)bcs[i].num_nodes;
    }

    PetscInt *rows = malloc(total_rows * sizeof(PetscInt));
    PetscScalar *values = malloc(total_rows * sizeof(PetscScalar));

    PetscInt idx = 0;
    for (size_t i = 0; i < num_bcs; i++) {
        printf("Applying BC on %s: %zu nodes, value = %.3f V\n",
               bcs[i].surface_name, bcs[i].num_nodes, bcs[i].value);

        for (size_t j = 0; j < bcs[i].num_nodes; j++) {
            rows[idx] = (PetscInt)bcs[i].node_indices[j];
            values[idx] = bcs[i].value;
            idx++;
        }
    }

    PetscErrorCode ierr;

    // Zero rows and set diagonal to 1.0
    ierr = MatZeroRows(K, total_rows, rows, 1.0, NULL, NULL); CHKERRQ(ierr);

    // Set RHS values
    ierr = VecSetValues(F, total_rows, rows, values, INSERT_VALUES); CHKERRQ(ierr);

    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

    free(rows);
    free(values);

    printf("Applied %zu boundary conditions (%d total DOFs constrained)\n",
           num_bcs, total_rows);

    return 0;
}
