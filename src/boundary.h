#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "mesh_io.h"
#include <petsc.h>

/**
 * Boundary condition data
 */
typedef struct {
    char *surface_name;
    int physical_id;
    double value;
    size_t num_nodes;
    int *node_indices;
} BoundaryCondition;

/**
 * Apply Dirichlet boundary conditions using row elimination method
 * Sets row i to identity, RHS to boundary value
 */
int boundary_apply_dirichlet(Mat K, Vec F, const BoundaryCondition *bc);

/**
 * Apply multiple Dirichlet boundary conditions
 */
int boundary_apply_dirichlet_all(Mat K, Vec F, const BoundaryCondition *bcs, size_t num_bcs);

#endif /* BOUNDARY_H */
