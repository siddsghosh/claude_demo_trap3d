#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "mesh_io.h"
#include <petsc.h>

/**
 * Compute shape function gradients for a tetrahedron
 * grads[4][3] = gradient of shape function i in direction j
 * Returns element volume
 */
double assembly_compute_gradients(const double *nodes, double grads[4][3]);

/**
 * Compute element stiffness matrix (4x4) for a tetrahedron
 * K_local[i][j] = integral of (grad N_i · grad N_j) dV
 */
void assembly_element_stiffness(const double *nodes, double K_local[4][4]);

/**
 * Compute element load vector (4x1) for constant source term
 * F_local[i] = f * V / 4, where V is element volume
 */
void assembly_element_load(const double *nodes, double f, double F_local[4]);

/**
 * Assemble global stiffness matrix and load vector
 * Uses PETSc Mat and Vec
 */
int assembly_global(const Mesh *mesh, double source_term, Mat *K, Vec *F);

#endif /* ASSEMBLY_H */
