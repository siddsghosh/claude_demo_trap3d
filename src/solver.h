#ifndef SOLVER_H
#define SOLVER_H

#include "mesh_io.h"
#include <petsc.h>

/**
 * Solve the linear system Kx = F using PETSc KSP
 * Returns solution vector x
 */
int solver_solve(Mat K, Vec F, Vec *x);

/**
 * Write solution to CSV file (x, y, z, phi)
 */
int solver_write_csv(const char *filename, const Mesh *mesh, Vec x);

/**
 * Write solution to VTK file for visualization
 */
int solver_write_vtk(const char *filename, const Mesh *mesh, Vec x);

#endif /* SOLVER_H */
