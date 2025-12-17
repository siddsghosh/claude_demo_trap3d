#ifndef MESH_IO_H
#define MESH_IO_H

#include <stddef.h>

/**
 * Mesh data structure for tetrahedral meshes
 */
typedef struct {
    // Node data
    size_t num_nodes;
    double *nodes;      // [num_nodes * 3] - x, y, z coordinates

    // Element data (tetrahedra only)
    size_t num_elements;
    int *elements;      // [num_elements * 4] - node indices (0-based)

    // Physical groups for boundary surfaces
    size_t num_physical_surfaces;
    int *physical_surface_ids;          // [num_physical_surfaces]
    char **physical_surface_names;      // [num_physical_surfaces]

    // Surface triangle data
    size_t *num_surface_triangles;      // [num_physical_surfaces]
    int **surface_triangles;            // [num_physical_surfaces][num_triangles * 3]
} Mesh;

/**
 * Read a Gmsh MSH 4.x format mesh file
 * Returns 0 on success, non-zero on failure
 */
int mesh_read_gmsh(const char *filename, Mesh *mesh);

/**
 * Free mesh memory
 */
void mesh_free(Mesh *mesh);

/**
 * Get nodes on a given physical surface
 * Returns array of node indices (caller must free)
 */
int* mesh_get_surface_nodes(const Mesh *mesh, int physical_id, size_t *num_nodes);

#endif /* MESH_IO_H */
