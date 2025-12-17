#include "mesh_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define GMSH_TET 4      // Tetrahedron element type
#define GMSH_TRI 2      // Triangle element type

/**
 * Read Gmsh MSH 4.x format mesh file
 */
int mesh_read_gmsh(const char *filename, Mesh *mesh) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open mesh file %s\n", filename);
        return 1;
    }

    // Initialize mesh
    memset(mesh, 0, sizeof(Mesh));

    char line[MAX_LINE_LENGTH];
    int in_section = 0;

    // Temporary storage
    size_t num_physical_names = 0;
    int *temp_phys_ids = NULL;
    char **temp_phys_names = NULL;

    // Entity to physical group mapping for surfaces
    int *entity_to_physical = calloc(10000, sizeof(int)); // Assume max 10000 entities

    while (fgets(line, sizeof(line), fp)) {
        // Check for section markers
        if (strncmp(line, "$MeshFormat", 11) == 0) {
            fgets(line, sizeof(line), fp);
            double version;
            int file_type, data_size;
            sscanf(line, "%lf %d %d", &version, &file_type, &data_size);
            if (version < 4.0) {
                fprintf(stderr, "Error: MSH format version %.1f not supported (need 4.x)\n", version);
                fclose(fp);
                return 1;
            }
            // Read $EndMeshFormat
            fgets(line, sizeof(line), fp);
        }
        else if (strncmp(line, "$PhysicalNames", 14) == 0) {
            fgets(line, sizeof(line), fp);
            sscanf(line, "%zu", &num_physical_names);

            temp_phys_ids = malloc(num_physical_names * sizeof(int));
            temp_phys_names = malloc(num_physical_names * sizeof(char*));

            for (size_t i = 0; i < num_physical_names; i++) {
                fgets(line, sizeof(line), fp);
                int dim, tag;
                char name[256];
                sscanf(line, "%d %d %s", &dim, &tag, name);

                // Remove quotes from name
                size_t len = strlen(name);
                if (len >= 2 && name[0] == '"' && name[len-1] == '"') {
                    name[len-1] = '\0';
                    memmove(name, name + 1, len - 1);
                }

                // Store only surface physical groups (dim = 2)
                if (dim == 2) {
                    temp_phys_ids[mesh->num_physical_surfaces] = tag;
                    temp_phys_names[mesh->num_physical_surfaces] = strdup(name);
                    mesh->num_physical_surfaces++;
                }
            }
            // Read $EndPhysicalNames
            fgets(line, sizeof(line), fp);

            // Allocate final arrays
            mesh->physical_surface_ids = malloc(mesh->num_physical_surfaces * sizeof(int));
            mesh->physical_surface_names = malloc(mesh->num_physical_surfaces * sizeof(char*));
            mesh->num_surface_triangles = calloc(mesh->num_physical_surfaces, sizeof(size_t));
            mesh->surface_triangles = calloc(mesh->num_physical_surfaces, sizeof(int*));

            memcpy(mesh->physical_surface_ids, temp_phys_ids,
                   mesh->num_physical_surfaces * sizeof(int));
            memcpy(mesh->physical_surface_names, temp_phys_names,
                   mesh->num_physical_surfaces * sizeof(char*));

            free(temp_phys_ids);
            free(temp_phys_names);
        }
        else if (strncmp(line, "$Entities", 9) == 0) {
            fgets(line, sizeof(line), fp);
            size_t num_points, num_curves, num_surfaces, num_volumes;
            sscanf(line, "%zu %zu %zu %zu", &num_points, &num_curves, &num_surfaces, &num_volumes);

            // Skip points
            for (size_t i = 0; i < num_points; i++) {
                fgets(line, sizeof(line), fp);
            }

            // Skip curves
            for (size_t i = 0; i < num_curves; i++) {
                fgets(line, sizeof(line), fp);
            }

            // Read surfaces and extract physical group tags
            for (size_t i = 0; i < num_surfaces; i++) {
                fgets(line, sizeof(line), fp);
                int entity_tag;
                double minx, miny, minz, maxx, maxy, maxz;
                size_t num_physical_tags;

                // Parse: tag minX minY minZ maxX maxY maxZ numPhysicalTags physicalTag(s) ...
                char *ptr = line;
                sscanf(ptr, "%d %lf %lf %lf %lf %lf %lf %zu",
                       &entity_tag, &minx, &miny, &minz, &maxx, &maxy, &maxz, &num_physical_tags);

                // Find the position after the bounding box values
                for (int j = 0; j < 8; j++) {
                    ptr = strchr(ptr, ' ');
                    if (ptr) ptr++;
                }

                // Read physical tag if present
                if (num_physical_tags > 0 && ptr) {
                    int physical_tag;
                    sscanf(ptr, "%d", &physical_tag);
                    entity_to_physical[entity_tag] = physical_tag;
                }
            }

            // Skip volumes
            for (size_t i = 0; i < num_volumes; i++) {
                fgets(line, sizeof(line), fp);
            }

            // Read $EndEntities
            fgets(line, sizeof(line), fp);
        }
        else if (strncmp(line, "$Nodes", 6) == 0) {
            fgets(line, sizeof(line), fp);
            size_t num_entity_blocks, num_nodes_total, min_tag, max_tag;
            sscanf(line, "%zu %zu %zu %zu", &num_entity_blocks, &num_nodes_total,
                   &min_tag, &max_tag);

            mesh->num_nodes = num_nodes_total;
            mesh->nodes = malloc(num_nodes_total * 3 * sizeof(double));

            size_t node_count = 0;

            for (size_t block = 0; block < num_entity_blocks; block++) {
                fgets(line, sizeof(line), fp);
                int entity_dim, entity_tag, parametric;
                size_t num_nodes_in_block;
                sscanf(line, "%d %d %d %zu", &entity_dim, &entity_tag,
                       &parametric, &num_nodes_in_block);

                // Read node tags
                size_t *node_tags = malloc(num_nodes_in_block * sizeof(size_t));
                for (size_t i = 0; i < num_nodes_in_block; i++) {
                    fgets(line, sizeof(line), fp);
                    sscanf(line, "%zu", &node_tags[i]);
                }

                // Read node coordinates
                for (size_t i = 0; i < num_nodes_in_block; i++) {
                    fgets(line, sizeof(line), fp);
                    size_t idx = node_tags[i] - 1; // Convert to 0-based
                    sscanf(line, "%lf %lf %lf",
                           &mesh->nodes[idx * 3 + 0],
                           &mesh->nodes[idx * 3 + 1],
                           &mesh->nodes[idx * 3 + 2]);
                    node_count++;
                }

                free(node_tags);
            }

            // Read $EndNodes
            fgets(line, sizeof(line), fp);
        }
        else if (strncmp(line, "$Elements", 9) == 0) {
            fgets(line, sizeof(line), fp);
            size_t num_entity_blocks, num_elements_total, min_tag, max_tag;
            sscanf(line, "%zu %zu %zu %zu", &num_entity_blocks, &num_elements_total,
                   &min_tag, &max_tag);

            // First pass: count tetrahedra and triangles per surface
            size_t num_tets = 0;
            size_t **tri_counts = calloc(mesh->num_physical_surfaces, sizeof(size_t*));

            long start_pos = ftell(fp);

            for (size_t block = 0; block < num_entity_blocks; block++) {
                fgets(line, sizeof(line), fp);
                int entity_dim, entity_tag, element_type;
                size_t num_elements_in_block;
                sscanf(line, "%d %d %d %zu", &entity_dim, &entity_tag,
                       &element_type, &num_elements_in_block);

                if (element_type == GMSH_TET) {
                    num_tets += num_elements_in_block;
                }
                else if (element_type == GMSH_TRI && entity_dim == 2) {
                    // Map entity tag to physical group
                    int physical_tag = entity_to_physical[entity_tag];
                    // Find which physical surface this belongs to
                    for (size_t i = 0; i < mesh->num_physical_surfaces; i++) {
                        if (mesh->physical_surface_ids[i] == physical_tag) {
                            mesh->num_surface_triangles[i] += num_elements_in_block;
                            break;
                        }
                    }
                }

                // Skip element lines
                for (size_t i = 0; i < num_elements_in_block; i++) {
                    fgets(line, sizeof(line), fp);
                }
            }

            // Allocate arrays
            mesh->num_elements = num_tets;
            mesh->elements = malloc(num_tets * 4 * sizeof(int));

            for (size_t i = 0; i < mesh->num_physical_surfaces; i++) {
                if (mesh->num_surface_triangles[i] > 0) {
                    mesh->surface_triangles[i] = malloc(mesh->num_surface_triangles[i] * 3 * sizeof(int));
                }
            }

            // Second pass: read elements
            fseek(fp, start_pos, SEEK_SET);

            size_t tet_idx = 0;
            size_t *tri_idx = calloc(mesh->num_physical_surfaces, sizeof(size_t));

            for (size_t block = 0; block < num_entity_blocks; block++) {
                fgets(line, sizeof(line), fp);
                int entity_dim, entity_tag, element_type;
                size_t num_elements_in_block;
                sscanf(line, "%d %d %d %zu", &entity_dim, &entity_tag,
                       &element_type, &num_elements_in_block);

                for (size_t i = 0; i < num_elements_in_block; i++) {
                    fgets(line, sizeof(line), fp);

                    if (element_type == GMSH_TET) {
                        size_t elem_tag;
                        int n0, n1, n2, n3;
                        sscanf(line, "%zu %d %d %d %d", &elem_tag, &n0, &n1, &n2, &n3);

                        // Store as 0-based indices
                        mesh->elements[tet_idx * 4 + 0] = n0 - 1;
                        mesh->elements[tet_idx * 4 + 1] = n1 - 1;
                        mesh->elements[tet_idx * 4 + 2] = n2 - 1;
                        mesh->elements[tet_idx * 4 + 3] = n3 - 1;
                        tet_idx++;
                    }
                    else if (element_type == GMSH_TRI && entity_dim == 2) {
                        // Map entity tag to physical group
                        int physical_tag = entity_to_physical[entity_tag];
                        // Find which physical surface this belongs to
                        for (size_t j = 0; j < mesh->num_physical_surfaces; j++) {
                            if (mesh->physical_surface_ids[j] == physical_tag) {
                                size_t elem_tag;
                                int n0, n1, n2;
                                sscanf(line, "%zu %d %d %d", &elem_tag, &n0, &n1, &n2);

                                size_t idx = tri_idx[j];
                                mesh->surface_triangles[j][idx * 3 + 0] = n0 - 1;
                                mesh->surface_triangles[j][idx * 3 + 1] = n1 - 1;
                                mesh->surface_triangles[j][idx * 3 + 2] = n2 - 1;
                                tri_idx[j]++;
                                break;
                            }
                        }
                    }
                }
            }

            free(tri_idx);

            // Read $EndElements
            fgets(line, sizeof(line), fp);
        }
    }

    fclose(fp);
    free(entity_to_physical);

    printf("Mesh loaded: %zu nodes, %zu tetrahedra, %zu physical surfaces\n",
           mesh->num_nodes, mesh->num_elements, mesh->num_physical_surfaces);

    return 0;
}

/**
 * Free mesh memory
 */
void mesh_free(Mesh *mesh) {
    if (mesh->nodes) free(mesh->nodes);
    if (mesh->elements) free(mesh->elements);

    if (mesh->physical_surface_ids) free(mesh->physical_surface_ids);

    if (mesh->physical_surface_names) {
        for (size_t i = 0; i < mesh->num_physical_surfaces; i++) {
            if (mesh->physical_surface_names[i]) {
                free(mesh->physical_surface_names[i]);
            }
        }
        free(mesh->physical_surface_names);
    }

    if (mesh->num_surface_triangles) free(mesh->num_surface_triangles);

    if (mesh->surface_triangles) {
        for (size_t i = 0; i < mesh->num_physical_surfaces; i++) {
            if (mesh->surface_triangles[i]) {
                free(mesh->surface_triangles[i]);
            }
        }
        free(mesh->surface_triangles);
    }

    memset(mesh, 0, sizeof(Mesh));
}

/**
 * Get nodes on a given physical surface
 */
int* mesh_get_surface_nodes(const Mesh *mesh, int physical_id, size_t *num_nodes) {
    // Find physical surface index
    size_t surf_idx = 0;
    int found = 0;
    for (size_t i = 0; i < mesh->num_physical_surfaces; i++) {
        if (mesh->physical_surface_ids[i] == physical_id) {
            surf_idx = i;
            found = 1;
            break;
        }
    }

    if (!found) {
        *num_nodes = 0;
        return NULL;
    }

    // Collect unique node indices from triangles
    size_t num_tris = mesh->num_surface_triangles[surf_idx];
    int *triangles = mesh->surface_triangles[surf_idx];

    // Use a simple marking array
    int *marked = calloc(mesh->num_nodes, sizeof(int));
    size_t count = 0;

    for (size_t i = 0; i < num_tris * 3; i++) {
        int node = triangles[i];
        if (!marked[node]) {
            marked[node] = 1;
            count++;
        }
    }

    // Allocate result array
    int *result = malloc(count * sizeof(int));
    size_t idx = 0;
    for (size_t i = 0; i < mesh->num_nodes; i++) {
        if (marked[i]) {
            result[idx++] = i;
        }
    }

    free(marked);

    *num_nodes = count;
    return result;
}
