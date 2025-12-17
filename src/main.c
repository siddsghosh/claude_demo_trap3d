#include "mesh_io.h"
#include "assembly.h"
#include "boundary.h"
#include "solver.h"
#include "yaml_reader.h"
#include <petsc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BCS 20

/**
 * Physical group IDs matching geometry.yaml
 */
#define PHYS_RF_CENTER      101
#define PHYS_DC_LEFT        102
#define PHYS_DC_RIGHT       103
#define PHYS_ENDCAP_PLUS_Y  104
#define PHYS_ENDCAP_MINUS_Y 105
#define PHYS_GROUND_PLANE   106
#define PHYS_BOX_X_MIN      107
#define PHYS_BOX_X_MAX      108
#define PHYS_BOX_Y_MIN      109
#define PHYS_BOX_Y_MAX      110
#define PHYS_BOX_Z_MAX      111

int main(int argc, char **argv) {
    PetscErrorCode ierr;

    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input.yaml>\n", argv[0]);
        PetscFinalize();
        return 1;
    }

    const char *input_file = argv[1];

    printf("======================================\n");
    printf("3D Poisson FEM Solver\n");
    printf("======================================\n\n");

    // Read input configuration
    char *mesh_file = NULL;
    char *output_prefix = NULL;
    double voltages[11]; // Up to 11 surfaces

    printf("Reading configuration from %s...\n", input_file);

    if (yaml_read_string(input_file, "mesh_file", &mesh_file) != 0 ||
        yaml_read_string(input_file, "output_prefix", &output_prefix) != 0) {
        fprintf(stderr, "Error: Failed to read input configuration\n");
        PetscFinalize();
        return 1;
    }

    // Read all voltage values
    if (yaml_read_double(input_file, "rf_center", &voltages[0]) != 0 ||
        yaml_read_double(input_file, "dc_left", &voltages[1]) != 0 ||
        yaml_read_double(input_file, "dc_right", &voltages[2]) != 0 ||
        yaml_read_double(input_file, "endcap_plus_y", &voltages[3]) != 0 ||
        yaml_read_double(input_file, "endcap_minus_y", &voltages[4]) != 0 ||
        yaml_read_double(input_file, "ground", &voltages[5]) != 0) {
        fprintf(stderr, "Error: Failed to read voltage values\n");
        free(mesh_file);
        free(output_prefix);
        PetscFinalize();
        return 1;
    }

    printf("  Mesh file: %s\n", mesh_file);
    printf("  Output prefix: %s\n", output_prefix);
    printf("  Voltages:\n");
    printf("    RF center: %.3f V\n", voltages[0]);
    printf("    DC left: %.3f V\n", voltages[1]);
    printf("    DC right: %.3f V\n", voltages[2]);
    printf("    Endcap +Y: %.3f V\n", voltages[3]);
    printf("    Endcap -Y: %.3f V\n", voltages[4]);
    printf("    Ground: %.3f V\n\n", voltages[5]);

    // Load mesh
    printf("Loading mesh...\n");
    Mesh mesh;
    if (mesh_read_gmsh(mesh_file, &mesh) != 0) {
        fprintf(stderr, "Error: Failed to load mesh\n");
        free(mesh_file);
        free(output_prefix);
        PetscFinalize();
        return 1;
    }
    printf("\n");

    // Assemble global system
    printf("Assembling FEM system...\n");
    Mat K;
    Vec F;
    double source_term = 0.0; // f = 0 for Laplace equation

    if (assembly_global(&mesh, source_term, &K, &F) != 0) {
        fprintf(stderr, "Error: Assembly failed\n");
        mesh_free(&mesh);
        free(mesh_file);
        free(output_prefix);
        PetscFinalize();
        return 1;
    }
    printf("\n");

    // Set up boundary conditions
    printf("Setting up boundary conditions...\n");
    BoundaryCondition bcs[MAX_BCS];
    size_t num_bcs = 0;

    // Map physical IDs to names and voltages
    typedef struct {
        int phys_id;
        const char *name;
        double voltage;
    } BCMapping;

    BCMapping bc_map[] = {
        {PHYS_RF_CENTER, "rf_center", voltages[0]},
        {PHYS_DC_LEFT, "dc_left", voltages[1]},
        {PHYS_DC_RIGHT, "dc_right", voltages[2]},
        {PHYS_ENDCAP_PLUS_Y, "endcap_plus_y", voltages[3]},
        {PHYS_ENDCAP_MINUS_Y, "endcap_minus_y", voltages[4]},
        {PHYS_GROUND_PLANE, "ground_plane", voltages[5]},
        {PHYS_BOX_X_MIN, "box_x_min", voltages[5]},
        {PHYS_BOX_X_MAX, "box_x_max", voltages[5]},
        {PHYS_BOX_Y_MIN, "box_y_min", voltages[5]},
        {PHYS_BOX_Y_MAX, "box_y_max", voltages[5]},
        {PHYS_BOX_Z_MAX, "box_z_max", voltages[5]}
    };

    size_t num_mappings = sizeof(bc_map) / sizeof(bc_map[0]);

    for (size_t i = 0; i < num_mappings; i++) {
        int *nodes = mesh_get_surface_nodes(&mesh, bc_map[i].phys_id,
                                            &bcs[num_bcs].num_nodes);
        if (nodes && bcs[num_bcs].num_nodes > 0) {
            bcs[num_bcs].surface_name = strdup(bc_map[i].name);
            bcs[num_bcs].physical_id = bc_map[i].phys_id;
            bcs[num_bcs].value = bc_map[i].voltage;
            bcs[num_bcs].node_indices = nodes;
            num_bcs++;
        }
    }

    if (num_bcs == 0) {
        fprintf(stderr, "Warning: No boundary conditions found\n");
    }

    // Apply boundary conditions
    if (boundary_apply_dirichlet_all(K, F, bcs, num_bcs) != 0) {
        fprintf(stderr, "Error: Failed to apply boundary conditions\n");
        mesh_free(&mesh);
        free(mesh_file);
        free(output_prefix);
        PetscFinalize();
        return 1;
    }
    printf("\n");

    // Solve system
    Vec x;
    if (solver_solve(K, F, &x) != 0) {
        fprintf(stderr, "Error: Solver failed\n");
        mesh_free(&mesh);
        free(mesh_file);
        free(output_prefix);
        PetscFinalize();
        return 1;
    }
    printf("\n");

    // Write output
    printf("Writing output files...\n");

    char csv_file[512];
    char vtk_file[512];
    snprintf(csv_file, sizeof(csv_file), "%s.csv", output_prefix);
    snprintf(vtk_file, sizeof(vtk_file), "%s.vtk", output_prefix);

    solver_write_csv(csv_file, &mesh, x);
    solver_write_vtk(vtk_file, &mesh, x);
    printf("\n");

    printf("======================================\n");
    printf("Simulation complete!\n");
    printf("======================================\n");

    // Cleanup
    for (size_t i = 0; i < num_bcs; i++) {
        free(bcs[i].surface_name);
        free(bcs[i].node_indices);
    }

    MatDestroy(&K);
    VecDestroy(&F);
    VecDestroy(&x);
    mesh_free(&mesh);
    free(mesh_file);
    free(output_prefix);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
