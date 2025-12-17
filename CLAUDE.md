# 3D Poisson FEM Solver

## Goal

Generate a C solver for the 3D Poisson equation using the Finite Element Method on tetrahedral meshes.

**Dependencies:** Gmsh (meshing), PETSc (linear algebra)

## Directory Structure

```
poisson_fem/
├── CLAUDE.md
├── specs/
│   ├── problem.yaml
│   ├── geometry.yaml
│   ├── solver.yaml
│   └── input.yaml
├── mesh/
│   ├── domain.geo
│   └── domain.msh
├── src/
│   ├── main.c
│   ├── mesh_io.h / mesh_io.c
│   ├── assembly.h / assembly.c
│   ├── boundary.h / boundary.c
│   ├── solver.h / solver.c
│   ├── yaml_reader.h / yaml_reader.c
│   └── Makefile
├── tests/
│   ├── test_mesh_io.c
│   ├── test_assembly.c
│   ├── test_boundary.c
│   ├── test_solver.c
│   ├── test_integration.c
│   ├── meshes/
│   │   ├── unit_cube.msh
│   │   └── single_tet.msh
│   └── Makefile
├── postprocessing/
│   └── extract_slice.py
└── output/
```

## Poisson Equation

**Strong form:**
```
-∇²φ = f    in Ω ⊂ ℝ³
```

**Weak form:**
```
∫_Ω ∇φ·∇v dΩ = ∫_Ω fv dΩ    for all v ∈ V
```

**Boundary conditions:**
- Dirichlet: φ = g (enforced by row elimination)
- Neumann: ∂φ/∂n = h (boundary integral)

## FEM Discretization

**P1 tetrahedra** with shape functions:
```
N₁ = 1 - ξ - η - ζ
N₂ = ξ
N₃ = η
N₄ = ζ
```

**Local stiffness (4×4):** Kᵉᵢⱼ = (∇Nᵢ · ∇Nⱼ) × Volume

**Local load:** Fᵉᵢ = f(centroid) × Volume / 4

**Dirichlet enforcement:** Set row to identity, RHS to boundary value

## Code Requirements

- C99 with PETSc
- Read `petsc_dir` and `gmsh_path` from specs/solver.yaml
- Read mesh file, output prefix, and voltages from specs/input.yaml
- Read Gmsh MSH 4.x format (3D tetrahedral mesh)
- Assemble sparse matrix via MatSetValues
- Solve with KSP (CG + GAMG recommended for 3D)
- Output: CSV (x,y,z,phi) and VTK

**Makefile:**
```makefile
include $(PETSC_DIR)/lib/petsc/conf/variables
include $(PETSC_DIR)/lib/petsc/conf/rules

CFLAGS += -Wall -Wextra -O2
solver: main.o mesh_io.o assembly.o boundary.o solver.o yaml_reader.o
	$(LINK.c) -o $@ $^ $(PETSC_LIB)
```

## Geometry Generation

Generate Gmsh .geo files from specs/geometry.yaml:

**Structure of geometry.yaml:**
- `factory`: OpenCASCADE or Built-in
- `units`: meters, mm, etc.
- `domain`: Main computational volume (box, cylinder, etc.)
- `electrodes`: Surface patches with parameters (rails, endcaps, pads)
- `mesh`: Characteristic length settings
- `physical_groups`: Named volumes and surfaces for BCs

**Generation rules:**
1. Create main domain volume (Box for ion trap)
2. Create electrode rectangles on z=0 plane
3. Use BooleanFragments to imprint patches onto domain surfaces
4. Assign physical groups for volume and each surface (electrodes, walls, ground)
5. Set mesh size parameters

**Meshing command:**
```
gmsh -3 mesh/domain.geo -o mesh/domain.msh
```

## Boundary Condition Mapping

The problem.yaml maps boundary names to keys in input.yaml:

| Surface | value_key | Description |
|---------|-----------|-------------|
| rf_center | rf_center | RF rail voltage |
| dc_left | dc_left | Left DC rail |
| dc_right | dc_right | Right DC rail |
| endcap_plus_y | endcap_plus_y | +Y endcap |
| endcap_minus_y | endcap_minus_y | -Y endcap |
| ground_plane, box_* | ground | Ground (0V) |

## Runtime Configuration (input.yaml)

All runtime parameters in one file:

```yaml
mesh_file: mesh/domain.msh
output_prefix: output/solution
rf_center: 100.0
dc_left: 5.0
dc_right: 5.0
endcap_plus_y: 10.0
endcap_minus_y: 10.0
ground: 0.0
```

**Usage:**
```bash
./solver specs/input.yaml
```

**Outputs generated:**
- `<output_prefix>.csv` — Node coordinates and potential (x,y,z,phi)
- `<output_prefix>.vtk` — VTK file for ParaView

## Postprocessing

- **VTK output** for ParaView (3D visualization)
- **extract_slice.py**: Extract 2D slice at specified z and plot contour
- **CSV output**: x,y,z,phi for each node

## Unit Tests

### Test Structure

Each module has a corresponding test file in `tests/`. Tests use simple assert macros and return 0 on success, non-zero on failure.

```c
#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)
#define ASSERT_NEAR(a, b, tol, msg) ASSERT(fabs((a)-(b)) < (tol), msg)
```

### test_mesh_io.c

| Test | Description |
|------|-------------|
| `test_read_single_tet` | Read mesh with 1 tetrahedron, verify 4 nodes, 1 element |
| `test_read_unit_cube` | Read unit cube mesh, verify node coordinates in [0,1]³ |
| `test_physical_groups` | Verify physical group IDs parsed correctly |
| `test_node_ordering` | Verify tetrahedral node ordering is consistent |
| `test_invalid_file` | Return error for non-existent file |

### test_assembly.c

| Test | Description |
|------|-------------|
| `test_shape_gradients` | Verify ∇N₁ + ∇N₂ + ∇N₃ + ∇N₄ = 0 |
| `test_element_volume` | Compute volume of known tetrahedron, compare to analytical |
| `test_stiffness_symmetry` | Verify Kᵉᵢⱼ = Kᵉⱼᵢ |
| `test_stiffness_rowsum` | Verify row sums of stiffness matrix = 0 (Laplacian property) |
| `test_reference_tet` | Compare Kᵉ for unit tet against hand-calculated values |

### test_boundary.c

| Test | Description |
|------|-------------|
| `test_dirichlet_single` | Apply Dirichlet to one node, verify row modified |
| `test_dirichlet_preserves_symmetry` | Verify K remains symmetric after BC application |
| `test_identify_boundary_nodes` | Correctly identify nodes on each physical surface |
| `test_value_lookup` | Verify voltage lookup from input.yaml keys |

### test_solver.c

| Test | Description |
|------|-------------|
| `test_petsc_init` | PETSc initializes without error |
| `test_matrix_assembly` | Assemble small system, verify non-zero pattern |
| `test_solve_identity` | Solve Ix = b, verify x = b |
| `test_solve_laplacian_1d` | 1D Laplacian with Dirichlet BCs, compare to analytical |

### test_integration.c

| Test | Description |
|------|-------------|
| `test_unit_cube_linear` | Cube with φ=0 at z=0, φ=1 at z=1; solution should be φ=z |
| `test_unit_cube_convergence` | Refine mesh, verify error decreases as O(h²) |
| `test_output_csv` | Verify CSV output is readable and contains correct columns |
| `test_output_vtk` | Verify VTK output is valid format |

### Test Meshes (tests/meshes/)

**single_tet.msh:** One tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)

**unit_cube.msh:** Unit cube [0,1]³ with physical groups:
- Surface 1: z=0 face
- Surface 2: z=1 face
- Volume 1: interior

### Running Tests

```makefile
# tests/Makefile
include $(PETSC_DIR)/lib/petsc/conf/variables

test: test_mesh_io test_assembly test_boundary test_solver test_integration
	./test_mesh_io
	./test_assembly
	./test_boundary
	./test_solver
	./test_integration
	@echo "All tests passed"

test_%: test_%.c ../src/*.o
	$(LINK.c) -o $@ $< ../src/*.o $(PETSC_LIB)
```

**Run all tests:**
```bash
cd tests && make test
```

### Validation Criteria

| Test Type | Pass Criteria |
|-----------|--------------|
| Unit tests | All assertions pass |
| Integration (linear solution) | L∞ error < 1e-10 |
| Integration (convergence) | Order ≥ 1.9 (expect 2.0 for P1) |
| Output format | Files readable by pandas/ParaView |
