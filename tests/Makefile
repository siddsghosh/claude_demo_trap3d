# Makefile for unit tests

# PETSc configuration
PETSC_DIR = /usr/lib/petscdir/petsc3.19
PETSC_ARCH = x86_64-linux-gnu-real
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables

# Compiler flags
CFLAGS += -Wall -Wextra -O2 -g -I../src

# Source files from src/
SRC_DIR = ../src
OBJS = $(SRC_DIR)/mesh_io.o $(SRC_DIR)/assembly.o $(SRC_DIR)/boundary.o $(SRC_DIR)/solver.o $(SRC_DIR)/yaml_reader.o

# Test executables
TESTS = test_mesh_io test_assembly test_boundary test_solver test_integration

# Default target: build and run all tests
all: $(TESTS)
	@echo "Running tests..."
	@./test_mesh_io
	@./test_assembly
	@./test_boundary
	@./test_solver
	@./test_integration
	@echo ""
	@echo "All tests completed!"

# Build individual tests
test_mesh_io: test_mesh_io.o $(SRC_DIR)/mesh_io.o
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

test_assembly: test_assembly.o $(SRC_DIR)/assembly.o $(SRC_DIR)/mesh_io.o
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

test_boundary: test_boundary.o $(SRC_DIR)/boundary.o $(SRC_DIR)/mesh_io.o
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

test_solver: test_solver.o $(SRC_DIR)/solver.o $(SRC_DIR)/mesh_io.o
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

test_integration: test_integration.o $(OBJS)
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

# Compile test files
%.o: %.c
	$(COMPILE.c) $(CFLAGS) $(PETSC_CC_INCLUDES) -o $@ $<

# Ensure src objects are built
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(MAKE) -C $(SRC_DIR) $(@F)

# Clean
clean:
	rm -f $(TESTS) *.o test_output.csv test_output.vtk

# Build only (no run)
build: $(TESTS)

.PHONY: all clean build
