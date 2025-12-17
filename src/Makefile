# Makefile for 3D Poisson FEM Solver

# PETSc configuration
PETSC_DIR = /usr/lib/petscdir/petsc3.19
PETSC_ARCH = x86_64-linux-gnu-real
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables
include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules

# Compiler flags
CFLAGS += -Wall -Wextra -O2 -g

# Target executable
TARGET = solver

# Source files
SRCS = main.c mesh_io.c assembly.c boundary.c solver.c yaml_reader.c
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Link target
$(TARGET): $(OBJS)
	$(LINK.c) -o $@ $^ $(PETSC_LIB)

# Compile C files
%.o: %.c
	$(COMPILE.c) $(CFLAGS) $(PETSC_CC_INCLUDES) -o $@ $<

# Dependencies
main.o: main.c mesh_io.h assembly.h boundary.h solver.h yaml_reader.h
mesh_io.o: mesh_io.c mesh_io.h
assembly.o: assembly.c assembly.h mesh_io.h
boundary.o: boundary.c boundary.h mesh_io.h
solver.o: solver.c solver.h mesh_io.h
yaml_reader.o: yaml_reader.c yaml_reader.h

# Clean
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
