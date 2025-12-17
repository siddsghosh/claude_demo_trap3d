// Unit cube [0,1]^3 with physical groups for testing
// Surface 1: z=0 face
// Surface 2: z=1 face
// Volume 1: interior

SetFactory("OpenCASCADE");

// Create unit cube
Box(1) = {0, 0, 0, 1, 1, 1};

// Physical groups
Physical Volume("INTERIOR", 1) = {1};

// Define physical surfaces
// Gmsh box surfaces:
// Surface 1: z=0 (bottom)
// Surface 2: x=1 (right)
// Surface 3: x=0 (left)
// Surface 4: y=0 (front)
// Surface 5: y=1 (back)
// Surface 6: z=1 (top)

Physical Surface("Z_MIN", 1) = {1};  // z=0
Physical Surface("Z_MAX", 2) = {6};  // z=1

// Mesh settings - relatively coarse for testing
Mesh.CharacteristicLengthMin = 0.2;
Mesh.CharacteristicLengthMax = 0.3;
