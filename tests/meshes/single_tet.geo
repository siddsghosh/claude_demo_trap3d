// Single tetrahedron test mesh
// Vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)

SetFactory("Built-in");

// Define points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {0, 0, 1, 1.0};

// Define lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

// Define surfaces
Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};  // Base triangle

Line Loop(2) = {1, 5, -4};
Plane Surface(2) = {2};

Line Loop(3) = {2, 6, -5};
Plane Surface(3) = {3};

Line Loop(4) = {3, 4, -6};
Plane Surface(4) = {4};

// Define volume
Surface Loop(1) = {1, 2, 3, 4};
Volume(1) = {1};

// Physical groups
Physical Volume("VOLUME", 1) = {1};
Physical Surface("SURFACE", 1) = {1, 2, 3, 4};

// Generate mesh with no subdivision
Mesh.CharacteristicLengthMin = 2.0;
Mesh.CharacteristicLengthMax = 2.0;
