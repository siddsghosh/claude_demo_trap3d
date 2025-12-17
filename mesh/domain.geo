// 3D Ion Trap Geometry - Generated from geometry.yaml
// Units: meters

SetFactory("OpenCASCADE");

// Main computational domain
Box(1) = {-0.0003, -0.0004, 0.0, 0.0006, 0.0008, 0.0003};

// Electrode patches on z=0 plane
// RF center rail
Rectangle(2) = {-3e-05, -0.00025, 0.0, 6e-05, 0.0005};

// DC left rail
Rectangle(3) = {-0.00011, -0.00025, 0.0, 6e-05, 0.0005};

// DC right rail
Rectangle(4) = {5.000000000000001e-05, -0.00025, 0.0, 6e-05, 0.0005};

// Endcap +y
Rectangle(5) = {-0.00022, 0.00025, 0.0, 0.00044, 0.00012};

// Endcap -y
Rectangle(6) = {-0.00022, -0.00037, 0.0, 0.00044, 0.00012};

// Imprint electrodes onto domain bottom surface
BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }

// Mesh sizing
Mesh.CharacteristicLengthMin = 8e-06;
Mesh.CharacteristicLengthMax = 3e-05;

// Physical groups
// Volume
Physical Volume("VACUUM", 1) = {1};

// Identify surfaces by geometric properties after BooleanFragments
// Tolerance for coordinate comparisons
tol = 1e-9;

// Get all surfaces in the model
surfaces[] = Surface "*";

// Arrays to store surface IDs
rf_center[] = {};
dc_left[] = {};
dc_right[] = {};
endcap_plus_y[] = {};
endcap_minus_y[] = {};
ground_plane[] = {};
box_x_min[] = {};
box_x_max[] = {};
box_y_min[] = {};
box_y_max[] = {};
box_z_max[] = {};

// Classify surfaces by bounding box and center of mass
For i In {0:#surfaces[]-1}
  s = surfaces[i];
  bbox[] = BoundingBox Surface{s};
  // bbox: xmin, ymin, zmin, xmax, ymax, zmax
  xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2];
  xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5];
  xcen = (xmin + xmax) / 2;
  ycen = (ymin + ymax) / 2;
  zcen = (zmin + zmax) / 2;

  // Check if surface is on z=0 plane
  If (Fabs(zmin) < tol && Fabs(zmax) < tol)
    // RF center rail
    If (xcen > -3.0001e-05 && xcen < 3.0001e-05 && ycen > -0.000249999 && ycen < 0.000249999)
      rf_center[] += s;
    // DC left rail
    ElseIf (xcen > -0.000110001 && xcen < -4.999900000000001e-05 && ycen > -0.000249999 && ycen < 0.000249999)
      dc_left[] += s;
    // DC right rail
    ElseIf (xcen > 4.999900000000001e-05 && xcen < 0.000110001 && ycen > -0.000249999 && ycen < 0.000249999)
      dc_right[] += s;
    // Endcap +y
    ElseIf (ycen > 0.000250001 && ycen < 0.000369999)
      endcap_plus_y[] += s;
    // Endcap -y
    ElseIf (ycen < -0.000250001 && ycen > -0.000369999)
      endcap_minus_y[] += s;
    // Ground plane (rest of z=0)
    Else
      ground_plane[] += s;
    EndIf
  // Check box walls
  ElseIf (Fabs(xmin - -0.0003) < tol && Fabs(xmax - -0.0003) < tol)
    box_x_min[] += s;
  ElseIf (Fabs(xmin - 0.0003) < tol && Fabs(xmax - 0.0003) < tol)
    box_x_max[] += s;
  ElseIf (Fabs(ymin - -0.0004) < tol && Fabs(ymax - -0.0004) < tol)
    box_y_min[] += s;
  ElseIf (Fabs(ymin - 0.0004) < tol && Fabs(ymax - 0.0004) < tol)
    box_y_max[] += s;
  ElseIf (Fabs(zmin - 0.0003) < tol && Fabs(zmax - 0.0003) < tol)
    box_z_max[] += s;
  EndIf
EndFor

// Assign physical groups
If (#rf_center[] > 0)
  Physical Surface("rf_center", 101) = {rf_center[]};
EndIf
If (#dc_left[] > 0)
  Physical Surface("dc_left", 102) = {dc_left[]};
EndIf
If (#dc_right[] > 0)
  Physical Surface("dc_right", 103) = {dc_right[]};
EndIf
If (#endcap_plus_y[] > 0)
  Physical Surface("endcap_plus_y", 104) = {endcap_plus_y[]};
EndIf
If (#endcap_minus_y[] > 0)
  Physical Surface("endcap_minus_y", 105) = {endcap_minus_y[]};
EndIf
If (#ground_plane[] > 0)
  Physical Surface("ground_plane", 106) = {ground_plane[]};
EndIf
If (#box_x_min[] > 0)
  Physical Surface("box_x_min", 107) = {box_x_min[]};
EndIf
If (#box_x_max[] > 0)
  Physical Surface("box_x_max", 108) = {box_x_max[]};
EndIf
If (#box_y_min[] > 0)
  Physical Surface("box_y_min", 109) = {box_y_min[]};
EndIf
If (#box_y_max[] > 0)
  Physical Surface("box_y_max", 110) = {box_y_max[]};
EndIf
If (#box_z_max[] > 0)
  Physical Surface("box_z_max", 111) = {box_z_max[]};
EndIf

// Generate 3D mesh
// Mesh 3;  // Commented out - mesh from command line instead
