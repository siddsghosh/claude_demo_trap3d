#!/usr/bin/env python3
"""
Generate Gmsh .geo file from geometry.yaml specification
"""

import yaml
import sys

def load_geometry(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def generate_geo(geom, output_file):
    """Generate Gmsh .geo script with physical groups"""

    # Extract domain parameters
    domain = geom['domain']
    x_min, x_max = float(domain['x_min']), float(domain['x_max'])
    y_min, y_max = float(domain['y_min']), float(domain['y_max'])
    z_min, z_max = float(domain['z_min']), float(domain['z_max'])

    # Extract electrode parameters
    elec = geom['electrodes']
    rail_length = float(elec['rail_length'])
    rf_width = float(elec['rf_rail_width'])
    dc_width = float(elec['dc_rail_width'])
    gap = float(elec['gap'])
    endcap_length = float(elec['endcap_length'])
    endcap_hw = float(elec['endcap_half_width'])

    # Calculate positions
    # RF rail centered at x=0
    rf_x_min = -rf_width / 2
    rf_x_max = rf_width / 2

    # DC rails offset by (rf_width/2 + gap + dc_width/2)
    dc_offset = rf_width / 2 + gap + dc_width / 2
    dc_left_x_min = -dc_offset - dc_width / 2
    dc_left_x_max = -dc_offset + dc_width / 2
    dc_right_x_min = dc_offset - dc_width / 2
    dc_right_x_max = dc_offset + dc_width / 2

    # Rail y extents (centered)
    rail_y_min = -rail_length / 2
    rail_y_max = rail_length / 2

    # Endcap positions (beyond rails)
    endcap_y_center = rail_length / 2 + endcap_length / 2
    endcap_y_plus_min = rail_length / 2
    endcap_y_plus_max = rail_length / 2 + endcap_length
    endcap_y_minus_min = -rail_length / 2 - endcap_length
    endcap_y_minus_max = -rail_length / 2

    # Mesh parameters
    mesh = geom['mesh']
    lc_min = float(mesh['characteristic_length_min'])
    lc_max = float(mesh['characteristic_length_max'])

    with open(output_file, 'w') as f:
        f.write('// 3D Ion Trap Geometry - Generated from geometry.yaml\n')
        f.write('// Units: meters\n\n')

        f.write('SetFactory("OpenCASCADE");\n\n')

        # Create main domain box
        f.write('// Main computational domain\n')
        f.write(f'Box(1) = {{{x_min}, {y_min}, {z_min}, '
                f'{x_max - x_min}, {y_max - y_min}, {z_max - z_min}}};\n\n')

        # Create electrode rectangles on z=0 plane
        f.write('// Electrode patches on z=0 plane\n')

        # RF center rail
        f.write(f'// RF center rail\n')
        f.write(f'Rectangle(2) = {{{rf_x_min}, {rail_y_min}, {z_min}, '
                f'{rf_width}, {rail_length}}};\n\n')

        # DC left rail
        f.write(f'// DC left rail\n')
        f.write(f'Rectangle(3) = {{{dc_left_x_min}, {rail_y_min}, {z_min}, '
                f'{dc_width}, {rail_length}}};\n\n')

        # DC right rail
        f.write(f'// DC right rail\n')
        f.write(f'Rectangle(4) = {{{dc_right_x_min}, {rail_y_min}, {z_min}, '
                f'{dc_width}, {rail_length}}};\n\n')

        # Endcap +y
        f.write(f'// Endcap +y\n')
        f.write(f'Rectangle(5) = {{-{endcap_hw}, {endcap_y_plus_min}, {z_min}, '
                f'{2 * endcap_hw}, {endcap_length}}};\n\n')

        # Endcap -y
        f.write(f'// Endcap -y\n')
        f.write(f'Rectangle(6) = {{-{endcap_hw}, {endcap_y_minus_min}, {z_min}, '
                f'{2 * endcap_hw}, {endcap_length}}};\n\n')

        # Use BooleanFragments to imprint electrodes onto domain
        f.write('// Imprint electrodes onto domain bottom surface\n')
        f.write('BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }\n\n')

        # Set mesh sizes
        f.write('// Mesh sizing\n')
        f.write('Mesh.CharacteristicLengthMin = {};\n'.format(lc_min))
        f.write('Mesh.CharacteristicLengthMax = {};\n\n'.format(lc_max))

        # Physical groups - assign by geometric criteria
        f.write('// Physical groups\n')
        f.write('// Volume\n')
        f.write('Physical Volume("VACUUM", 1) = {1};\n\n')

        f.write('// Identify surfaces by geometric properties after BooleanFragments\n')
        f.write('// Tolerance for coordinate comparisons\n')
        f.write('tol = 1e-9;\n\n')

        # Write surface identification logic
        f.write('// Get all surfaces in the model\n')
        f.write('surfaces[] = Surface "*";\n\n')

        f.write('// Arrays to store surface IDs\n')
        f.write('rf_center[] = {};\n')
        f.write('dc_left[] = {};\n')
        f.write('dc_right[] = {};\n')
        f.write('endcap_plus_y[] = {};\n')
        f.write('endcap_minus_y[] = {};\n')
        f.write('ground_plane[] = {};\n')
        f.write('box_x_min[] = {};\n')
        f.write('box_x_max[] = {};\n')
        f.write('box_y_min[] = {};\n')
        f.write('box_y_max[] = {};\n')
        f.write('box_z_max[] = {};\n\n')

        f.write('// Classify surfaces by bounding box and center of mass\n')
        f.write('For i In {0:#surfaces[]-1}\n')
        f.write('  s = surfaces[i];\n')
        f.write('  bbox[] = BoundingBox Surface{s};\n')
        f.write('  // bbox: xmin, ymin, zmin, xmax, ymax, zmax\n')
        f.write('  xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2];\n')
        f.write('  xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5];\n')
        f.write('  xcen = (xmin + xmax) / 2;\n')
        f.write('  ycen = (ymin + ymax) / 2;\n')
        f.write('  zcen = (zmin + zmax) / 2;\n\n')

        f.write('  // Check if surface is on z=0 plane\n')
        f.write('  If (Fabs(zmin) < tol && Fabs(zmax) < tol)\n')

        # RF center rail: z=0, x in [-3e-5, 3e-5], y in [-2.5e-4, 2.5e-4]
        f.write(f'    // RF center rail\n')
        f.write(f'    If (xcen > {rf_x_min - 1e-9} && xcen < {rf_x_max + 1e-9} && ycen > {rail_y_min + 1e-9} && ycen < {rail_y_max - 1e-9})\n')
        f.write(f'      rf_center[] += s;\n')

        # DC left rail
        f.write(f'    // DC left rail\n')
        f.write(f'    ElseIf (xcen > {dc_left_x_min - 1e-9} && xcen < {dc_left_x_max + 1e-9} && ycen > {rail_y_min + 1e-9} && ycen < {rail_y_max - 1e-9})\n')
        f.write(f'      dc_left[] += s;\n')

        # DC right rail
        f.write(f'    // DC right rail\n')
        f.write(f'    ElseIf (xcen > {dc_right_x_min - 1e-9} && xcen < {dc_right_x_max + 1e-9} && ycen > {rail_y_min + 1e-9} && ycen < {rail_y_max - 1e-9})\n')
        f.write(f'      dc_right[] += s;\n')

        # Endcap +y
        f.write(f'    // Endcap +y\n')
        f.write(f'    ElseIf (ycen > {endcap_y_plus_min + 1e-9} && ycen < {endcap_y_plus_max - 1e-9})\n')
        f.write(f'      endcap_plus_y[] += s;\n')

        # Endcap -y
        f.write(f'    // Endcap -y\n')
        f.write(f'    ElseIf (ycen < {endcap_y_minus_max - 1e-9} && ycen > {endcap_y_minus_min + 1e-9})\n')
        f.write(f'      endcap_minus_y[] += s;\n')

        # Ground plane (remaining z=0 surfaces)
        f.write(f'    // Ground plane (rest of z=0)\n')
        f.write(f'    Else\n')
        f.write(f'      ground_plane[] += s;\n')
        f.write(f'    EndIf\n')

        # Box walls
        f.write(f'  // Check box walls\n')
        f.write(f'  ElseIf (Fabs(xmin - {x_min}) < tol && Fabs(xmax - {x_min}) < tol)\n')
        f.write(f'    box_x_min[] += s;\n')
        f.write(f'  ElseIf (Fabs(xmin - {x_max}) < tol && Fabs(xmax - {x_max}) < tol)\n')
        f.write(f'    box_x_max[] += s;\n')
        f.write(f'  ElseIf (Fabs(ymin - {y_min}) < tol && Fabs(ymax - {y_min}) < tol)\n')
        f.write(f'    box_y_min[] += s;\n')
        f.write(f'  ElseIf (Fabs(ymin - {y_max}) < tol && Fabs(ymax - {y_max}) < tol)\n')
        f.write(f'    box_y_max[] += s;\n')
        f.write(f'  ElseIf (Fabs(zmin - {z_max}) < tol && Fabs(zmax - {z_max}) < tol)\n')
        f.write(f'    box_z_max[] += s;\n')
        f.write(f'  EndIf\n')
        f.write('EndFor\n\n')

        # Assign physical groups
        f.write('// Assign physical groups\n')
        f.write('If (#rf_center[] > 0)\n')
        f.write('  Physical Surface("rf_center", 101) = {rf_center[]};\n')
        f.write('EndIf\n')
        f.write('If (#dc_left[] > 0)\n')
        f.write('  Physical Surface("dc_left", 102) = {dc_left[]};\n')
        f.write('EndIf\n')
        f.write('If (#dc_right[] > 0)\n')
        f.write('  Physical Surface("dc_right", 103) = {dc_right[]};\n')
        f.write('EndIf\n')
        f.write('If (#endcap_plus_y[] > 0)\n')
        f.write('  Physical Surface("endcap_plus_y", 104) = {endcap_plus_y[]};\n')
        f.write('EndIf\n')
        f.write('If (#endcap_minus_y[] > 0)\n')
        f.write('  Physical Surface("endcap_minus_y", 105) = {endcap_minus_y[]};\n')
        f.write('EndIf\n')
        f.write('If (#ground_plane[] > 0)\n')
        f.write('  Physical Surface("ground_plane", 106) = {ground_plane[]};\n')
        f.write('EndIf\n')
        f.write('If (#box_x_min[] > 0)\n')
        f.write('  Physical Surface("box_x_min", 107) = {box_x_min[]};\n')
        f.write('EndIf\n')
        f.write('If (#box_x_max[] > 0)\n')
        f.write('  Physical Surface("box_x_max", 108) = {box_x_max[]};\n')
        f.write('EndIf\n')
        f.write('If (#box_y_min[] > 0)\n')
        f.write('  Physical Surface("box_y_min", 109) = {box_y_min[]};\n')
        f.write('EndIf\n')
        f.write('If (#box_y_max[] > 0)\n')
        f.write('  Physical Surface("box_y_max", 110) = {box_y_max[]};\n')
        f.write('EndIf\n')
        f.write('If (#box_z_max[] > 0)\n')
        f.write('  Physical Surface("box_z_max", 111) = {box_z_max[]};\n')
        f.write('EndIf\n\n')

        f.write('// Generate 3D mesh\n')
        f.write('Mesh 3;\n')

def main():
    if len(sys.argv) < 2:
        yaml_file = 'specs/geometry.yaml'
    else:
        yaml_file = sys.argv[1]

    output_file = 'mesh/domain.geo'

    geom = load_geometry(yaml_file)
    generate_geo(geom, output_file)

    print(f'Generated {output_file} from {yaml_file}')
    print(f'Note: Physical surface IDs need manual assignment after BooleanFragments')
    print(f'Run: gmsh {output_file} -0 to see surface numbering in GUI')

if __name__ == '__main__':
    main()
