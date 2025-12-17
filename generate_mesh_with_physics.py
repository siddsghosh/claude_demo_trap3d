#!/usr/bin/env python3
"""
Generate mesh with physical groups using Gmsh Python API
This ensures proper physical group assignment after BooleanFragments
"""

import gmsh
import sys

def main():
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    # Set OpenCASCADE factory
    gmsh.model.add("ion_trap")
    factory = gmsh.model.occ

    # Domain parameters (meters)
    x_min, x_max = -300e-6, 300e-6
    y_min, y_max = -400e-6, 400e-6
    z_min, z_max = 0.0, 300e-6

    # Create main box
    box = factory.addBox(x_min, y_min, z_min,
                         x_max - x_min, y_max - y_min, z_max - z_min)

    # Electrode parameters
    rf_width = 60e-6
    dc_width = 60e-6
    gap = 20e-6
    rail_length = 500e-6
    endcap_length = 120e-6
    endcap_hw = 220e-6

    # Calculate positions
    rf_x_min = -rf_width / 2
    rf_x_max = rf_width / 2

    dc_offset = rf_width / 2 + gap + dc_width / 2
    dc_left_x_min = -dc_offset - dc_width / 2
    dc_left_x_max = -dc_offset + dc_width / 2
    dc_right_x_min = dc_offset - dc_width / 2
    dc_right_x_max = dc_offset + dc_width / 2

    rail_y_min = -rail_length / 2
    rail_y_max = rail_length / 2

    endcap_y_plus_min = rail_length / 2
    endcap_y_plus_max = rail_length / 2 + endcap_length
    endcap_y_minus_min = -rail_length / 2 - endcap_length
    endcap_y_minus_max = -rail_length / 2

    # Create electrode rectangles on z=0
    rf_center = factory.addRectangle(rf_x_min, rail_y_min, z_min, rf_width, rail_length)
    dc_left = factory.addRectangle(dc_left_x_min, rail_y_min, z_min, dc_width, rail_length)
    dc_right = factory.addRectangle(dc_right_x_min, rail_y_min, z_min, dc_width, rail_length)
    endcap_plus_y = factory.addRectangle(-endcap_hw, endcap_y_plus_min, z_min, 2*endcap_hw, endcap_length)
    endcap_minus_y = factory.addRectangle(-endcap_hw, endcap_y_minus_min, z_min, 2*endcap_hw, endcap_length)

    # Synchronize before Boolean operations
    factory.synchronize()

    # Perform BooleanFragments to imprint electrodes
    out_dimtags, out_dimtags_map = factory.fragment(
        [(3, box)],  # Object (3D volume)
        [(2, rf_center), (2, dc_left), (2, dc_right), (2, endcap_plus_y), (2, endcap_minus_y)]  # Tools (2D surfaces)
    )

    # Synchronize after Boolean operations
    factory.synchronize()

    # Get all surfaces
    all_surfaces = gmsh.model.getEntities(dim=2)

    print(f"Found {len(all_surfaces)} surfaces after BooleanFragments")

    # Classify surfaces by their bounding box
    tol = 1e-6  # 1 micron tolerance

    rf_center_surfs = []
    dc_left_surfs = []
    dc_right_surfs = []
    endcap_plus_y_surfs = []
    endcap_minus_y_surfs = []
    ground_plane_surfs = []
    box_x_min_surfs = []
    box_x_max_surfs = []
    box_y_min_surfs = []
    box_y_max_surfs = []
    box_z_max_surfs = []

    for dim, tag in all_surfaces:
        bbox = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        xcen = (xmin + xmax) / 2
        ycen = (ymin + ymax) / 2
        zcen = (zmin + zmax) / 2

        print(f"Surface {tag}: bbox=({xmin:.3e}, {ymin:.3e}, {zmin:.3e}) to ({xmax:.3e}, {ymax:.3e}, {zmax:.3e}), center=({xcen:.3e}, {ycen:.3e}, {zcen:.3e})")

        # Check if on z=0 plane
        if abs(zmin) < tol and abs(zmax) < tol:
            # Calculate surface dimensions
            x_span = xmax - xmin
            y_span = ymax - ymin

            # RF center rail (should be approximately rf_width x rail_length)
            if (rf_x_min - tol < xcen < rf_x_max + tol and
                rail_y_min + tol < ycen < rail_y_max - tol and
                x_span < 2 * rf_width and y_span < 2 * rail_length):
                rf_center_surfs.append(tag)
                print(f"  Surface {tag}: rf_center")
            # DC left rail
            elif (dc_left_x_min - tol < xcen < dc_left_x_max + tol and
                  rail_y_min + tol < ycen < rail_y_max - tol and
                  x_span < 2 * dc_width and y_span < 2 * rail_length):
                dc_left_surfs.append(tag)
                print(f"  Surface {tag}: dc_left")
            # DC right rail
            elif (dc_right_x_min - tol < xcen < dc_right_x_max + tol and
                  rail_y_min + tol < ycen < rail_y_max - tol and
                  x_span < 2 * dc_width and y_span < 2 * rail_length):
                dc_right_surfs.append(tag)
                print(f"  Surface {tag}: dc_right")
            # Endcap +y
            elif (endcap_y_plus_min + tol < ycen < endcap_y_plus_max - tol and
                  x_span < 4 * endcap_hw and y_span < 2 * endcap_length):
                endcap_plus_y_surfs.append(tag)
                print(f"  Surface {tag}: endcap_plus_y")
            # Endcap -y
            elif (endcap_y_minus_min + tol < ycen < endcap_y_minus_max - tol and
                  x_span < 4 * endcap_hw and y_span < 2 * endcap_length):
                endcap_minus_y_surfs.append(tag)
                print(f"  Surface {tag}: endcap_minus_y")
            # Ground plane (rest of z=0)
            else:
                ground_plane_surfs.append(tag)
                print(f"  Surface {tag}: ground_plane")
        # Box walls
        elif abs(xmin - x_min) < tol and abs(xmax - x_min) < tol:
            box_x_min_surfs.append(tag)
            print(f"  Surface {tag}: box_x_min")
        elif abs(xmin - x_max) < tol and abs(xmax - x_max) < tol:
            box_x_max_surfs.append(tag)
            print(f"  Surface {tag}: box_x_max")
        elif abs(ymin - y_min) < tol and abs(ymax - y_min) < tol:
            box_y_min_surfs.append(tag)
            print(f"  Surface {tag}: box_y_min")
        elif abs(ymin - y_max) < tol and abs(ymax - y_max) < tol:
            box_y_max_surfs.append(tag)
            print(f"  Surface {tag}: box_y_max")
        elif abs(zmin - z_max) < tol and abs(zmax - z_max) < tol:
            box_z_max_surfs.append(tag)
            print(f"  Surface {tag}: box_z_max")

    # Create physical groups
    if rf_center_surfs:
        gmsh.model.addPhysicalGroup(2, rf_center_surfs, 101)
        gmsh.model.setPhysicalName(2, 101, "rf_center")

    if dc_left_surfs:
        gmsh.model.addPhysicalGroup(2, dc_left_surfs, 102)
        gmsh.model.setPhysicalName(2, 102, "dc_left")

    if dc_right_surfs:
        gmsh.model.addPhysicalGroup(2, dc_right_surfs, 103)
        gmsh.model.setPhysicalName(2, 103, "dc_right")

    if endcap_plus_y_surfs:
        gmsh.model.addPhysicalGroup(2, endcap_plus_y_surfs, 104)
        gmsh.model.setPhysicalName(2, 104, "endcap_plus_y")

    if endcap_minus_y_surfs:
        gmsh.model.addPhysicalGroup(2, endcap_minus_y_surfs, 105)
        gmsh.model.setPhysicalName(2, 105, "endcap_minus_y")

    if ground_plane_surfs:
        gmsh.model.addPhysicalGroup(2, ground_plane_surfs, 106)
        gmsh.model.setPhysicalName(2, 106, "ground_plane")

    if box_x_min_surfs:
        gmsh.model.addPhysicalGroup(2, box_x_min_surfs, 107)
        gmsh.model.setPhysicalName(2, 107, "box_x_min")

    if box_x_max_surfs:
        gmsh.model.addPhysicalGroup(2, box_x_max_surfs, 108)
        gmsh.model.setPhysicalName(2, 108, "box_x_max")

    if box_y_min_surfs:
        gmsh.model.addPhysicalGroup(2, box_y_min_surfs, 109)
        gmsh.model.setPhysicalName(2, 109, "box_y_min")

    if box_y_max_surfs:
        gmsh.model.addPhysicalGroup(2, box_y_max_surfs, 110)
        gmsh.model.setPhysicalName(2, 110, "box_y_max")

    if box_z_max_surfs:
        gmsh.model.addPhysicalGroup(2, box_z_max_surfs, 111)
        gmsh.model.setPhysicalName(2, 111, "box_z_max")

    # Add physical volume
    volumes = gmsh.model.getEntities(dim=3)
    if volumes:
        vol_tags = [tag for dim, tag in volumes]
        gmsh.model.addPhysicalGroup(3, vol_tags, 1)
        gmsh.model.setPhysicalName(3, 1, "VACUUM")

    # Set mesh sizes
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 8e-6)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 30e-6)

    # Generate 3D mesh
    print("\nGenerating mesh...")
    gmsh.model.mesh.generate(3)

    # Write mesh
    output_file = "mesh/domain.msh"
    gmsh.write(output_file)
    print(f"\nMesh written to {output_file}")

    # Get mesh statistics
    nodes = gmsh.model.mesh.getNodes()
    elements = gmsh.model.mesh.getElements(dim=3)
    print(f"Mesh statistics:")
    print(f"  Nodes: {len(nodes[0])}")
    print(f"  Tetrahedra: {len(elements[2][0])}")

    # Cleanup
    gmsh.finalize()

if __name__ == "__main__":
    main()
