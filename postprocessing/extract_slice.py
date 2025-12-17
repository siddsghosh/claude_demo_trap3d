#!/usr/bin/env python3
"""
Extract 2D slice from 3D FEM solution VTK file and plot contours

Usage:
    python extract_slice.py <solution.vtk> <z_value> [--output output.png]
    python extract_slice.py output/solution.vtk 50e-6 -o slice_z50um.png
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.interpolate import griddata

def load_vtk_solution(filename):
    """
    Load FEM solution from VTK file
    Returns pyvista UnstructuredGrid object
    """
    print(f"Loading VTK file: {filename}")
    mesh = pv.read(filename)
    print(f"  Mesh type: {type(mesh)}")
    print(f"  Number of points: {mesh.n_points}")
    print(f"  Number of cells: {mesh.n_cells}")

    # List available arrays
    print(f"  Point data arrays: {mesh.point_data.keys()}")

    return mesh

def extract_slice_from_vtk(mesh, z_value, tolerance=1e-6):
    """
    Extract horizontal slice at z = z_value from VTK mesh
    Returns arrays of x, y, phi values at the slice
    """
    print(f"\nExtracting slice at z = {z_value:.3e} m (tolerance = {tolerance:.3e} m)...")

    # Get point coordinates
    points = mesh.points
    z_coords = points[:, 2]

    # Find points within tolerance of z_value
    mask = np.abs(z_coords - z_value) < tolerance

    if not np.any(mask):
        # Try to find closest z value
        z_min, z_max = z_coords.min(), z_coords.max()
        print(f"Warning: No points found at z = {z_value:.3e}")
        print(f"Z range in mesh: [{z_min:.3e}, {z_max:.3e}]")

        # Find closest z value
        closest_z = z_coords[np.argmin(np.abs(z_coords - z_value))]
        print(f"Using closest z = {closest_z:.3e} instead")
        mask = np.abs(z_coords - closest_z) < tolerance

    # Extract slice data
    x_slice = points[mask, 0]
    y_slice = points[mask, 1]
    z_slice = points[mask, 2]

    # Get potential values (assuming field is named 'phi')
    if 'phi' in mesh.point_data:
        phi_slice = mesh.point_data['phi'][mask]
    else:
        # Try other common names
        available_arrays = list(mesh.point_data.keys())
        if len(available_arrays) > 0:
            field_name = available_arrays[0]
            print(f"  Using field '{field_name}' (phi not found)")
            phi_slice = mesh.point_data[field_name][mask]
        else:
            raise ValueError("No scalar field found in VTK file")

    print(f"  Found {len(x_slice)} points in slice")
    print(f"  Mean z coordinate: {z_slice.mean():.3e} m")

    return x_slice, y_slice, phi_slice

def interpolate_to_grid(x, y, phi, resolution=200):
    """
    Interpolate unstructured slice data onto regular grid
    """
    print(f"\nInterpolating to {resolution}×{resolution} regular grid...")

    # Determine bounds
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    print(f"  X range: [{x_min:.3e}, {x_max:.3e}] m")
    print(f"  Y range: [{y_min:.3e}, {y_max:.3e}] m")

    # Add 5% padding
    dx = (x_max - x_min) * 0.05
    dy = (y_max - y_min) * 0.05

    # Create regular grid
    xi = np.linspace(x_min - dx, x_max + dx, resolution)
    yi = np.linspace(y_min - dy, y_max + dy, resolution)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolate using linear method
    Zi = griddata((x, y), phi, (Xi, Yi), method='linear')

    return Xi, Yi, Zi

def plot_potential_slice(Xi, Yi, Zi, z_value, output_file=None,
                         show_points=False, x=None, y=None):
    """
    Create contour plot of potential distribution
    """
    print("\nGenerating contour plot...")

    fig, ax = plt.subplots(figsize=(12, 9))

    # Determine contour levels based on data range
    phi_min = np.nanmin(Zi)
    phi_max = np.nanmax(Zi)

    print(f"  Potential range: [{phi_min:.3f}, {phi_max:.3f}] V")

    # Create contour levels
    if phi_max - phi_min > 0.1:
        num_levels = 20
    else:
        num_levels = 10

    # Filled contours
    contourf = ax.contourf(Xi * 1e6, Yi * 1e6, Zi, levels=num_levels,
                           cmap='RdYlBu_r', extend='both')

    # Contour lines
    contour_lines = ax.contour(Xi * 1e6, Yi * 1e6, Zi, levels=num_levels,
                               colors='black', linewidths=0.5, alpha=0.4)

    # Add contour labels
    ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.1f V')

    # Colorbar
    cbar = plt.colorbar(contourf, ax=ax, pad=0.02)
    cbar.set_label('Potential φ (V)', fontsize=14, weight='bold')
    cbar.ax.tick_params(labelsize=11)

    # Optionally show data points
    if show_points and x is not None and y is not None:
        ax.scatter(x * 1e6, y * 1e6, c='white', s=2, alpha=0.3,
                  edgecolors='black', linewidths=0.3, label='Mesh nodes')
        ax.legend(fontsize=10)

    # Labels and title
    ax.set_xlabel('x (μm)', fontsize=14, weight='bold')
    ax.set_ylabel('y (μm)', fontsize=14, weight='bold')
    ax.set_title(f'Potential Distribution at z = {z_value*1e6:.1f} μm',
                fontsize=16, weight='bold', pad=20)

    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax.tick_params(labelsize=11)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        print(f"  Saved plot to: {output_file}")
    else:
        plt.show()

    plt.close()

def print_statistics(phi):
    """
    Print statistics about the potential distribution
    """
    print("\n" + "="*50)
    print("SLICE STATISTICS")
    print("="*50)
    print(f"  Minimum potential:     {phi.min():10.6f} V")
    print(f"  Maximum potential:     {phi.max():10.6f} V")
    print(f"  Mean potential:        {phi.mean():10.6f} V")
    print(f"  Std deviation:         {phi.std():10.6f} V")
    print(f"  Potential range:       {phi.max() - phi.min():10.6f} V")
    print("="*50)

def main():
    parser = argparse.ArgumentParser(
        description='Extract 2D horizontal slice from 3D VTK solution and plot contours',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract slice at z=50 microns and display
  python extract_slice.py output/solution.vtk 50e-6

  # Extract slice and save to file
  python extract_slice.py output/solution.vtk 50e-6 -o slice_50um.png

  # Show mesh points and statistics
  python extract_slice.py output/solution.vtk 100e-6 -o slice.png --show-points --stats
        """
    )

    parser.add_argument('vtk_file', help='Input VTK file (output from FEM solver)')
    parser.add_argument('z_value', type=float,
                       help='Z coordinate for slice extraction (meters)')
    parser.add_argument('--output', '-o',
                       help='Output image file (PNG, PDF, SVG, etc.)')
    parser.add_argument('--tolerance', '-t', type=float, default=1e-6,
                       help='Tolerance for z matching in meters (default: 1e-6)')
    parser.add_argument('--resolution', '-r', type=int, default=200,
                       help='Grid resolution for interpolation (default: 200)')
    parser.add_argument('--show-points', action='store_true',
                       help='Show mesh nodes on plot')
    parser.add_argument('--stats', action='store_true',
                       help='Print detailed statistics')

    args = parser.parse_args()

    # Load VTK mesh
    mesh = load_vtk_solution(args.vtk_file)

    # Extract slice at specified z
    x, y, phi = extract_slice_from_vtk(mesh, args.z_value, tolerance=args.tolerance)

    # Print statistics if requested
    if args.stats:
        print_statistics(phi)

    # Interpolate to regular grid
    Xi, Yi, Zi = interpolate_to_grid(x, y, phi, resolution=args.resolution)

    # Create plot
    plot_potential_slice(Xi, Yi, Zi, args.z_value,
                        output_file=args.output,
                        show_points=args.show_points,
                        x=x, y=y)

    print("\n✓ Done!")

if __name__ == '__main__':
    main()
