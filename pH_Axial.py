import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import argparse
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d


def get_data_from_vtk(vtkfilename):
    """
        extract data for vtkfile
        retruns name and data=[x,y,z]
    """
    # Load the VTK file
    grid = pv.read(vtkfilename)

    # Check the type
    print(type(grid))  # should show pyvista.core.pointset.StructuredGrid

    # Access the point coordinates (x, y, z)
    points = grid.points
    print("Points shape:", points.shape)  # (N, 3)

    # Access the grid dimensions
    dims = grid.dimensions
    print("Grid dimensions:", dims)

    # Access point data (scalars, vectors, etc.)
    print("Available point data arrays:", grid.point_data.keys())

    # Access cell data
    print("Available cell data arrays:", grid.cell_data.keys())

    for name in grid.cell_data:
        data = grid.cell_data[name]
        print(name, data.shape)

    # Check cell data only contain one array
    if len(list(grid.cell_data.keys())) > 1:
        raise ValueError("vtk file contain more then one array")

    # Reshape scalars to 3D array matching dimensions
    nx, ny, nz = dims
    data_3d = data.reshape((nz - 1, ny - 1, nx - 1))  # note the order depends on VTK storage

    print("3D scalar field shape:", data_3d.shape)

    return (name, data_3d)

def make_av_axial_2D(data_3d, i):
    i += 1
    dr = 0.5  # nm
    dz = 0.5  # nm

    Nx, Ny, Nz = data_3d.shape
    print("Nx:", Nx,"Ny:",Ny,"Nz:",Nz)

    # Index-based coordinates → physical units
    x = (np.arange(Nx) - (Nx - 1) / 2) * dr
    y = (np.arange(Ny) - (Ny - 1) / 2) * dr
    z_vals = np.arange(Nz) * dz
    print(x ,y )

    # Expand to 3D
    x = x[:, None, None]
    y = y[None, :, None]
    z = np.broadcast_to(z_vals[None, None, :], data_3d.shape)
    # Radius (FIXED)
    x0, y0 = 0.0, 0.0
    r2d = np.sqrt(x ** 2 + y ** 2)
    r = np.broadcast_to(r2d, data_3d.shape)
    r_max = (Nx/2)*dr

    # Flatten (true NumPy arrays)
    f_flat = np.asarray(data_3d, dtype=float).ravel()
    r_flat = np.asarray(r, dtype=float).ravel()
    z_flat = np.asarray(z, dtype=float).ravel()

    mask=r_flat <=r_max

    f_flat_limited = f_flat[mask]
    r_flat_limited = r_flat[mask]
    z_flat_limited = z_flat[mask]

    # Bin + average
    nr = Nx / 2
    nz = Nz
    avg, r_edges, z_edges, _ = binned_statistic_2d(
        r_flat_limited, z_flat_limited, f_flat_limited,
        statistic="mean",
        bins=[nr, nz]
    )


    # Plot
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    r_pos = r_centers
    r_neg = -r_centers[1:][::-1]

    avg_pos = avg
    avg_neg = avg[1:, :][::-1, :]

    # Combine
    r_full = np.concatenate((r_neg, r_pos))
    avg_full = np.vstack((avg_neg, avg_pos))



    R, Z = np.meshgrid(r_full, z_centers, indexing="ij")

    out = np.column_stack([
        R.ravel(),
        Z.ravel(),
        avg_full.ravel()
    ])
    #plot_axial_average(z_centers, r_full, avg_full)
    np.savetxt(
        f"axisymmetric_avg_pH_{i:03}.csv",
        out,
        delimiter=",",
        header="r [nm],z [nm],pH",
        comments=""
    )

def plot_axial_average(z_centers, r_full, avg_full):
    plt.pcolormesh(
        z_centers,
        r_full,
        avg_full,
        shading="auto"
    )
    plt.xlabel("z [nm]")
    plt.ylabel("r [nm]")
    plt.colorbar(label="pH")
    plt.show()


def get_pHvalues():
    sysnames = gb.glob("sys*.dat")
    sysnames.sort()
    pH = []
    for fname in sysnames:
        file = open(fname, 'r')
        lines=file.readlines()
        file.close()
        pHval=[]
        for line in lines:
            word=line.split()
            if word[0]=="pH" :
                pHval.append(float(word[2]))
        pH.append(pHval[0])
    return(pH)


def get_average():

        avpolnames = gb.glob("avpol.*.vtk")
        avpolnames.sort()
        avAvpol = []
        for file in avpolnames:
            average = get_data_from_vtk(file)
            avAvpol.append(float(average))
        #average_poly , count = read_vtk_file("avpol.001.vtk")
        print(avAvpol)
        return(avAvpol)


def main():
    pHnames = gb.glob("ph.*.vtk")
    pHnames.sort()
    pH = get_pHvalues()
    for i in range(len(pH)):
        name , data = get_data_from_vtk(pHnames[i])
        make_av_axial_2D(data, i)




if __name__ == "__main__":
    main()