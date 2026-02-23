import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import sys


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
    if len(list(grid.cell_data.keys()))>1:
        raise ValueError("vtk file contain more then one array")
 
    # Reshape scalars to 3D array matching dimensions
    nx, ny, nz = dims
    data_3d = data.reshape((nz-1, ny-1, nx-1))   # note the order depends on VTK storage

    print("3D scalar field shape:", data_3d.shape)

    return(name, data_3d)



def main():
    print(len(sys.argv))

    if len(sys.argv) != 3:
        print("Usage: python convert_vtk_to_numpy_vista.py <input_vtk_file>  npy or txt")
        sys.exit(1)

    vtkfilename=sys.argv[1]
    name, data = get_data_from_vtk(vtkfilename)
    
    match sys.argv[2]:
        case "npy":
            npyfilename=vtkfilename.replace(".vtk",".npy")
            # Save the array to a .npy file
            np.save(npyfilename, data)
        case "txt":
            # Reshape the 3D array to 2D
            reshaped_array = data.reshape(data.shape[0], -1)
            txtfilename=vtkfilename.replace(".vtk",".txt")
            # Save the reshaped array to a text file
            np.savetxt(txtfilename, reshaped_array, delimiter=',')

    # To load and reshape back:
    # loaded_reshaped_array = np.loadtxt('output.txt', delimiter=',')
    # loaded_3d_array = loaded_reshaped_array.reshape(data.shape)


if __name__ == "__main__":
    main()





