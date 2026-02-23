import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import argparse as argp
import vtk as pyvtk

def convert_vtk_to_ascii(input_file, output_file):
    """
    Converts a VTK file (e.g., .vtp, .vtu, .vtk) to its ASCII representation.

    Args:
        input_file (str): Path to the input VTK file.
        output_file (str): Path for the output ASCII VTK file.
    """
    # Determine the appropriate reader based on the file extension
    if input_file.endswith(".vtp"):
        reader = pyvtk.vtkXMLPolyDataReader()
        writer = pyvtk.vtkXMLPolyDataWriter()
    elif input_file.endswith(".vtu"):
        reader = pyvtk.vtkXMLUnstructuredGridReader()
        writer = pyvtk.vtkXMLUnstructuredGridWriter()
    elif input_file.endswith(".vtk"):
        reader = pyvtk.vtkGenericDataObjectReader() # For legacy .vtk files
        writer = pyvtk.vtkDataSetWriter() # For legacy .vtk files
    else:
        print(f"Error: Unsupported VTK file type for input file: {input_file}")
        return

    reader.SetFileName(input_file)
    reader.Update()

    writer.SetFileName(output_file)
    writer.SetInputConnection(reader.GetOutputPort())
    # writer.SetDataModeToAscii() # Crucial for setting ASCII output
    writer.SetFileTypeToASCII()
    writer.Write()

    #print(f"Successfully converted '{input_file}' to ASCII format in '{output_file}'")

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


def get_dims(sysname):
    """
        obtains x-, y nd z dimension of value from file named 
        input: *str
                sysname  
        output: *int 
                nx, ny, nz  
    """
    fname=sysname
    file=open(fname,'r')
    lines=file.readlines()
    file.close()
    for line in lines:
        word=line.split()
        if word[0]=="nx" :
            nx=int(word[2])
        if word[0]=="ny" :
            ny=int(word[2])
        if word[0]=="nz" :
            nz=int(word[2])

    return(nx,ny,nz)

def get_data_from_dat(datfilename):
    """
        extract data for datfile 
        returns data=[x,y,z]
    """
    # Load the dat file
    data = np.loadtxt(datfilename)
    #
    sysname="system."+datfilename.split(".",1)[1]
    nx, ny, nz = get_dims(sysname)
    # Reshape scalars to 3D array matching dimensions
    data3d=data.reshape(nx,ny,nz)
    data3d.shape
    return(data3d)




def  make_plot_xslice(data_3d,array_name,nslice):
    """ 
        make a contour plot from slice in x-direction of data
    """
    #nx, ny, nz =data_3d.shape
    x_slice = data_3d[nslice, :, :]
    plt.imshow(x_slice, origin="lower", cmap="viridis")
    plt.title(f"x-slice of {array_name} at x={nslice}")
    plt.colorbar(label=array_name)
    plt.show()


def  make_plot_yslice(data_3d,array_name,nslice):
    """ 
        make a contour plot from slice in y-direction of data
    """
    #nx, ny, nz =data_3d.shape
    x_slice = data_3d[:, nslice, :]
    plt.imshow(x_slice, origin="lower", cmap="viridis")
    plt.title(f"y-slice of {array_name} at y={nslice}")
    plt.colorbar(label=array_name)
    plt.show()


def  make_plot_zslice(data_3d,array_name,nslice):
    """ 
        make a contour plot from slice in z-direction of data
    """
    #nx, ny, nz =data_3d.shape
    x_slice = data_3d[ :, :,nslice ]
    plt.imshow(x_slice, origin="lower", cmap="viridis")
    plt.title(f"z-slice of {array_name} at z={nslice}")
    plt.colorbar(label=array_name)
    plt.show()


def parse_cmdline():   
    ap = argp.ArgumentParser(
        description="Computes and plots slice of 3D data."
    )
    ap.add_argument("filename", type=str, help="Input filename")
    ap.add_argument("--filetype", type=str, default="vtk", help="Type of input file: vtk or dat" )
    ap.add_argument("--plane",type=str, default="z",help="Select plane of slice (default z)")
    ap.add_argument("--coordplane",type=int, default="-1",help="Select coordinate plane of slice (default middle)")
    args = ap.parse_args()

    if args.plane not in ["z","x","y"]:
        raise ValueError("plane : x, y, or z")
    if args.filetype not in ["vtk","dat"]:
        raise ValueError("filetype: vtk","dat")

    return(args)	


def determine_value_coordplane(args,dims):
    """
        Determine coordinate (int value ) of plane slice based in args and dims data
        return: nslice
    """
    nx, ny, nz = dims
    nslice = args.coordplane
    if args.coordplane== -1  :

        match args.plane:
            case "x":
                nslice = nx // 2
            case "y":
                nslice = ny // 2    
            case "z":
                nslice = nz // 2  
    is_in_range = True
    match args.plane:
        case "x":
            is_in_range = ((0<=nslice) and ( nslice <= nx-1))
        case "y":
            is_in_range = ((0<=nslice) and ( nslice <= ny-1))
        case "z":
            is_in_range = ((0<=nslice) and ( nslice <= nz-1))
    if  not is_in_range:
        raise ValueError(f"coordinate plane {args.coordplane} outside range {shape}")  
    return(nslice)    
   

def main():

    args = parse_cmdline()

    match args.filetype:
        case "vtk":
            vtkfilename=args.filename
            name, data = get_data_from_vtk(vtkfilename)
        case "dat":
            datfilename=args.filename
            data = get_data_from_dat(datfilename)
            name=""  
            
    
    nslice = determine_value_coordplane(args,data.shape)

    match args.plane:
        case "x":
            make_plot_xslice(data,name,nslice)
        case "y":
            make_plot_yslice(data,name,nslice)
        case "z":
            make_plot_zslice(data,name,nslice)
            
            
        


if __name__ == "__main__":
    main()





