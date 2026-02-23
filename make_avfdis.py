import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import argparse  

def get_data_from_vtk(vtkfilename):
    """
        extract data for vtkfile 
        retruns name and data=[x,y,z]
    """
    # Load the VTK file
    grid = pv.read(vtkfilename)

    # Access the point coordinates (x, y, z)
    points = grid.points

    # Access the grid dimensions
    dims = grid.dimensions

    # Access cell data
    for name in grid.cell_data:
        data = grid.cell_data[name]

    # Check cell data only contain one array
    if len(list(grid.cell_data.keys()))>1:
        raise ValueError("vtk file contain more then one array")

    # Reshape scalars to 3D array matching dimensions
    nx, ny, nz = dims
    data_3d = data.reshape((nz-1, ny-1, nx-1))   # note the order depends on VTK storage

    return(name, data_3d)


def get_pHvalues():
    """
        obtains pH values for system files returned in a sorted array
    """
    sysnames=gb.glob("sys*.dat")
    sysnames.sort()
    pH=[]
    for fname in sysnames:
        file=open(fname,'r')
        lines=file.readlines()
        file.close()
        pHval=[]
        for line in lines:
            word=line.split()
            if word[0]=="pH" :
                pHval.append(float(word[2]))
        pH.append(pHval[0])
    return(pH)



def make_avfdis():
    """
        Computes average disociation values returned in a sorted array
    """
    # get avpol*.vtk and frdis*.vtk file
    frdisnames=gb.glob("frdis*.vtk")
    frdisnames.sort()
    avpolnames=gb.glob("avpol*.vtk")
    avpolnames.sort()


    # compute avfdis 
    avfdis=[]
    for i in range(len(avpolnames)):
        name, fdis= get_data_from_vtk(frdisnames[i])
        name, xpol = get_data_from_vtk(avpolnames[i])
        frac= np.sum(fdis*xpol)/np.sum(xpol)
        avfdis.append(frac)
    return(avfdis)


def make_plot(x,y,xlabel,ylabel):
    """ 
       plots array x vs array y 
    """	
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot()
    ax.plot(x,y)
    ax.set_xlabel(xlabel) 
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    plt.show()


def parse_cmdline():   
    ap = argparse.ArgumentParser(description="Computes spatial average degree of dissociation.")
    ap.add_argument("--out", default="pHvsavfdis.dat", help="Output filename")
    ap.add_argument("--plot", action="store_true", help="Show a xy plot of data (matplotlib)")
    args = ap.parse_args()
    return(args)	


def main():
    args = parse_cmdline()	
    pH = get_pHvalues()
    avfdis = make_avfdis()
    # output
    data = [(x,y) for (x,y) in zip(pH,avfdis)]
    np.savetxt(args.out,data)
    if args.plot:
        make_plot(pH,avfdis,"pH","avfdis")		

if __name__ == "__main__":
    main()





