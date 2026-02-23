import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import argparse  as argp


Na=6.02214076e+23
vsol=0.03
conversion_from_l_to_nm3=(0.1**3)/(1e-27) # 1 dm^3 =1l to nm^3

def volume(radius):
    return( (4.0*np.pi/3.0)*(radius**3.0))

def convert_volfraction_numdensity(x,vol):
    rho=[ elem/vol for elem in x]
    return(rho)

def convert_volfraction_concentration(x,vol):
    con=[elem/((Na/conversion_from_l_to_nm3)*vol) for elem in x]
    return(con)

def convert_numdensity_concentration(rho):
    con=[elem/((Na/conversion_from_l_to_nm3)) for elem in rho]
    return(con)

class moleclist:
    def __init__(self,sol,Na,Cl,Hplus,OHmin, automatic):
        self.sol =sol
        self.Na =Na
        self.Cl =Cl
        self.Hplus =Hplus
        self.OHmin =OHmin
       
        
radius=moleclist(0,0,0,0,0,0)
radius.Na=0.20          # radius of Na+ in nm
radius.Cl=0.20           # radius of K+ in nm

volion=moleclist(0,0,0,0,0,0)
volion.Na=volume(radius.Na)/vsol
volion.Cl=volume(radius.Cl)/vsol
volion.Hplus=1.0
volion.OHmin=1.0

def write_data_file(data,outfname):
    fname=open(outfname, 'w')   
    for elem in data:
        strtmp=str(elem[0])+" "+str(elem[1])+"\n"
        fname.write(strtmp)       
    fname.close()  

def density_ion(xsol,pot,expmu,vol_ion,z_ion):
    xion=[]
    for i in range(len(xsol)):
        xion.append(expmu*(xsol[i]**vol_ion)*np.exp(-pot[i]*z_ion))
    return(xion)
        

def init_xbulk(pHbulk,cNaCl): 
    """ 
        computes volume fraction of comppmnets of reservoir of given pH and cNaCl concentration  
    """
    pHbulk=float(pHbulk)
    cNaCl=float(cNaCl)
    xbulk=moleclist(0,0,0,0,0,0)
    pKw=14.0
    # .. initializations of input dependent variables
    cHplus = (10.0)**(-pHbulk) #concentration H+ in bulk
    pOHbulk = pKw -pHbulk       
    cOHmin  = (10.0)**(-pOHbulk) # concentration OH- in bulk
    xbulk.Hplus = (cHplus*Na/(1.0e24))*(vsol) # volume fraction H+ in bulk vH+=vsol
    xbulk.OHmin = (cOHmin*Na/(1.0e24))*(vsol) # volume fraction OH- in bulk vOH-=vsol
    # NaCl in solution 
    xNaClsalt = (cNaCl*Na/(1.0e24))*((volion.Na+volion.Cl)*vsol) # volume fraction NaCl salt in mol/l
    if(pHbulk<=7):       # pH<= 7
        xbulk.Na=xNaClsalt*volion.Na/(volion.Na+volion.Cl)  
        xbulk.Cl=xNaClsalt*volion.Cl/(volion.Na+volion.Cl) +(xbulk.Hplus -xbulk.OHmin)*volion.Cl  # NaCl+ HCl
    else:                      # pH >7 
        xbulk.Na=xNaClsalt*volion.Na/(volion.Na+volion.Cl) +(xbulk.OHmin -xbulk.Hplus)*volion.Na  # NaCl+ NaOH  
        xbulk.Cl=xNaClsalt*volion.Cl/(volion.Na+volion.Cl)  

    xbulk.sol=1.0 -xbulk.Hplus-xbulk.OHmin-xbulk.Cl-xbulk.Na
    return(xbulk)

def init_expmu(xbulk):
    """ 
        computes expmu prefactor for computation of number density 
    """
    expmu =moleclist(0,0,0,0,0,0)
    expmu.Na    = xbulk.Na   /(xbulk.sol**volion.Na) 
    expmu.Cl    = xbulk.Cl   /(xbulk.sol**volion.Cl)
    expmu.Hplus = xbulk.Hplus/xbulk.sol  
    expmu.OHmin = xbulk.OHmin/xbulk.sol 
    return(expmu)

def make_xion(xsol,pot,expmu,vol_ion): 
    """ 
        computes volume fraction of ions 
    """
    xion=moleclist(0,[],[],[],[],[])
    xion.Na=density_ion(xsol,pot,expmu.Na,vol_ion.Na,1)
    xion.Cl=density_ion(xsol,pot,expmu.Cl,vol_ion.Cl,-1)
    xion.Hplus=density_ion(xsol,pot,expmu.Hplus,vol_ion.Hplus,1)
    xion.OHmin=density_ion(xsol,pot,expmu.OHmin,vol_ion.OHmin,-1) 
    return(xion)

def make_concentration_ion(xion,vol_ion):
    """ 
        computes concentration of ions in mol/l 
    """
    cion=moleclist(0,[],[],[],[],[])
    cion.Na=convert_volfraction_concentration(xion.Na,vol_ion.Na*vsol)
    cion.Cl=convert_volfraction_concentration(xion.Cl,vol_ion.Cl*vsol)
    cion.Hplus=convert_volfraction_concentration(xion.Hplus,vol_ion.Hplus*vsol)
    cion.OHmin=convert_volfraction_concentration(xion.OHmin,vol_ion.OHmin*vsol)
    return(cion)

def get_data_from_vtk(vtkfilename):
    """
        extract data for vtkfile 
        retuns name of data and data=[x,y,z]
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


def get_pH_csalt_values():
    """
        obtains pH values for system files returned in a sorted array
    """
    sysnames=gb.glob("sys*.dat")
    sysnames.sort() # important to allign value with other output file
    pH=[]
    csalt=[]
    for fname in sysnames:
        file=open(fname,'r')
        lines=file.readlines()
        file.close()
        pHval=[]
        csaltval=[]
        for line in lines:
            word=line.split()
            if word[0]=="pH" :
                pHval.append(float(word[2]))
            if word[0]=="csalt" :
                csaltval.append(float(word[2]))    
        pH.append(pHval[0])
        csalt.append(csaltval[0])
    return(pH,csalt)

def get_delta_value():
    """
        obtains delta value for systemfile : assume at least on system file present
    """
    sysnames=gb.glob("sys*.dat")
    sysnames.sort() # important to allign value with other output file
    for fname in [sysnames[0]]:
        file=open(fname,'r')
        lines=file.readlines()
        file.close()
        pHval=[]
        csaltval=[]
        for line in lines:
            word=line.split()
            if word[0]=="delta" :
                delta=float(word[2])
    return(delta)



def write_numpy_to_vtkfile(data,fnamevtk,delta): 
    """
        converts numpy 3d data array into a vtkdatafile named fnamevtk
        delta is spacing of volume cell
    """

    dimx,dimy,dimz = data.shape
    f=open(fnamevtk,'w')
    s="# vtk DataFile Version 2.0\n"
    s=s+"title\n"
    s=s+"ASCII\n"
    s=s+"DATASET STRUCTURED_GRID\n"
    s=s+"DIMENSIONS "+str(dimz+1)+" "+str(dimy+1)+" "+str(dimx+1)+"\n"
    s=s+"POINTS "+str((dimx+1)*(dimy+1)*(dimz+1))+" float\n"
    f.write(s)
    for i in range(0,dimx+1):
        for j in range(0,dimy+1):
            for k in range(0,dimz+1):
                s=str(i*delta)+"   "+str(j*delta)+"   "+str(k*delta)+"\n"
                f.write(s)
    s="CELL_DATA "+str(dimx*dimy*dimz)+"\n"
    s=s+"SCALARS var float 1\n"
    s=s+"LOOKUP_TABLE default\n"
    f.write(s)
    for i in range(dimx):
        for j in range(dimy):
            for k in range(dimz):    
                s = str(data[i][j][k])+"\n"
                f.write(s)
    f.close()          


def make_pHdistribution():
    """
        Computes pH form xsol and potential found in vtk files
    """
    # get avpol*.vtk and frdis*.vtk file
    xsolnames=gb.glob("avsol*.vtk")
    xsolnames.sort()
    potnames=gb.glob("poten*.vtk")
    potnames.sort()
   
    # get pH and salt concentration  from  volume fraction 
    pHval, csaltval = get_pH_csalt_values()
    # get delta
    delta = get_delta_value()

    for i in range(len(pHval)):
        pH=pHval[i]
        cNaCl=csaltval[i]
        # Calculate local volume fractions and density of ions
        xbulk = init_xbulk(pH,cNaCl)
        expmu = init_expmu(xbulk)
        # get pi and psi 
        xsolname, xsol = get_data_from_vtk(xsolnames[i])
        potname, pot = get_data_from_vtk(potnames[i])
        # make ion volume fraction and concentration 
        xion = make_xion(xsol,pot,expmu,volion)
        cion = make_concentration_ion(xion,volion)
        # convert Hplus concentration to log10 scale
        pH = np.array(-np.log10(cion.Hplus))
        # write to vtk file 
        fnamevtk=xsolnames[i].replace("avsol","pH",1)
        write_numpy_to_vtkfile(pH,fnamevtk,delta)


def  make_plot_xslice(data_3d,array_name):
    """
        make a contour plot from slice in x-direction of data
    """
    nx, ny, nz =data_3d.shape
    x_slice = data_3d[nx // 2, :, :]
    print(f"Middle x-slice shape: {x_slice.shape}")
    plt.imshow(x_slice, origin="lower", cmap="viridis")
    plt.title(f"Middle x-slice of {array_name}")
    plt.colorbar(label=array_name)
    plt.show()



def main():
    make_pHdistribution()



if __name__ == "__main__":
    main()





