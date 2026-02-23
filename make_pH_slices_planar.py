import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import argparse  as argp


def linearIndexFromCoordinate(nx,ny,nz,x,y,z):
    idx = x + nx*y + nx*ny*z  -(nx*(1+ny)) 
    return(idx)

def coordinateFromLinearIndex(nx,ny,idx):
    idxtmp=idx
    x = ((idxtmp-1)%nx)+1
    idxtmp =int((idxtmp-1)/nx)+1
    y = ((idxtmp-1)%ny)+1
    idxtmp = int((idxtmp-1)/ny)+1
    z = idxtmp
    return([int(x),int(y),int(z)])

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
    def __init__(self,sol,Na,Cl,K,Ca,Mg,Hplus,OHmin,Rb, automatic):
        self.sol =sol
        self.Na =Na
        self.Cl =Cl
        self.K =K
        self.Ca =Ca
        self.Mg =Mg
        self.Hplus =Hplus
        self.OHmin =OHmin
        self.Rb =Rb
             

radius=moleclist(0,0,0,0,0,0,0,0,0,0)
radius.Na=0.102          # radius of Na+ in nm
radius.K=0.138           # radius of K+ in nm
radius.Cl=0.181          # radius of Cl- in nm
radius.Ca=0.106          # radius of Ca2+ in nm
radius.Rb=0.152          # radius of Rb+ in nm 
radius.Mg=0.072          # radius of Mg++ in nm
      


volion=moleclist(0,0,0,0,0,0,0,0,0,0)
volion.Na=volume(radius.Na)/vsol
volion.K=volume(radius.K)/vsol
volion.Cl=volume(radius.Cl)/vsol
volion.Ca=volume(radius.Ca)/vsol
volion.Rb=volume(radius.Rb)/vsol
volion.Mg=volume(radius.Mg)/vsol
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
        

def init_xbulk(pHbulk,cNaCl,cKCl,cRbCl,cMgCl2,cCaCl2):
    pHbulk=float(pHbulk)
    cNaCl=float(cNaCl)
    cKCl=float(cKCl)
    cRbCl=float(cRbCl)
    cMgCl2=float(cMgCl2)
    cCaCl2=float(cCaCl2)
    xbulk=moleclist(0,0,0,0,0,0,0,0,0,0)
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

    # RbCl in solution 
    xRbClsalt = (cRbCl*Na/(1.0e24))*((volion.Rb+volion.Cl)*vsol)
    xbulk.Rb = xRbClsalt*volion.Rb/(volion.Rb+volion.Cl)  
    xbulk.Cl = xbulk.Cl+xRbClsalt*volion.Cl/(volion.Rb+volion.Cl)  
    # KCl in solution 
    xKClsalt = (cKCl*Na/(1.0e24))*((volion.K+volion.Cl)*vsol) # volume fraction KCl salt in mol/l
    xbulk.K  = xKClsalt*volion.K/(volion.K+volion.Cl)  
    xbulk.Cl = xbulk.Cl+xKClsalt*volion.Cl/(volion.K+volion.Cl)  
    # CaCl2 in solution 
    xCaCl2salt = (cCaCl2*Na/(1.0e24))*((volion.Ca+2.0*volion.Cl)*vsol) # volume fraction CaCl2 in mol/l
    xbulk.Ca=xCaCl2salt*volion.Ca/(volion.Ca+2.0*volion.Cl)
    xbulk.Cl=xbulk.Cl+ xCaCl2salt*2.0*volion.Cl/(volion.Ca+2.0*volion.Cl)
    # MgCl2 in solution 
    xMgCl2salt = (cMgCl2*Na/(1.0e24))*((volion.Mg+2.0*volion.Cl)*vsol) # volume fraction MgCl2 in mol/l
    xbulk.Mg=xMgCl2salt*volion.Mg/(volion.Mg+2.0*volion.Cl)
    xbulk.Cl=xbulk.Cl+ xMgCl2salt*2.0*volion.Cl/(volion.Mg+2.0*volion.Cl)

    xbulk.sol=1.0 -xbulk.Hplus-xbulk.OHmin-xbulk.Cl-xbulk.Na-xbulk.K-xbulk.Ca-xbulk.Rb-xbulk.Mg
    
    return(xbulk)



def init_expmu(xbulk):
    """ 
        computes volume fraction of comppmnets of reservoir of given pH and cNaCl concentration  
    """
    expmu =moleclist(0,0,0,0,0,0,0,0,0,0)
    expmu.Na    = xbulk.Na   /(xbulk.sol**volion.Na) 
    expmu.Cl    = xbulk.Cl   /(xbulk.sol**volion.Cl)
    expmu.K     = xbulk.K    /(xbulk.sol**volion.K)
    expmu.Ca    = xbulk.Ca   /(xbulk.sol**volion.Ca) 
    expmu.Mg    = xbulk.Mg   /(xbulk.sol**volion.Mg) 
    expmu.Rb    = xbulk.Rb   /(xbulk.sol**volion.Rb)
    expmu.Hplus = xbulk.Hplus/xbulk.sol  
    expmu.OHmin = xbulk.OHmin/xbulk.sol 
    return(expmu)



def make_xion(xsol,pot,expmu,vol_ion): 
    """ 
        computes volume fraction of ions 
    """
    xion=moleclist(0,[],[],[],[],[],[],[],[],[])
    xion.Na=density_ion(xsol,pot,expmu.Na,vol_ion.Na,1)
    xion.K =density_ion(xsol,pot,expmu.K,vol_ion.K,1)
    xion.Cl=density_ion(xsol,pot,expmu.Cl,vol_ion.Cl,-1)
    xion.Mg=density_ion(xsol,pot,expmu.Mg,vol_ion.Mg,2)
    xion.Ca=density_ion(xsol,pot,expmu.Ca,vol_ion.Ca,2)
    xion.Hplus=density_ion(xsol,pot,expmu.Hplus,vol_ion.Hplus,1)
    xion.OHmin=density_ion(xsol,pot,expmu.OHmin,vol_ion.OHmin,-1) 
    return(xion)
    

def make_concentration_ion(xion,vol_ion):
    """ 
        computes concentration of ions in mol/l 
    """
    cion=moleclist(0,[],[],[],[],[],[],[],[],[])
    cion.Na=convert_volfraction_concentration(xion.Na,volion.Na*vsol)
    cion.K=convert_volfraction_concentration(xion.K,volion.K*vsol)
    cion.Cl=convert_volfraction_concentration(xion.Cl,volion.Cl*vsol)
    cion.Mg=convert_volfraction_concentration(xion.Mg,volion.Mg*vsol)
    cion.Ca=convert_volfraction_concentration(xion.Ca,volion.Ca*vsol)
    cion.Hplus=convert_volfraction_concentration(xion.Hplus,volion.Hplus*vsol)
    cion.OHmin=convert_volfraction_concentration(xion.OHmin,volion.OHmin*vsol)
    return(cion)


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


def get_pH_csalt_value(sysname):
    """
        obtains pH values for system files returned in a sorted array
    """
   
    fname=sysname
    file=open(fname,'r')
    lines=file.readlines()
    file.close()

    pHval=[]
    cNaClval=[]
    cKClval =[]
    cMgCl2val = []  

    for line in lines:
        word=line.split()
        if word[0]=="pHbulk":
            pHval.append(float(word[2]))
        if word[0]=="cNaCl" :
            cNaClval.append(float(word[2]))    
        if word[0]=="cKCl" :
            cKClval.append(float(word[2]))    
        if word[0]=="cMgCl2" :
            cMgCl2val.append(float(word[2]))    
    pH=pHval[0]
    cNaCl=cNaClval[0]
    cKCl=cKClval[0]
    cMgCl2=cMgCl2val[0]

    return(pH,cNaCl,cKCl,cMgCl2)

def get_delta_value(sysname):
    """
        obtains delta value for systemfile : assume at least on system file present
    """

    file=open(sysname,'r')
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

def write_numpy_to_datfile(data,fnamedat,delta): 
    """
        converts numpy 3d data array into a dat datafile named fnamedat
        delta is spacing of volume cell
    """
    print(data.shape)
    n=np.prod(data.shape)
    print(n)
    data1d=data.reshape(n)
    np.savetxt(fnamedat,data1d)



def make_pHdistribution_vtk(systemfile):
    """
        Computes pH form xsol and potential found in vtk files
    """
    # get avpol*.vtk and frdis*.vtk file
    xsolname=systemfile.replace("system","xsol")
    xsolname=xsolname.replace("dat","vtk")
    potname=xsolname.replace("xsol","potential")
    # get pH and salt concentration  from  volume fraction 
    pH, cNaCl, cKCl, cMgCl2 = get_pH_csalt_values()
    # get delta
    delta = get_delta_value()

    # Calculate local volume fractions and density of ions
    xbulk=init_xbulk(pH,cNaCl,cKCl,0,cMgCl2,0)
    expmu=init_expmu(xbulk)
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
    return (pH)



def make_pHdistribution_dat(systemfile):
    """
        Computes pH form xsol and potential found in vtk files
    """
    # get pH and salt concentration from volume fraction 
    pHval, cNaCl,cKCl,cMgCl2 = get_pH_csalt_value(systemfile)

    # get avsol*.vtk and frdis*.vtk file
    xsolname=systemfile.replace("system","xsol")
    potname=xsolname.replace("xsol","potential")
   
    # get delta
    delta = get_delta_value(systemfile)
   
    # Calculate local volume fractions and density of ions
    xbulk=init_xbulk(pHval,cNaCl,cKCl,0,cMgCl2,0)
    expmu=init_expmu(xbulk)
    # get pi and psi 
    xsol = get_data_from_dat(xsolname)
    pot = get_data_from_dat(potname)
    # make ion volume fraction and concentration 
    xion = make_xion(xsol,pot,expmu,volion)
    cion = make_concentration_ion(xion,volion)
    # convert Hplus concentration to log10 scale
    pH = np.array(-np.log10(cion.Hplus))
    # write to dat file 
    fnamedat=xsolname.replace("xsol","pH",1)
    write_numpy_to_datfile(pH,fnamedat,delta)
    return (pH)



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
        description="Computes pH and plots slice of 3D data."
    )
    ap.add_argument("filename", type=str, help="Input filename")
    ap.add_argument("--filetype", type=str, default="dat", help="Type of input file: vtk or dat" )
    ap.add_argument("--plane",type=str, default="z",help="Select plane of slice (default z)")
    ap.add_argument("--coordplane",type=int, default="-1",help="Select coordinate plane of slice (default middle)")
    args = ap.parse_args()

    if args.plane not in ["z","x","y"]:
        raise ValueError("plane : x, y, or z")
    if args.filetype not in ["vtk","dat"]:
        raise ValueError("filetype: vtk","dat")
    if not ( "system" in args.filename ) :
        raise ValueError("file not a system file")

    return(args)	


def determine_value_coordplane(args,shape):
    """
    """
    nx, ny, nz = shape  
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
            is_in_range = ((0<=nslice) and ( nslice  <= ny-1))
        case "z":
            is_in_range = ((0<=nslice) and ( nslice <= nz-1))
    if  not is_in_range:
        raise ValueError("coordinate plane {args.coordplane} outside range {shape}")  
    return(nslice)    
   


def main():

    args = parse_cmdline()

    print(args.filetype)

    match args.filetype:
        case "vtk":
            vtkfilename=args.filename
            data = make_pHdistribution_vtk(args.filename)
            name ="pH"
        case "dat":
            datfilename=args.filename
            data = make_pHdistribution_dat(args.filename)
            name="pH"  
            
            
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






