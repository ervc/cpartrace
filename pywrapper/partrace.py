from typing import Union, Any

import numpy as np
import numpy.typing as npt

from . import constants as const

def grainLabel(grainsize: float) -> str:
    if grainsize >= 1:
        return f'{grainsize:.0f} cm'
    elif grainsize >= 0.1:
        return f'{grainsize*10:.0f} mm'
    elif grainsize >= 1.e-4:
        return f'{grainsize*1.e4:.0f} '+r'$\mu$m'
    else:
        return f'{grainsize*1.e4:g} '+r'$\mu$m'

class ModelParams:
    def __init__(self,directory: str) -> None:
        """Wrapper for parameters for the partrace run

        Args:
            directory (str): output directory to look in
        """
        self.directory = directory
        self.paramfile = directory+'/inputs.in'

        self.params = self.read_params()

    def read_params(self) -> dict[str, Any]:
        params = {}
        typefuncs = {
            "Strings" : str,
            "Doubles" : float,
            "Integers" : int,
            "Booleans" : bool
        }
        typefunc: type = str
        with open(self.paramfile,"r") as f:
            for line in f:
                if line.strip() in list(typefuncs.keys()):
                    typefunc = typefuncs[line.strip()]
                else:
                    try:
                        key,val = line.split()
                    except ValueError as e:
                        print(line.split())
                        raise e
                    params[key.strip()] = typefunc(val.strip())
        return params
    
    def get_param(self,key):
        return self.params[key.upper()]

    def __getitem__(self,key):
        return self.get_param(key)


class ParticleOutput:
    def __init__(self,outputdir: str, partnumber: int) -> None:
        """Wrapper to help read the particle trajectory outputs. Reads the file
        outputdir/particle{partnumber}.txt

        Args:
            outputdir (str): Output directory
            partnumber (int): number of particle
        """
        params = ModelParams(outputdir)
        t0: float = float(params["t0"])
        tf: float = float(params["tf"])
        dtout: float = float(params["dtout"])
        nt: int = int( (tf-t0)/dtout )+1
        self.fname = outputdir+f'/particle{partnumber}.txt'
        self.x = np.ones(nt)*np.nan
        self.y = np.ones(nt)*np.nan
        self.z = np.ones(nt)*np.nan
        self.vx = np.ones(nt)*np.nan
        self.vy = np.ones(nt)*np.nan
        self.vz = np.ones(nt)*np.nan
        self.times = np.ones(nt)*np.nan
        self.lvls = np.ones(nt)*np.nan

        self.read_partfile()

    def new_read_partfile(self) -> None:
        ts = []
        xs = []
        ys = []
        zs = []
        vxs = []
        vys = []
        vzs = []
        lvls = []
        with open(self.fname,"r") as f:
            for line in f:
                if line=='': continue
                t,x,y,z,vx,vy,vz,lvl = map(float,line.split())
                ts.append(t)
                xs.append(x)
                ys.append(y)
                zs.append(z)
                vxs.append(vx)
                vys.append(vy)
                vzs.append(vz)
                lvls.append(lvl)
        self.times = np.array(ts)
        self.x = np.array(xs)
        self.y = np.array(ys)
        self.z = np.array(zs)
        self.vx = np.array(vxs)
        self.vy = np.array(vys)
        self.vz = np.array(vzs)
        self.lvl = np.array(lvls)
        return

    def read_partfile(self) -> None:
        with open(self.fname,"r") as f:
            i: int = 0
            for line in f:
                if line=='': continue
                (self.times[i],self.x[i],self.y[i],self.z[i],
                 self.vx[i],self.vy[i],self.vz[i],self.lvls[i]) = map(float,line.split())
                i+=1

class ParticleArray:
    def __init__(self,outputdir: str, N_parts: Union[int, list[int]]) -> None:
        """Helper wrapper for a collection of particles from an output array

        Args:
            outputdir (str): Output directory
            N_parts (int | list[int]): List of particle numbers to read in. If
                an integer is supplied, then N_parts = list(range(N_parts))
        """
        params = ModelParams(outputdir)
        t0: float = float(params["t0"])
        tf: float = float(params["tf"])
        dtout: float = float(params["dtout"])
        nt: int = int( (tf-t0)/dtout )+1 #plus one for the initial position
        self.outputdir = outputdir
        N_parts_: list[int] = [0]
        if type(N_parts) == int:
            N_parts_ = list(range(N_parts))
        elif type(N_parts) == list:
            N_parts_ = N_parts
        else:
            raise TypeError("N_parts must be an int or list of ints")
        self.N_parts = N_parts_
        self.len = len(self.N_parts)
        shape = (self.len,nt)
        self.x = np.ones(shape)*np.nan
        self.y = np.ones(shape)*np.nan
        self.z = np.ones(shape)*np.nan
        self.vx = np.ones(shape)*np.nan
        self.vy = np.ones(shape)*np.nan
        self.vz = np.ones(shape)*np.nan
        self.times = np.ones(shape)*np.nan
        self.lvls = np.ones(shape)*np.nan
        for i,n in enumerate(self.N_parts):
            print(i,n,end='\r',flush=True)
            part = ParticleOutput(self.outputdir,n)
            self.x[i] = part.x
            self.y[i] = part.y
            self.z[i] = part.z
            self.vx[i] = part.vx
            self.vy[i] = part.vy
            self.vz[i] = part.vz
            self.times[i] = part.times
            self.lvls[i] = part.lvls
        print('Done reading in particles')

class AllParts():
    def __init__(self,outputdir: str) -> None:
        self.outputdir = outputdir
        self.fname = outputdir+f'/allparts.txt'

        self.read_allparts()

    def read_allparts(self):
        times = []
        starts = []
        ends = []
        statuses = []
        Ntot = 0
        with open(self.fname,'r') as f:
            for n,line in enumerate(f):
                if n==0: continue # header
                try:
                    t,x0,y0,z0,x,y,z,status = line.split()
                except ValueError:
                    continue
                times.append(float(t))
                starts.append(list(map(float,(x0,y0,z0))))
                ends.append(list(map(float,(x,y,z))))
                statuses.append(int(status))
                Ntot+=1
        self.times = np.array(times)
        self.starts = np.array(starts)
        self.ends = np.array(ends)
        self.statuses = np.array(statuses)
        self.Ntot = Ntot

    def get_oob_times(self):
        oobtimes = []
        for s,t in zip(self.statuses, self.times):
            if s==2: oobtimes.append(t)
        return np.array(oobtimes)
                
    def get_accretion_times(self):
        acctimes = []
        for s,t in zip(self.statuses, self.times):
            if s==3: acctimes.append(t)
        return np.array(acctimes)
    
    def get_cross_times(self):
        crosstimes = []
        for s,t in zip(self.statuses, self.times):
            if s==4: crosstimes.append(t)
        return np.array(crosstimes)


def make_biggrid(smallarr: npt.NDArray, flip: str='') -> npt.NDArray:
    """Take an array of the half the disk and extend to the full disk.
    Extends across the 0th axis (theta) of the array.

    Args:
        smallarr (NDArray): half disk array.
        flip (str, optional): How to extend the other side of the disk. 
            "THETA": extends as PI-smallarr.
            "REFLECT": extends as -smallarr
            "" : (default) extends as smallarr

    Returns:
        NDArray: Array extended to the full disk
    """
    nz,ny,nx = smallarr.shape
    bigarr = np.zeros((2*nz,ny,nx))
    bigarr[:nz] = smallarr
    if flip.upper()=="THETA":
        bigarr[nz:] = np.pi-smallarr[::-1]
    elif flip.upper()=="REFLECT":
        bigarr[nz:] = -smallarr[::-1]
    else:
        bigarr[nz:] = smallarr[::-1]
    return bigarr

def get_closest_arg(arr: npt.NDArray, val: float) -> int:
    """Return the index of the closest value in the array
    
    Args:
        arr (NDArray): 1D array
        val (float): value to find
        
    Returns:
        int: Index of closest value
    """
    return np.argmin(np.abs(arr-val))

class ResidenceTimes:
    def __init__(self,outputdir: str) -> None:
        """Helper to read the residence time file from CPartrace

        Args:
            outputdir (str): Output directory
        """
        params = ModelParams(outputdir)
        nz = 32
        ny = 256
        nx = 512
        self.resfile = outputdir+'/residenceTimes.dat'
        resout = np.fromfile(self.resfile)
        self.nparts = resout[0]
        self.restimes = resout[1:].reshape((2*nz,ny,nx))

class Velocities:
    def __init__(self,outputdir: str) -> None:
        """Helper to read the velocities from CPartrace

        Args:
            outputdir (str): Output directory
        """
        nz = 32
        ny = 256
        nx = 2048
        self.velfile = outputdir+'/velocities.dat'
        velout = np.fromfile(self.velfile)
        self.nparts = velout[0]
        self.allvels = velout[1:].reshape((2*nz,ny,nx,3))
        self.vx = self.allvels[:,:,:,0]/self.nparts
        self.vy = self.allvels[:,:,:,1]/self.nparts
        self.vz = self.allvels[:,:,:,2]/self.nparts

class PartTemps:
    def __init__(self,outputdir: str) -> None:
        """Helper to read the binary particle temperatures from interpolation
        
        Args:
            outputdir (str): Output directory
        """
        self.outputdir = outputdir
        self.parttemps = self.read_parttemps()
    
    def read_parttemps(self) -> np.ndarray:
        tempfilename = self.outputdir+"/ctemp.bdat"
        with open(tempfilename, "rb") as f:
            alldata = f.read()
        from struct import unpack
        npart = unpack("i", alldata[:4])[0]
        ntime = unpack("i", alldata[4:8])[0]
        parttemps = np.array(unpack("d"*npart*ntime, alldata[8:]))
        return parttemps.reshape((npart,ntime))
    
class PartLocations:
    def __init__(self, outputdir: str) -> None:
        """Helper to read the zipped particle locations
        
        Args:
            outputdir (str): Output directory
        """
        self.outputdir = outputdir
        self.partlocs = self.read_partlocs()
        self.x = self.partlocs[:,:,0]
        self.y = self.partlocs[:,:,1]
        self.z = self.partlocs[:,:,2]
        self.r = np.sqrt(self.x**2 + self.y**2)
        self.phi = np.arctan2(self.y, self.x)
        # this is a hack but I will fix this later
        self.times = np.linspace(0,1.e6,int(1e5+1))*const.YR

    def read_partlocs(self) -> np.ndarray:
        locfilename = self.outputdir+"/allpos.bdat"
        with open(locfilename,"rb") as f:
            alldata = f.read()
        npart = np.frombuffer(alldata, dtype='i', count=1, offset=0)[0]
        ntime = np.frombuffer(alldata, dtype='i', count=1, offset=4)[0]
        partlocs = np.frombuffer(alldata, dtype='d', count=npart*ntime*3, offset=8)
        partlocs = partlocs.reshape((npart,ntime,3))
        return partlocs

class Crossings:
    def __init__(self,outputdir: str) -> None:
        """Helper to read in the particle crossing times and
        locations from CPartrace
        
        Args:
            outputdir (str): Output directory
        """
        self.crossfile = outputdir+'/partCrossings.txt'
        self.params = ModelParams(outputdir)
        crosstimes = []
        crossx     = []
        crossy     = []
        crossz     = []
        try:
            with open(self.crossfile,'r') as f:
                for line in f:
                    t,x,y,z,*_ = list(map(float,line.split()))
                    crosstimes.append(t)
                    crossx.append(x)
                    crossy.append(y)
                    crossz.append(z)
        except FileNotFoundError:
            pass
        self.crosstimes = np.array(crosstimes)
        self.crossx = np.array(crossx)
        self.crossy = np.array(crossy)
        self.crossz = np.array(crossz)

        self.crossr = np.sqrt(self.crossx**2
                              + self.crossy**2
                              + self.crossz**2)
        self.crossphi = np.arctan2(self.crossy,self.crossx)

    from matplotlib.lines import Line2D
    from matplotlib.axes import Axes
    def plot_cdf(self,ax: Axes, fraction=False, *args,**kwargs) -> list[Line2D]:
        """
        plot the (e)cdf of particles crossed as a function of time
        """
        if len(self.crosstimes)==0:
            X = np.array([0,self.params['tf']])
            Y = np.array([0,0])
        else:
            X = np.sort(self.crosstimes)
            Y = np.arange(1,len(X)+1)
            X = np.concatenate(([0],X,[self.params['tf']]))
            Y = np.concatenate(([0],Y,[Y[-1]]))
        if fraction:
            Y = Y/1000
        return ax.plot(X/const.YR,Y,*args,**kwargs)


import ctypes as ct

def main(infile: str):
    libpartrace = ct.CDLL("src/libpartrace.so")
    libpartrace.run_partrace(ct.c_char_p(str.encode(infile)))

def read_equilibrium_velocity(eqfile: str) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    eqvels = np.fromfile(eqfile).reshape((32,256,2048,3))
    eqvx = eqvels[:,:,:,0]
    eqvy = eqvels[:,:,:,1]
    eqvz = eqvels[:,:,:,2]
    return eqvx,eqvy,eqvz

def solve_equilibrium_velocity_grid(fargodir: str, nout: int|str,
        outfile: str, grainsize: float) -> None:
    libgrainvels = ct.CDLL("src/libgrainvelocities.so")
    res = libgrainvels.py_solve_grid(
        ct.create_string_buffer(fargodir.encode(),100),
        ct.create_string_buffer(str(nout).encode(),5),
        ct.create_string_buffer(outfile.encode(),100),
        ct.c_double(grainsize)
    )
    if res!=0:
        print("something went wrong")

def get_equilibrium_velocity(fargodir: str, nout: int|str, x: float, y: float, z: float, s: float
                             ) -> tuple[float, float, float]:
    libgrainvels = ct.CDLL("src/libgrainvelocities.so")
    pywrapper = libgrainvels.py_get_equilibrium_velocity
    vx = ct.c_double(0.0)
    vy = ct.c_double(0.0)
    vz = ct.c_double(0.0)
    nout = str(nout)
    # return ct.c_wchar_p(fargodir).value
    pywrapper(
        ct.c_char_p(str.encode(fargodir)),
        ct.c_char_p(str.encode(nout)),
        ct.c_double(x),
        ct.c_double(y),
        ct.c_double(z),
        ct.c_double(s),
        ct.byref(vx),
        ct.byref(vy),
        ct.byref(vz)
    )
    return vx.value, vy.value, vz.value

if __name__ == '__main__':
    main("inputs/example.in")