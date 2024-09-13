from typing import Union, Any

import numpy as np
import numpy.typing as npt

from . import constants as const

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
            "Strings:" : str,
            "Doubles:" : float,
            "Integers:" : int,
            "Booleans:" : bool
        }
        typefunc: type = str
        with open(self.paramfile,"r") as f:
            for line in f:
                if line.strip() in list(typefuncs.keys()):
                    typefunc = typefuncs[line.strip()]
                else:
                    key,val = line.split(':')
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
        t0: float = params["t0"]
        tf: float = params["tf"]
        dtout: float = 1*const.YR # params["dtout"]
        nt: int = int( (tf-t0)/dtout )+1+1 #plus one just in case
        self.fname = outputdir+f'/particle{partnumber}.txt'
        self.x = np.ones(nt)*np.nan
        self.y = np.ones(nt)*np.nan
        self.z = np.ones(nt)*np.nan
        self.vx = np.ones(nt)*np.nan
        self.vy = np.ones(nt)*np.nan
        self.vz = np.ones(nt)*np.nan
        self.times = np.ones(nt)*np.nan

        self.read_partfile()

    def read_partfile(self) -> None:
        with open(self.fname,"r") as f:
            i: int = 0
            for line in f:
                if line=='': continue
                (self.times[i],self.x[i],self.y[i],self.z[i],
                 self.vx[i],self.vy[i],self.vz[i]) = map(float,line.split())
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
        t0: float = params["t0"]
        tf: float = params["tf"]
        dtout: float = 1*const.YR # params["dtout"]
        nt: int = int( (tf-t0)/dtout )+1+1 #plus one for the initial positiona and an extra just in case
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
        for i,n in enumerate(self.N_parts):
            part = ParticleOutput(self.outputdir,n)
            self.x[i] = part.x
            self.y[i] = part.y
            self.z[i] = part.z
            self.vx[i] = part.vx
            self.vy[i] = part.vy
            self.vz[i] = part.vz
            self.times[i] = part.times



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

class ResidenceTimes:
    def __init__(self,outputdir: str) -> None:
        """Helper to read the residence time file from CPartrace

        Args:
            outputdir (str): Output directory
        """
        nz = 32
        ny = 256
        nx = 2048
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
    def plot_cdf(self,ax: Axes, fraction=False, *args,**kwargs) -> Line2D:
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