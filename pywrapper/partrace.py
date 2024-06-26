from typing import Union, Any

import numpy as np
import numpy.typing as npt

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
        dtout: float = params["dtout"]
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
        dtout: float = params["dtout"]
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
        self.vx = self.allvels[:,:,:,0]
        self.vy = self.allvels[:,:,:,1]
        self.vz = self.allvels[:,:,:,2]
