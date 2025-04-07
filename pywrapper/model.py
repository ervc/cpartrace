import numpy as np
import numpy.typing as npt

from . import constants as const
from .interpolation import interp3d

from .partrace import ModelParams

class Model():
    def __init__(self, fargodir: str, nout: str|int, rescale=True):
        self.fargodir = fargodir
        self.nout = str(nout)
        self.rescale = rescale

        self.read_domain()
        self.fargovars = self.read_varfile()

        self.dusttemp: None | np.ndarray = None

        ### TODO: read in model params

    @classmethod
    def from_partracedir(cls, directory: str, rescale=True):
        params = ModelParams(directory)
        return cls(params['fargodir'], params['nout'], rescale)
    
    def read_varfile(self):
        variables = {}
        with open(self.fargodir+"/variables.par", "r") as f:
            for line in f:
                k,v = line.split()
                variables[k] = v
        return variables

    def get_omega(self, x: float, y: float, z: float) -> float:
        r = np.sqrt(x*x + y*y + z*z)
        return np.sqrt(const.G*const.MSUN/r)/r

    def get_scaleheight(self,x: float,y: float,z: float) -> float:
        ### TODO: make this general
        r = np.sqrt(x*x + y*y + z*z)
        h0 = float(self.fargovars['ASPECTRATIO'])
        flaring = float(self.fargovars['FLARINGINDEX'])
        return h0*r*(r/const.R0)**(flaring)
    
    def get_soundspeed(self,x: float, y: float,z: float) -> float:
        OM = self.get_omega(x,y,z)
        H  = self.get_scaleheight(x,y,z)
        return H*OM
        
    def get_planet_envelope(self, mplan, sma=None):
        ### TODO: make this general
        """
        double hillRadius = sma * pow(model->planetmass/3/MSUN,1.0/3.0);
        double soundspeed = get_soundspeed(model,sma);
        double bondiRadius = 2*G*model->planetmass/soundspeed/soundspeed;
        // planet enevlope is min(hillRadius/4, bondiRadius)
        model->planetEnvelope = (hillRadius/4. < bondiRadius) ? hillRadius/4 : bondiRadius;
        """
        if sma is None:
            sma = const.R0
        m_pl = mplan*const.MEARTH
        hill = sma * (m_pl/3/const.MSUN)**(1/3)
        cs = self.get_soundspeed(sma,0,0)
        bondi = 2*const.G*m_pl/cs/cs
        return min(hill/4, bondi)
        


    def read_domain(self):
        self.phi_edges, self.phi_centers = self.read_domfile(
            self.fargodir+'/domain_x.dat', ghostcells=0, scale=1.
            )
        self.r_edges, self.r_centers = self.read_domfile(
            self.fargodir+'/domain_y.dat', ghostcells=0, scale=1.
            )
        self.theta_edges, self.theta_centers = self.read_domfile(
            self.fargodir+'/domain_z.dat', ghostcells=0, scale=1.
            )
        self.nx = len(self.phi_centers)
        self.ny = len(self.r_centers)
        self.nz = len(self.theta_centers)
        self.shape = (self.nz,self.ny,self.nx)

    def read_domfile(self, filename: str, ghostcells: int = 0, scale: float = 1) -> tuple[npt.NDArray, npt.NDArray]:
        xedges = []
        with open(filename,'r') as f:
            for line in f:
                xedges.append(float(line))
        if ghostcells > 0:
            edges = np.array(xedges[ghostcells:-ghostcells])
        else:
            edges = np.array(xedges)
        edges = edges*scale
        centers = (edges[1:] + edges[:-1])/2
        return edges, centers
    
    def get_spheregrid(self) -> list[npt.NDArray]:
        T = self.theta_centers
        R = self.r_centers
        P = self.phi_centers
        return np.meshgrid(T,R,P,indexing='ij')
    
    def get_cartgrid(self) -> list[npt.NDArray]:
        tt,rr,pp = self.get_spheregrid()
        xx = rr*np.cos(pp)*np.sin(tt)
        yy = rr*np.sin(pp)*np.sin(tt)
        zz = rr*np.cos(tt)
        return [xx,yy,zz]
    
    def read_state(self, state: str) -> npt.NDArray:
        filename = self.fargodir+f'/{state}{self.nout}.dat'
        scale = 1.
        if (state=='gasdens') and (self.rescale):
            scale = const.MASS/const.LEN/const.LEN/const.LEN
        elif ('gasv' in state) and (self.rescale):
            scale = const.LEN/const.TIME
        arr = np.fromfile(filename).reshape(self.shape)
        return arr*scale
    
    def get_rhogrid(self) -> npt.NDArray:
        return self.read_state('gasdens')
    
    def get_rho(self, x: float, y: float, z: float) -> float:
        rhogrid = self.get_rhogrid()
        domain = (self.phi_centers, self.r_centers, self.theta_centers)
        return interp3d(rhogrid, domain, (x,y,z))
    
    def get_Stokes(self, x: float, y: float, z: float, s: float, rho_s: float=2.0) -> float:
        rho = self.get_rho(x,y,z)
        cs = self.get_soundspeed(x,y,z)
        return rho_s*s/rho/cs * self.get_omega(x,y,z)
    
    def get_spherevelgrid(self) -> list[npt.NDArray]:
        gasvphi = self.read_state('gasvx')
        gasvr = self.read_state('gasvy')
        gasvtheta = self.read_state('gasvz')
        return [gasvphi, gasvr, gasvtheta]
    
    def get_cartvelgrid(self) -> list[npt.NDArray]:
        gasvphi, gasvr, gasvtheta = self.get_spherevelgrid()
        tt,rr,pp = self.get_spheregrid()

        gasvx = (gasvr*np.cos(pp)*np.sin(tt)
                + -gasvphi*np.sin(pp)*np.sin(tt)
                + gasvtheta*np.cos(pp)*np.cos(tt))
        gasvy = (gasvr*np.sin(pp)*np.sin(tt)
                + gasvphi*np.cos(pp)*np.sin(tt)
                + gasvtheta*np.sin(pp)*np.cos(tt))
        gasvz = (gasvr*np.cos(tt) + -gasvtheta*np.sin(tt))

        return [gasvx,gasvy,gasvz]
    
    def get_cartvel(self, x: float, y: float, z: float) -> tuple[float, float, float]:
        gasvx, gasvy, gasvz = self.get_cartvelgrid()
        domain = (self.phi_centers, self.r_centers, self.theta_centers)
        pos = (x,y,z)
        gvx = interp3d(gasvx,domain,pos)
        gvy = interp3d(gasvy,domain,pos)
        gvz = interp3d(gasvz,domain,pos)
        return gvx,gvy,gvz
    
    def get_dusttempgrid(self) -> npt.NDArray:
        if self.dusttemp is None:
            self.dusttemp = self.read_state('dusttemp')
        return self.dusttemp
    
    def get_dusttemp(self, x: float, y: float, z: float) -> float:
        dusttemp = self.get_dusttempgrid()
        domain = (self.phi_centers, self.r_centers, self.theta_centers)
        pos = (x,y,z)
        return interp3d(dusttemp, domain, pos)