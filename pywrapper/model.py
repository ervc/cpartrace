import numpy as np
import numpy.typing as npt

from . import constants as const
from .interpolation import interp3d

class Model():
    def __init__(self, fargodir: str, nout: str|int):
        self.fargodir = fargodir
        self.nout = str(nout)

        self.read_domain()

        ### TODO: read in model params

    def get_scaleheight(self,x: float,y: float,z: float):
        ### TODO: make this general
        r = np.sqrt(x*x + y*y + z*z)
        return 0.05*r*(r/const.R0)**(1/4)
        
    def get_planet_envelope(self):
        ### TODO: make this general
        """
        double hillRadius = sma * pow(model->planetmass/3/MSUN,1.0/3.0);
        double soundspeed = get_soundspeed(model,sma);
        double bondiRadius = 2*G*model->planetmass/soundspeed/soundspeed;
        // planet enevlope is min(hillRadius/4, bondiRadius)
        model->planetEnvelope = (hillRadius/4. < bondiRadius) ? hillRadius/4 : bondiRadius;
        """
        sma = const.R0
        m_pl = 300*const.MEARTH
        hill = sma * (m_pl/3/const.MSUN)**(1/3)
        cs = self.get_scaleheight(sma,0,0)
        bondi = 2*const.G*m_pl/cs/cs
        return min(hill/4, bondi)
        


    def read_domain(self):
        self.phi_edges, self.phi_centers = self.read_domfile(
            self.fargodir+'/domain_x.dat', ghostcells=0, scale=1.
            )
        self.r_edges, self.r_centers = self.read_domfile(
            self.fargodir+'/domain_y.dat', ghostcells=3, scale=const.LEN
            )
        self.theta_edges, self.theta_centers = self.read_domfile(
            self.fargodir+'/domain_z.dat', ghostcells=3, scale=1.
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
        if state=='gasdens':
            scale = const.MASS/const.LEN/const.LEN/const.LEN
        elif 'gasv' in state:
            scale = const.LEN/const.TIME
        arr = np.fromfile(filename).reshape(self.shape)
        return arr*scale
    
    def get_rhogrid(self) -> npt.NDArray:
        return self.read_state('gasdens')
    
    def get_rho(self, x: float, y: float, z: float) -> float:
        rhogrid = self.get_rhogrid()
        domain = (self.phi_centers, self.r_centers, self.theta_centers)
        return interp3d(rhogrid, domain, (x,y,z))
    
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
