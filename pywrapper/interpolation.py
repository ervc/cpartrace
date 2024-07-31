import numpy as np
import numpy.typing as npt

class InterpolationError(Warning):
    pass
class OutOfBoundsError(Warning):
    pass

def interp3d(arr: npt.NDArray, domain: tuple[npt.NDArray,...], coords: tuple[float,...]) -> float:
    phicenters, rcenters, thetacenters = domain
    nx = len(phicenters)
    ny = len(rcenters)
    nz = len(thetacenters)
    domshape = (nz,ny,nx)
    arrshape = arr.shape
    if domshape != arrshape:
        raise InterpolationError("Interpolation array and domain sizes do not match")
    
    x,y,z = coords
    if z<0:
        z=-z
    
    phi = np.arctan2(y,x)
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    assert theta <= np.pi/2

    if (r<rcenters[0] or r>rcenters[-1]):
        raise OutOfBoundsError("R is out of bounds on interpolation")
    if theta < thetacenters[0]:
        raise OutOfBoundsError("Theta is out of bounds")
    
    if phi < phicenters[0]:
        i = nx-1
        phi += 2*np.pi
    else:
        for i in range(nx):
            if phicenters[i] > phi:
                break
        i-=1
    for j in range(ny):
        if rcenters[j] > r:
            break
    j-=1
    if theta > thetacenters[-1]:
        k = nz-1
    for k in range(nz):
        if thetacenters[k] < theta:
            break
    k-=1

    if i==nx-1:
        ip = 0
    else:
        ip = i+1
    jp = j+1
    if k!=nz-1:
        # full 3d interpolation
        kp = k+1
        c000 = arr[k ,j ,i ]
        c100 = arr[k ,j ,ip]
        c010 = arr[k ,jp,i ]
        c110 = arr[k ,jp,ip]
        c001 = arr[kp,j ,i ]
        c101 = arr[kp,j ,ip]
        c011 = arr[kp,jp,i ]
        c111 = arr[kp,jp,ip]
        x0 = phicenters[i ]
        x1 = phicenters[ip]
        if ip==0:
            x1+=2*np.pi
        y0 = rcenters[j]
        y1 = rcenters[jp]
        z0 = thetacenters[k]
        z1 = thetacenters[kp]

        xd = (phi-x0)/(x1-x0)
        yd = (r-y0)/(y1-y0)
        zd = (theta-z0)/(z1-z0)

        c00 = c000*(1-xd)+c100*xd
        c01 = c001*(1-xd)+c101*xd
        c10 = c010*(1-xd)+c110*xd
        c11 = c011*(1-xd)+c111*xd

        c0 = c00*(1-yd)+c10*yd
        c1 = c01*(1-yd)+c11*yd

        return c0*(1-zd)+c1*zd
    else:
        # 2d midplane interpolation
        c00 = arr[k,j ,i ]
        c10 = arr[k,j ,ip]
        c01 = arr[k,jp,i ]
        c11 = arr[k,jp,ip]

        x0 = phicenters[i ]
        x1 = phicenters[ip]
        if ip==0:
            x1+=2*np.pi
        y0 = rcenters[j]
        y1 = rcenters[jp]

        xd = (phi-x0)/(x1-x0)
        yd = (r-y0)/(y1-y0)

        c0 = c00*(1-xd)+c10*xd
        c1 = c01*(1-xd)+c11*xd

        return c0*(1-yd)+c1*yd


