import ctypes as ct
libpartrace = ct.CDLL("src/libpartrace.so")

class C_DOMAIN(ct.Structure):
    _fields_ = [
        ("cfargodir", ct.POINTER(ct.c_char)),
        ("nx", ct.c_size_t),
        ("ny", ct.c_size_t),
        ("nz", ct.c_size_t),
        ("phiCenters",   ct.POINTER(ct.c_double)),
        ("rCenters",     ct.POINTER(ct.c_double)),
        ("thetaCenters", ct.POINTER(ct.c_double)),
        ("phiEdges",     ct.POINTER(ct.c_double)),
        ("rEdges",       ct.POINTER(ct.c_double)),
        ("thetaEdges",   ct.POINTER(ct.c_double))
    ]

    def __init__(self, fargodir: str, nx: int, ny: int, nz: int) -> None:
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.fargodir = fargodir
        self.cfargodir = ct.create_string_buffer(str.encode(fargodir),100)

        self.phiCenters = (ct.c_double * nx)()
        self.phiEdges   = (ct.c_double * (nx+1))()
        self.rCenters = (ct.c_double * ny)()
        self.rEdges   = (ct.c_double * (ny+1))()
        self.thetaCenters = (ct.c_double * nz)()
        self.thetaEdges   = (ct.c_double * (nz+1))()

        xfile = fargodir+'/domain_x.dat'
        yfile = fargodir+'/domain_y.dat'
        zfile = fargodir+'/domain_z.dat'

        print('Reading x')
        libpartrace.read_domfile(
            self.phiEdges,
            self.phiCenters,
            ct.c_size_t(nx),
            ct.c_char_p(str.encode(xfile)),
            ct.c_int(0),ct.c_int(0)
            )
        print('Done')
        print(self.phiCenters[0])
        print('Reading y')
        libpartrace.read_domfile(
            self.rEdges,
            self.rCenters,
            ny,ct.c_char_p(str.encode(yfile)),3,1
            )
        print('Reading z')
        libpartrace.read_domfile(
            self.thetaEdges,
            self.thetaCenters,
            nz,ct.c_char_p(str.encode(zfile)),3,0
            )


if __name__ == "__main__":
    domain = C_DOMAIN("/Users/ericvc/fargo/outputs/alpha3_mplan0",2048,256,32)
    print(domain.fargodir)
    print(domain.rCenters[:32])