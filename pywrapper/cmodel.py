import ctypes as ct
libpartrace = ct.CDLL("src/libpartrace.so")

from cmesh import C_MESH
from cdomain import C_DOMAIN

class C_MODEL(ct.Structure):
    _fields_ = [
        ("fargodir", ct.POINTER(ct.c_char)),
        ("nout", ct.POINTER(ct.c_char)),
        ("nx", ct.c_size_t),
        ("ny", ct.c_size_t),
        ("nz", ct.c_size_t),
        ("alpha", ct.c_double),
        ("aspect", ct.c_double),
        ("flaring", ct.c_double),
        ("planetmass", ct.c_double),
        ("omegaframe", ct.c_double),
        ("planetpos", ct.c_double * 3),
        ("sunpos", ct.c_double * 3),
        ("planetEnvelope", ct.c_double),
        ("gasdens", ct.POINTER(C_MESH)),
        ("domain", ct.POINTER(C_DOMAIN)),
        ("gasvx", ct.POINTER(C_MESH)),
        ("gasvy", ct.POINTER(C_MESH)),
        ("gasvz", ct.POINTER(C_MESH)),
        ("drhodx", ct.POINTER(C_MESH)),
        ("drhody", ct.POINTER(C_MESH)),
        ("drhodz", ct.POINTER(C_MESH)),
    ]

    def __init__(self, fargodir: str, nout: str|int, nx: int, ny: int, nz: int) -> None:
        self.fargodir = ct.create_string_buffer(fargodir.encode(),100)
        self.nout = ct.create_string_buffer(str(nout).encode(),5)
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # read in data
        rhofile = fargodir+f"/gasdens{nout}.dat"
        vphifile = fargodir+f"/gasdens{nout}.dat"
        vrfile = fargodir+f"/gasdens{nout}.dat"
        vthetafile = fargodir+f"/gasdens{nout}.dat"
        varfile = fargodir+f"/variables.par"

        self.gasdens = self.init_MeshField_fromFile(rhofile,1)
        gasvphi     = self.init_MeshField_fromFile(vphifile,2)
        gasvr       = self.init_MeshField_fromFile(vrfile,2)
        gasvtheta   = self.init_MeshField_fromFile(vthetafile,2)

        self.domain = self.init_Domain()

        print("Making cartvels...")
        self.make_cartvels(gasvphi,gasvr,gasvtheta)

        print("Getting Gradrho...")
        self.init_gradrho()

    def __repr__(self) -> str:
        return f"{self.fargodir = }\n{self.nout = }"


    def init_MeshField_fromFile(self, file: str, rescale: int):
        libpartrace.init_MeshField_fromFile.restype = ct.POINTER(C_MESH)
        return libpartrace.init_MeshField_fromFile(
            ct.c_char_p(str.encode(file)),
            ct.c_size_t(self.nx), ct.c_size_t(self.ny), ct.c_size_t(self.nz),
            ct.c_int(rescale)
        )
    
    def init_Domain(self):
        libpartrace.init_Domain.restype = ct.POINTER(C_DOMAIN)
        return libpartrace.init_Domain(
            self.fargodir,
            ct.c_size_t(self.nx), ct.c_size_t(self.ny), ct.c_size_t(self.nz)
        )
    
    def make_cartvels(self, gasvphi, gasvr, gasvtheta):
        libpartrace.make_cartvels(
            ct.byref(self),
            gasvphi, gasvr, gasvtheta
        )

    def init_gradrho(self):
        libpartrace.init_gradrho(
            ct.byref(self)
        )


if __name__ == '__main__':
    model = C_MODEL("/Users/ericvc/fargo/outputs/alpha3_mplan0","50",2048,256,32)
    print(model)