import ctypes as ct
libpartrace = ct.CDLL("src/libpartrace.so")

class C_MESH(ct.Structure):
    _fields_ = [
        ("nx", ct.c_size_t),
        ("ny", ct.c_size_t),
        ("nz", ct.c_size_t),
        ("data", ct.POINTER(ct.c_double))
    ]

    def __init__(self, nx: int, ny: int, nz: int, fname: str=''):
        self.nx = ct.c_size_t(nx)
        self.ny = ct.c_size_t(ny)
        self.nz = ct.c_size_t(nz)
        self.size = nx*ny*nz
        self.data = (ct.c_double * self.size)()
        if fname:
            print("reading file")
            rescale: int = 0
            if 'gasdens' in fname:
                rescale = 1 # density scaling
            elif 'gasv' in fname:
                rescale = 2 # velocity scaling
            self.read_datfile(fname,rescale)
            print("Read the file")
    
    def read_datfile(self,fname,rescale=0):
        libpartrace.read_datfile(
            ct.byref(self),
            ct.c_char_p(str.encode(fname)),
            ct.c_int(rescale)
        )

    def get_idx(self, k: int, j: int, i: int) -> int:
        libpartrace.get_idx.restype = int
        idx = libpartrace.get_idx(
            ct.byref(self),
            ct.c_size_t(k),
            ct.c_size_t(j),
            ct.c_size_t(i)
        )
        return idx
    
    def get_data(self, k: int, j: int, i: int) -> float:
        idx = self.get_idx(k,j,i)
        return self.data[idx]

if __name__ == '__main__':
    mesh = C_MESH(2048,256,32,
        "/Users/ericvc/fargo/outputs/alpha3_mplan0/gasdens50.dat")
    print('created mesh')
    rho = mesh.get_data(0,0,0)
    print(rho)