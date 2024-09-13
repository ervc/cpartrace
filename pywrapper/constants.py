import numpy as np

MSUN = 1.9891e33
MEARTH = 5.97e27
AU = 1.49597871e13
G = 6.674e-8
YR = 3.1557600e7
RJUP = 6.995e9

R0 = ( 5.2*AU )

LEN = ( R0 )
MASS = ( MSUN )
TIME = ( np.sqrt(R0*R0*R0/G/MSUN) )