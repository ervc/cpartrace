import numpy as np

MH = 1.66053907e-24
BK = 1.38e-16 # Boltzmann constant [cgs]
MSUN = 1.9891e33
MEARTH = 5.97e27
AU = 1.49597871e13
G = 6.674e-8
YR = 3.1557600e7
RJUP = 6.995e9

R0 = ( 75*AU )

LEN = ( R0 )
MASS = ( MSUN )
TIME = ( np.sqrt(R0*R0*R0/G/MSUN) )