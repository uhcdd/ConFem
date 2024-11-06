from numpy import array
from numpy.linalg import eigh, norm #, det #, eigvalsh

dEps = [1.5200e-03, -5.1549e-04, -5.7045e-04,  0.0000e+00,  0.0000e+00,  4.0000e-05]
dEps = [1.4177e-03, -5.2116e-04, -5.4405e-04,  0.0000e+00,  0.0000e+00,  5.8410e-06]
dEps = [3.7766,   0.9582,   0.0230,   0.0000,   0.0000,   0.4982 ]

dEps = [2.4985e-03,  9.5226e-04, -2.3924e-03,  0.0000e+00,  0.0000e+00,  1.9970e-03]
dEps = [2.4593e-03,  8.9584e-04, -2.3684e-03,  0.0000e+00,  0.0000e+00,  1.9520e-03]
dEps = [1.8395,   2.3416,  -0.0051,   0.0000,   0.0000,   0.6569]
xx = array([ [dEps[0], dEps[5], dEps[4]] , [dEps[5], dEps[1], dEps[3]] , [dEps[4], dEps[3], dEps[2]] ])
PdEps, gg_ = eigh( xx)
print xx
print PdEps
print gg_
