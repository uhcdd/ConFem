"""
=======================================
Contour plot of irregularly spaced data
=======================================

"""

import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.ndimage.filters import gaussian_filter

def ReadData(fileName):
    def rL():
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split()
        return z3
        
    ff = open(fileName,'r')
    x, y, axB, ayB, axT, ayT = [], [], [], [], [], []
    z3 = rL()
    while len(z3):
        if z3[0] in ['**','\n']:
            pass
        else:
            x   += [float(z3[0])]   # x coordinate
            y   += [float(z3[1])]   # y coordinate
            axB += [float(z3[3])]   # reinf x Bottom
            ayB += [float(z3[4])]   # reinf y Bottom
            axT += [float(z3[6])]   # reinf x Bottom
            ayT += [float(z3[7])]   # reinf y Bottom
        z3 = rL()
    return x,y,axB, ayB, axT, ayT

fileName = "E7-01.reinforcement.txt"
x, y, axB, ayB, axT, ayT = ReadData(fileName)
asB = [ax+ay for  ax, ay in zip(axB, ayB)]
asT = [ax+ay for  ax, ay in zip(axT, ayT)]
# smooth whole stuff
triang = tri.Triangulation(x, y)
refiner = tri.UniformTriRefiner(triang)
tri_refi, asB_refi = refiner.refine_field(asB, subdiv=3)

P0 = plt.figure()
p0 = P0.add_subplot(111)
if True:
    p0.tricontour(x, y, asB, linewidths=0.5, colors='k')
    cntr2 = p0.tricontourf(x, y, asB, cmap="Purples")
else:
    p0.tricontour( tri_refi, asB_refi, linewidths=0.5, colors='k')
    cntr2 = p0.tricontourf( tri_refi, asB_refi, cmap="Purples")
P0.colorbar(cntr2, ax=p0)
p0.axis('equal')
p0.plot(x, y, 'ko', ms=1)
#
P1 = plt.figure()
p1 = P1.add_subplot(111)
p1.tricontour(x, y, asT, linewidths=0.5, colors='k')
cntr2 = p1.tricontourf(x, y, asT, cmap="Purples")
P1.colorbar(cntr2, ax=p1)
p1.axis('equal')
p1.plot(x, y, 'ko', ms=1)
#
plt.show()
