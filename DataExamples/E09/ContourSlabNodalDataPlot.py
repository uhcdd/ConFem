"""
=======================================
Contour plot of irregularly spaced data
=======================================

Comparison of a contour plot of irregularly spaced data interpolated
on a regular grid versus a tricontour plot for an unstructured triangular grid.

Since `~.axes.Axes.contour` and `~.axes.Axes.contourf` expect the data to live
on a regular grid, plotting a contour plot of irregularly spaced data requires
different methods. The two options are:

* Interpolate the data to a regular grid first. This can be done with on-board
  means, e.g. via `~.tri.LinearTriInterpolator` or using external functionality
  e.g. via `scipy.interpolate.griddata`. Then plot the interpolated data with
  the usual `~.axes.Axes.contour`.
* Directly use `~.axes.Axes.tricontour` or `~.axes.Axes.tricontourf` which will
  perform a triangulation internally.

This example shows both methods in action.
"""

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

def ReadData(fileName):
    ff = open(fileName,'r')
    x, y, z, R = [], [], [], []
    z1 = ff.readline()
    z1 = ff.readline()
    z2 = z1.strip()
    z3 = z2.split()
    while z1!="":
        x += [float(z3[1])]   # x coordinate
        y += [float(z3[2])]   # y coordinate
        z += [float(z3[4])]   # deflection
        R += [float(z3[10])]  # reaction forces
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split()
    return x,y,z, R

fileName = "E7-04.nodeout.txt"
x, y, z, R = ReadData(fileName)
P0 = plt.figure()
p0 = P0.add_subplot(111)
P1 = plt.figure()
p1 = P1.add_subplot(111)

# ----------
# Tricontour
# ----------
# Directly supply the unordered, irregularly spaced coordinates
# to tricontour.

p0.tricontour(x, y, z, linewidths=0.5, colors='k')
cntr2 = p0.tricontourf(x, y, z, cmap="bone")
P0.colorbar(cntr2, ax=p0)
p0.axis('equal')
p0.plot(x, y, 'ko', ms=1)

p1.axis('equal')
area = [5.0e4*r**2 for r in R]
plt.scatter(x, y, s=area, c='lightgrey')
for x_, y_, R_ in zip(x, y, R):
    if abs(R_)>1.0e-3: p1.text(x_,y_,"%5.3f"%(R_),
                               color=('green' if R_>0. else 'red'),
                               fontsize='large')
print('total reactions ',sum(R))
plt.show()
