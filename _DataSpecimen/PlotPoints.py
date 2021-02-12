import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def PlotF( P0, ff, fx, fy, col):
    xL1, xL2, yL = [0.], [0.], [0.]
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        xL1 += [fx*float(z2[0])]
        yL  += [fy*float(z2[1])]
        z1 = ff.readline()
        z2 = z1.split()
    P0.plot(xL1,yL, color=col) #,linewidth=4)

def Plot0(P0,fileName, fileN2 ):
    LL = 50.
    fx = 1.0e3/LL # -1.0e3
    fy = 2.0 # -4.0
    ff = open(fileName,'r')
    PlotF( P0, ff, fx, fy, "magenta" )
    ff.close()
    if fileN2 <> None:
        ff = open(fileN2,'r')
        PlotF( P0, ff, fx, fy, "blue" )
        ff.close()
#    xL1, xL2, yL = [0.], [0.], [0.]
#    z1 = ff.readline()
#    z2 = z1.split()
#    while z1<>"":
#        xL1 += [fx*float(z2[0])]
#        yL  += [fy*float(z2[1])]
#        z1 = ff.readline()
#        z2 = z1.split()
#    ff.close()
#    P0.plot(xL1,yL, color="magenta") #,linewidth=4)

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
P3, p3 = DefPlot("X")
Plot0(p3,'UniaxTensionBar2D.timeout.txt', 'UniaxTensionBar2D.timeout - Kopie.txt')
P3.autofmt_xdate()
plt.show()
