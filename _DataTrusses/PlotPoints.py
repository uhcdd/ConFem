import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def Plot0(P0,fileName,col,ind0,ind1):
    ff = open(fileName,'r')
    fx = -1.0
    fy = -1.0
    xL, yL = [0.], [0.]
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        print 'ZZ', z2
        xL += [fx*(float(z2[ind0]))]
        yL += [fy*float(z2[ind1])]
        z1 = ff.readline()
        z2 = z1.split()
    P0.plot(xL,yL, color=col) #,linewidth=4)
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
P3, p3 = DefPlot("P_u")
#Plot0(p3,'E8-01c.timeout.txt','green')
Plot0(p3,'staebe3D2.timeout.txt','green', 0, 1)
Plot0(p3,'abaqus.rpt.txt','red', 6, 3)
P3.autofmt_xdate()
plt.show()
