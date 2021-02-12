import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def Plot0(P0,fileName,col):
    ff = open(fileName,'r')
    fx = 1.0
    fy = 1.0
    xL, yL = [0.], [0.]
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        xL += [fx*(float(z2[0]))]
        yL += [fy*float(z2[1])]
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
Plot0(p3,'E2-01.timeout.txt','green')
P3.autofmt_xdate()
plt.show()
