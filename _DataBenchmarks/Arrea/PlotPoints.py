import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def Plot0(P0,P1,fileName):
    ff = open(fileName,'r')
    fx = -1.0e5
    fy = -1.0e3
    xL1, xL2, yL = [], [], []
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        xL1 += [fx*(float(z2[4])-float(z2[2]))]
        xL2 += [-fx*(float(z2[8])-float(z2[6]))] 
        yL  += [fy*float(z2[1])]
        z1 = ff.readline()
        z2 = z1.split()
    P0.plot(xL1,yL, color="magenta") #,linewidth=4)
    P1.plot(xL2,yL, color="magenta")
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
P3, p3 = DefPlot("P - CMSD")
P4, p4 = DefPlot("P - CMOD")
Plot0(p3,p4,'arrea2D.timeout.txt')
P3.autofmt_xdate()
plt.show()
