import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def Plot0(P0,P1,fileName):
    ff = open(fileName,'r')
    fx = 1.0
    fy = 1.0
    xL, yL1, yL2, yL3, yL4, yL5, yL6, yL7 = [], [], [], [], [], [], [], []
    yM1, yM2, yM3 = [], [], []
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        xL += [fx*float(z2[0])]
        yL1 += [fy*float(z2[2])]
        yL2 += [fy*float(z2[4])]
        yL3 += [fy*float(z2[6])]
        yL4 += [fy*float(z2[8])]
        yL5 += [fy*float(z2[10])]
        yL6 += [fy*float(z2[12])]
        yL7 += [fy*float(z2[14])]
        y1 = fy*(float(z2[4]) +float(z2[6]) +float(z2[8]))/3.
        y2 = fy*(float(z2[10])+float(z2[12])+float(z2[14]))/3.
        yM1 += [y1]
        yM2 += [y2]
        yM3 += [y1-y2]
        z1 = ff.readline()
        z2 = z1.split()
#    P0.plot(xL,yL1, color="magenta") #,linewidth=4)
    P0.plot(xL,yL2, color="red") #,linewidth=4)
    P0.plot(xL,yL3, color="blue") #,linewidth=4)
    P0.plot(xL,yL4, color="green") #,linewidth=4)
    P0.plot(xL,yL5,'--', color="red") #,linewidth=4)
    P0.plot(xL,yL6,'--', color="blue") #,linewidth=4)
    P0.plot(xL,yL7,'--', color="green") #,linewidth=4)
    P0.plot(xL,yM1,'o', color="red") #,linewidth=4)
    P0.plot(xL,yM2,'o', color="blue") #,linewidth=4)
    P0.plot(xL,yM3,'o', color="green") #,linewidth=4)
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
P3, p3 = DefPlot("L_D")
P4, p4 = DefPlot("P_u")
Plot0(p3,p4,'Nooru-1550-5.timeout.txt')
P3.autofmt_xdate()
plt.show()
