import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan,mean
ZeroD = 1.e-9

def Plot0(fileName, p0, col):
    ff = open(fileName,'r')
    TT, uu, rB, pp = [], [0.], [0.], []
    z1 = ff.readline()
    z2 = z1.strip()
    z3 = z2.split()
    while z1<>"":
        uu_, rB_ = [], []
        uu += [float(z3[0])]
        rB += [float(z3[1])]
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split()
    p0.plot(uu,rB, color=col) #,linewidth=4)
#    p0.plot(xx,xL, color='blue') #,linewidth=4)
#    p0.plot([xx[0],xx[-1]],[mean(xL),mean(xL)],color='black')
#    p0.set_ylim(0,1.05*max(xL))
#    p1.plot(xx,rx_bot, color='red') #,linewidth=4)
#    p1.plot(xx,rx_top, color='blue') #,linewidth=4)
#    p1.plot([xx[0],xx[-1]],[mean(rx_top),mean(rx_top)],color='black')
    ff.close()
#    print fileName, counter, sum(xL), sum(yL), sum(rx_bot), sum(rx_top)

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
P3, p3 = DefPlot("R_y")
Plot0('E8-01du02.timeout - Kopie (4).txt', p3, 'indianred') # elastic lt
Plot0('E8-01du02.timeout - Kopie (3).txt', p3, 'magenta')  # isaodamage
Plot0('E8-01du02.timeout - Kopie.txt', p3, 'red')      # microlane 
Plot0('E8-01du02.timeout.txt', p3, 'blue')
P3.autofmt_xdate()
plt.show()
