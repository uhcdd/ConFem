import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan, mean
ZeroD = 1.e-9

# Lines linestyle '-' or 'solid', '--' or 'dashed', '-.' or 'dashdot', ':' or 'dotted'
def Plot0(P0,fileName,col,ind0,ind1,fx,fy,LineS, Lab):
    ff = open(fileName,'r')
    xL, yL = [0.], [0.]
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        ll = len(z2)-2
        sum = 0.
        for i in range(ll/2):
            sum += fy*(float(z2[2*i+1]))
        if (ll % 2) <> 0: raise NameError("input line inconsistent")
#        print 'ZZ', z2, sum
        xL += [fx*(float(z2[ind0]))]
        yL += [sum]
        z1 = ff.readline()
        z2 = z1.split()
    P0.plot(xL,yL, color=col, linestyle=LineS, label=Lab) #,linewidth=4)
    ff.close()

def Plot0S(fileName, p0, f1, f2, col):
    ff = open(fileName,'r')
    TT, uu, rB, pp = [], [0.], [0.], []
    z1 = ff.readline()
    z2 = z1.strip()
    z3 = z2.split(',')
    while z1<>"":
        uu_, rB_ = [], []
        for i in xrange(int(len(z3)/2)):
            uu_ += [f1*float(z3[2*i])]
            rB_ += [f2*float(z3[2*i+1])]
        uu += [mean(uu_)]
        rB += [sum(rB_)]
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split(',')
    print 'XXX', uu, rB
    p0.plot(uu,rB, color=col) #,linewidth=4)
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
