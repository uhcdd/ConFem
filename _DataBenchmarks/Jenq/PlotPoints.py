import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

method = 2

if method==2:
    def plo(fileName):
        ff = open(fileName,'r')
        fx = 1.0
        fy = 1.0
        xL = []
        yL = []
        z1 = ff.readline()
        z2 = z1.split()
        while z1<>"":
            xL += [fx*(float(z2[0]))]
            yL += [fy*float(z2[1])]
#            print(z2)
            z1 = ff.readline()
            z2 = z1.split()
#        plt.axis('equal')
        plt.plot(xL,yL)
#        plt.autofmt_xdate()
        ff.close()
    def plo1(fileName):
        ff = open(fileName,'r')
        ss = 0.08
        z1 = ff.readline()
        z2 = z1.split()
        while z1<>"":
            xL = [float(z2[1]),float(z2[1])+ss*float(z2[4])]
            yL = [float(z2[2]),float(z2[2])+ss*float(z2[5])]
            print xL,yL
            plt.plot(xL,yL,color='red')
            z1 = ff.readline()
            z2 = z1.split()
#        plt.autofmt_xdate()
        ff.close()
else:
    def plo(fileName,fileName2):
        ff = open(fileName,'r')
        f2 = open(fileName2,'r')
        fx = 100.0
        fy = 1.0
        xL, xL_, yL, yL_ = [], [], [], []
        z1 = ff.readline()
        z2 = z1.split()
        z1_= f2.readline()
        z2_= z1_.split()
        while z1<>"":
#            xL += [(float(z2[1])+fx*float(z2[4]))]
#            yL += [fy*float(z2[5])]
#            xL_ += [(float(z2_[1])+fx*float(z2_[4]))]
#            yL_ += [fy*float(z2_[5])]
            xL += [float(z2[0])]
            xL_+= [float(z2_[0])]
            yL += [fy*float(z2[1])]
            yL_+= [fy*float(z2_[1])]
            z1 = ff.readline()
            z2 = z1.split()
            z1_= f2.readline()
            z2_= z1_.split()
        plt.plot(xL,yL,linewidth=4)
        plt.plot(xL,yL_,linewidth=4)
        plt.xticks(fontsize='x-large')
        plt.yticks(fontsize='x-large')
        plt.grid()

def Plot1(fileName,P0,i0,iX,i1,i2,i3, P1,k1,k2,k3):
    ff = open(fileName,'r')
    z1 = ff.readline()
    z2 = z1.split()
    x_ = 0.
    xP, y1P, y2P, y3P, z1P, z1P_ = [], [], [], [], [], []
    while z1<>"":
        x, xx, y1, y2, y3 = float(z2[i0]), float(z2[iX]), float(z2[i1]), float(z2[i2]), float(z2[i3])
        if y2<-0.01: x=xx    # quick and dirty
        if P1<>None: xx, yy, xy = float(z2[k1]), float(z2[k2]), float(z2[k3])
        if x <> x_:              # read same values of x before, last one relevant
            xP += [x_]    
            y1P += [y1_]    
            y2P += [y2_]    
            y3P += [y3_]
            if P1<>None:
                dd, phi, phi_ = xx_-yy_, 0., 0.
                if abs(dd)>ZeroD:
                    c2phi = 0.5*(xx_-yy_)/sqrt((0.5*(xx_-yy_))**2+xy_**2)
                    phi = 0.5*arccos(c2phi)*180./pi
                    phi_= 0.5*arctan(2.*xy_/(xx_-yy_))*180./pi
                z1P += [phi]
                z1P_+= [phi_]
        x_, y1_, y2_, y3_ = x, y1, y2, y3
        if P1<>None: xx_, yy_, xy_ = xx, yy, xy
        z1 = ff.readline()
        z2 = z1.split()
    ff.close()
    P0.plot(xP,y1P)
    P0.plot(xP,y2P)
    P0.plot(xP,y3P)
    if P1<>None:
        P1.plot(xP,z1P)
        P1.plot(xP,z1P_)

#plo('UniaxTens.timeout.txt')
#plo('UniaxTens.timeout_.txt')
#plo('UniaxTens.timeout.k15.txt')
#Plot1('WillamsTest.stuff.txt',0,6,7,8)
#exit

def Plot0(P0,fileName):
    ff = open(fileName,'r')
    fx = -1.0
    fy = -1.0
    xL, yL = [], []
    z1 = ff.readline()
    z2 = z1.split()
    while z1<>"":
        xL += [fx*float(z2[0])]
        yL += [fy*float(z2[1])]
        z1 = ff.readline()
        z2 = z1.split()
    P0.plot(xL,yL) #,linewidth=4)
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
#    P0.autofmt_xdate()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
#P0 = plt.figure()
#P0.autofmt_xdate()
#p0 = P0.add_subplot(111)
#p0.set_title('strains',fontsize='x-large')
#p0.tick_params(axis='x', labelsize='large')
#p0.tick_params(axis='y', labelsize='large')
#p0.grid()
#Plot1('WillamsTest.stuff.txt',p0,0,1,2,4, p1, 0,1,4)
#p1.set_title('principal strain orientation',fontsize='x-large')
#p3.set_title('principal stress direction',fontsize='x-large')

P3, p3 = DefPlot("L_D")
Plot0(p3,'Jenq-L3-12122007-Final.timeout.txt')
Plot0(p3,'Jenq-L3-12122007-Final.timeout - Kopie.txt')

P3.autofmt_xdate()

#p2 = DefPlot("stresses")
#Plot1('TenCom2P.stuff.txt',p2,1,2,6,7,8, None, 5,6,7)
#Plot1('TenCom2P.stuff_.txt',p2,1,2,6,7,8, None, 5,6,7)
#Plot1('WillamsTest.stuff_.txt',p2,0,5,6,7, p3, 6,7,8)
plt.show()
