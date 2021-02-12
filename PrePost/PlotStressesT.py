import matplotlib.pyplot as plt
from numpy import sqrt, arccos, pi, arctan
ZeroD = 1.e-9

def Plot0(P0,P1,P2,P3,fileName,sym,col,fy):
    ff = open(fileName,'r')
    fx = 1.
    tt, s1,s2, sigxx,sigyy,sigxy, phi,phiEps = [0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]
    Cr1Tn1, Cr1Ts1, Cr1Tn2, Cr1Ts2, Cr2Tn1, Cr2Ts1, Cr2Tn2, Cr2Ts2 = [0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]
    Cr1Wn1, Cr1Ws1, Cr1Wn2, Cr1Ws2, Cr2Wn1, Cr2Ws1, Cr2Wn2, Cr2Ws2 = [0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]
    z1 = ff.readline()
    z3 = z1.strip()
    z2 = z3.split()
    while z1<>"":
#        print z2
        tt += [fx*float(z2[0])]
        t = fx*float(z2[0])
        exx,eyy,exy, sxx,syy,sxy = float(z2[3]),float(z2[4]),float(z2[6]), float(z2[7]),float(z2[8]),float(z2[10]) 
        # stresses
        sigxx += [sxx]
        sigyy += [syy]
        sigxy += [sxy]
        # principal stresses
        a = 0.5*(sxx-syy)
        b = sqrt(a**2+sxy**2)
        if b>ZeroD:
            c = 0.5*(sxx+syy)
            s1 += [c + b]
            s2 += [c - b]
            if sxy<0.: si = -1.
            else:      si = 1.
            phi_ = 0.5*si*arccos(a/b)
        else:
            s1   = [sxx]
            s2   = [syy]
            phi_ = 0.
        if s1<s2: phi_ = phi_+0.5*pi
        phi += [phi_/pi*180.]
        # principal strain direction
        a = 0.5*(exx-eyy)
        b = sqrt(a**2+exy**2)
        if b>ZeroD:
            c = 0.5*(sxx+syy)
            e1 = c + b
            e2 = c - b
            if exy<0.: si = -1.
            else:      si = 1.
            phi_ = 0.5*si*arccos(a/b)
        else: phi_ = 0.
        if e1<e2: phi_ = phi_+0.5*pi
        phiEps += [phi_/pi*180.]
        # crack tractions
        if len(z2)>13:
            Cr1Ws1 += [float(z2[11])]
            Cr1Wn1 += [float(z2[12])]
            Cr1Ws2 += [float(z2[13])]
            Cr1Wn2 += [float(z2[14])]
            Cr1Ts1 += [float(z2[15])]
            Cr1Tn1 += [float(z2[16])]
            Cr1Ts2 += [float(z2[17])]
            Cr1Tn2 += [float(z2[18])]
            Cr2Ws1 += [float(z2[19])]
            Cr2Wn1 += [float(z2[20])]
            Cr2Ws2 += [float(z2[21])]
            Cr2Wn2 += [float(z2[22])]
            Cr2Ts1 += [float(z2[23])]
            Cr2Tn1 += [float(z2[24])]
            Cr2Ts2 += [float(z2[25])]
            Cr2Tn2 += [float(z2[26])]
        else:
            Cr1Ws1 += [0.]
            Cr1Wn1 += [0.]
            Cr1Ws2 += [0.]
            Cr1Wn2 += [0.]
            Cr1Ts1 += [0.]
            Cr1Tn1 += [0.]
            Cr1Ts2 += [0.]
            Cr1Tn2 += [0.]
            Cr2Ws1 += [0.]
            Cr2Wn1 += [0.]
            Cr2Ws2 += [0.]
            Cr2Wn2 += [0.]
            Cr2Ts1 += [0.]
            Cr2Tn1 += [0.]
            Cr2Ts2 += [0.]
            Cr2Tn2 += [0.]
        # next line
        z1 = ff.readline()
        z3 = z1.strip()
        z2 = z3.split()
    P0.plot(tt,sigxx,sym,color=col,label='sxx') #,linewidth=4)
    P0.plot(tt,sigyy,sym,color='blue',label='syy') #,linewidth=4)
    P0.plot(tt,sigxy,sym,color='green',label='sxy') #,linewidth=4)
    P0.legend(loc='upper left')
    P1.plot(tt,phi,label='phi_s')
    P1.plot(tt,phiEps,color='red',label='phi_e')
    P1.legend(loc='lower left')
    P2.plot(tt,s1,label='s1')
    P2.plot(tt,s2,label='s2')
    P2.legend(loc='lower left')
    print len(tt), len(Cr1Tn1), len(Cr1Tn2)
    P3.plot(tt,Cr1Tn1,color='green',label='1tn1')
    P3.plot(tt,Cr1Tn2,color='seagreen',label='1tn2')
    P3.plot(tt,Cr2Tn1,color='firebrick',label='2tn1')
    P3.plot(tt,Cr2Tn2,color='red',label='2tn1')
    P3.legend(loc='upper left')
    P3t = P3.twinx()
    P3t.plot(tt,Cr1Wn1,color='blue',label='1wn1')
    P3t.plot(tt,Cr1Wn2,color='navy',label='1wn2')
    P3t.plot(tt,Cr2Wn1,color='magenta',label='2wn1')
    P3t.plot(tt,Cr2Wn2,color='mediumorchid',label='2wn2')
    P3t.legend(loc='lower right')
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.set_xlabel('loading time')
    p0.grid()
    return P0, p0 

def Assemble(Case):
    P3, p3 = DefPlot(Case+" stress components")
    P4, p4 = DefPlot(Case+" principal directions")
    P5, p5 = DefPlot(Case+" principal stresses")
    P6, p6 = DefPlot(Case+" crack traction/width")
    return P3,p3, P4,p4, P5,p5, P6,p6 

