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
        tt += [fx*float(z2[0])]
        t = fx*float(z2[0])
        exx,eyy,exy, sxx,syy,sxy = float(z2[3]),float(z2[4]),float(z2[6]), float(z2[7]),float(z2[8]),float(z2[10]) 
        # stresses
        sigxx += [sxx]
        sigyy += [syy]
        sigxy += [sxy]
        # principal stresses
        ds = 0.5*(sxx-syy)
        ss  = 0.5*(sxx+syy)
        if abs(ds)>ZeroD and abs(sxy)>ZeroD:
            deno = sqrt(ds**2+sxy**2)
            c2p = ds/deno
            s1 += [ss + deno]
            s2 += [ss - deno]
        else:
            c2p = 0.
            s1 += [sxx]
            s2 += [syy]
        phi += [0.5*arccos(c2p)/pi*180.]
        # principal strain
        ds = 0.5*(exx-eyy)
        if abs(ds)>ZeroD and abs(exy)>ZeroD:
            c2p = ds/sqrt(ds**2+(0.5*exy)**2)
        else:
            c2p = 0.
        phiEps += [0.5*arccos(c2p)/pi*180.]
        # crack tractions
        if len(z2)>11:
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
    P0.plot(tt,sigxx,sym,color=col) #,linewidth=4)
    P0.plot(tt,sigyy,sym,color='blue') #,linewidth=4)
    P0.plot(tt,sigxy,sym,color='green') #,linewidth=4)
    P1.plot(tt,phi)
    P1.plot(tt,phiEps,color='red')
    P2.plot(tt,s1)
    P2.plot(tt,s2)
    print len(tt), len(Cr1Tn1), len(Cr1Tn2)
    P3.plot(tt,Cr1Tn1)
    P3.plot(tt,Cr1Tn2)
    P3.plot(tt,Cr2Tn1)
    P3.plot(tt,Cr2Tn2)
    P3t = P3.twinx()
    P3t.plot(tt,Cr1Wn1,color='red')
    P3t.plot(tt,Cr1Wn2,color='red')
    P3t.plot(tt,Cr2Wn1,color='violet')
    P3t.plot(tt,Cr2Wn2,color='violet')
    ff.close()

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 
    
Case = ' '
P3, p3 = DefPlot(Case+" stress components")
P4, p4 = DefPlot(Case+" principal stress directions")
P5, p5 = DefPlot(Case+" principal stresses")
P6, p6 = DefPlot(Case+" crack tractions")

Plot0(p3,p4,p5,p6,'LSP-2079-7_6_DAMMK.elemout.6_0.txt','','red',200.)
P3.autofmt_xdate()
plt.show()
