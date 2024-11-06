# ConPlaD -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from math import *
from time import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from os import path as pth
from numpy import degrees

import ConFemBasics
from ConFemBasics import *
from ConFemMat import *
import ConFemElem
import ConFemSteps
from ConFemInOut import *
from ConLinAlg import *

def Reinforcement2D( ElList, NodeList,NoIndToCMInd, ReinfDes, Sca, KeyElset,KeyMat, ff, PlotFlag, Type, Time):
    # data from binary database
    def MSIG(pp):
        if pp[0]<0:
#            phi = atan(pp[2]/pp[1])                                         # direction of principal compressive stress
            if abs(pp[1])>ZeroD: phi = atan(pp[2]/pp[1])
            else:                phi = -(0.5*pi-atan(pp[1]/pp[2]))
            sig = pp[0]                                                     # value of principal compressive stress
        elif pp[3]<0:
            if abs(pp[4])>ZeroD: phi = atan(pp[5]/pp[4])
            else:                phi = -(0.5*pi-atan(pp[4]/pp[5]))
            sig = pp[3]
        else:
            phi = None
            sig = None
        return phi, sig
    def KIXC( sigxy, rhox, xx):                                             # NR iteration matrix for reinforcement in x-direction sigsy prescribed, sigsx unknown 
        sigc=xx[0]
        phi = xx[2]
        KI = array([[-2*sigc**2*cos(phi)/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2),0,
(-sin(phi)**2+cos(phi)**2)/sin(phi)/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2)*sigc],
[2*cos(phi)*sigc**2*(cos(phi)**2+sin(phi)**2)/rhox/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2),1/rhox,
 -cos(phi)*(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc)/rhox/sin(phi)/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2)],
[sin(phi)/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2)*sigc,0,
 -sigxy/sigc/sin(phi)/(-2*sigxy*cos(phi)-sin(phi)**3*sigc+sin(phi)*sigc*cos(phi)**2)]])
        return KI
    def KIYC( sigxy, rhoy, xx):                                             # NR iteration matrix for reinforcement in x-direction sigsx prescribed, sigsy unknown 
        sigc=xx[0]
        phi = xx[2]
        KI = array([[2*sigc**2*sin(phi)/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc),
-(sin(phi)**2-cos(phi)**2)*sigc/cos(phi)/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc),0],
 [-2*sigc**2*sin(phi)*(cos(phi)**2+sin(phi)**2)/rhoy/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc),
  sin(phi)*(2*sigxy*cos(phi)+sin(phi)**3*sigc-sin(phi)*sigc*cos(phi)**2)/cos(phi)/rhoy/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc),1/rhoy],
 [cos(phi)/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc)*sigc,
  -sigxy/sigc/cos(phi)/(2*sigxy*sin(phi)-cos(phi)*sigc*sin(phi)**2+cos(phi)**3*sigc),0]])
        return KI
    def SLABRE(mx, my, mxy, fy, fcd, dd, ff):                               # bending reinforcement for slabs
        Tol = 1.e-6
        pp = PrinC( mx, my, mxy)                                            # principal stresses
        if pp[0]<=0 and pp[3]<=0: return 0,0,0,0, [0,0,0,0,0,0],0
        chi = 0.95
        kk = 0.8
        fc_ = kk*chi*fc
        if mxy>0: phic=-pi/4                                                # prescription of compressive principal stress direction
        else:     phic= pi/4
        mc = -2*abs(mxy)                                                    # default value for concrete moment
        indi = -2*mc/fc_/dd**2                                            # indicator for concrete demand
        if indi>0.99: NameError ("Reinforcement2D: slab concrete demand")
        zz = 0.5*dd * ( 1 + sqrt(1-indi) )                                  # internal lever arm
        as1 = (mx -mc*cos(phic)**2)/(zz*fy)                                 # reinforcement x - direction
        as2 = (my -mc*sin(phic)**2)/(zz*fy)                                 # reinforcement y - direction
        if as1<0 and as2<0: raise NameError ("Reinforcement2D: unexpected condition")
        if as1<0 or as2<0:                                                  # not a reasonable reinforcement value, iteration for concrete strut direction necessary
            if abs(mxy)<0.01:                                               # improve iteration starter
                phic = atan(pp[5]/pp[4])
                mc = pp[3]
            if as1<0: xx = array([mc,as2,zz,phic])                          # iteration initial value
            else:     xx = array([mc,as1,zz,phic])
            N1 = 21
            for k in range(N1):                                             # Newton Raphson iteration 
                indi = -2*xx[0]/fc_/dd**2                                   # indicator for concrete demand
                if indi>0.99: NameError ("Reinforcement2D: slab concrete demand")
                if as1<0:
                    RR = array([cos(xx[3])*sin(xx[3])*xx[0]-mxy,
                                -mx+xx[0]*cos(xx[3])**2,
                                fy*xx[2]*xx[1]-my+xx[0]*sin(xx[3])**2,
                                xx[2]-0.5*dd*(1+sqrt(1-indi))])
                    KM = array([[cos(xx[3])*sin(xx[3]),0,0,-xx[0]*sin(xx[3])**2+xx[0]*cos(xx[3])**2],
                                [cos(xx[3])**2,0,0,-2*cos(xx[3])*sin(xx[3])*xx[0]],
                                [sin(xx[3])**2,fy*xx[2], fy*xx[1],2*cos(xx[3])*sin(xx[3])*xx[0]],
                                [-1/(dd*sqrt(1-indi)*fc_),0,1,0]])
                else:
                    RR = array([cos(xx[3])*sin(xx[3])*xx[0]-mxy,
                                fy*xx[2]*xx[1]-mx+xx[0]*cos(xx[3])**2,
                                -my+xx[0]*sin(xx[3])**2,
                                xx[2]-0.5*dd*(1+sqrt(1-indi))])
                    KM = array([[cos(xx[3])*sin(xx[3]),0,0,-xx[0]*sin(xx[3])**2+xx[0]*cos(xx[3])**2],
                                [cos(xx[3])**2,fy*xx[2],fy*xx[1],-2*cos(xx[3])*sin(xx[3])*xx[0]],
                                [sin(xx[3])**2,0,fy*as2,2*cos(xx[3])*sin(xx[3])*xx[0]],
                                [-1/(dd*sqrt(1-indi)*fc_),0,1,0]]) 
                if norm(RR)<Tol: break                # found solution
                xn = xx - solve(KM,RR) #,overwrite_a=True,overwrite_b=True)
                xx[:] = xn[:]
            if k<N1-1:
                mc = xx[0]                                                  # concrete moment
                if as1<0: 
                    as1 = 0
                    as2 = xx[1]                                             # reinforcement
                else:
                    as1 = xx[1]
                    as2 = 0
                zz = xx[2]                                                  # internal lever arm
                phic = xx[3]                                                # concrete moment direction
            else: raise NameError ("Reinforcement2D: no convergence") # no convergence reached
#            print('XXX', k, mx, my, mxy,'__',mc, as1, as2, zz, phic,norm(RR), file=ff)
        x = 2*(dd-zz)/kk                                # compression zone height
        if (fcd*chi*kk*x + mc/zz) > 100*Tol: raise NameError ("Reinforcement2D: moment control",(fcd*chi*kk*x + mc/zz)/(-mc/zz) )
        return phic, zz, as1, as2, pp, x
    def SLABSH(mx, my, mxy, vx, vy, as1B, as2B, as1T, as2T, zz, fy, fcd, dd): # shear calc for slabs
        pp = PrinC( mx, my, mxy)                                            # principal stresses
        phi1 = atan(pp[2]/pp[1])                                            # pp1, pp2 indicate x, y-direction of larger principal moment m1 which is pp[0] --> phi1 direction angle of m1
        phi2 = phi1+0.5*pi
        vphi1 = vx*cos(phi1)+vy*sin(phi1)
        vphi2 = vx*cos(phi2)+vy*sin(phi2)
#        print('XXX',degrees(phi1),'__',vphi1,vphi2,'__',vx,vy) for control purposes
        return phi1, vphi1, vphi2
    def RHOITER( stri, sig, phi, N1):
        xx = array([sig,0,phi])                                             # iteration initial value
        rho = 1
        for k in range(N1):                                                 # Newton Raphson iteration for sigsy --> xx[1]
            if   stri=='rhoy': RR =array([cos(xx[2])*sin(xx[2])-sigxy/xx[0],     -sigx+xx[0]*cos(xx[2])**2,xx[1]-sigy+xx[0]*sin(xx[2])**2])
            elif stri=='rhox': RR =array([cos(xx[2])*sin(xx[2])-sigxy/xx[0],xx[1]-sigx+xx[0]*cos(xx[2])**2,     -sigy+xx[0]*sin(xx[2])**2])
            if norm(RR)<1.e-4: break    # found solution
            if   stri=='rhoy': KI = KIYC( sigxy, rho, xx)
            elif stri=='rhox': KI = KIXC( sigxy, rho, xx)
            xn = xx - dot(KI,RR)
            xx = xn.copy() #copy(xn)
        if k<N1-1:
            sigc = xx[0]
            rho = xx[1]/fy
            phic = xx[2]
        else: raise NameError ("Reinforcement2D: no convergence") # no convergence reached
        return sigc, rho, phic, k
    # end of defs

    fy =  ReinfDes[0]                                                       # admissible reinforcement tensile stress
    fc = -ReinfDes[1]                                                       # admissible concrete compressive stress
    minRhoX = ReinfDes[2]/100.                                              # minimum reinforcement ratio absolute from percentage   
    minRhoY = ReinfDes[3]/100.                                              #                         / structural height
    hh      = ReinfDes[4]                                                   # geometric height
    dd      = ReinfDes[5]                                                   # structural height
    specWeSt= ReinfDes[6]                                                   # specific weight reinforcing steel
#    N1 = 10                                                                 # max number of iterations for NR
    ASx, ASy, Vol = 0, 0, 0                                                 # sum reinforcement x direction, sum reinforcement y direction, concrete volume 
    if  Type.upper()=="PLATE":  
        ff.write('             x       y     sig_x   sig_y  sig_xy     phi_c  sig_c     rho_x  rho_y      as_x    as_y\n')
        PlotTitles = [KeyElset+" plate reinforcement time "+'%.2f'%Time,KeyElset+" plate concrete struts time "+'%.2f'%Time,KeyElset+" plate reinforcement ties time "+'%.2f'%Time]
    elif Type.upper()=="SLAB": 
        ff.write('**                                  bottom                       top\n**     x       y     phicB    as_x    as_y     phicT    as_x    as_y     phi_1     v_1     v_2\n')
        PlotTitles = ["slab lower bending reinforcement time "+'%.2f'%Time,"slab upper bending reinforcement time "+'%.2f'%Time,"slab shear time "+'%.2f'%Time]
    else: raise NameError ("Reinforcement2D: unknown type")
#    P0 = plt.figure().add_subplot(111,title='plate reinforcement / slab lower bending reinforcement: ')
    P0 = plt.figure()
    p0 = P0.add_subplot(111,title=PlotTitles[0])
    p0.axis('equal')
    p0.grid()
    P1 = plt.figure()
    p1 = P1.add_subplot(111,title=PlotTitles[1])
    p1.axis('equal')
    p1.grid()
    P2 = plt.figure()
    p2 = P2.add_subplot(111,title=PlotTitles[2])
    p2.axis('equal')
    p2.grid()

    NoPoint, minSigC, maxRhoX, maxRhoY, PointmSC,PointmRX,PointmRY, xIL, yIL, vL = 0, 0., 0., 0., None, None, None, [], [], []
    Col_ = { "PLATEpos": '-r', "PLATEneg": '-g', "SLABpos": '-g', "SLABneg": '-r'}      
    for i, Elem in enumerate(ElList):
        if Elem.Type not in ['CPE4','CPE3','CPS4','CPS3','SB3'] or Elem.MatN != KeyMat: continue # Elem.MatN: name of material
        # element geometry
        if Elem.Type=='CPE4' or Elem.Type=='CPS4':
            xN = [Elem.X0,Elem.X1,Elem.X2,Elem.X3,Elem.X0]
            yN = [Elem.Y0,Elem.Y1,Elem.Y2,Elem.Y3,Elem.Y0]
            p0.plot(xN,yN, 'b--')
            p1.plot(xN,yN, 'b--')
            p2.plot(xN,yN, 'b--')
            xN = [Elem.X0,Elem.X1,Elem.X2,Elem.X3]
            yN = [Elem.Y0,Elem.Y1,Elem.Y2,Elem.Y3]
        elif Elem.Type=='CPE3' or Elem.Type=='CPS3' or Elem.Type=='SB3':
            xN = [Elem.X0,Elem.X1,Elem.X2,Elem.X0]
            yN = [Elem.Y0,Elem.Y1,Elem.Y2,Elem.Y0]
            p0.plot(xN,yN, 'b--')
            p1.plot(xN,yN, 'b--')
            p2.plot(xN,yN, 'b--',linewidth=0.5)
            xN = [Elem.X0,Elem.X1,Elem.X2]
            yN = [Elem.Y0,Elem.Y1,Elem.Y2]
#        xC = dot( Elem.FormX(0,0,0), xN)                                    # center coordinates
#        yC = dot( Elem.FormX(0,0,0), yN)
#        P0.text(xC,yC,("%i "%(Elem.Label)),ha='left',va='top',color='black', fontsize='medium') # plot element label
        Thickness = Elem.Geom[1,0]
        # integration point loop
        for j in range(Elem.nIntLi):                                        
            r = SamplePoints[Elem.IntT,Elem.nInt-1,j][0]
            s = SamplePoints[Elem.IntT,Elem.nInt-1,j][1]
            t = SamplePoints[Elem.IntT,Elem.nInt-1,j][2]
            f = Elem.Geom[1,0]*Elem.Geom[0,0]*SampleWeight[Elem.IntT,Elem.nInt-1,j]*Elem.JacoD(r,s, t) # weighting factor
                                                                            # Geom[0,0] = 1.# dummy for Area / Jacobi determinant used instead
                                                                            # Geom[1,0]                  # thickness
            xI = dot( Elem.FormX(r,s,0), xN)                                # global integration point coordinate
            yI = dot( Elem.FormX(r,s,0), yN)                                # global integration point coordinate
            if PlotFlag:
                xIL += [xI]
                yIL += [yI]
#           sigx, sigy, sigxy = Elem.Data[j,3], Elem.Data[j,4], Elem.Data[j,5]  # retrieve stresses/moments / loading
            if   Elem.Type=='CPE4' or Elem.Type=='CPE3' or Elem.Type=='CPS4' or Elem.Type=='CPS3': sigx, sigy, sigxy = Elem.Data[j,4], Elem.Data[j,5], Elem.Data[j,7]
            elif Elem.Type=='SB3':                                                                 sigx, sigy, sigxy = Elem.Data[j,3], Elem.Data[j,4], Elem.Data[j,5]
            else: print('Element type not supported')
            pp = PrinC( sigx, sigy, sigxy)                                  # principal stresses, moments
#            scaleP = Sca[KeyElset][0]
            if PlotFlag:
                scaleP = Sca[KeyElset][0]
                SS = pp[0]                                                  # principal stress plot
                if SS>=0:                
                    p0.plot(                     [xI-SS*scaleP*pp[1],xI+SS*scaleP*pp[1]],[yI-SS*scaleP*pp[2],yI+SS*scaleP*pp[2]], Col_[Type.upper()+'pos'])
                    if Elem.Type=='SB3': p1.plot([xI-SS*scaleP*pp[1],xI+SS*scaleP*pp[1]],[yI-SS*scaleP*pp[2],yI+SS*scaleP*pp[2]], Col_[Type.upper()+'pos'])
                else:                    
                    p0.plot(                     [xI-SS*scaleP*pp[1],xI+SS*scaleP*pp[1]],[yI-SS*scaleP*pp[2],yI+SS*scaleP*pp[2]], Col_[Type.upper()+'neg'])
                    if Elem.Type=='SB3': p1.plot([xI-SS*scaleP*pp[1],xI+SS*scaleP*pp[1]],[yI-SS*scaleP*pp[2],yI+SS*scaleP*pp[2]], Col_[Type.upper()+'neg'])
                SS = pp[3]
                if SS>=0: 
                    p0.plot(                     [xI-SS*scaleP*pp[2],xI+SS*scaleP*pp[2]],[yI+SS*scaleP*pp[1],yI-SS*scaleP*pp[1]], Col_[Type.upper()+'pos'])
                    if Elem.Type=='SB3': p1.plot([xI-SS*scaleP*pp[2],xI+SS*scaleP*pp[2]],[yI+SS*scaleP*pp[1],yI-SS*scaleP*pp[1]], Col_[Type.upper()+'pos'])
                else:     
                    p0.plot(                     [xI-SS*scaleP*pp[2],xI+SS*scaleP*pp[2]],[yI+SS*scaleP*pp[1],yI-SS*scaleP*pp[1]], Col_[Type.upper()+'neg'])
                    if Elem.Type=='SB3': p1.plot([xI-SS*scaleP*pp[2],xI+SS*scaleP*pp[2]],[yI+SS*scaleP*pp[1],yI-SS*scaleP*pp[1]], Col_[Type.upper()+'neg'])
            if Elem.Type=='CPE4' or Elem.Type=='CPE3' or Elem.Type=='CPS4' or Elem.Type=='CPS3':
                NoPoint +=1                                                 # counter for all points
                rhox, rhoy, percx, percy, asx, asy, phic = 0., 0., 0., 0., 0., 0., 0.
                if pp[0]<0 and pp[3]<0: 
                    sigc = min(pp[0],pp[3])                                 # compressive principal stresses only
                    if pp[0]<pp[3]: sigc, phic = pp[0], atan(pp[2]/pp[1])
                    else:           sigc, phic = pp[3], atan(pp[5]/pp[4])
                    rhox, rhoy = minRhoX, minRhoY                           # minimum reinforcement 
                else:                                                       # all other stress states
                    phi, sig = MSIG(pp)                                     # direction and value of compressive principal stress
                    if sigxy>0: phic=-pi/4                                  # prescription of compressive principal stress direction
                    else:       phic= pi/4
                    sigc = sigxy/(cos(phic)*sin(phic))                      # corresponding concrete stress
                    rhox = (sigx - sigc*cos(phic)**2)/fy                    # reinforcement ratio from equilibrium condition
                    rhoy = (sigy - sigc*sin(phic)**2)/fy
                    if rhoy>=0 and rhox<0:                                  # undesired case, iterate for better solution with rhox=0
                        rhox = 0
                        sigc, rhoy, phic, k = RHOITER( 'rhoy', sig, phi, 10)
                        ff.write('%8.4f%8.4f rhoy iteration %3i\n'%(xI,yI,k))
#                        ff.write('%8.4f%8.4f  %8.3f%8.3f%8.3f  %10.7f  %8.5f%8.5f\n'%(xI,yI,sigx,sigy,sigxy,phic,100*rhox,100*rhoy))
                    elif rhoy<0 and rhox>=0:                                # undesired case, iterate for better solution with rhoy=0
                        rhoy = 0
                        sigc, rhox, phic, k = RHOITER( 'rhox', sig, phi, 10)
                        ff.write('%8.4f%8.4f rhox iteration %3i\n'%(xI,yI,k))
#                        ff.write('%8.4f%8.4f  %8.3f%8.3f%8.3f  %10.7f  %8.5f%8.5f\n'%(xI,yI,sigx,sigy,sigxy,phic,100*rhox,100*rhoy))
                    elif rhoy<0 and rhox<0: raise NameError ("Reinforcement2D: an exceptional case")
                    rxy = sigc*cos(phic)*sin(phic) - sigxy                  # equilibrium control
                    rxx = fy*rhox - sigx + sigc*cos(phic)**2                #
                    ryy = fy*rhoy - sigy + sigc*sin(phic)**2                #
                    RR = sqrt(rxy**2+rxx**2+ryy**2)
                    if   RR >1.e-3: print((i,j,'no equilibrium',rxx,ryy,rxy,RR))
                    if rhox < minRhoX: rhox = minRhoX 
                    if rhoy < minRhoY: rhoy = minRhoY
                if PlotFlag: # and Elem.Label in [48,61,91,118]: 
                    scc, scr = 0.7, 100*2.
                    p1.plot([xI-sigc*scc*scaleP*cos(phic),xI+sigc*scc*scaleP*cos(phic)],[yI-sigc*scc*scaleP*sin(phic),yI+sigc*scc*scaleP*sin(phic)],'g-')
                    p2.plot([xI-rhox*scr*scaleP*1.       ,xI+rhox*scr*scaleP*1.       ],[yI-0.                       ,yI+0.                       ],'r-')
                    p2.plot([xI-0.                       ,xI+0.                       ],[yI-rhoy*scr*scaleP*1.       ,yI+rhoy*scr*scaleP*1.       ],'r-')
                   
                if sigc < minSigC: minSigC, PointmSC = sigc, NoPoint
                if rhox > maxRhoX: maxRhoX, PointmRX = rhox, NoPoint
                if rhoy > maxRhoY: maxRhoY, PointmRY = rhoy, NoPoint
                percx, percy = 100*rhox, 100*rhoy                           # percentage
                asx, asy     = rhox*Thickness*1.0e4, rhoy*Thickness*1.0e4   # cm**2/m
#                if Elem.Label in [49,62,92,123]:
#                    print('XXX',Elem.Label,f"  {xI:.2f} {yI:.2f} __ {sigx:.2f} {sigy:.2f} {sigxy:.2f} __ {pp[0]:.2f} {pp[3]:.2f} __ {100*rhox:.2f} {100*rhoy:.2f} {degrees(phic):.1f} {sigc:.1f}")
                if PlotFlag and rhox>minRhoX: p0.text(xI,yI,format(100*rhox,".2f"),ha='center',va='bottom',color='blue',fontsize='small') # x-large, medium xx-large, small
                if PlotFlag and rhoy>minRhoY: p0.text(xI,yI,format(100*rhoy,".2f"),ha='center',va='top',color='blue', rotation=90,fontsize='small')
                ASx = ASx + rhox*f                                          # total reinforcement volume in x direction
                ASy = ASy + rhoy*f                                          # total reinforcement volume in y direction
                Vol = Vol + f                                               # total volume
                ff.write('%6i%8.4f%8.4f  %8.3f%8.3f%8.3f  %8.5f%7.2f   %7.4f%7.4f  %8.1f%8.1f'%(NoPoint,xI,yI,sigx,sigy,sigxy,phic,sigc,percx,percy,asx,asy))
                if sigc<fc: 
                    ff.write(' - concrete stress %8.4f'%(sigc))
                    print('allowed concrete stress exceeded',xI,yI,sigc)
                ff.write('\n')

            elif Elem.Type=='SB3':
                # reinforcement bottom side,                    admiss reinfor stress, conc stress, minimum reinf ratio
                phicB,zzB,asB1,asB2, ppB,xB =SLABRE( sigx, sigy, sigxy, fy, -fc, dd, ff)   
                phicT,zzT,asT1,asT2, ppT,xT =SLABRE(-sigx,-sigy,-sigxy, fy, -fc, dd, ff)   # reinforcement top side
                phi1, vphi1, vphi2 = SLABSH(sigx,sigy,sigxy,Elem.Data[j,6],Elem.Data[j,7],asB1,asB2,asT1,asT2, 0.5*(zzB+zzT), ReinfDes[0],ReinfDes[1], dd) # shear ratio, B/T index, req. as stirrup
                if PlotFlag:
                    as1B = "%.2f"%(10000*asB1)                              # plot lower reinforcement
                    as2B = "%.2f"%(10000*asB2)
                    if asB1>0: p0.text(xI,yI,as1B,ha='center',va='bottom',color='black')#,fontsize='small')
                    if asB2>0: p0.text(xI,yI,as2B,ha='center',va='top',color='black', rotation=90)#,fontsize='small')
                    as1T = "%.2f"%(10000*asT1)                              # plot upper reinforcement
                    as2T = "%.2f"%(10000*asT2)
                    if asT1>0: p1.text(xI,yI,as1T,ha='center',va='bottom',color='black')#,fontsize='small')
                    if asT2>0: p1.text(xI,yI,as2T,ha='center',va='top',color='black', rotation=90)#,fontsize='small')
#                    p2.text(xI,yI,"%.1f"%(1000*(abs(vphi1)+abs(vphi2))))
                    vL += [abs(vphi1)+abs(vphi2)]
                ff.write('%8.4f%8.4f  %8.2f%8.5f%8.5f  %8.2f%8.5f%8.5f  %8.2f%8.4f%8.4f'%(xI,yI,degrees(phicB),10000*asB1,10000*asB2,degrees(phicT),10000*asT1,10000*asT2,degrees(phi1),vphi1,vphi2)) # 1e4 from m2 to cm2
#                if xT:   ff.write('%8.4f%8.4f%8.3f_%8.5f%8.5f_%8.5f%8.5f'%(ppB[0],ppB[3],degrees(arctan(ppB[2]/ppB[1])),zzB/dd,xB/dd,zzT/dd,xT/dd)) # for control purposes
#                if Elem.Label==13: ff.write('%8.4f%8.4f%8.3f_%8.5f%8.5f_%8.5f%8.5f'%(ppB[0],ppB[3],degrees(arctan(ppB[2]/ppB[1])),zzB/dd,xB/dd,zzT/dd,xT/dd)) # for control purposes
                ff.write('\n')
                ASx = ASx + (asB1+asT1)*f                                   # total reinforcement volume in x direction
                ASy = ASy + (asB2+asT2)*f                                   # total reinforcement volume in y direction
                Vol = Vol + f*hh
            else: print('Element type not supported')
        # end of integration point loop
    # end of element loop
    # finish data writing and shear plotting
    if Vol<ZeroD: raise NameError("ConPlad::Reinforcement2D: no volume, check material types with KeyMat")
    if  Type.upper()=="PLATE":    ff.write('\nfy             %8.4f\nfc             %8.4f\nmin rhox       %8.4f\nmin rhoy       %8.4f\nsteel spec w%8.1f\n\n'\
             'asx volume     %8.4f\nasy volume     %8.4f\nconcrete vol %8.2f\nmass x / vol %8.2f\nmass y / vol %8.2f\nmin sig_c   %8.1f\npoint      %8i\nmax rho_x   %8.2f\npoint      %8i\nmax rho_y   %8.2f\npoint      %8i\n'\
                %(fy,fc,100*minRhoX,100*minRhoY,specWeSt,ASx,ASy,Vol,
                  ASx*specWeSt/Vol,ASy*specWeSt/Vol,
                  minSigC,PointmSC,
                  100*maxRhoX,PointmRX,100*maxRhoY,PointmRY))
    elif Type.upper()=="SLAB": 
        ff.write('\n** fy             %8.4f\n** fc             %8.4f\n** min rhox       %8.4f\n** min rhoy       %8.4f\n** steel spec w%8.1f\n\n** asx volume     %8.4f\n** asy volume     %8.4f\n** concrete vol %8.2f\n** mass x / vol %8.2f\n** mass y / vol %8.2f\n'
                %(fy,fc,100*minRhoX,100*minRhoY,specWeSt,ASx,ASy,Vol, ASx*specWeSt/Vol,ASy*specWeSt/Vol))
        if PlotFlag:
            p2.tricontour(xIL, yIL, vL, linewidths=0.5, colors='k')
            cntr2 = p2.tricontourf(xIL, yIL, vL, cmap="Purples")
            P2.colorbar(cntr2, ax=p2)

    return ASx, ASy, Vol

class ConPlaD:
    def __init__(self):
        pass
    def Run(self, Name, KeyElset, PloF, Type):
        KeyElset = KeyElset.strip()
#        ScaleStress, ScD =  [], {KeyElset: [1.0]}
        print("ConPlaD: ", Name, KeyElset)
        stime = process_time()
        # collect required data
        if pth.isfile(Name+".pkl"):                                         # read restart file if there is one
            fd = open(Name+'.pkl', 'rb')                                    # has to be in sync with pickle.dumo in ConFem, ConSimFem                                   
            NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);FlElasticLT=pickle.load(fd);\
                VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
                VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
                TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);Step=pickle.load(fd);                Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);\
                NoLabToNoInd=pickle.load(fd);NoIndToCMInd=pickle.load(fd);ContinuumNodes=pickle.load(fd);CoNoToNoLi=pickle.load(fd);SecDic=pickle.load(fd);LinAlgFlag=pickle.load(fd);ResultTypes=pickle.load(fd);Header=pickle.load(fd);
            fd.close()
        else: raise NameError ("ConPlaD: .pkl file missing -- maybe wrong name definition")
        if KeyElset in SecDic:  KeyMat = SecDic[KeyElset].Mat
        else:                   raise NameError("ConPlaD: unknown elset", KeyElset)
        if pth.isfile(Name+".opt.txt"):                                     # read options, mandatory for reinforcement design
            f4=open( Name+".opt.txt", 'r')
            WrNodes, LineS, ReDes, MaxType, _ = ReadOptionsFile(f4, NodeList,NoLabToNoInd,NoIndToCMInd)
            f4.close()
            if not (len(ReDes)==7): raise NameError("ConPlaD: reinforcement design parameters missing") 
        else: 
            raise NameError ("ConPlaD: options file missing")
        # read plot options file if there is any -- for superposed scaling factors only
        if PloF:
            if pth.isfile(Name+".plt.txt"):
                f1=open( Name+".plt.txt", 'r')
                from ConFemInOut import ReadPlotOptionsFile
                ScaleDis, ScaleDis2D, PE2DFlag, PE3DFlag, PlotTimes, Post1DFlag, ScaleStress, PostNodeFl, ShellL, ScaleShellL, _ = ReadPlotOptionsFile(f1, SecDic) # to modify plot scaling factors
                f1.close()
            # this is quick and dirty to get data to derive scaling factors
            from ConFemPost import PostScales, ReadResultsFromFile
            ff  = open(Name+".elemout_"+".txt",'r')
            ffN = open(Name+".nodeout_"+".txt",'r')
            EndFlag = False
            while not EndFlag:
                EndFlag, _, Time_, _, ElemResults,_,ElResults = ReadResultsFromFile( MatList, ff, ffN) # reads until next time marker, will presumably work only for elemenout_, nodeout_
                if Time_ == Time: break
            if Time_ != Time: raise NameError("ConPlaD: inconsistency in time -- last computed time unequal to last output time")
            ff.close()
            ffN.close()
            Sc = PostScales( SecDic, ElList,ElResults, ScaleStress)
        else:
            Sc = []
        # reinforcement design
        f4=open( Name+".reinforcement.txt", 'w')
        ASx, ASy, Vol = Reinforcement2D( ElList, NodeList,NoIndToCMInd, ReDes, Sc, KeyElset,KeyMat, f4, PloF, Type, Time)
        f4.close()
        print(ASx,ASy,Vol,ASx*ReDes[6]/Vol,ASy*ReDes[6]/Vol)
        print(process_time()-stime)
        # post
        if PloF:
            plt.show()
            return 0
        else:
            import hashlib
            mmm = hashlib.md5()
            fp= open( Name+".reinforcement.txt", "r")
            while True:
                data= fp.read(65536)
                if not data: break
                mmm.update(data.encode())
            fp.close()
            RC = mmm.hexdigest()
            print(RC)
            return RC

if __name__ == "__main__":
    PlotFlag = True
#    numpy.seterr(all='raise')
    Name, KeyElset, Type ="../DataExamples/E08/E8-01",' EL1', "plate"                 # ConPlad Plate E4_01 linear elastic, reinforcement design        11s
#    Name, KeyElset, Type ="../DataExamples/E09/E9-01",'PROP1', "slab"                         # input data name
    ConPlaD_ = ConPlaD()
    RC = ConPlaD_.Run(Name, KeyElset, PlotFlag, Type)
