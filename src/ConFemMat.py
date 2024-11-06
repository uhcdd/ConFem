# ConFemMat -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
import matplotlib.pyplot as plt
from numpy import array, sin, exp, log, cos, dot, pi, double, zeros, sqrt, outer, transpose, tan, tanh, log, linspace, sign, arctan, arcsin, ones
from numpy.linalg import eigh, norm, det, solve #, eigvalsh
from numpy import fabs
from scipy import integrate
from bisect import bisect_left, bisect_right
from time import strftime#import sys 

from ConFemBasics import ZeroD, PrinCLT_, I21Points, I21Weights, SampleWeightRCShell, Echo # StrNum0
from ConFemRandom import RandomField_Routines
try:
    from ConFemMatC import *
    ConFemMatCFlag = True    
except ImportError:
    ConFemMatCFlag = False    

def SpringSpline( s1, F1, s2, F2, s0, d0, eps):
    if eps<0: si=-1
    else:        si= 1
    s = abs(eps)
    if s<= s0:
        ssig=si*d0*s
        dsig=d0
    elif s<= s1:
        t5 = s*s
        t8 = s0*s0
        t9 = d0*t8
        t11 = s0*F1
        t14 = s0*s1*d0
        t17 = 3.0*s1*F1
        t18 = s1*s1
        t19 = t18*d0
        t20 = 2.0*t19
        ssig = (-(d0*s0+s1*d0-2.0*F1)*t5*s
                   +(2.0*t9-3.0*t11+2.0*t14-t17+t20)*t5
                   -s1*(t14+t19+4.0*t9-6.0*t11)*s
                   +t8*(t20+t11-t17)) / (t8*s0+3.0*s0*t18-t18*s1-3.0*s1*t8)
        ssig = si*ssig
        dsig = (-3.0*(d0*s0+s1*d0-2.0*F1)*t5+2.0*(2.0*t9-3.0*t11+2.0*t14-3.0*s1*F1+2.0*t19)*s-s1*(t14+t19+4.0*t9-6.0*t11))/(t8*s0+3.0*s0*t18-t18*s1-3.0*s1*t8)
    elif s<=s2:
        t1  = -F1+F2
        t20 = s1*s1
        t21 = t20*s1
        t14 = s2*s2
        t15 = t14*s2
        t17 = t14*s1
        t2  = s*s
        ssig = (-2.0*t1*t2*s+3.0*t1*(s2+s1)*t2-6.0*t1*s2*s1*s+t15*F1-3.0*t17*F1-F2*t21+3.0*F2*s2*t20)/(-t21-3.0*t17+t15+3.0*s2*t20)
        ssig = si*ssig
        dsig = 6.0*(-t1*t2+t1*(s2+s1)*s-t1*s2*s1)/(-t20*s1-3.0*s1*t14+t14*s2+3.0*s2*t20)
    else:
        ssig = si*F2
        dsig = 0.
    return ssig, dsig

def EigenJacobiSym(A,N):                                        # from applied numerical methods 4.8
    itmax = 10 # 10
    Eps1,Eps2,Eps3 = 1.0e-10, 1.0e-10, 1.0e-7
    # initializations
    T       = zeros((N,N), dtype=float)
    AIK     = zeros((N), dtype = float)
    EIGEN   = zeros((N), dtype = float)
    sigma1, OffDsq = 0., 0.
    for i in range(N):
        T[i,i] = 1.0
        sigma1 = sigma1 + A[i,i]**2
        for j in range(i+1,N): OffDsq = OffDsq + A[i,j]**2      # considers symmetry
    S = 2.0*OffDsq + sigma1                                     # this should not change during transformations
    if sqrt(S)<Eps1: 
        eV = zeros((N), dtype=float)
        eV[0] = 1.
        return 0., eV
    # iteration loop
    for iter in range(itmax):
        for i in range(N-1):
            for j in range(i+1,N):
#                if abs(A[i,j])<= Eps2: break    #    this is an error uhc190612. If, e.g. A[0,1] is zero, and A[0,2] not the loop is left to early
                Q = abs(A[j,j]-A[i,i])                          # diagonal i, j difference
                # compute sine and cosine of rotation angle
                if Q>=Eps1:                                     # see carnahan: applied numerical methods, p. 251
                    siQ = Q/(A[i,i]-A[j,j])
                    P   = 2.*A[i,j] * siQ
                    alpha = P/Q
                    CSA = sqrt( 0.5*(1. + 1./sqrt(1.+alpha**2)) ) # this is obviously larger than sqrt( 1/2 )
                    SNA = sign(P)/( 2.*CSA*sqrt(1. + 1./alpha**2) )
                else:
                    CSA = 1./sqrt(2.)
                    SNA = CSA
                # update columns i and j of T --> T_old * T
                for k in range(N):
                    HOLDKI = T[k,i]
                    T[k,i] = HOLDKI*CSA + T[k,j]*SNA
                    T[k,j] = HOLDKI*SNA - T[k,j]*CSA
                # rows i and j --> transpose(T) * A
                for k in range(i,N):
                    if k<=j:
                        AIK[k] = A[i,k] 
                        A[i,k] = AIK[k]*CSA + A[k,j]*SNA
                        if k==j: A[j,k] = AIK[k]*SNA - A[j,k]*CSA
                    else:
                        HOLDKI = A[i,k] 
                        A[i,k] = HOLDKI*CSA + A[j,k]*SNA
                        A[j,k] = HOLDKI*SNA - A[j,k]*CSA
                
                AIK[j] = SNA*AIK[i] - CSA*AIK[j]
                # columns i and j --> A * T
                for k in range(j+1):
                    if k>i:
                        A[k,j] = SNA*AIK[k] - CSA*A[k,j]
                    else:
                        HOLDKI = A[k,i]
                        A[k,i] = HOLDKI*CSA + A[k,j]*SNA
                        A[k,j] = HOLDKI*SNA - A[k,j]*CSA
#                A[i,j] = 0.
        # find sigma2 for transformed A and test for convergence
#        print('ZZZ', iter, A[0,0],A[0,1],A[0,2],A[1,0],A[1,1],A[1,2],A[2,0],A[2,1],A[2,2])
        sigma2 = 0.
        for i in range(N):
            EIGEN[i] = A[i,i]
            sigma2 += EIGEN[i]**2
        if (1.-sigma1/sigma2) < Eps3 and abs(sigma1-S)<Eps3:
            # convergence
            val, ind = -1.0e3, -1
            for i in range(N):                                  # find largest eigenvalue
                if EIGEN[i]>val:
                    ind = i
                    val = EIGEN[i] 
            return EIGEN[ind], T[:,ind]
        sigma1 = sigma2
    # no convergence
    raise NameError(":EigenJacobiSym 1")
def LinearRegression( x, f):
    nx = len(x)
    nf = len(f)
    if nx!=nf: raise NameError("linear regression error 1")
    X = ones((nx,2), dtype=float)
    for i in range(nx): 
        X[i,0] = x[i]
    XX = dot(transpose(X),X)
    Xf = dot(transpose(X),f)
    a = solve(XX,Xf)
    return a
def QuadraticRegression( x, f):
    nx = len(x)
    nf = len(f)
    if nx!=nf: raise NameError("quadratic regression error 1")
    X = ones((nx,3), dtype=float)
    for i in range(nx): 
        X[i,0] = x[i]**2
        X[i,1] = x[i]
    XX = dot(transpose(X),X)
    Xf = dot(transpose(X),f)
    a = solve(XX,Xf)
    return a

class Material(object):
    Density   = 0
    Symmetric = False
    RType     = None
    Used      = False                               # flag whether this material is actually used (referenced to by elements), may be later set to True
    Update    = False
    Updat2    = False
    Conc      = None                                # used as handle if actual material is a wrapper, e.g. WraRCShell
    
    ## temporarily
#    Sig1, Sig2, Sig3, Eps1, DEps, CC, C1 = [0], [0], [0], [0], [0], [0], [0]
#    Incr, Itera = [0], [0]
#    Sig1_, Sig2_, Sig3_, Eps1_, DEps_, CC_, C1_ = 0, 0, 0, 0, 0, 0, 0 #Sig1_, Sig1 otherwise used
    ##

    def __init__(self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal, TypeVal):
        self.Symmetric = SymmetricVal               # flag (True/False) for symmetry of material matrices
        self.RType     = RTypeVal                   # string for type of regularization (None: no regularization)
        self.Update    = UpdateVal                  # flag whether state variables have update method
        self.Updat2    = Updat2Val                  # flag whether state variables have update method connected to special conditions
        self.StateVar  = StateVarVal                # number of state variables (None: no state variables)
        self.NData     = NDataVal                   # number of data items for output                NData overridden in ComFemInOut::DataInput
        self.Conc      = None                       # used if actual material is a wrapper, e.g. WraRCShell
        self.CrBwN     = 50                         # crack band regularization: number of intervals for scaling factor interpolation
        self.Type      = TypeVal                    #
        self.matbStiff = 1.0                        # used for beams / MisesBeam2D
    def Mass(self, Elem):                           # mass matrix
        if Elem.dim in [10,12]:                                             # bernoulli beam
            val = self.Density*Elem.Geom[1,3]                               # beam mass per unit length from element integration point
            MatS = array([[val,0],[0,val]])                                 # beam mass matrix
        elif Elem.dim in [11,13]:                                           # timoshenko beam
            val = self.Density*Elem.Geom[1,3]                               # beam mass per unit length from element integration point
            MatS = array([[val,0,0],[0,val,0],[0,0,val]])                   # beam mass matrix -- rotational mass set to translational mass for numerical aspects
#            MatS = array([[val,0,0],[0,val,0],[0,0,0]])                   # beam mass matrix -- rotational mass set to translational mass for numerical aspects
        elif Elem.dim in [2, 4, 96, 97]:                                    # CP*, CAX, axisym bond element, 2D bond element
            val = self.Density                                              # plate (quad) mass per unit thickness -- thickness is considered in IntForces
            MatS = array([[val,0],[0,val]])
        elif Elem.dim in [5]:                                               # TAX2, TAX3, TAX2I
            val = self.Density #   * Elem.Geom[1,0] #is in system integration
            MatS = array([[val, 0], [0, val]])
        elif Elem.dim==1:
            val = self.Density # *Elem.Geom[1,0] is in system integration      #
            if Elem.Type in ["T3D2","T3D3","T3D2I","T3D3I","Bond3D2","Bond3D3"]: MatS = array([[val,0,0],[0,val,0],[0,0,val]])
            else:                                                            MatS = array([[val,0],[0,val]])
        elif Elem.dim==21:                      
            val = self.Density                                              # thickness is considered in IntForces
            MatS = array([[val,0,0],            # !!! remains to be validated
                          [0,val,0],
                          [0,0,val]])
        elif Elem.dim in [3,95]:                                            # 3D continuum, 3D bond element
            val = self.Density
            MatS = array([[val,0,0],[0,val,0],[0,0,val]])
        else: raise NameError("ConFemMaterial::Material.Mass: mass not yet defined for this element type", Elem.dim)
        return MatS
    def ViscExten3D(self, Dt, eta, Dps, Elem, ipI, sI):
        VepsOld = array([Elem.StateVar[ipI,sI],Elem.StateVar[ipI,sI+1],Elem.StateVar[ipI,sI+2],Elem.StateVar[ipI,sI+3],Elem.StateVar[ipI,sI+4],Elem.StateVar[ipI,sI+5]]) # strain rate of previous time step 
        if norm(VepsOld)<ZeroD or dot(Dps,VepsOld)<0.: 
            if Dt>ZeroD: VepsOld = Dps/Dt
            else:        VepsOld = zeros((6), dtype=float)
        if Dt>ZeroD and norm(Dps)>ZeroD: Veps = 2.*Dps/Dt - VepsOld             # actual strain rate
        else:                            Veps = VepsOld
        if Dt>ZeroD: zz = 2.*self.eta/Dt
        else:        zz = 0.
        Elem.StateVarN[ipI,sI]   = Veps[0]
        Elem.StateVarN[ipI,sI+1] = Veps[1]
        Elem.StateVarN[ipI,sI+2] = Veps[2]
        Elem.StateVarN[ipI,sI+3] = Veps[3]
        Elem.StateVarN[ipI,sI+4] = Veps[4]
        Elem.StateVarN[ipI,sI+5] = Veps[5]
        return zz, Veps
    def SpeCrEnergySample(self, lam, eps_ct, PlotF,P0,epsEndP, UniaxTension): # crack energy uniaxial tension for unit volume (specific crack energy)
                                                                            # lambda is scaling factor
        def IntEnd( UniaxTension ):                                         # stress strain in strain softening range with upper limit
            xP, yP = [], []
            tol = 0.001*self.fct
            de = 0.01e-3
            ee = eps_ct                                                     # start of integration: with strength strain tension
            x = 1.*self.fct
            while x>tol:                                                    # to determine a suitable end for integration of specific crack energy
                x = UniaxTension( lam, ee)                                  # stress
                ee += de
                xP += [ee]                                                  # collect strain
                yP += [x]                                                   # collect stress
            return ee, xP, yP                                               # return last strain value
        # end def
        def IntCrE( UniaxTension ):
            nI = 100                                                        # number of integration samples
            sia = zeros((nI+1),dtype=double)                                # stress
            ee = eps_ct
            de = (epsEnd-ee)/nI
            for i in range(nI+1):                                           # generate samples
                sia[i] = UniaxTension( lam, ee)
                ee += de
            return integrate.trapz( sia, x=None, dx=de)                     # integrate samples
        # end def
        epsEnd,xP,yP = IntEnd( UniaxTension )
        CrE = IntCrE( UniaxTension )                                        # uses epsEnd
        if PlotF:
            xP_, yP_ = [], []
            if epsEndP!=None:
                for i, e in enumerate(xP):
                    if e <= epsEndP:
                        xP_ += [e]
                        yP_ += [yP[i]]
                    else:
                        break
            else:
                xP_, yP_ = xP, yP
            P0.plot(xP_,yP_,label=str(lam)[0:6])
            P0.grid()
        return CrE                                                          # unit volume - specific crack energy
    def SpeCrEnergyData(self, type, la1, la2, n, eps_ct, UniaxTensionScaled):
                                                                            # scaled uniaxial tension: samples of scaling factor depending on specific crack energy
        qq = pow(la2/la1,1./n)                                              # factor of geometric row for sampling points
        CrX, CrY = zeros((n+1),dtype=float), zeros((n+1),dtype=float)
        for i in range(n+1):
            arg = la1*qq**i                                                 # actual scaling factor
            CrX[n-i] = self.SpeCrEnergySample( arg, eps_ct, False,None,None, UniaxTensionScaled )/self.SpecCrEn  # scaled specific crack enerty related to specific crack energy
                                                                            # SpeCrEnergy-->specific crack energy scaled by arg; SpecCrEn-->unscaled specific crack energy (specific is per unit volume) 
                                                                            # should correspond to (crack energy) / (scaled crack band width)
            CrY[n-i] = arg                                                  # current scaling factor
        # following looks weird - compare ConFemElem::CrBScaleType
        if False:
            X = self.RegPar/self.SpecCrEn                                       # physical crack band width
            CrL = X/CrX                                                         # list of scaled element length
            sep = self.ElCharLenBound*X # self.RegPar/self.SpecCrEn             # separates the approaches for scaling factors
            a   =  LinearRegression(CrL,CrY)
            b   =  QuadraticRegression(CrL,CrY)
            if self.Type in ['ISODAMAGE']:
                if   type == "SmallEl": self.SmallElReg = [min(CrL),max(CrL),sep,a[0],a[1]]
                elif type == "LargeEl": self.LargeElReg = [min(CrL),max(CrL),sep,b[0],b[1],b[2]]
                else: raise NameError("ConFemMat::SpeCrEnergyData: unknown regularization scaling type")
            if False:                                                            # plot and write data and crack band regularization data
                PP = plt.figure()
                P0 = PP.add_subplot(111,title='L_c -> lam')
                P0.plot(CrL,CrY)
                P0.grid()
                P0.plot([sep,sep],[0,max(CrY)])
                xr = linspace(min(CrL),max(CrL),10)
                y0 = ones((len(xr)), dtype=float)
                yr = a[0]*xr + a[1]*y0
                P0.plot(xr,yr)
                yr = b[0]*xr*xr + b[1]*xr + b[2]*y0
                P0.plot(xr,yr)
                exportTime = strftime("%Y-%m-%d_%H-%M-%S")
                if self.Type in ['ISODAMAGE']:
                    ff = open("../LogFiles/IsoDam."+exportTime+".txt",'w')
                    ff.write("Young's modulus      %14.6e\nPoission's ratio     %14.6e\ncompressive strength %14.6e\ntensile strength     %14.6e\nstrain com. strength %14.6e\nstrain tens. strength%14.6e\n"\
                             %(self.Emod,self.nu,self.fc,self.fct,self.epsC,self.eps_ct))
                    ff.write("crack energy         %14.6e\nspec. crack energy   %14.6e\ncrack band width     %14.6e\n"\
                             %(self.RegPar,self.SpecCrEn,self.RegPar/self.SpecCrEn))
                    ff.write("regularization\nsmallest el. length  %14.6e\nlargest el. length   %14.6e\nlength separator     %14.6e\n"\
                             %(min(CrL),max(CrL),sep))
                    ff.write("lin. regression lin. %14.6e\nlin. regression cons.%14.6e\n"%(a[0],a[1]))
                    ff.write("quad. reg. quad.     %14.6e\nquad. reg. lin.      %14.6e\nlquad. reg. cons.    %14.6e"%(b[0],b[1],b[2]))
                    ff.close()
        return CrX, CrY
    def Rankine(self, ConFemMatCFlag, Elem, vec ):
        if ConFemMatCFlag:
            laMax= zeros((1),dtype=float)
            eVec = zeros((3),dtype=float)
            if eigJacobiSymWrapper( vec, laMax, eVec) != 100:
                raise NameError("ConFemMat::ElasticLT.Sig: Eigensystem iteration failed: %12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e"%(vec[0,0], vec[1,1], vec[2,2], vec[1,2], vec[0,2], vec[0,1]))
        else:
            laMax, eVec = EigenJacobiSym( array( [[vec[0],vec[5],vec[4]], [vec[5],vec[1],vec[3]], [vec[4],vec[3],vec[2]]] ), 3)
        if Elem.dim==2 and abs(eVec[2])>ZeroD:  laMax=0.                    # lateral extension to plane, no lateral cracks allowed
        if (eVec[0]+eVec[1]+eVec[2]<0.):                                    # for scalar product with (1,1,1) -- directions indicate the same plane (approximately for different IPs) but might have opposite directions. Has to be avoided 
            eVec[0] = -eVec[0]
            eVec[1] = -eVec[1]
            eVec[2] = -eVec[2]
        return laMax, eVec
    def RankineUpd(self, Elem, ipI, sig, eVec, laMax, offset):
        tt = [sig[0]*eVec[0]+sig[5]*eVec[1]+sig[4]*eVec[2],                 # crack traction in global system -- might have nonzero local shear components due to viscous contributions  
              sig[5]*eVec[0]+sig[1]*eVec[1]+sig[3]*eVec[2], 
              sig[4]*eVec[0]+sig[3]*eVec[1]+sig[2]*eVec[2]]
        if norm(tt) < ZeroD : laMax=0.                                      # no crack allowed in case of orthogonal sig and eVec e.g. uniaxial compr. 3D - Ahmad
        Elem.StateVarN[ipI,6+offset]  = laMax
        Elem.StateVarN[ipI,7+offset]  = eVec[0]                             # eigenvector of largest principal strain / stress
        Elem.StateVarN[ipI,8+offset]  = eVec[1]
        Elem.StateVarN[ipI,9+offset]  = eVec[2]
        Elem.StateVarN[ipI,10+offset] = tt[0]                               # stress vector derived from stress tensor and above eigenvector 
        Elem.StateVarN[ipI,11+offset] = tt[1]
        Elem.StateVarN[ipI,12+offset] = tt[2]
    def CrackBandScale(self, elLabel, Lch_, ScaleType):                     # called by ConFemElem:CrBScaleType, Lch_ comes from Elem.Lch_
        if self.RType ==2:                                                  # find scaling factor for band width regularization 
            x = self.bw/Lch_                                                # physical crack band width related to char element length -> scaled g_f / unscaled g_f
            # scaling approach 1                                            #       |_  1/Lch_ * G_f/g_f,unscaled = 1/Lch_ * g_f,scaled * Lch_ / g_f,unscaled
            terminate = False
            if ScaleType == 1:
                CrX, CrY = self.CrX, self.CrY                               # support points CrX: bw=Lchar, CrY: scaling factor - ascending order - smaller value for larger element
                i = bisect_left(CrX, x)
                if i>0 and i<(self.CrBwN+1):                                # BrBwN number of intervals for scaling factor interpolation
                    CrBwS = CrY[i-1] + (x-CrX[i-1])/(CrX[i]-CrX[i-1])*(CrY[i]-CrY[i-1]) # scaling factor by linear interpolation
                else:                                             # Lch_ falls below admissible element length self.bw/CrX
                    print(f"ConFemMat::Material.CrackBandScale: el {elLabel:d}, Lchar {Lch_:.2f} "
                          f"req.x=g_f,scaled/g_f,unscaled={x:.2f}: falls out of Type 1 regularization range [{CrX[0]:.2f}..{CrX[-1]:.2f}]")
                    terminate = True
            # scaling approach 2
            elif ScaleType == 2:
                CrX, CrY = list(reversed(self.CrX2)), list(reversed(self.CrY2)) # support points for tensile softening scaling factor approach 2 for small elements -- see remarks above
                i = bisect_left(CrX, x)
                if i>0 and i<(self.CrBwN+1):                                # to kick out larger elements with scaling approach 2
                    CrBwS = CrY[i-1] + (x-CrX[i-1])/(CrX[i]-CrX[i-1])*(CrY[i]-CrY[i-1]) # scaling factor by linear interpolation
                else:
                    print(f"ConFemMat::Material.CrackBandScale: el {elLabel:d}, Lchar {Lch_:.2f}"
                          f" req.x=g_f,scaled/g_f,unscaled {x:.2f}: falls out of Type 2 regularization range [{CrX[0]:.2f}..{CrX[-1]:.2f}]")
                    terminate = True
            if terminate: raise NameError("terminated with crack band scaling")
        else:
            CrBwS = 1.0
        return CrBwS

class Template(Material):                              # elastoplastic Mises
    def __init__(self, PropMat):
#    def __init__(self,         SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,      False,     1,           8)
#        self.Symmetric = True                       # flag for symmetry of material matrices
#        self.Update = True                          # has specific update procedure for update of state variables
#        self.Updat2 = False                         # no 2 stage update procedure
#        self.StateVar = 1                           # number of state variables per integration point
#        self.NData = 8                              # number of data items

        self.Emod = PropMat[0]                      # Young's modulus
        self.nu = PropMat[1]                        # Poissons's ratio
        self.sigY = PropMat[2]                      # uniaxial yield stress
        self.sigU = PropMat[3]                      # strength
        self.epsU = PropMat[4]                      # limit strain
        self.alphaT = PropMat[5]                    # thermal expansion coefficient
        self.Density = PropMat[6]                   # specific mass
        self.epsY = self.sigY/self.Emod             # uniaxial yield strain  ????
        self.Etan = (self.sigU-self.sigY)/(self.epsU-self.epsY)# tangential / hardening modulus
        self.H = self.Etan/(1-self.Etan/self.Emod)  # hardening modulus
        
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        sigy = max(self.sigY,Elem.StateVar[ipI][0]) # current state parameter - current uniaxial yield stress
        Elem.StateVarN[ipI][0] = 0                  # current values of state variables have to initialized again
        if Elem.dim==2 or Elem.dim==3:           # plane stress/strain biaxial
            nu = self.nu                            # Poisson's ratio
            mu = self.Emod/(2.*(1.+nu))             # shear modulus
            C0 = self.Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),nu/(1-nu),0,0,0],
                                                           [nu/(1-nu),1,nu/(1-nu),0,0,0],
                                                           [nu/(1-nu),nu/(1-nu),1,0,0,0],
                                                           [0,0,0,(1-2*nu)/(2*(1-nu)),0,0],
                                                           [0,0,0,0,(1-2*nu)/(2*(1-nu)),0],
                                                           [0,0,0,0,0,(1-2*nu)/(2*(1-nu))]]) # triaxial isotropic elasticity 
            if Elem.dim==2:                         # plate plane stress / plane strain
                Sig = array( [Elem.DataP[ipI,3],Elem.DataP[ipI,4],Elem.DataP[ipI,6],0.,0.,Elem.DataP[ipI,5]] ) # stress of previous increment
                dEps = array([Dps[0],Dps[1],0.,0.,0.,Dps[2]]) # total strain increment
            elif Elem.dim==3:
                Sig = array([Elem.DataP[ipI,0],Elem.DataP[ipI,1],Elem.DataP[ipI,2],Elem.DataP[ipI,3],Elem.DataP[ipI,4],Elem.DataP[ipI,5]] ) # stress of previous increment
                dEps = array([Dps[0],Dps[1],Dps[2],Dps[3],Dps[4],Dps[5]]) # total strain increment
            Fn = sqrt(3.*(((Sig[0]-Sig[1])**2+(Sig[0]-Sig[2])**2+(Sig[1]-Sig[2])**2)/6.+Sig[3]**2+Sig[4]**2+Sig[5]**2))-sigy # distance to yield surface of previous step
            Eflag = False
            dEpp = zeros((6),dtype=float)           # initial value plastic strain increment
            dSig = zeros((6),dtype=float)           # initial value plastic strain increment
            dLam = 0.                               # initial value plastic multiplier increment
            ni = 20                                 # iteration limit
            for i in range(ni):
                if Elem.PlSt: dEps[2]=-( C0[2,0]*(dEps[0]-dEpp[0])
                                        +C0[2,1]*(dEps[1]-dEpp[1])
                                        +C0[2,3]*(dEps[3]-dEpp[3])
                                        +C0[2,4]*(dEps[4]-dEpp[4])
                                        +C0[2,5]*(dEps[5]-dEpp[5]))/C0[2,2]+dEpp[2] # lateral strain for plane stress
                dSig = dot(C0,dEps-dEpp)
                SigN = Sig + dSig
                J2 = ((SigN[0]-SigN[1])**2+(SigN[0]-SigN[2])**2+(SigN[1]-SigN[2])**2)/6.+SigN[3]**2+SigN[4]**2+SigN[5]**2 # 2nd stress deviator invariant
                if sqrt(3.*J2)-sigy<1.e-9:          # elastic loading, unloading or reloading
                    Eflag = True
                    break
                sm = (SigN[0]+SigN[1]+SigN[2])/3.   # 1st stress invariant / mean stress predictor stress
                rr = sqrt(3./(4.*J2))*array([SigN[0]-sm,SigN[1]-sm,SigN[2]-sm,SigN[3],SigN[4],SigN[5]]) # yield gradient predictor stress
                dL = (Fn + dot(rr,dSig)+rr[3]*dSig[3]+rr[4]*dSig[4]+rr[5]*dSig[5] - self.H*dLam)/(3.*mu+self.H) # with Voigt notation correction
                if dL<1.e-9: break
                dLam = dLam + dL                    # update plastic multiplier
                dEpp = dLam*array([rr[0],rr[1],rr[2],2.*rr[3],2.*rr[4],2.*rr[5]]) # plastic strain incremen with Voigt notation correction
            if i>=ni-1:  
                print(elI, ipI, i, ni, Sig, dEps, sigy, dLam, dEpp, dSig)  
                raise NameError ("ConFemMaterials::Mises.Sig: no convergence")
            Sig = SigN
            if Eflag:                               # elastic loading or unloading / reloading
                if Elem.dim==2:                     # plane stress / strain
                    if Elem.PlSt: MatM = self.Emod/(1-nu**2)*array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]]) # plane stress
                    else:         MatM = self.Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),0],[nu/(1-nu),1,0],[0,0,(1-2*nu)/(2*(1-nu))]]) # plane strain
                else: MatM = C0
            else:                                   # plastic loading
                Elem.StateVarN[ipI][0]= sqrt(3.*J2) # new equivalent yield limit 
                A = 1./(self.H + 3.*mu)             # --> 2.*mu*dot(rr,rr) --> dot(xx,rr) --> dot(C0,rr)
                CC = C0 -  A*4.*mu**2*outer(rr,rr)  # --> A*outer(xx,xx)            # tangential material stiffness
                cD = 1./CC[2,2]
                if Elem.dim==2:
                    if Elem.PlSt: MatM=array([[CC[0,0]-CC[0,2]*CC[2,0]*cD,CC[0,1]-CC[0,2]*CC[2,1]*cD,CC[0,5]-CC[0,2]*CC[2,5]*cD],
                                              [CC[1,0]-CC[1,2]*CC[2,0]*cD,CC[1,1]-CC[1,2]*CC[2,1]*cD,CC[1,5]-CC[1,2]*CC[2,5]*cD], 
                                              [CC[5,0]-CC[5,2]*CC[2,0]*cD,CC[5,1]-CC[5,2]*CC[2,1]*cD,CC[5,5]-CC[5,2]*CC[2,5]*cD]])
                    else:         MatM=array([[CC[0,0],CC[0,1],CC[0,5]],[CC[1,0],CC[1,1],CC[1,5]],[CC[5,0],CC[5,1],CC[5,5]]])
                elif Elem.dim==3:
                    MatM=array([[CC[0,0],CC[0,1],CC[0,2],CC[0,3],CC[0,4],CC[0,5]],
                                [CC[1,0],CC[1,1],CC[1,2],CC[1,3],CC[1,4],CC[1,5]],
                                [CC[2,0],CC[2,1],CC[2,2],CC[2,3],CC[2,4],CC[2,5]],
                                [CC[3,0],CC[3,1],CC[3,2],CC[3,3],CC[3,4],CC[3,5]],
                                [CC[4,0],CC[4,1],CC[4,2],CC[4,3],CC[4,4],CC[4,5]],
                                [CC[5,0],CC[5,1],CC[5,2],CC[5,3],CC[5,4],CC[5,5]]])
            if Elem.dim==2:
                sig = array([Sig[0],Sig[1],Sig[5]])
                return sig, MatM, [Eps[0], Eps[1], Eps[2], Sig[0], Sig[1], Sig[5], Sig[2]] # data
            else:                                   # 3D, shell
                sig = array([Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5]])
                return sig, MatM, [Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5], 0.] # data
        else: raise NameError ("ConFemMaterials::Mises.Sig: not implemented for this element type")
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):    # loop over integration points 
            if Elem.StateVarN[j,0]>Elem.StateVar[j,0]: Elem.StateVar[j,0] = Elem.StateVarN[j,0]
        return False


class PolyBeam(Material):  # elastoplastic Mises
    def __init__(self, Label, eMin, eMax, PropE, PropS):
        #    def __init__(self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,      False,     1,           8, "POLYBEAM")
        #        self.Symmetric = True                       # flag for symmetry of material matrices
        #        self.Update = True                          # has specific update procedure for update of state variables
        #        self.Updat2 = False                         # no 2 stage update procedure
        #        self.StateVar = 1                           # number of state variables per integration point
        #        self.NData = 8                              # number of data items
        self.Label = Label
        self.epsCU = eMin  # limit for compressive strain (signed)
        self.epsTU = eMax  # limit for tensile strain  (signed)
        if len(PropE) != len(PropS): raise NameError("CrazyMat 1")
        self.epsL = PropE
        self.sigL = PropS
    def SigEps(self, eps):
        for i, e in enumerate(self.epsL):
            if i > 0:
                Emod = (self.sigL[i] - self.sigL[i - 1]) / (self.epsL[i] - self.epsL[i - 1])
                sig = self.sigL[i - 1] + Emod * (eps - self.epsL[i - 1])
                if eps <= e: break
        return sig, Emod
    def IntCrossSec(self, epsRef, kap, z1, z2, Elem, nI=50):
        dz = (z2-z1)/nI                                     # z1 local lower, z2 local upper; cross sectional height / number of integration points
        nn = zeros((nI+1), dtype=float)
        mm = zeros((nI+1), dtype=float)
        nde= zeros((nI+1), dtype=float)
        ndk= zeros((nI+1), dtype=float)
        mdk= zeros((nI+1), dtype=float)
        bb = Elem.Geom[1,1]                                     # width of cross section, may eventually be overridden
        z = z1                                              # cross section / strain coordinate, initial value
        for i in range(nI+1):                               # numerical integration of concrete contribution
            eps  = epsRef - z*kap                           # fiber strain
            sig, dsig = self.SigEps( eps)                    # uniaxial stress strain definition
            sig_ = bb*sig                                   # consider width
            dsig_= bb*dsig
            nn[i]  =    sig_                                # normal force contribution
            mm[i]  = -z*sig_                                # moment
            nde[i] =       dsig_                            # normal force grad eps
            ndk[i] = -z   *dsig_                            # normal force grad kappa
            mdk[i] =  z**2*dsig_                            # moment grad kappa
            z = z + dz
        dz_ = fabs(dz)
        NN  = integrate.trapz( nn, x=None, dx=dz_)          # numerical integration normal force
        MM  = integrate.trapz( mm, x=None, dx=dz_)          # numerical integration moment
        NDE = integrate.trapz( nde, x=None, dx=dz_)         # numerical normal force grad eps
        NDK1= integrate.trapz( ndk, x=None, dx=dz_)         # numerical normal force grad kappa - same as moment grad eps ?
        MDK = integrate.trapz( mdk, x=None, dx=dz_)         # numerical moment grad kappa
        return NN, MM, NDE, NDK1, MDK
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        if Elem.CrSecType=="POLYLINE": z1, z2, =  Elem.Geom[1,1], Elem.Geom[1,2]
        else:                          z1, z2, = -Elem.Geom[1,2]/2., Elem.Geom[1,2]/2. # coordinate of bottom fibre/lower edge, top fibre/upper edge
        epsR = Eps[0]
        kap  = Eps[1]
        if fabs(kap)>ZeroD:
            zCu = (epsR - self.epsCU) / kap                     # position of compressive strain limitation
            zTu = (epsR - self.epsTU) / kap                     # position of tensile strain limitation
        else:
            zCu, zTu = z2, z1
        if zTu < z1: zTu = z1                               # tensile limit line below cross section - upward bending compression
        if zTu > z2: zTu = z2                               # tensile limit line above cross section - downward bending compression
        if zCu < z1: zCu = z1                               # compressive limit below cross ssection - downward bending tension
        if zCu > z2: zCu = z2                               # compressive limit above cross section  - upward bending tension
        # concrete / bulk contribution
        NN, MM, NDE, NDK1, MDK = self.IntCrossSec(epsR, kap, zTu, zCu, Elem)
        NDK2 = NDK1
        if Elem.dim==10:                                            # Bernoulli beam
            MatM = array( [[NDE,NDK1],[NDK2,MDK]] )                 # material tangential stiffness (Script Eq. (3.21))
            sig = array( [NN,MM] )                                  # stress vector
            return sig, MatM, [Eps[0],Eps[1],sig[0],sig[1]]
        else:
            raise NameError("ConFemMaterials::PolyBeam: not implemented for this element type",Elem.Type)

    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):  # loop over integration points
            if Elem.StateVarN[j, 0] > Elem.StateVar[j, 0]: Elem.StateVar[j, 0] = Elem.StateVarN[j, 0]
        return False

class Elastic(Material):                                                    # linear elastic
    def __init__(self, PropMat):
        self.PhaseField = False
        self.Visco = False
        if PropMat[3] > 0.:                                                 # phase field
#                            (self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
            Material.__init__(self, True,         4,        True,      False,     8,           8, 'ELASTIC_PHASEFIELD') # state variables 2 + 6 for strain rate 
            self.Gf  = PropMat[3]                                           # fracture toughness
            self.charLen = PropMat[4]                                       # char. length
            self.eta = PropMat[5]                                           # artificial viscosity
            self.PhaseField = True
        elif PropMat[6] > 0.:                                               # uniaxial visco-elasticity
            self.phi = PropMat[6]                                           # creep number
            self.zeta= 1./PropMat[7]                                           # viscosity parameter
            self.Emod= PropMat[0]
            self.alpha = 0.5                            # parameter for numerical integration with trapezoidal rule
            Material.__init__(self, True,         None,     True,      False,     8,           8, 'ELASTIC_VISCO')
            self.Visco = True
        else:
            Material.__init__(self, True,         None,     False,     False,     None,        8, 'ELASTIC')
            self.PhaseField = False
        self.PropMat = PropMat
        self.alphaT  = PropMat[2]                                           # thermal expansion coefficient
        self.Density = PropMat[8]                                           # specific mass
        self.Dam = False                                                    # flag for damage
    def C3(self, Emod, nu):
        ff = Emod / ( ( 1. + nu ) * ( 1. - 2.*nu ) )
        return ff*array([[1.-nu,nu,nu,0,0,0],[nu,1.-nu,nu,0,0,0],[nu,nu,1.-nu,0,0,0],[0,0,0,(1.-2.*nu)/2.,0,0],[0,0,0,0,(1.-2.*nu)/2.,0],[0,0,0,0,0,(1.-2.*nu)/2.]])
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        Emod = self.PropMat[0]
        nu = self.PropMat[1]
        #
        if self.PhaseField:
            Gf = self.Gf
            charLen = self.charLen
            if Elem.dim==1:
                scaleE = 1. # Emod #1.0e3  # scaling != requires un-symmetric approach
                du_old = Elem.StateVar[ipI,1]
                d      = EpsR[1]                                            # damage variable in integration point
                du     = Eps[0]                                             # strain in integration point
                if abs(du)<abs(du_old): du = du_old                         # in case of unloading keeps damage
                CC   = array([ [ (1-d)*(1-d)*Emod, 0. ] , [ 0., scaleE*Gf*charLen] ])
                CR   = array([ [ 0., -2.*Emod*(1.-d)*du ] , [ scaleE*(-2.*Emod*(1.-d)*du), scaleE*(Emod*du*du + Gf/charLen) ] ])
                sig  = dot(CC,Eps)                                          # Eps should hold derivatives ?of field variables which are currently eps and derivative of d 
                sigR = array([ 0.,  scaleE*(Gf/charLen*d - Emod*(1-d)*du*du) ])      # 
                Elem.StateVarN[ipI,0] = d                                   # presumably not used
                Elem.StateVarN[ipI,1] = du
                # viscous regularization
                if self.eta>0.:
                    Dps_ = array([Dps[0],-nu*Dps[0],-nu*Dps[0],0,0,0])     # voigt notation
                    zz, Veps = self.ViscExten3D( Dt, self.eta, Dps_, Elem, ipI, 2)
                else:
                    zz = 0.
                    Veps = zeros((1),dtype=float)
                sigV = self.eta*Veps[0]
                sig[0]  += sigV
                CC[0,0] += zz
                # 
                return sig, CC, sigR, CR, [Eps[0], sig[0], EpsR[1]]
            else: raise NameError("ConFemMaterials::Elastic.PhaseField: not implemented for this element type")
        #
        elif self.Visco:
            if Elem.dim==1:
                depT = array([self.alphaT*dTmp,0.])                         # increment of imposed strain uniaxial !!!
                epsT = array([self.alphaT*Temp,0.])                         # imposed strain uniaxial !!!
                V = Dt*self.zeta*self.Emod                                  # auxiliary values for creep and relaxation depending on Dt
                W = self.zeta*(1+self.phi)
                W_hat = 1+self.alpha*Dt * W
                W_I = 1/W_hat
                W_bar = W_I*(1-(1-self.alpha)*Dt*W) - 1
                C_bar = W_I*self.Emod
                V_bar = W_I*V
                cc = C_bar + self.alpha*V_bar
                epsP = array( [Elem.DataP[ipI,0], 0.] )
                sigP = array( [Elem.DataP[ipI,1], 0.] )
                MatM = array([[cc,0],[0,0]] )                               # material tangential stiffness
                siV = W_bar*sigP + V_bar*epsP                               # viscoelastic part of stress
                sig = sigP + dot(MatM,(Dps-depT)) + siV                     # total stress
                Elem.StateVarN[ipI,1] = sig[0]
                Elem.StateVarN[ipI,4] = Eps[0]
                return sig, MatM, [ (Eps[0]-epsT[0]), sig[0]]               # ! returns stress inducing strain
            else: raise NameError("ConFemMaterials::Elastic.Visco: not implemented for this element type")
        #
        else:
            if Elem.dim==1 or Elem.dim==99:                                 # uniaxial / spring
                MatM = array([[Emod,0],[0,0]])                              # material stiffness uniaxial
                sig = [Elem.DataP[ipI,1],0] + dot(MatM,Dps)                 # stress incrementally
#                sig = dot(MatM,Eps)
                return sig, MatM, [Eps[0], sig[0]]
            elif Elem.dim==2:
                if Elem.PlSt: MatM = Emod/(1-nu**2)               *array([[1,nu,0],       [nu,1,0],       [0,0,(1-nu)/2]           ]) # plane stress
                else:         MatM = Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),0],[nu/(1-nu),1,0],[0,0,(1-2*nu)/(2*(1-nu))]]) # plane strain
                sig = dot(MatM,Eps)                     # stress
                return sig, MatM, [Eps[0], Eps[1], 0., Eps[2], sig[0], sig[1], 0., sig[2]]
            elif Elem.dim==3: 
                MatM = self.C3( Emod, nu)                                   # triaxial isotropic elasticity
                sig = dot(MatM,Eps)                                         # stress
                return sig, MatM, [Eps[0],Eps[1],Eps[2],Eps[3],Eps[4],Eps[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
            elif Elem.dim==4:                                               # CAX 2D axisymmetric
                MatM = Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),nu/(1-nu),0],
                                                            [nu/(1-nu),1,nu/(1-nu),0],
                                                            [nu/(1-nu),nu/(1-nu),1,0],
                                                            [0,0,0,(1-2*nu)/(2*(1-nu))]])
                sig = dot(MatM, Eps)
                return sig, MatM, [Eps[0],Eps[1],Eps[2],Eps[3], sig[0],sig[1],sig[2],sig[3]]
            elif Elem.dim==5:                                               # 2D bar axisymmetric, should also work for 2D bond elements
                MatM = array([[Emod, 0.], [0., Emod]])
                sig = dot(MatM, Eps)
                return sig, MatM, [Eps[0], Eps[1], sig[0], sig[1]]
            elif Elem.dim in [10,12]:                                       # bernoulli beam
                hh =  Elem.Geom[1,2]                                        # height
                z1, z2 = Elem.zLow, Elem.zUpp                               # lower, upper coordinate
                AA, SS, JJ = Elem.Geom[1,3], Elem.Geom[1,4], Elem.Geom[1,5] 
                Tr = 0.5*(Temp[0]+Temp[1])                                  # temperature of reference / middle axis
                Tg = (Temp[0]-Temp[1])/hh                                   # temperature gradient
                Tps = self.alphaT*array([Tr,Tg])                            # temperature strain
                MatM = array([[Emod*AA,-Emod*SS],[-Emod*SS,Emod*JJ]])       # beam tangential stiffness
                sig = dot(MatM,Eps-Tps)
                eps1= Eps[0] - z1*Eps[1]                                    # lower strain
                eps2= Eps[0] - z2*Eps[1]                                    # upper strain
                return sig, MatM, [Eps[0],Eps[1],sig[0],sig[1],Emod*eps1,Emod*eps2]
            elif Elem.dim in [11,13]:                                       # timoshenko beam, axisymmetric timoshenko beam
                AA = Elem.Geom[1,3]                                         # cross sectional area
                JJ = Elem.Geom[1,5]                                         # moment of inertia
                alpha = 0.8                                                 # factor for shear stiffness
    #            Tr = 0.5*(Temp[0]+Temp[1])              # temperature of reference / middle axis
    #            Tg = (Temp[0]-Temp[1])/hh               # temperature gradient
    #            Tps = self.alphaT*array([Tr,Tg])        # temperature strain
#                MatM = array([[Emod*bb*hh,0,0],[0,Emod*bb*hh**3/12.,0],[0,0,bb*hh*alpha*0.5*Emod/(1+nu)]])# beam tangential stiffness
                MatM = array([[Emod*AA,0,0],[0,Emod*JJ,0],[0,0,AA*alpha*0.5*Emod/(1+nu)]])# beam tangential stiffness
                sig = dot(MatM,Eps)                     #
                return sig, MatM, [Eps[0],Eps[1],Eps[2],sig[0],sig[1],sig[2]]
            elif Elem.dim==20:                                              # Kirchhoff slab
                hh =  Elem.Geom[1,1]                                        # height
#                KK = Emod*hh**3/(12*(1-nu))             # slab stiffness
                KK = Emod*hh**3/(12.*(1.-nu**2))                            # slab stiffness
#                MatM = array([[KK,nu*KK,0],[nu*KK,KK,0],[0,0,(1-nu)*KK]])# slab stiffness. (1-nu)*KK] should be strictly divided by 2 and B_xy doubled, see SB3.FormB.
                MatM = array([[KK,nu*KK,0],[nu*KK,KK,0],[0,0,0.5*(1-nu)*KK]])
                sig = dot(MatM,Eps)                 #
                return sig, MatM, [Eps[0],Eps[1],Eps[2],sig[0],sig[1],sig[2]]
            elif Elem.dim==21:                                              # Continuum based shell
                MatM = array([[Emod/(1-nu**2), nu*Emod/(1-nu**2), 0., 0., 0., 0.],
                              [nu*Emod/(1-nu**2), Emod/(1-nu**2), 0., 0., 0., 0.],
                              [0., 0., 0., 0., 0., 0.],
                              [0., 0., 0., Emod/(2+2*nu), 0., 0.],
                              [0., 0., 0., 0., Emod/(2+2*nu), 0],
                              [0., 0., 0., 0., 0., Emod/(2+2*nu)]]) 
                sig = dot(MatM,Eps)                     #
                if Elem.ShellRCFlag: return sig, MatM, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps[0],Eps[1],Eps[5], sig[0],sig[1],sig[5], 0., 0.]
                else:                return sig, MatM, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]] # used to calculate internal forces --> DataOutStress, WriteElemData
            elif Elem.dim in [97,96]:                                       # 2D bond elements
                MatM = array([[Emod, 0.], [0., nu*Emod]])                   # nu is not poisson's ratio here but multoplier for lateral stiffness
                sig = dot(MatM, Eps)  #
                return sig, MatM, [Eps[0], Eps[1], sig[0], sig[1]]
            elif Elem.dim==98:                                              # 2D spring with 3 dofs -- looks like spring only, see also ElasticOrtho
                MatM = array([[0.,0.,0.],[0.,0.,0.],[0.,0.,Emod]])
                sig = dot(MatM,Eps)                     #
                return sig, MatM, [sig[0],sig[1],sig[2],Eps[0],Eps[1],Eps[2]]
            else: raise NameError ("ConFemMaterials::Elastic.Sig: not implemented for this element type")
    def UpdateStateVar(self, Elem, ff):                                     # for phase field
        if self.PhaseField:
            for j in range(Elem.StateVar.shape[0]):
                if Elem.StateVarN[j,1]>Elem.StateVar[j,1]:                  # damage measures
                    Elem.StateVar[j,0] = Elem.StateVarN[j,0]
                    Elem.StateVar[j,1] = Elem.StateVarN[j,1]

class ElasticR1D(Material):                                                 # linear elastic axisymmetric in radial direction
    def __init__(self, PropMat):
        Material.__init__(self, True, None, False, False, None, 8, 'ELASTIC1DR')
        self.PropMat = PropMat
        self.alphaT = PropMat[2]                                            # thermal expansion coefficient
        self.Density = PropMat[3]                                           # specific mass
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        Emod = self.PropMat[0]
        if Elem.dim == 5:                                                   # 2D bar axisymmetric
            MatM = array([[Emod, 0.], [0., 0.]])
            sig = dot(MatM, Eps)
            return sig, MatM, [Eps[0], Eps[1], sig[0], sig[1]]
        else: raise NameError ("ConFemMaterials::ElasticR1D.Sig: not implemented for this element type")
class ElasticC1D(Material):                                                 # linear elastic axisymmetric in circumferential direction
    def __init__(self, PropMat):
        Material.__init__(self, True, None, False, False, None, 8, 'ELASTICC1D')
        self.PropMat = PropMat
        self.alphaT = PropMat[2]                                            # thermal expansion coefficient
        self.Density = PropMat[3]                                           # specific mass
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        Emod = self.PropMat[0]
        if Elem.dim in [5]:                                                 # 2D bar axisymmetric
            MatM = array([[0., 0.], [0., Emod]])
            sig = dot(MatM, Eps)
            return sig, MatM, [Eps[0], Eps[1], sig[0], sig[1]]
        else:
            raise NameError ("ConFemMaterials::ElasticR1D.Sig: not implemented for this element type")

class ElasticOrtho(Material):                           # linear orthotropic elastic
    def __init__(self, PropMat):
        Material.__init__(self, True, None, False, False, None, 6, "ElasticOrtho")
        self.E1   = PropMat[0]
        self.E2   = PropMat[1]
        self.E3   = PropMat[2]
        self.nu12 = PropMat[3]
        self.nu13 = PropMat[4]
        self.nu23 = PropMat[5]
        self.G1   = PropMat[6]
        self.G2   = PropMat[7]
        self.G3   = PropMat[8]
        self.Density   = PropMat[9]                 # specific density
        if self.E1<ZeroD or self.E2<ZeroD: raise NameError ("ConFemMat::ElasticOrtho.__init__: E1 and E2 should be larger than zero")
        self.nu21 = self.nu12*self.E2/self.E1           # to have symmetry
        self.nu31 = self.nu13*self.E3/self.E1
        self.nu32 = self.nu23*self.E3/self.E2
        self.ff = -1./(-1+self.nu32*self.nu23+self.nu21*self.nu12+self.nu21*self.nu32*self.nu13+self.nu31*self.nu12*self.nu23+self.nu31*self.nu13) 
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if Elem.dim==3:
            MatM = self.ff*array([[(1.-self.nu32*self.nu23)*self.E1,(self.nu12+self.nu32*self.nu13)*self.E2,(self.nu12*self.nu23+self.nu13)*self.E3,0,0,0],
                                  [(self.nu21+self.nu31*self.nu23)*self.E1,(1.-self.nu31*self.nu13)*self.E2,(self.nu23+self.nu21*self.nu13)*self.E3,0,0,0],
                                  [(self.nu21*self.nu32+self.nu31)*self.E1,(self.nu32+self.nu31*self.nu12)*self.E2,(1.-self.nu21*self.nu12)*self.E3,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,0,0],
                                  [0,0,0,0,0,0]])
            MatM[3,3], MatM[4,4], MatM[5,5] = self.G1, self.G2, self.G3
            sig = dot(MatM,Eps)                     # stress
            return sig, MatM, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
#            return sig, MatM, [Eps[0],Eps[1],Eps[2],Eps[3],Eps[4],Eps[5]]
        elif Elem.dim==98:
            MatM = array([[self.E1,0.,0.],[0.,self.E2,0.],[0.,0.,self.G3]])
            sig = dot(MatM,Eps)                     #
            return sig, MatM, [sig[0],sig[1],sig[2],Eps[0],Eps[1],Eps[2]]
        else: raise NameError ("ConFemMaterials::ElasticOrtho.Sig: not implemented for this element type")

class ElasticLT(Material):                          # linear elastic with limited tensile strength -- see also ConFemNotes
    def __init__(self, PropMat):
#                        (self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, False,        None,     True,      True,      15,          10, 'ELASTICLT')
#        self.Symmetric = False                       # flag for symmetry of material matrices
#        self.Update = True                          # has specific update procedure
#        self.Updat2 = True                          # two stage update procedure
#        self.StateVar = 15                          # number of state variables per integration point
                                                    # [0] state, [1] max. prin stress,
                                                    # [2] largest crack width reached ever, [3] crack traction related to [2], 
                                                    # [4] max prin strain, [5] direction corresponding to 1st crack upon initiation
                                                    # [6] state2, [7] largest crack width reached ever 2, [8] crack traction related to [7] 2
                                                    # [9] principal strain 0 of previous step, [10] principal strain 1 of previous step, 
                                                    # [11] crack width 0 prev. step, [12] crack width 1 prev. step --  [13] crack width velocity 0 prev. step, [14] crack width velocity 1 prev. step 
#        self.NData = 8                              # number of data items
        self.PropMat = PropMat                      # for wrappers, e.g. RCSHELL
        self.Emod = PropMat[0]
        self.nu = PropMat[1]
        self.fct = PropMat[2]                       # uniaxial tensile strength
        self.Gf = PropMat[3]                        # fracture energy
        self.bw = PropMat[4]                        # crack bandwidth
        self.epsct = PropMat[2]/PropMat[0]          # uniaxial strain at tensile strength
        if self.fct>ZeroD: self.wcr = 2*self.Gf/self.fct           # critical crack width, may be reset by ConfemRandom
        else:              self.wcr=0
#            if self.epsct*self.bw/self.wcr > 0.01: raise NameError("ElasticLT.init: incoherent material parameters")
        self.phi = 0. #PropMat[5]                       # creep number
        self.zeta = 0. #PropMat[6]                      # viscosity parameter
        self.cCrit = 0 # PropMat[7]                     # type of criterion for tensile failure 0: stress, 1: strain
        self.CrackVisc = PropMat[6] # PropMat[8]                 # crack viscosity
        self.rho = PropMat[5] #PropMat[9]                       # related residual strength with full crack
        self.Density = PropMat[7] # PropMat[10]                  # specific mass
        if self.cCrit == 0: self.cVal = 0. #self.fct   
        else:               self.cVal = self.epsct  # may be reset by ConfemRandom
        self.alpha = 0.5                            # parameter for numerical integration with trapezoidal rule
        self.alphaT = 1.e-5                         # thermal expansion coefficient -- should currently not be used for *ELASTICLT

        self.sigM = 0.                              # determine max tensile stress/strain for attached element set
        self.Elem = None                            # " -> corresponding element index
        self.IPoint = None                          # " -> corresponding integration point index
        self.sigM2 = 0.                             # determine max tensile stress/strain for attached element set
        self.Elem2 = None                           # " -> corresponding element index
        self.IPoint2 = None                         # " -> corresponding integration point index
        self.NoCracks=0                             # Number of cracks
        
        alpha        = 10.                                                  # this is overridden for dim=1 with prescribing bw and GF
        self.epscu   = alpha*self.epsct
        self.epsDelt = (alpha-1.)*self.epsct
        self.bw_     = 2.*self.Gf/(self.fct*self.epsDelt)                   # crack band width from Book2nd (7.56) with prescription of alpha and Gf 
        self.wcr_    = 2.*self.Gf/self.fct * alpha/(alpha-1.)
        
#    @profile
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        def State1C(i, ww ):                                # single crack - loading
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _MatML[0,0] = -self.fct/ds1*(-etaH*epsu+1+etaH*epst)
            _MatML[0,1] = -self.fct/ds1*(etaH*epst*nu-etaH*epsu*nu-etaH*epst*xi*nu+nu+etaH*epsu*xi*nu-xi*nu)
            _MatML[1,0] = -Emod/ds1*(-nu*etaH*epsu*epst+nu*etaH*epst2+nu*epst)
            _MatML[1,1] = -Emod/ds1*(-xi*epsu+etaH*epsu*xi*epst-etaH*epst2*xi+etaH*epst2-etaH*epsu*epst+epst)
            _sigL[0] = -self.fct/ds1*(xi*etaH*eps_c11_k*epsu+epsL[0]-xi*epsu+etaH*epsu*xi*nu*epsL[1]+etaH*epst*epsL[0]+etaH*epst*nu*epsL[1]+xi*etaB*deps_c11_k*epsu-etaH*epsu*nu*epsL[1]-xi*etaB*deps_c11_k*epst-etaH*epsu*epsL[0]+nu*epsL[1]-etaH*epst*xi*nu*epsL[1]-xi*nu*epsL[1]-xi*etaH*eps_c11_k*epst) 
            _sigL[1] = -Emod/ds1*(epsL[1]*epst+nu*epsL[0]*epst-nu*epsu*xi*epst-epsL[1]*xi*epsu+epsL[1]*etaH*epst2-epsL[1]*etaH*epsu*epst+nu*xi*etaH*eps_c11_k*epsu*epst+nu*etaH*epst2*epsL[0]+nu*xi*etaB*deps_c11_k*epsu*epst-nu*xi*etaB*deps_c11_k*epst2-nu*etaH*epsu*epsL[0]*epst-nu*xi*etaH*eps_c11_k*epst2+epsL[1]*etaH*epsu*xi*epst-epsL[1]*etaH*epst2*xi) 
            Elem.StateVarN[ipI,i+2] = ww
            Elem.StateVarN[ipI,i+3] = _sigL[i]
            return 1, ww, _sigL, _MatML
        def State2C( ww ):                                  # single crack - unloading
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _MatML[0,0] = -self.fct/ds2 * (-tP1*self.bw_-self.fct*wP1*etaH) 
            _MatML[0,1] = -self.fct/ds2 * (-self.fct*wP1*etaH*nu+self.fct*wP1*etaH*xi*nu+self.bw_*tP1*xi*nu-nu*tP1*self.bw_)
            _MatML[1,0] = -Emod/ds2     * (-nu*self.bw_*tP1*epst-nu*self.fct*wP1*etaH*epst) 
            _MatML[1,1] = -Emod/ds2     * (-xi*self.fct*wP1-tP1*self.bw_*epst+self.bw_*tP1*xi*epst+self.fct*wP1*etaH*xi*epst-self.fct*wP1*etaH*epst) 
            _sigL[0] = -self.fct/ds2 * (-self.bw_*epsL[0]*tP1+xi*self.fct*wP1*etaH*eps_c11_k+xi*self.fct*wP1*etaB*deps_c11_k-self.fct*wP1*etaH*nu*epsL[1]-self.fct*wP1*etaH*epsL[0]+self.fct*wP1*etaH*xi*nu*epsL[1]+self.bw_*tP1*xi*nu*epsL[1]-tP1*self.bw_*nu*epsL[1])
            _sigL[1] = -Emod/ds2     * (-nu*self.bw_*epsL[0]*tP1*epst+nu*xi*self.fct*wP1*etaH*eps_c11_k*epst+nu*xi*self.fct*wP1*etaB*deps_c11_k*epst-nu*self.fct*wP1*etaH*epsL[0]*epst-epsL[1]*self.fct*xi*wP1-epsL[1]*tP1*self.bw_*epst+epsL[1]*self.bw_*tP1*xi*epst+epsL[1]*self.fct*wP1*etaH*xi*epst-epsL[1]*self.fct*wP1*etaH*epst)
            return 2, ww, _sigL, _MatML
        def State3C():                                      # single crack - closure
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            xxx = Emod/(1-nu**2)
            _MatML[0,0] = xxx
            _MatML[0,1] = xxx*nu
            _MatML[1,0] = xxx*nu
            _MatML[1,1] = xxx 
            _sigL[0] = _MatML[0,0]*epsL[0] + _MatML[0,1]*epsL[1]
            _sigL[1] = _MatML[1,0]*epsL[0] + _MatML[1,1]*epsL[1]
            return 3, 0., _sigL, _MatML
        def State4C( ww ):                                  # single cracking break through after softening
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _sigL[0] = -self.fct/ds4 * (xi*etaH*eps_c11_k+xi*etaB*deps_c11_k-etaH*nu*epsL[1]-epsL[0]*etaH+etaH*xi*nu*epsL[1]-xi*rho)
            _MatML[0,0] = self.fct*etaH                   /ds4 #(xi+etaH*epst-etaH*nu2*epst-etaH*xi*epst+etaH*xi*nu2*epst)
            _MatML[0,1] = self.fct*(etaH*nu-etaH*xi*nu)   /ds4 #(xi+etaH*epst-etaH*nu2*epst-etaH*xi*epst+etaH*xi*nu2*epst)
            _sigL[1] = -Emod/ds4     * (nu*xi*etaH*eps_c11_k*epst+nu*xi*etaB*deps_c11_k*epst-nu*epsL[0]*etaH*epst-epsL[1]*xi-epsL[1]*etaH*epst+epsL[1]*etaH*xi*epst-nu*xi*rho*epst)
            _MatML[1,0] = self.fct*nu*etaH                /ds4 # (xi+etaH*epst-etaH*nu2*epst-etaH*xi*epst+etaH*xi*nu2*epst)
            _MatML[1,1] = Emod*(xi+etaH*epst-etaH*xi*epst)/ds4 # (xi+etaH*epst-etaH*nu2*epst-etaH*xi*epst+etaH*xi*nu2*epst)
            Elem.StateVarN[ipI,2] = ww
            return  4, ww, _sigL, _MatML
        # dual cracking - assumes that first crack is in state 4 
        def State1C2( ww ):                                 # dual cracking - loading - 1st crack with full crack
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _sigL[0]    =  self.fct/dd4 * (xi*etaH*eps_c11_k+xi*etaB*deps_c11_k-etaH*epsL[0]-xi*rho)
            _sigL[1]    = -self.fct/dd1 * (xi*epsu-xi*etaH*eps_c21_k*epsu+xi*etaH*eps_c21_k*epst-xi*etaB*deps_c21_k*epsu+xi*etaB*deps_c21_k*epst-epsL[1]-etaH*epst*epsL[1]+etaH*epsu*epsL[1])
            _MatML[0,0] = -self.fct/dd4 *  etaH             # !!!!! dd4 presumably negative
            _MatML[1,1] = -self.fct/dd1 * (-1-etaH*epst+etaH*epsu)
            Elem.StateVarN[ipI,7] = ww
            Elem.StateVarN[ipI,8] = _sigL[1]
            return 1, ww, _sigL, _MatML
        def State2C2( ww ):                                 # dual cracking - unloading
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _sigL[0]    =  self.fct/dd4 * (xi*etaH*eps_c11_k+xi*etaB*deps_c11_k-etaH*epsL[0]-xi*rho)
            _sigL[1]    =  self.fct/dd2 * (-self.bw_*tP2*epsL[1]+xi*self.fct*wP2*etaH*eps_c11_k+xi*self.fct*wP2*etaB*deps_c11_k-self.fct*wP2*etaH*epsL[1])
            _MatML[0,0] = -self.fct/dd4 *  etaH
            _MatML[1,1] =  self.fct/dd2 * (-tP2*self.bw_-self.fct*wP2*etaH)
            return 2, ww, _sigL, _MatML
        def State3C2():                                     # dual cracking - crack closure
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _sigL[0]    =  self.fct/dd4 * (xi*etaH*eps_c11_k+xi*etaB*deps_c11_k-etaH*epsL[0]-xi*rho)
            _sigL[1]    =  MatML[1,1]*epsL[1]
            _MatML[0,0] = -self.fct/dd4 *  etaH
            _MatML[1,1] =  Emod 
            return 3, 0., _sigL, _MatML
        def State4C2( ww ):                                 # dual cracking - full crack
            _MatML = zeros((3,3), dtype=float)
            _sigL = zeros((3), dtype=float)
            _sigL[0]    =  self.fct/dd4 * (xi*etaH*eps_c11_k+xi*etaB*deps_c11_k-etaH*epsL[0]-xi*rho)
            _sigL[1]    =  self.fct/dd4 * (xi*etaH*eps_c21_k+xi*etaB*deps_c21_k-etaH*epsL[1]-xi*rho)
            _MatML[0,0] = -self.fct/dd4 *  etaH
            _MatML[1,1] =  _MatML[0,0]
            Elem.StateVarN[ipI,7] = ww
            return  4, ww, _sigL, _MatML

        self.fct=RandomField_Routines().get_property(Elem.Set,Elem.Label,self.fct)[0]
        self.wcr,self.cVal=RandomField_Routines().more_properties(self.fct,self.Gf,self.cCrit,self.epsct,self.wcr,self.cVal,Elem.Set,Elem.Label)
                                                                            # cVal presumably random variable around 0 for tensile strength around fct
        depT = array([self.alphaT*dTmp,0.])                                 # increment of imposed strain uniaxial !!!
        epsT = array([self.alphaT*Temp,0.])                                 # imposed strain uniaxial !!!
        if CalcType==0:
            if self.cCrit==0: sip = Elem.StateVarN[ipI,1]-self.fct          # retrieve largest principal stress distance to tensile strength, regards randomness
            else:             sip = Elem.StateVarN[ipI,4]                   # retrieve largest principal strain  -- this should be checked again uhc 22-01-25
            if Elem.StateVarN[ipI,0]==0:                                    # check for 1st crack
                if sip>self.sigM:                                           # determine point with largest principal stress/strain
                    self.sigM = sip
                    self.Elem = elI
                    self.IPoint = ipI
            elif self.cCrit==0:                                             # check for 2nd crack in case 1st crack is there / stress based only !!!
#                if sip>(self.sigM2-self.fct):                        # determine point with largest principal stress
                if sip>(self.sigM2-0.0*self.fct):         # 0.5               # determine point with largest principal stress
                    self.sigM2 = sip
                    self.Elem2 = elI
                    self.IPoint2 = ipI
            return [], [], []
        else:
            # first crack check CalcType 1
            if self.Elem==elI and self.IPoint==ipI and self.sigM>self.cVal and Elem.StateVarN[ipI,0]==0:
                Elem.StateVarN[ipI,0] = 1                                   # first crack starts
                self.NoCracks = self.NoCracks+1
                Echo(f"El {Elem.Label:d} ip {ipI:d} ElasticLT element cracks  fct {self.fct:f}, delSig {self.sigM:f}, no cracks {self.NoCracks:d}", ff)
                self.sigM = 0.
                self.Elem = None
                self.IPoint = None
            # second crack check
            if Elem.StateVarN[ipI,0]==4 and Elem.StateVarN[ipI,6]==0 and self.Elem2==elI and self.IPoint2==ipI and self.sigM2>self.cVal:
                Elem.StateVarN[ipI,6] = 1                                   # second crack starts
                self.NoCracks = self.NoCracks+1
                Echo(f"El {Elem.Label:d} ip {ipI:d} ElasticLT element cracks twice fct {self.fct:f}, delSig {self.sigM2:f}, no cracks {self.NoCracks:d}, {Elem.StateVarN[ipI,1]:f}", ff)
                self.sigM2 = 0.
                self.Elem2 = None
                self.IPoint2 = None
            # 
            if Elem.dim==1:                                                 # uniaxial -- visco elastic has gone to *ELASTIB option VISCO
                eps = Eps[0]
                eps_ct = self.epsct
                xi = self.bw / Elem.Lch_                                    # bw different to bw_ used below
                eps_cu = (2.*self.Gf*self.Emod/(self.fct**2) / (xi*Elem.Lch_) +1) * eps_ct
                State = int(Elem.StateVarN[ipI,0])                          # current state of actual time step
                if State>0:
                    if self.Gf > 0.:
                        alpha = 2*self.Gf*self.Emod/self.fct**2/self.bw + 1 # factor for eps_cu - with prescription of bw and Gf
                        self.epscu   = alpha*self.epsct
                        self.epsDelt = (alpha-1.)*self.epsct
                        self.bw_     = 2.*self.Gf/(self.fct*self.epsDelt)   # crack band width from Book2nd (7.56)
                        self.wcr_    = 2.*self.Gf/self.fct * alpha/(alpha-1.)
                        Emod = self.Emod
                        rho = self.rho                 # related residual strength with full crack
                        TOL = -1.e-15                                       # 
                        wcr  = self.wcr_    - rho*(self.wcr_ - self.bw_*self.epsct)                                # effective critical crack width
                        xi   = self.bw_/Elem.Lch_                           # related crack band width
                        etaB = self.CrackVisc/self.fct                      # related crack viscosity
                        if Dt>ZeroD: etaH = 2.*etaB/Dt                                   # related crack viscosity 2
                        else:        etaH = 0.
                        epsu = self.epscu                                   # failure strain
                        epst = self.epsct                                   # tensile strength strain
                        epst2= epst*epst
                        # state of previous time step
                        StateP1 = Elem.StateVar[ipI,0]                      # final state of last time step 1st crack
                        wP1 = Elem.StateVar[ipI,2]                          # largest crack ever reached of last time step first crack
                        tP1 = Elem.StateVar[ipI,3]                          # corresponding crack traction 
                        epsLP = array([Elem.StateVar[ipI,9],Elem.StateVar[ipI,10],0.]) # local strain previous step
                        # predictors crack width 1st crack
                        nu = 0.
                        nu2 = 0.
                        ds1 = xi*epsu+nu2*epst-etaH*epst2+etaH*epsu*epst-etaH*epsu*xi*epst+etaH*epst2*xi-xi*nu2*epst+etaH*epsu*xi*nu2*epst+etaH*epst2*nu2-etaH*epsu*nu2*epst-etaH*epst2*xi*nu2-epst
                        ds1_= xi*epsu                                                                                               -epst
                        ds2 = self.fct*wP1*xi+tP1*self.bw_*epst-tP1*self.bw_*nu2*epst-self.bw_*tP1*xi*epst+self.bw_*tP1*xi*nu2*epst-self.fct*wP1*etaH*xi*epst+self.fct*wP1*etaH*xi*nu2*epst+self.fct*wP1*etaH*epst-self.fct*wP1*etaH*nu2*epst
                        ds4 = xi+etaH*epst-nu2*etaH*epst-xi*etaH*epst+nu2*xi*etaH*epst
                        eps_c11_k = Elem.StateVar[ipI,11]/self.bw_          # 1st crack strain of previous time step for viscous contributions
                        deps_c11_k= Elem.StateVar[ipI,13]/self.bw_          # 1st crack strain velocity of previous time step for viscous contributions
                        epsL0 = eps
                        epsL1 = 0.
                        epsL = array([eps,0.])
                        ww1_1 = self.bw_ /ds1*(-epsu*epst-eps*epst+eps*epsu-nu2*etaH*eps_c11_k*epsu*epst-xi*etaH*eps_c11_k*epsu*epst-xi*etaB*deps_c11_k*epsu*epst+nu2*xi*etaH*eps_c11_k*epsu*epst-epsu*xi*nu2*epst+nu2*xi*etaB*deps_c11_k*epsu*epst-nu*epsL1*epst+xi*nu*epsL1*epst+nu*epsL1*epsu+nu2*etaH*eps_c11_k*epst2-etaH*eps_c11_k*epst2-nu*epsL1*xi*epsu-etaB*deps_c11_k*epst2+nu2*etaB*deps_c11_k*epst2+xi*etaH*eps_c11_k*epst2+xi*etaB*deps_c11_k*epst2-nu2*xi*etaB*deps_c11_k*epst2-nu2*xi*etaH*eps_c11_k*epst2+epsu*nu2*epst+xi*epsu*epst-nu2*etaB*deps_c11_k*epsu*epst+etaH*eps_c11_k*epsu*epst+etaB*deps_c11_k*epsu*epst)
                        if abs(ds2)>ZeroD: ww1_2 = self.bw_*wP1*self.fct/ds2 *(eps_c11_k*etaH*epst+deps_c11_k*etaB*epst-nu2*etaB*deps_c11_k*epst-nu2*etaH*eps_c11_k*epst+nu*epsL1+nu2*xi*etaB*deps_c11_k*epst+nu2*xi*etaH*eps_c11_k*epst-nu*epsL1*xi-xi*etaH*eps_c11_k*epst-xi*etaB*deps_c11_k*epst+epsL0)
                        else:              ww1_2 = 0.
                        ww1_4 = self.bw_/ds4 *(eps_c11_k*etaH*epst+deps_c11_k*etaB*epst-nu2*etaH*eps_c11_k*epst-nu2*etaB*deps_c11_k*epst+nu*epsL1-xi*etaH*eps_c11_k*epst-xi*etaB*deps_c11_k*epst+nu2*xi*etaH*eps_c11_k*epst+nu2*xi*etaB*deps_c11_k*epst-nu*epsL1*xi+epsL0+rho*epst*(nu2-1+xi-nu2*xi))
#                        
                        if ww1_1>wcr:
                            State_,ww1,sigL,MatML = State4C(  ww1_4)        # starting state 4 loading beyond critical crack width   
#                        elif eps>epsLP[0]+TOL and ww1_1>wP1:
                        elif  ww1_1>wP1:              
                            State_,ww1,sigL,MatML = State1C(0,ww1_1)        # state 1 loading 
                        elif ww1_2>0:                                     
                            State_,ww1,sigL,MatML = State2C(  ww1_2)        # state 2 unloading
                        elif ww1_2<=0:                                    
                            State_,ww1,sigL,MatML = State3C()               # state 3 crack closure
#
                        Elem.StateVarN[ipI,0]  = State_
                        Elem.StateVarN[ipI,1]  = max(sigL[0],sigL[1])       # used for both 1st and 2nd crack in CalcType 0    
                        Elem.StateVarN[ipI,11] = ww1                        # current crack width 1
                        if Dt>ZeroD: Elem.StateVarN[ipI,13] = 2.*(ww1-Elem.StateVar[ipI,11])/Dt - Elem.StateVar[ipI,13]  # current crack width velocity 1 
                        else:        Elem.StateVarN[ipI,13] = 0.
                        sig = [sigL[0],0.]
                        MatM= [[MatML[0,0],0.],[0.,0.]]
                        sig_, ww = sig[0], ww1
#                         
                        Eps[0] = eps
                    else:
                        MatM = array([[0,0],[0,0]])                             # material tangential stiffness
                        Eps[0] = 0
                        sig = array([0,0], dtype=double)
                else:                                                       # uncracked state
                    cc = self.Emod #C_bar #+ self.alpha*V_bar
                    epsP = array( [Elem.DataP[ipI,0], 0.] )
                    sigP = array( [Elem.DataP[ipI,1], 0.] )
                    MatM = array([[cc,0],[0,0]] )                           # material tangential stiffness
                    sig = sigP + dot(MatM,(Dps-depT)) # + siV               # total stress
                    Elem.StateVarN[ipI,1] = sig[0]
                    Elem.StateVarN[ipI,4] = Eps[0]
                ##
                return sig, MatM, [ (Eps[0]-epsT[0]), sig[0]]               # ! returns stress inducing strain
            #
            elif Elem.dim==2 or Elem.dim==21:
                if not ConFemMatCFlag: #False: #False: #False: #True: # flag for C-version
#                if True:
                    if Elem.dim==2:
                        if not Elem.PlSt: raise NameError("ConFemMaterials::ElasticLT.sig: ElasticLT not yet defined for plane strain")
                        Eps_ = Eps
                    elif Elem.dim==21: 
                        Eps_ = array([Eps[0],Eps[1],Eps[5]])
                    nu = self.nu                                            # Poisson's ratio
                    Emod = self.Emod
                    ww1, ww2 = 0., 0.                                       # initial value crack width
                    pep,phe,pep_ = PrinCLT_( Eps_[0], Eps_[1], 0.5*Eps_[2]) # principal strains, largest value, corr. direction, lower value
                    Elem.StateVarN[ipI,9] = pep                             # 1st larger principal strain - updated in UpdateStat2Var
                    Elem.StateVarN[ipI,10]= pep_                            # 2nd principal strain - updated in UpdateStat2Var
                    State = int(Elem.StateVarN[ipI,0])                      # current state of actual time step
#                    if State == 0: Elem.StateVarN[ipI,5] = phe              # save for updated potential 1st crack direction -- presumably not needed currently  uhc 22-01-25
                    # cracked state
                    if State > 0:
                        rho = self.rho # 0.005 #0.02 # 0.01 #0.                # related residual strength with full crack
                        TOL = -1.e-15                                       # 
                        State2 = int(Elem.StateVarN[ipI,6])                 # current state 2 of actual time step
                        if State2 > 0: nu=0.
                        # transformation to local system
                        Trans=array([[cos(phe)**2,sin(phe)**2, cos(phe)*sin(phe)],
                                     [sin(phe)**2,cos(phe)**2,-cos(phe)*sin(phe)],
                                     [-2*cos(phe)*sin(phe),2*cos(phe)*sin(phe),cos(phe)**2-sin(phe)**2]])# transformation matrix for strains from global to local
                        epsL  = dot(Trans,Eps_)                             # local strain
                        # state of previous time step
                        StateP1 = Elem.StateVar[ipI,0]                      # final state of last time step 1st crack
                        wP1 = Elem.StateVar[ipI,2]                          # largest crack ever reached of last time step first crack
                        tP1 = Elem.StateVar[ipI,3]                          # corresponding crack traction 
                        StateP2 = Elem.StateVar[ipI,6]                      # final state of last time step 2nd crack
                        wP2 = Elem.StateVar[ipI,7]                          # largest crack ever reached of last time step first crack
                        tP2 = Elem.StateVar[ipI,8]                          # corresponding crack traction
                        epsLP = array([Elem.StateVar[ipI,9],Elem.StateVar[ipI,10],0.]) # local strain previous step
                        # auxiliary values
                        wcr  = self.wcr_    - rho*(self.wcr_ - self.bw_*self.epsct)                                # effective critical crack width
                        xi   = self.bw_/Elem.Lch_                           # related crack band width
                        etaB = self.CrackVisc/self.fct                      # related crack viscosity
                        etaH = 2.*etaB/Dt                                   # related crack viscosity 2
                        epsu = self.epscu                                   # failure strain
                        epst = self.epsct                                   # tensile strength strain
                        nu2  = nu*nu
                        epst2= epst*epst
                        # auxiliary values - 1st crack
                        ds1 = xi*epsu+nu2*epst-etaH*epst2+etaH*epsu*epst-etaH*epsu*xi*epst+etaH*epst2*xi-xi*nu2*epst+etaH*epsu*xi*nu2*epst+etaH*epst2*nu2-etaH*epsu*nu2*epst-etaH*epst2*xi*nu2-epst
                        ds2 = self.fct*wP1*xi+tP1*self.bw_*epst-tP1*self.bw_*nu2*epst-self.bw_*tP1*xi*epst+self.bw_*tP1*xi*nu2*epst-self.fct*wP1*etaH*xi*epst+self.fct*wP1*etaH*xi*nu2*epst+self.fct*wP1*etaH*epst-self.fct*wP1*etaH*nu2*epst
                        ds4 = xi+etaH*epst-nu2*etaH*epst-xi*etaH*epst+nu2*xi*etaH*epst
                        eps_c11_k = Elem.StateVar[ipI,11]/self.bw_          # 1st crack strain of previous time step for viscous contributions
                        deps_c11_k= Elem.StateVar[ipI,13]/self.bw_          # 1st crack strain velocity of previous time step for viscous contributions
                        # predictors crack width 1st crack
                        ww1_1 = self.bw_ /ds1*(-epsu*epst-epsL[0]*epst+epsL[0]*epsu-nu2*etaH*eps_c11_k*epsu*epst-xi*etaH*eps_c11_k*epsu*epst-xi*etaB*deps_c11_k*epsu*epst+nu2*xi*etaH*eps_c11_k*epsu*epst-epsu*xi*nu2*epst+nu2*xi*etaB*deps_c11_k*epsu*epst-nu*epsL[1]*epst+xi*nu*epsL[1]*epst+nu*epsL[1]*epsu+nu2*etaH*eps_c11_k*epst2-etaH*eps_c11_k*epst2-nu*epsL[1]*xi*epsu-etaB*deps_c11_k*epst2+nu2*etaB*deps_c11_k*epst2+xi*etaH*eps_c11_k*epst2+xi*etaB*deps_c11_k*epst2-nu2*xi*etaB*deps_c11_k*epst2-nu2*xi*etaH*eps_c11_k*epst2+epsu*nu2*epst+xi*epsu*epst-nu2*etaB*deps_c11_k*epsu*epst+etaH*eps_c11_k*epsu*epst+etaB*deps_c11_k*epsu*epst)
                        if abs(ds2)>ZeroD: ww1_2 = self.bw_*wP1*self.fct/ds2 *(eps_c11_k*etaH*epst+deps_c11_k*etaB*epst-nu2*etaB*deps_c11_k*epst-nu2*etaH*eps_c11_k*epst+nu*epsL[1]+nu2*xi*etaB*deps_c11_k*epst+nu2*xi*etaH*eps_c11_k*epst-nu*epsL[1]*xi-xi*etaH*eps_c11_k*epst-xi*etaB*deps_c11_k*epst+epsL[0])
                        else:              ww1_2 = 0.
                        ww1_4 = self.bw_/ds4 *(eps_c11_k*etaH*epst+deps_c11_k*etaB*epst-nu2*etaH*eps_c11_k*epst-nu2*etaB*deps_c11_k*epst+nu*epsL[1]-xi*etaH*eps_c11_k*epst-xi*etaB*deps_c11_k*epst+nu2*xi*etaH*eps_c11_k*epst+nu2*xi*etaB*deps_c11_k*epst-nu*epsL[1]*xi+epsL[0]+rho*epst*(nu2-1+xi-nu2*xi))
                        #
                        State_ = None                                       # initial value current states
                        State2_ = 0
                        # actual state of cracking depending on previous state and crack width predictors 
                        if StateP1>=4:                                      # state 4 or 5 - maybe single or dual cracking as direction 1 decoupled from 2 for the following (nu = 0)
                            if ww1_4>0:                                     # state 4 open crack loading or unloading 
                                State_,                                               ww1,sigL,MatML = State4C(ww1_4)          # 
                                if State2>0:                                # dual cracking crack 2
                                    # auxiliary values -  2nd crack
                                    dd1 = -xi*epsu-etaH*epsu*epst+epst+etaH*epst2-etaH*epst2*xi+etaH*epsu*xi*epst
                                    dd2 = -tP2*self.bw_*epst+self.bw_*tP2*xi*epst-self.fct*wP2*xi+self.fct*wP2*etaH*xi*epst-self.fct*wP2*etaH*epst 
                                    dd4 = -xi+etaH*epst*xi-etaH*epst
                                    eps_c21_k = Elem.StateVar[ipI,12]/self.bw_          # 2nd crack strain of previous time step for viscous contributions
                                    deps_c21_k= Elem.StateVar[ipI,14]/self.bw_          # 2nd crack strain velocity of previous time step for viscous contributions
                                    # predictors for crack width 2nd crack
                                    ww2_1 = -self.bw_/dd1*(-epsu*epst+etaH*eps_c21_k*epsu*epst-etaH*eps_c21_k*epst2+etaB*deps_c21_k*epsu*epst-etaB*deps_c21_k*epst2+epsu*xi*epst-xi*etaH*eps_c21_k*epsu*epst+xi*etaH*eps_c21_k*epst2-xi*etaB*deps_c21_k*epsu*epst+xi*etaB*deps_c21_k*epst2-epst*epsL[1]+epsu*epsL[1])
                                    if abs(dd2)>ZeroD: ww2_2 = self.bw_*wP2*self.fct/dd2 *(-etaH*eps_c21_k*epst-etaB*deps_c21_k*epst+xi*etaH*eps_c21_k*epst+xi*etaB*deps_c21_k*epst-epsL[1])
                                    else:              ww2_2 = 0.
                                    ww2_4 = self.bw_/dd4 *(-etaH*eps_c21_k*epst-etaB*deps_c21_k*epst+xi*etaH*eps_c21_k*epst+xi*etaB*deps_c21_k*epst-epsL[1])
                                    #
                                    if StateP2>=4:                          # state 4 or 5
                                        if   ww2_4>0:                         State2_,ww2,sigL,MatML = State4C2(ww2_4)         # state 4 open crack loading or unloading
                                        else:                                 State2_,ww2,sigL,MatML = State3C2(); State2_= 5  # state 5 crack closure
                                    else:                                   # crack 2 state 1 or 2 or 3 
                                        if   ww2_1>wcr:                       State2_,ww2,sigL,MatML = State4C2(ww2_4)         # new state 4
                                        elif epsL[1]>epsLP[1]+TOL and ww2_1>wP2:  State2_,ww2,sigL,MatML = State1C2(ww2_1)     # state 1 loading
                                        elif ww2_2 >0:                        State2_,ww2,sigL,MatML = State2C2(ww2_2)         # state 2 unloading
                                        elif ww2_2<=0:                        State2_,ww2,sigL,MatML = State3C2()              # state 3 crack closure
                                        else: raise NameError ("ConFemMaterials::ElasticLT.sig: crack exception 2")
                                else:                                         State2_ = 0
                            else:                                             State_,ww1,sigL,MatML = State3C(); State_= 5    # state 5 - crack closure
                        else:
                            if ww1_1>wcr:                                     State_,ww1,sigL,MatML = State4C(  ww1_4)        # starting state 4 loading beyond critical crack width   
                            elif epsL[0]>epsLP[0]+TOL and ww1_1>wP1:          State_,ww1,sigL,MatML = State1C(0,ww1_1)        # state 1 loading 
                            elif ww1_2>0:                                     State_,ww1,sigL,MatML = State2C(  ww1_2)        # state 2 unloading
                            elif ww1_2<=0:                                    State_,ww1,sigL,MatML = State3C()               # state 3 crack closure
                            else:
                                print(elI,ipI,Eps_,epsL,wcr,State,wP1,ww1_1,ww1_2,ww1_4)
                                if ff!=None: print(elI,ipI,Eps_,epsL,wcr,State,wP1,ww1_1,ww1_2,ww1_4, file=ff)
                                raise NameError ("ConFemMaterials::ElasticLT.sig: crack exception")
                        # finish 
                        Elem.StateVarN[ipI,0]  = State_
                        Elem.StateVarN[ipI,1]  = max(sigL[0],sigL[1])       # used for both 1st and 2nd crack in CalcType 0    
                        Elem.StateVarN[ipI,6]  = State2_
                        Elem.StateVarN[ipI,11] = ww1                        # current crack width 1
                        Elem.StateVarN[ipI,12] = ww2                        # current crack width 2
                        Elem.StateVarN[ipI,13] = 2.*(ww1-Elem.StateVar[ipI,11])/Dt - Elem.StateVar[ipI,13]  # current crack width velocity 1 
                        Elem.StateVarN[ipI,14] = 2.*(ww2-Elem.StateVar[ipI,12])/Dt - Elem.StateVar[ipI,14]  # current crack width velocity 2
                        #
                        MatM  = dot(Trans.transpose(),dot(MatML,Trans))     # transformation of tangent stiffness into global system
                        sig = dot(Trans.transpose(),sigL)                   # transformation of stress into global system

                        print('XXX',Elem.Label,ipI,'_',rho,wcr,State_,ww1,Eps,sig)
                        
                        if Elem.dim==21:
                            cc = 0.5*Emod/(2.+0.*nu)                        # zero poisson's ratio assumed for out of plane shear
                            sig_ = array([sig[0],sig[1],0.,cc*Eps[3],cc*Eps[4],sig[2]])
                            MatM_ = array([[MatM[0,0], MatM[0,1], 0., 0., 0., MatM[0,2]],
                                           [MatM[1,0], MatM[1,1], 0., 0., 0., MatM[1,2]],
                                           [0., 0., 0., 0., 0., 0.],
                                           [0., 0., 0., cc, 0., 0.],
                                           [0., 0., 0., 0., cc, 0],
                                           [MatM[2,0], MatM[2,1], 0., 0., 0., MatM[2,2]]])
                    # uncracked state 
                    else:
                        if Elem.dim==2:
                            xxx = Emod/(1-nu**2)
                            MatM = array([[xxx,xxx*nu,0],[xxx*nu,xxx,0],[0,0,0.5*xxx*(1-nu)]]) # isotropic linear elastic plane stress
                            sig = dot(MatM,Eps_) # stress
                        elif Elem.dim==21:
                            MatM_ = array([[Emod/(1-nu**2), nu*Emod/(1-nu**2), 0., 0., 0., 0.],
                                           [nu*Emod/(1-nu**2), Emod/(1-nu**2), 0., 0., 0., 0.],
                                           [0., 0., 0., 0., 0., 0.],
                                           [0., 0., 0., Emod/(2+2*nu), 0., 0.],
                                           [0., 0., 0., 0., Emod/(2+2*nu), 0],
                                           [0., 0., 0., 0., 0., Emod/(2+2*nu)]]) 
                            sig_ = dot(MatM_,Eps)
                            sig = array([sig_[0],sig_[1],sig_[5]])          #
                        pig, _, _ = PrinCLT_( sig[0], sig[1], sig[2])       # principal stresses, larger value, corr. direction, lower value
                        Elem.StateVarN[ipI,1] = pig                         # larger 1st principal stress for dual cracking control
                    if Elem.dim==2:
                        return sig, MatM, [ Eps_[0], Eps_[1], 0., Eps_[2], sig[0], sig[1], 0., sig[2], ww1, ww2]
                    elif Elem.dim==21:
#                        return sig_, MatM_, [Eps_[0],Eps_[1],Eps_[2],sig[0],sig[1],sig[2],ww1,ww2,sig_[0],sig_[1],sig_[2],sig_[3],sig_[4],sig_[5]]
#                        return sig_, MatM_, [sig_[1],sig_[2],sig_[3],sig_[4],sig_[5],Eps_[0],Eps_[1],Eps_[2],sig[0],sig[1],sig[2],ww1,ww2,sig_[0]] # strange ???
                        return sig_, MatM_, [sig_[0],sig_[1],sig_[2],sig_[3],sig_[4],sig_[5],Eps_[0],Eps_[1],Eps_[2],sig[0],sig[1],sig[2],ww1,ww2]
                # C-Version
                else:
                    DataOut = zeros((1),dtype=float)  # should have an eye on this whether if conforms to C-code
                    sig  = zeros((6),dtype=float)
                    MatM_ = zeros((36),dtype=float)
                    ww   = zeros((2),dtype=float)   
                    rc = ElasticLTC2(Elem.dim, Elem.PlSt, Elem.StateVar[ipI], Elem.StateVarN[ipI], 
                                    Eps, sig, MatM_, ww, self.nu, self.Emod, Dps, Dt,
                                    self.rho, self.wcr_, self.bw_, self.epsct, Elem.Lch_, self.CrackVisc, self.fct, self.epscu,
                                    DataOut)                    
                    if rc!=110: raise NameError("ConFemMaterials::ElasticLT:sig:ElasticLTC1 RC "+str(rc))
                    if Elem.dim==2:
                        MatM = array([[MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]])
                        return [sig[0], sig[1], sig[2]], MatM, [ Eps[0], Eps[1], 0., Eps[2], sig[0], sig[1], 0., sig[2], ww[0], ww[1]]
                    elif Elem.dim==21:
                        MatM = array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5]],\
                                      [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],\
                                      [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],\
                                      [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],\
                                      [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],\
                                      [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]])
#                        return sig, MatM, [Eps[0],Eps[1],Eps[5],sig[0],sig[1],sig[5],ww[0],ww[1],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
                        return sig, MatM, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],Eps[0],Eps[1],Eps[5],sig[0],sig[1],sig[5],ww[0],ww[1]]
    def UpdateStateVar(self, Elem, ff):
        if Elem.Type=='SH4': 
            if   Elem.nInt==2 :nn = min(Elem.StateVar.shape[0],16) # loop over integration and integration sub points excluding reinforcement with SH4 elements
            elif Elem.nInt==5 :nn = min(Elem.StateVar.shape[0],20)
        else:                  nn = Elem.StateVar.shape[0]    # loop over integration points
        SFlag = False   
        for j in range(nn):
            if (Elem.StateVar[j][0]!=Elem.StateVarN[j][0]): 
                Echo(f"El {Elem.Label:d} ip {j:d} ElasticLT state change {Elem.StateVar[j][0]:.0f} -> {Elem.StateVarN[j][0]:.0f}", ff)
                SFlag = True
            if (Elem.StateVar[j][6]!=Elem.StateVarN[j][6]): 
                print('El', Elem.Label, j, 'ElasticLT state 2 change', Elem.StateVar[j][6],'->',Elem.StateVarN[j][6])
                if ff!=None: print('El', Elem.Label, j, 'ElasticLT state 2 change', Elem.StateVar[j][6],'->',Elem.StateVarN[j][6], file=ff)
                Echo(f"El {Elem.Label:d} ip {j:d} ElasticLT state change 2 {Elem.StateVar[j][6]:.0f} -> {Elem.StateVarN[j][6]:.0f}", ff)
                SFlag = True
            for k in (0,1,4,5,6): Elem.StateVar[j,k] = Elem.StateVarN[j,k] # [0] state, [1] max. prin stress, [4] max prin strain, [5] direction corresponding to 1st crack upon initiation, [6] state2
            if Elem.StateVarN[j,2]>Elem.StateVar[j,2]: 
                Elem.StateVar[j,2] = Elem.StateVarN[j,2]    # maximum crack width ever reached
                Elem.StateVar[j,3] = Elem.StateVarN[j,3]    # related crack traction
            if Elem.StateVarN[j,7]>Elem.StateVar[j,7]: 
                Elem.StateVar[j,7] = Elem.StateVarN[j,7]    # maximum crack width ever reached 2 
                Elem.StateVar[j,8] = Elem.StateVarN[j,8]    # related crack traction 2
        return SFlag
    def UpdateStat2Var(self, Elem, ff, SFlag, TimeTargetActiveFlag):
        if Elem.Type=='SH4': 
            if   Elem.nInt==2 :nn = min(Elem.StateVar.shape[0],16) # loop over integration and integration sub points excluding reinforcement with SH4 elements
            elif Elem.nInt==5 :nn = min(Elem.StateVar.shape[0],20)
        else:                  nn = Elem.StateVar.shape[0]  # loop over integration points
        if TimeTargetActiveFlag or SFlag:                                   # LoFl: flag for loading increment, SFlag
            for j in range(nn):
                for k in (9,11,13):  Elem.StateVar[j,k] = Elem.StateVarN[j,k]
                for k in (10,12,14): Elem.StateVar[j,k] = Elem.StateVarN[j,k]
    
class ElasticSDA(Material):                                                          # linear elastic with limited tensile strength
    def __init__(self, PropMat, RedInt, RotCrack,PrinStrain,ShearRetFac, S2ndCrack,S2ndCrackA ):
        Material.__init__(self, True, 3, True, False, 13, 8, 'ELASTICLT_SDA')
#        self.Symmetric = True                              # flag for symmetry of material matrices
#        self.RType                                         # type of regularization by default 3: SDA with crack tractions
#        self.Update = True                                 # has specific update procedure for update of state variables
#        Update2
#        self.StateVar = 13                                 # [0-5] strain rate of previous time step (voigt?),
                                                            # 1D [6] float flag for failure 
                                                            # 2D 
                                                            # [6] largest principal strain, [7-9] largest principal strain direction, [10-12] stress vector in largest principal strain direction
# obsolete  # --> two integration points along discontinuity: [13-14] velocity of discontinuity normal (2D) !!!: number of integration points along discontinuity should not exceed number of bulk integration points
            #                                                 [15,16]                           tangent
            # crack width normal values values with two integration points are currently assigned to integration point 0 overruling the previous line
        self.Emod= PropMat[0]
        self.nu  = PropMat[1]
        self.fct = PropMat[2] 
        self.RegPar= PropMat[3]                                 # crack energy for  SDA RegType 3
        self.eta   = PropMat[4]                                 # artificial viscosity
        self.etaCr = PropMat[5]                                 # artificial viscosity for crack
        self.Density= PropMat[6]                                # specific mass
        self.eps_ct = self.fct/self.Emod
        self.RedInt = RedInt                                    # reduced integration for CPS4 using this material in case of SDA, see CaeFemElement:ElementC2D_SDA()
        self.RotCrack = RotCrack                                # rotating crack for CPS4 using this material
        self.PrinStrains = PrinStrain                           # SDA rankine criterion with principal strains / stresses
        self.S2ndCrack = S2ndCrack                              # allows for second crack
        if self.S2ndCrack:  self.S2ndCrackA = cos(S2ndCrackA*pi/180.)         # cos(deviation 2nd crack  from 1st crack in rad) should be smaller than this value, i.e. should be more or less perpendicular
        if ShearRetFac > 0: self.ShearRetFac = ShearRetFac
        else:               self.ShearRetFac = 1.0e-6
        self.RankineFactor = 1.0
        nu = self.nu
        ff = self.Emod / ( ( 1. + nu ) * ( 1. - 2.*nu ) )
        self.C3_ = ff*array([[1.-nu,nu,nu,0,0,0],[nu,1.-nu,nu,0,0,0],[nu,nu,1.-nu,0,0,0],[0,0,0,(1.-2.*nu)/2.,0,0],[0,0,0,0,(1.-2.*nu)/2.,0],[0,0,0,0,0,(1.-2.*nu)/2.]]) # linear elastic stiffness matrix
        
    def Sig(self,  ff, CalcType, Dt, elI, ipI, Elem, Dps_, Eps_, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        Emod = self.Emod
        nu = self.nu
        # C-Version
#        if ConFemMatCFlag:
        if False:   # C-Version not yet implemented ElasticLTC2 is for limited tensile strength with smeared crack
            DataOut = zeros((1),dtype=float)
            if Elem.dim==1:
                sig_  = zeros((2),dtype=float)
                MatM_ = zeros((4),dtype=float)
                PlSt_ = False
            elif Elem.dim==2:
                sig_  = zeros((3),dtype=float)
                MatM_ = zeros((9),dtype=float)
                PlSt_ = Elem.PlSt
            elif Elem.dim==3:
                sig_  = zeros((6),dtype=float)
                MatM_ = zeros((36),dtype=float)
                PlSt_ = False
            rc = ElasticLTC2( Elem.dim, PlSt_, self.PrinStrains, Elem.StateVar[ipI], Elem.StateVarN[ipI], Eps_, sig_, MatM_, nu, Emod, Dps_, self.eta, Dt, DataOut)
###############################
#            if Elem.Label==1: print('ZZZ', ipI, DataOut[0:4] ) #,'\n',DataOut[8:14],'\n',DataOut[14:20],'\n',DataOut[20:26],'\n',DataOut[26:32],'\n',DataOut[32:38])
#            if Elem.Label==38 and ipI==0:       print( 'ZZZ', Elem.Label, ipI, Elem.StateVar[ipI]) #  array([ [MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]])
###############################                
            if rc>110:
                raise NameError("ConFemMaterials::ElasticLTC2:MatC: RC "+str(rc)) 
            if Elem.dim==1:
                return [sig_[0],0], array([[MatM_[0],MatM_[1]],[MatM_[2],MatM_[3]]]), [Eps_[0], sig_[0]]
            elif Elem.dim==2:
                return ([sig_[0],sig_[1],sig_[2]]), array([ [MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]]), [Eps_[0],Eps_[1],0.,Eps_[2],sig_[0],sig_[1],0.,sig_[2]]
            elif Elem.dim==3:
                return [sig_[0],sig_[1],sig_[2],sig_[3],sig_[4],sig_[5]],\
                array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
                       [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
                       [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
                       [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
                       [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
                       [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
                       [sig_[0],sig_[1],sig_[2],sig_[3],sig_[4],sig_[5]]
        else:
            CC = zeros((6,6), dtype=float)
            CC[:] = self.C3_[:]                                     # triaxial tangential material stiffness
            if Elem.dim==1:
                if Elem.StateVar[ipI,6]>0.:  
                    MatM = array([[0.01*Emod,0],[0,0]])             # failure -- softening currently not yet implemented -- anyway currently overridden by element erosion, see Update  
                else:                        
                    MatM = array([[     Emod,0],[0,0]])
                sig = dot(MatM,Eps_)
                # viscous contribution
                if Elem.StateVar[ipI,6]>0. and self.eta>0.:
                    Dps__ = array([Dps_[0],0.,0.,0.,0.,0.]) 
                    zz, Veps = self.ViscExten3D( Dt, self.eta, Dps__, Elem, ipI, 0)
                    MatM[0,0] = MatM[0,0] + zz
                    sig[0] = sig[0] + self.eta*Veps[0] 
                return sig, MatM, [Eps_[0], sig[0]]
            else:
                if Elem.dim==2:
                    if Elem.PlSt:
                        alpha = self.nu/(1-self.nu)
                        Eps = array([[Eps_[0],0.5*Eps_[2],0],[0.5*Eps_[2],Eps_[1],0],[0,0,-alpha*(Eps_[0]+Eps_[1])]])# plane stress --> strain tensor 
                        Eps__ = array([Eps_[0],Eps_[1],-self.nu*(Eps_[0]+Eps_[1])/(1-self.nu),0.,0.,Eps_[2]]) # plane stress --> strain Voigt notation
                        Dps__ = array([Dps_[0],Dps_[1],-self.nu*(Dps_[0]+Dps_[1])/(1-self.nu),0.,0.,Dps_[2]]) # plane stress --> strain Voigt notation
                    else:         
                        Eps = array([[Eps_[0],0.5*Eps_[2],0],[0.5*Eps_[2],Eps_[1],0],[0,0,0]]) # plane strain --> strain tensor
                        Eps__ = array([Eps_[0],Eps_[1],0.,0.,0.,Eps_[2]]) # plane strain --> strain Voigt notation
                        Dps__ = array([Dps_[0],Dps_[1],0.,0.,0.,Dps_[2]]) # plane strain --> strain Voigt notation
                elif Elem.dim==3: 
                    Eps = array([[Eps_[0],0.5*Eps_[5],0.5*Eps_[4]],[0.5*Eps_[5],Eps_[1],0.5*Eps_[3]],[0.5*Eps_[4],0.5*Eps_[3],Eps_[2]]]) # triaxial strain tensor
                    Eps__ = array([Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5]])                            # triaxial strain Voigt notation
                    Dps__ = array([Dps_[0],Dps_[1],Dps_[2],Dps_[3],Dps_[4],Dps_[5]])                            # triaxial strain Voigt notation
                else: raise NameError("CaeFemMaterials::ElasticLT.Sig: element type not implemented for this material")
                # stress without viscous contribution
                sig  = dot(self.C3_,Eps__)
                # principal strains (EvFlag = True), stresses
                if False:
                    if ConFemMatCFlag:
                        laMax= zeros((1),dtype=float)
                        eVec = zeros((3),dtype=float)
                        if self.PrinStrains:
                            if eigJacobiSymWrapper(array([Eps[0,0], Eps[1,1], Eps[2,2], Eps[1,2], Eps[0,2], Eps[0,1]]), laMax, eVec) != 100:
                                raise NameError("CaeFemMaterials::ElasticLT.Sig: Eigensystem iteration failed: %12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e"%(Eps[0,0],Eps[1,1],Eps[2,2],Eps[1,2],Eps[0,2],Eps[0,1]))
                        else:
                            if eigJacobiSymWrapper(array([sig[0], sig[1], sig[2], sig[3], sig[4], sig[5]]), laMax, eVec) != 100:
                                raise NameError("CaeFemMaterials::ElasticLT.Sig: stress Eigensystem iteration failed: %12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e"%(sig[0], sig[1], sig[2], sig[3], sig[4], sig[5]))
                    else:
                        if self.PrinStrains: laMax, eVec = EigenJacobiSym(Eps, 3)
                        else:                laMax, eVec = EigenJacobiSym(array([[sig[0],sig[5],sig[4]], [sig[5],sig[1],sig[3]], [sig[4],sig[3],sig[2]]]), 3)
                    if Elem.dim==2 and abs(eVec[2])>ZeroD:  laMax=0.    # lateral extension to plane, no lateral cracks allowed
                    if (eVec[0]+eVec[1]+eVec[2]<0.):                    # for scalar product with (1,1,1) -- directions indicate the same plane (approximately for different IPs) but might have opposite directions. Has to be avoided 
                        eVec[0] = -eVec[0]
                        eVec[1] = -eVec[1]
                        eVec[2] = -eVec[2]
                    
                else:
                    if self.PrinStrains: laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([Eps[0,0], Eps[1,1], Eps[2,2], Eps[1,2], Eps[0,2], Eps[0,1]]) )
                    else:                laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([sig[0],   sig[1],   sig[2],   sig[3],   sig[4],   sig[5]]) )
                # viscous contribution   
                if self.eta>0.:
                    zz, Veps = self.ViscExten3D( Dt, self.eta, Dps__, Elem, ipI, 0)
                    for k in range(6): CC[k,k] = CC[k,k] + zz       # viscous extension for tangential material stiffness
                    sig  = sig  + self.eta*Veps                     # triaxial stress with viscous extension
                # direction & traction
                if False: #True:
                    tt = [sig[0]*eVec[0]+sig[5]*eVec[1]+sig[4]*eVec[2], # crack traction in global system -- might have nonzero local shear components due to viscous contributions  
                          sig[5]*eVec[0]+sig[1]*eVec[1]+sig[3]*eVec[2], 
                          sig[4]*eVec[0]+sig[3]*eVec[1]+sig[2]*eVec[2]]
                    if norm(tt) < ZeroD : laMax=0.       # no crack allowed incase of orthog. sig and eVec e.g. uniaxial compr. 3D - Ahmad
                    Elem.StateVarN[ipI,6]  = laMax
                    Elem.StateVarN[ipI,7]  = eVec[0]                    # eigenvector of largest principal strain
                    Elem.StateVarN[ipI,8]  = eVec[1]
                    Elem.StateVarN[ipI,9]  = eVec[2]
                    Elem.StateVarN[ipI,10] = tt[0]                      # stress vector derived from stress tensor and above eigenvector 
                    Elem.StateVarN[ipI,11] = tt[1]
                    Elem.StateVarN[ipI,12] = tt[2]
                else:
                    self.RankineUpd( Elem, ipI, sig, eVec, laMax, 0)
                # more values may be updated in ViscExtenCr2D; sorry for that inconsistency (uhc)
                if Elem.dim==2:
                    if Elem.PlSt: CC_= array([[CC[0,0]-CC[0,2]*CC[2,0]/CC[2,2], CC[0,1]-CC[0,2]*CC[2,1]/CC[2,2], CC[0,5]-CC[0,2]*CC[2,5]/CC[2,2]],
                                              [CC[1,0]-CC[1,2]*CC[2,0]/CC[2,2], CC[1,1]-CC[1,2]*CC[2,1]/CC[2,2], CC[1,5]-CC[1,2]*CC[2,5]/CC[2,2]], 
                                              [CC[5,0]-CC[5,2]*CC[2,0]/CC[2,2], CC[5,1]-CC[5,2]*CC[2,1]/CC[2,2], CC[5,5]-CC[5,2]*CC[2,5]/CC[2,2]]])
                    else: 
                        CC_= array([[CC[0,0],CC[0,1],CC[0,5]],[CC[1,0],CC[1,1],CC[1,5]],[CC[5,0],CC[5,1],CC[5,5]]])
                    return [sig[0],sig[1],sig[5]], CC_ , [Eps__[0],Eps__[1],0.,Eps__[5],sig[0],sig[1],0.,sig[5]]
                elif Elem.dim==3:
                    return sig, CC, [Eps__[0],Eps__[1],Eps__[2],Eps__[3],Eps__[4],Eps__[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]

    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):
            if Elem.dim==1 and Elem.Active:                                 # unaxial
                sig0 = Elem.Data[j,1]
                if sig0>self.fct: 
                    Elem.StateVarN[j,6]  = 1.0                              # failure condition
                    print("1D ElasticLT failure",  Elem.Label, j, sig0, self.fct, file=ff)
                    Elem.Active = False                                     # element erosion

            Elem.StateVar[j,0] = Elem.StateVarN[j,0]                        # strain velocities - used in ViscExten3D
            Elem.StateVar[j,1] = Elem.StateVarN[j,1]
            Elem.StateVar[j,2] = Elem.StateVarN[j,2]
            Elem.StateVar[j,3] = Elem.StateVarN[j,3]
            Elem.StateVar[j,4] = Elem.StateVarN[j,4]
            Elem.StateVar[j,5] = Elem.StateVarN[j,5]
            
            Elem.StateVar[j,6] = Elem.StateVarN[j,6]                        # crack data
            Elem.StateVar[j,7] = Elem.StateVarN[j,7]
            Elem.StateVar[j,8] = Elem.StateVarN[j,8]
            Elem.StateVar[j,9] = Elem.StateVarN[j,9]
            Elem.StateVar[j,10]= Elem.StateVarN[j,10]
            Elem.StateVar[j,11]= Elem.StateVarN[j,11]
            Elem.StateVar[j,12]= Elem.StateVarN[j,12]
#            Elem.StateVar[j,13]= Elem.StateVarN[j,13]
#            Elem.StateVar[j,14]= Elem.StateVarN[j,14]
#            Elem.StateVar[j,15]= Elem.StateVarN[j,15]
#            Elem.StateVar[j,16]= Elem.StateVarN[j,16]
        return False

class IsoDamage(Material):                                                  # isotropic damage
    def __init__(self, PropMat, RotCrack,PrinStrain,ShearRetFac, S2ndCrack,S2ndCrackA):
#                        (self, SymmetricVal, RTypeVal,   UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, False,        PropMat[8], True,     False,      16,           10, 'ISODAMAGE') # same as in ComFemInOut::DataInput for newly created SDA elements
#        Material.__init__(self, False,        PropMat[8], True,     False,      16,           9, 'ISODAMAGE')  # NData overridden in ComFemInOut::DataInput
#        self.Symmetric = False                                             # flag for symmetry of material matrices
#        self.RType = PropMat[8]                                            # type of regularization 0: without, 1: gradient 2: crack band 3: SDA
#        self.Update = True                                                 # has specific update procedure for update of state variables
#        self.StateVar = 16                                                 # number of state variables [0] damage [1] equivalent damage strain, [2-7]strain rate of previous time step (voigt?), 
                                                                            # [8] crack energy, [9] largest principal strain/sresss, [10-12] largest principal val direction, [13-15] traction in largest principal strain direction
                                                                            # [16-21] strain increment of current equilibrium iteration
                                                                            # 2 + 6 + 1 + 7 + 6
        self.PropMat= PropMat                                               # for wrappers, e.g. RCSHELL
        self.Emod= PropMat[0]
        self.nu  = PropMat[1]
        self.fc  = PropMat[2]
        self.fct = PropMat[3]
        self.LiTy= PropMat[4]                                               # type of limit function
        self.beta= PropMat[5]
        self.ga  = PropMat[6]
        self.a3  = PropMat[7]
        self.edt,self.ed,self.gd,self.kapStrength, self.cc0,self.cc1,self.cc2,self.cc3 = self.EvalPar_( self.Emod,self.nu, self.fc,self.fct, self.beta,self.ga,self.a3)
        # [8] set with material
        self.RegPar = PropMat[9]                                            # regularization parameter;  char. length for gradient regularization RegType 1
                                                                            #                            crack energy for crack band approach RegType 2
        self.eta = PropMat[10]                                              # artificial viscosity
        self.etaCr = PropMat[11]                                            # artificial viscosity for SDA crack
        self.Density = PropMat[13]  # specific mass
        self.DStrength = 1.0 - exp(-pow((self.kapStrength - self.edt) / self.ed, self.gd))  # scalar damage
        dDestr   = 1.-1.e-6
        self.kapUlt   = exp(log(-log(1.-dDestr))/self.gd)*self.ed+self.edt
        # following for regularization
        # uses material data of *MATERIAL section  -- random field data are not considered --> random field data should correspond to *MATERIAL data
        self.ElCharLenBound = 0.278 # PropMat[12]                                   # boundary (related to crack band width) to distinguish between scaling types if Rtype==2-crack band - might be 0 in that case type 1 is currently used
        a, b = self.cc0*(1.+self.nu)**2/3., self.cc1*(1.+self.nu)/sqrt(3.) + self.cc2 + self.cc3*(1.-2.*self.nu) #auxiliary values
        self.alpha = 0.5*b + sqrt(0.25*b**2+a)                              # scaling factor for tension - required
#        eps_ct = self.kapStrength / self.alpha
        self.eps_ct = self.kapStrength / self.alpha
        self.SpecCrEn = self.SpeCrEnergySample( None, self.eps_ct, False,None,None,                      self.UniaxTensionUnit)
                                                                            # from Material - specific crack energy - is unscaled her - presumably used by element init
        self.gam2    = 350.                                                 # parameter for scaling of equivalent damage strain
        self.lamLow  = 0.01
        self.lamHigh = 100.
        self.lam2Low = 0.01
        self.lam2High= 0.5
        if self.RType==2:
            self.bw = self.RegPar/self.SpecCrEn                             # crack energy (RegPar) is crack band width times volume specific crack energy (factor 1 in preceding)
            self.CrX, self.CrY   = self.SpeCrEnergyData( "LargeEl",self.lamLow,self.lamHigh,  self.CrBwN, self.eps_ct, self.UniaxTensionScaled) # arrays to determine scaling factor for given characteristic element length
            self.CrX2, self.CrY2 = self.SpeCrEnergyData( "SmallEl",self.lam2Low,self.lam2High,self.CrBwN, self.eps_ct, self.UniaxTensionScaled2) # arrays to determine scaling factor for given characteristic element length
        else:
            self.bw = 0.
            self.CrX, self.CrY, self.CrX2, self.CrY2 = None, None, None, None
        # regularization finished
        self.svrTol = 0.01                                                  # minimum reference value for viscous stresses - for control computations only
        nu = self.nu
        ff = self.Emod / ( (1.+nu) * (1.-2.*nu) )
        self.C3_ = ff*array([[1.-nu,nu,nu,0,0,0],[nu,1.-nu,nu,0,0,0],[nu,nu,1.-nu,0,0,0],[0,0,0,(1.-2.*nu)/2.,0,0],[0,0,0,0,(1.-2.*nu)/2.,0],[0,0,0,0,0,(1.-2.*nu)/2.]]) # linear elastic stiffness matrix
        self.RotCrack = RotCrack                                            # rotating crack for CPS4 using this material      in case of SDA 
        self.PrinStrains = PrinStrain                                       # SDA rankine criterion with principal strains / stresses 
        self.S2ndCrack = S2ndCrack                                          # allows for second crack
        if self.S2ndCrack: self.S2ndCrackA = cos(S2ndCrackA*pi/180.)        # cos(deviation 2nd crack  from 1st crack in rad) should be smaller than this value, i.e. should be more or less perpendicular
        if ShearRetFac > 0: self.ShearRetFac = ShearRetFac
        else:               self.ShearRetFac = 1.0e-6
        self.RankineFactor = 0.85 #1.0                                      # for SDA, PrinStrains = False i.e. stress only
    def EvalPar_(self, Emod,nu, fC, fT, beta,ga,a3):
        ll   = exp(-1./2.)
        epsC = fC/Emod/ll                                                   # strain with compressive strength
        gd = 2.0
        edt= epsC*(2.*log(fC/Emod/epsC)+1.)                                 # gd in here?
        ed = 2.*sqrt(-log(fC/Emod/epsC))*epsC
#        kapStrength = 0.5*(edt+sqrt(edt**2+2.*ed**2))                       # equivalent strain for compressive strenght -- same as epsC ?
#        print('XXX',epsC, kapStrength)
        alpha = fT / fC
        A = -3./2.*(alpha*ga*beta-beta*a3*alpha+beta*alpha-ga*alpha+ga*beta)/alpha/(-ga+2.*a3*ga-a3**2-ga**2+a3-ga*alpha+ga*beta)/beta
        B = -1./3.*(-ga*beta**2*alpha+2*alpha*beta**2-2*a3*beta**2*alpha+3*ga*beta**2+4*ga*alpha**2*beta+3*ga**2*beta*alpha+2*a3*ga*beta-ga**2*beta-6*a3*ga*beta*alpha-3*beta*alpha+3*a3**2*beta*alpha+beta*a3-4*ga*beta-a3**2*beta+alpha**2*beta-alpha**2*beta*a3+4*a3*ga*alpha+ga*alpha-2*ga**2*alpha-2*a3**2*alpha-3*ga*alpha**2+2*a3*alpha)/alpha/beta*sqrt(6.)/(ga*alpha-2*a3*ga+ga+ga**2-a3+a3**2-ga*beta)
        C =  1./2.*(-ga*beta**2-alpha*beta**2+4*a3*ga*beta*alpha+2*beta*alpha-2*a3*ga*alpha-2*a3**2*beta*alpha+a3*beta**2*alpha+ga*beta**2*alpha-2*ga**2*beta*alpha-2*ga*alpha**2*beta-ga*alpha+a3**2*alpha+ga**2*alpha-a3*alpha+ga*alpha**2+2*ga*beta)/alpha/beta*sqrt(6.)/(-ga+2*a3*ga-a3**2-ga**2+a3-ga*alpha+ga*beta)
        D = -1./3.*(-a3*beta**2*alpha+alpha*beta**2+ga*beta**2*alpha+a3**2*beta+alpha**2*beta*a3-beta*a3-alpha**2*beta+ga*beta-2*a3*ga*beta-ga*alpha**2*beta+ga**2*beta+2*a3*ga*alpha+a3*alpha-ga*alpha-a3**2*alpha-ga**2*alpha)*sqrt(3.)/alpha/(-ga+2*a3*ga-a3**2-ga**2+a3-ga*alpha+ga*beta)/beta
        cc0 = 2.*A/(1.+nu)**2
        cc1 = C*sqrt(2.)/(1+nu)
        cc2 = B*sqrt(1.5)/(1+nu)
        cc3 = D/sqrt(3.)/(1-2*nu)-B/sqrt(6.)/(1+nu)
        return edt,ed,gd,epsC, cc0,cc1,cc2,cc3
    # for regularization
    def UniaxTensionScaled(self, gam1, eps):                                # uniaxial tension with scaling type 1 to be applied for eps>eps_ct
        x = self.alpha*eps - self.kapStrength                               # tensile strain delta after tensile strength strain
        kap_ =gam1*x + (1.-gam1)/self.gam2 * (1.-exp(-self.gam2*x)) + self.kapStrength
        return self.Emod*exp(-pow(( kap_ - self.edt)/self.ed ,self.gd)) * eps
    def UniaxTensionScaled2(self, beta, eps):                               # uniaxial tension with scaling type 2 to be applied for eps>eps_ct
        kap = self.alpha*eps
        kap_=  (1-beta)*self.kapStrength * log( ( kap-beta*self.kapStrength )/( (1-beta)*self.kapStrength) ) + self.kapStrength
        return self.Emod*exp(-pow(( kap_ - self.edt)/self.ed ,self.gd)) * eps
    def UniaxTensionUnit(self, dummy, eps):                                 # uniaxial tension without scaling
        kap = self.alpha*eps
        return self.Emod*exp(-pow(( kap - self.edt)/self.ed ,self.gd)) * eps
    # equivalent strain from multiaxial strain
    def EquivStrain1(self, I1,J2,J2S, Eps,EpsD, cc0,cc1,cc2,cc3):
        la,v = eigh(Eps)                                                    # principal values and eigenvectors of strain tensor
        ila = 0                                                             # largest principal value
        if la[1]>la[0]: ila=1                       
        if la[2]>la[1] and la[2]>la[0]: ila=2
        laMax = la[ila]
#        xxx = 0.25*pow((self.cc1*J2S+self.cc2*laMax+self.cc3*I1),2.0) + self.cc0*J2
        xxx = 0.25*pow((      cc1*J2S+     cc2*laMax+     cc3*I1),2.0) +      cc0*J2
        if xxx<0.: raise NameError("ConFemMaterials::IsoDamage.EquivStrain1: inconsistent state")
#        kap_ = 0.5*(self.cc1*J2S+self.cc2*laMax+self.cc3*I1)+(sqrt(xxx))    # equivalent damage strain
        kap_ = 0.5*(      cc1*J2S+     cc2*laMax+     cc3*I1)+(sqrt(xxx))    # equivalent damage strain
        if J2S>ZeroD:
            eins = array([1.,1.,1.,0.,0.,0.])                               # auxiliary vector
            nd_ = array([v[0,ila]*v[0,ila],v[1,ila]*v[1,ila],v[2,ila]*v[2,ila],v[1,ila]*v[2,ila],v[0,ila]*v[2,ila],v[0,ila]*v[1,ila]]) # corresponding to Voigt notation
#            nd  = (self.cc0+self.cc1*kap_/(2.0*J2S))*EpsD + self.cc2*kap_*nd_ +self.cc3*kap_*eins # gradient of damage function
            nd = (       cc0 +    cc1*kap_/(2.0*J2S))*EpsD +      cc2*kap_*nd_ +     cc3*kap_*eins # gradient of damage function
#            Hd = -( self.cc1*J2S + self.cc2*laMax + self.cc3*I1 -2*kap_ )   # H_d: dF/dkappa local
            Hd = -(       cc1*J2S +      cc2*laMax +      cc3*I1 -2*kap_ )   # H_d: dF/dkappa local
        else:
            nd = zeros((6))
            Hd = 1
        return kap_, nd, Hd, laMax
#    def EquivStrain2(self, I1, J2):
#        kap_ = self.k0 * I1 + sqrt( (self.k1 * I1)**2 + self.k2*J2 )
#        Hd = 1.                             #
#        xx = sqrt((self.k1*I1)**2+self.k2*J2)
#        if xx>ZeroD: nd = (self.k0 + self.k1**2*I1/xx ) * eins + 0.5*self.k2/xx*EpsD #
#        else:        nd = zeros((6),dtype=float)
#        return kap_, nd, Hd

#    @profile
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps_, Eps_, dTmp, Temp, EpsR):  # Dps_ should be strain increment related to time step increment - not equilibrium iteration element
        if CalcType == 0: return [], [], []
        Emod = self.Emod
        nu = self.nu
        r = self.fct
        AlgFlag = False #True                                                     # flag for algorithmic modulus modification not yet implemented in C-version
        if Elem.RandomData != []:
            cc0 = Elem.RandomData[ipI,0]
            cc1 = Elem.RandomData[ipI,1]
            cc2 = Elem.RandomData[ipI,2]
            cc3 = Elem.RandomData[ipI,3]
            r   = Elem.RandomData[ipI,4]
        else:
            cc0 = self.cc0
            cc1 = self.cc1
            cc2 = self.cc2
            cc3 = self.cc3
                                                                            # take care for additional StateVar in __init__ and UpdateStateVar#
        if not ConFemMatCFlag: # False: # flag for not C-version
#        if True:
            D = Elem.StateVar[ipI,0]                                        # D of last converged load / time increment
            kapOld = Elem.StateVar[ipI,1]                                   # kappa of last converged load / time increment
            if Elem.dim==1:
                Eps = array([[Eps_[0],0,0],[0,-self.nu*Eps_[0],0],[0,0,-self.nu*Eps_[0]]]) # tensor notation
                Eps__ = array([Eps[0,0],Eps[1,1],Eps[2,2],0,0,0])           # voigt notation
                Dps__ = array([Dps_[0],-self.nu*Dps_[0],-self.nu*Dps_[0],0,0,0]) # voigt notation
            elif Elem.dim==2:
                if Elem.PlSt:
                    alpha = self.nu/(1-self.nu)
                    Eps = array([[Eps_[0],0.5*Eps_[2],0],[0.5*Eps_[2],Eps_[1],0],[0,0,-alpha*(Eps_[0]+Eps_[1])]])# plane stress --> strain tensor 
                    Eps__ = array([Eps_[0],Eps_[1],-self.nu*(Eps_[0]+Eps_[1])/(1-self.nu),0.,0.,Eps_[2]]) # plane stress --> strain Voigt notation
                    Dps__ = array([Dps_[0],Dps_[1],-self.nu*(Dps_[0]+Dps_[1])/(1-self.nu),0.,0.,Dps_[2]]) # plane stress --> strain Voigt notation
                else:         
                    Eps = array([[Eps_[0],0.5*Eps_[2],0],[0.5*Eps_[2],Eps_[1],0],[0,0,0]]) # plane strain --> strain tensor
                    Eps__ = array([Eps_[0],Eps_[1],0.,0.,0.,Eps_[2]])       # plane strain --> strain Voigt notation
                    Dps__ = array([Dps_[0],Dps_[1],0.,0.,0.,Dps_[2]])       # plane strain --> strain Voigt notation
            elif Elem.dim == 4:
                Eps = array([[    Eps_[0],0.5*Eps_[3],0],
                             [0.5*Eps_[3],    Eps_[1],0],
                             [0,              0,      Eps_[2]]])  # plane strain --> strain tensor
                Eps__ = array([Eps_[0], Eps_[1], Eps_[2], 0., 0., Eps_[3]])  # plane strain --> strain Voigt notation
                Dps__ = array([Dps_[0], Dps_[1], Dps_[2], 0., 0., Dps_[3]])  # plane strain --> strain Voigt notation
            elif Elem.dim==3:
                Eps = array([[Eps_[0],0.5*Eps_[5],0.5*Eps_[4]],[0.5*Eps_[5],Eps_[1],0.5*Eps_[3]],[0.5*Eps_[4],0.5*Eps_[3],Eps_[2]]]) # triaxial strain tensor
                Eps__ = array([Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5]])  # triaxial strain Voigt notation
                Dps__ = array([Dps_[0],Dps_[1],Dps_[2],Dps_[3],Dps_[4],Dps_[5]])  # triaxial strain Voigt notation
            elif Elem.dim==21:                                              # continuum based shell
                xxx = -self.nu/(1-self.nu)*(Eps_[0]+Eps_[1])
                Eps = array([[Eps_[0],0.5*Eps_[5],0.5*Eps_[4]],[0.5*Eps_[5],Eps_[1],0.5*Eps_[3]],[0.5*Eps_[4],0.5*Eps_[3],xxx]]) # triaxial strain tensor
                Eps__ = array([Eps_[0],Eps_[1],xxx,Eps_[3],Eps_[4],Eps_[5]])# triaxial strain Voigt notation
                xxx = -self.nu/(1-self.nu)*(Dps_[0]+Dps_[1])
                Dps__ = array([Dps_[0],Dps_[1],xxx,Dps_[3],Dps_[4],Dps_[5]])# triaxial strain Voigt notation
            else: raise NameError("ConFemMaterials::IsoDamage.Sig: element type not implemented for this material")

            I1=Eps[0,0]+Eps[1,1]+Eps[2,2]                                   # 1st invariant strain tensor
            pp=I1/3.                                                        # volumetric strain
            EpsD = array([Eps[0,0]-pp,Eps[1,1]-pp,Eps[2,2]-pp,Eps[1,2],Eps[0,2],Eps[0,1]]) # deviatoric strain tensor components
            J2 =0.5*(EpsD[0]*EpsD[0]+EpsD[1]*EpsD[1]+EpsD[2]*EpsD[2])+EpsD[3]*EpsD[3]+EpsD[4]*EpsD[4]+EpsD[5]*EpsD[5]# 2nd invariant of deviator
            J2S=sqrt(J2)
            if self.LiTy==1:
                kap_, nd, Hd, laMax = self.EquivStrain1( I1,J2,J2S, Eps,EpsD, cc0,cc1,cc2,cc3)
            else:
                raise NameError("ConFemMaterials::IsoDamage.Sig: unknown type of limit function")
            
            if self.RType==1:                                               # gradient damage
                if Elem.dim==  1: kap = EpsR[1]
                elif Elem.dim==2: kap = EpsR[3]
                elif Elem.dim==3: kap = EpsR[6]
                dkk = 1.
            elif self.RType==2:                                             # crack band
                if kap_> self.kapStrength:
                    beta = Elem.CrBwS                                       # element specific scaling factor from Mat.CCrackBandScale applied for each element depending on characteristic length
                    if Elem.ScaleType==1:
                        kap = beta*(kap_-self.kapStrength) + (1.-beta)/self.gam2 * (1-exp(-self.gam2*(kap_-self.kapStrength))) + self.kapStrength # scaled damage strain
                        dkk = beta                         + (1.-beta)           *    exp(-self.gam2*(kap_-self.kapStrength)) # scaling factor for tangential material stiffness
                    else:
                        kap = (1-beta)*self.kapStrength * log( ( kap_-beta*self.kapStrength )/( (1-beta)*self.kapStrength) ) + self.kapStrength
                        dkk = (1-beta)*self.kapStrength/(kap_-beta*self.kapStrength)
                else: 
                    kap = kap_
                    dkk = 1. 
            else:                                                           # no regularization or SDA
                kap = kap_
                dkk = 1.
            if kap>self.kapUlt: kap=self.kapUlt                             # Damage should not be zero to avoid numerical singularity. This constrains D to dDestr
            #
            sig0  = dot(self.C3_,Eps__)                         # for what is this needed - maybe later for gradient damage
            #
            if kap>self.edt and kap>=kapOld and J2S>ZeroD:                  # case loading with nonzero strain deviator
            #  if kap>self.edt and kap>=(kapOld+kapZero) and J2S>ZeroD:  # case loading with nonzero strain deviator -- in CaeFem - but same ConFemMatC
                D = 1.0 - exp(-pow(( kap-self.edt)/self.ed ,self.gd))       # scalar damage
                hdI = pow(( kap-self.edt)/self.ed ,self.gd)*self.gd/(kap-self.edt)*exp(-pow((kap-self.edt)/self.ed,self.gd)) * dkk # dD/dkappa
                if not AlgFlag: 
                    sig0 = dot(self.C3_,Eps__)
                else:                                                       # considering algorithmic modulus
                    delDps = Dps__ - Elem.StateVar[ipI,16:22]        
                    sig0 = dot(self.C3_,(Eps__ + delDps))                   # uhc for algorithmic modulus Dps__ should have same 3D Voigt format as Eps__
                if self.RType==1: CD = zeros((6,6))         
                else:             CD = hdI/Hd * outer(sig0,nd)              # tangential material stiffness loading (voigt for Eps__ should be correct as C3_*Eps__ is a stress)
            else: 
                CD = zeros((6,6),dtype=float)                               # case of unloading or zero strain deviator
                nd = zeros((6))
                Hd = 1
                hdI = 0
            #
#            if Elem.Label==131: print('XXX',Elem.Label,ipI,self.kapStrength,kap_,kap,'_',D)
            #
            if AlgFlag:
                Elem.StateVar[ipI,16:22]  = Dps__                               # uhc for algorithmic modulus 22-01-09
                delD = D - Elem.StateVarN[ipI,0]                                # uhc for algorithmic modulus 22-01-09
            Elem.StateVarN[ipI,0] = D                                       # store damage of actual iteration
            Elem.StateVarN[ipI,1] = kap                                     # store equivalent damage strain of actual iteration
            if self.eta>0.:                                                 # viscous regularization
                zz, Veps = self.ViscExten3D( Dt, self.eta, Dps__, Elem, ipI, 2)
            else:
                zz = 0.
                Veps = zeros((6),dtype=float)
            #
            if not AlgFlag: CC = (1-D)     *self.C3_ - CD                                        # triaxial tangential material stiffness
            else:           CC = (1-D-delD)*self.C3_ - CD
            #
            for k in range(6): CC[k,k] = CC[k,k] + zz                       # viscous extension for tangential material stiffness
            sigV = self.eta*Veps
            sig = (1-D)*dot(self.C3_,Eps__)
            #
            if self.RType==3:                                               # before viscous contribution
                if self.PrinStrains: laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([Eps[0,0], Eps[1,1], Eps[2,2], Eps[1,2], Eps[0,2], Eps[0,1]]) )
                else:                laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([sig[0],   sig[1],   sig[2],   sig[3],   sig[4],   sig[5]]) )
            #
            svs = 0.
            for i in range(6):                                              # to record size of viscous stress 
                if abs(sig[i])>self.svrTol: xxx = sigV[i]/sig[i]
                else:                       xxx = 0.
                if abs(xxx)>abs(svs): svs = xxx         
            sig = sig  + sigV                                               # triaxial stress with viscous extension
            #
            if self.RType==3: self.RankineUpd( Elem, ipI, sig, eVec,laMax,3)# needed for explicit SDA creation, last value offset to 6 in StateVAr: damage, damage strain, crack energy on index 8
            #                              
            if kap>self.kapStrength and kap>kapOld: Elem.StateVarN[ipI,8] = dot(sig,Dps__) # Crack energy increment for step
            else:                                   Elem.StateVarN[ipI,8] = 0.
            if Elem.dim==1:
                if self.RType==1:
                    ccc = 0.5*self.RegPar**2
                    fact = 0.5*self.Emod
                    nd[0]=nd[0]-self.nu*nd[1]-self.nu*nd[2]
                    CR = array([[0,-self.Emod*Eps_[0]*hdI],[-fact*nd[0]/Hd,fact]])
                    return [sig[0],ccc*Eps_[1]*fact], array([[CC[0,0],0],[0,ccc*fact]]), [0,(EpsR[1]-kap_)*fact], CR, [Eps_[0], sig[0], D]
                else:
                    return [sig[0],0], array([[CC[0,0],0],[0,0]]), [Eps_[0], sig[0],D]
            elif Elem.dim==2:
                if self.RType==1:                                           # gradient damage
                    ccc = 0.5*self.RegPar**2
                    fact = 0.5*self.Emod
                    if Elem.PlSt:                                           # plane stress
                        CC_= array([[CC[0,0]-CC[0,2]*CC[2,0]/CC[2,2],CC[0,1]-CC[0,2]*CC[2,1]/CC[2,2],CC[0,5]-CC[0,2]*CC[2,5]/CC[2,2], 0, 0],
                                    [CC[1,0]-CC[1,2]*CC[2,0]/CC[2,2],CC[1,1]-CC[1,2]*CC[2,1]/CC[2,2],CC[1,5]-CC[1,2]*CC[2,5]/CC[2,2], 0, 0], 
                                    [CC[5,0]-CC[5,2]*CC[2,0]/CC[2,2],CC[5,1]-CC[5,2]*CC[2,1]/CC[2,2],CC[5,5]-CC[5,2]*CC[2,5]/CC[2,2], 0, 0],
                                    [                              0,                              0,                              0, ccc*fact,0],
                                    [                              0,                              0,                              0, 0, ccc*fact]])
                        CR = array([[                           0,                           0,             0, -sig0[0]*hdI],
                                    [                           0,                           0,             0, -sig0[1]*hdI],
                                    [                           0,                           0,             0, -sig0[5]*hdI],
                                    [-fact*(nd[0]-alpha*nd[2])/Hd,-fact*(nd[1]-alpha*nd[2])/Hd,-fact*nd[5]/Hd, fact]])
                    else:                                                   # plane strain
                        CC_= array([[CC[0,0],CC[0,1],CC[0,5], 0, 0],
                                    [CC[1,0],CC[1,1],CC[1,5], 0, 0], 
                                    [CC[5,0],CC[5,1],CC[5,5], 0, 0],
                                    [      0,      0,      0, ccc*fact, 0],
                                    [      0,      0,      0, 0, ccc*fact]])
                        CR = array([[             0,             0,             0, -sig0[0]*hdI],
                                    [             0,             0,             0, -sig0[1]*hdI],
                                    [             0,             0,             0, -sig0[5]*hdI],
                                    [-fact*nd[0]/Hd,-fact*nd[1]/Hd,-fact*nd[5]/Hd, fact]])
                    return [sig[0],sig[1],sig[5],ccc*Eps_[3]*fact,ccc*Eps_[4]*fact], CC_, [0,0,0,(EpsR[3]-kap_)*fact], CR, [Eps_[0], Eps_[1], Eps_[2], sig[0], sig[1], sig[2]]                            
                else:
                    if Elem.PlSt: 
                        CC_= array([[CC[0,0]-CC[0,2]*CC[2,0]/CC[2,2],CC[0,1]-CC[0,2]*CC[2,1]/CC[2,2],CC[0,5]-CC[0,2]*CC[2,5]/CC[2,2]],
                                              [CC[1,0]-CC[1,2]*CC[2,0]/CC[2,2],CC[1,1]-CC[1,2]*CC[2,1]/CC[2,2],CC[1,5]-CC[1,2]*CC[2,5]/CC[2,2]], 
                                              [CC[5,0]-CC[5,2]*CC[2,0]/CC[2,2],CC[5,1]-CC[5,2]*CC[2,1]/CC[2,2],CC[5,5]-CC[5,2]*CC[2,5]/CC[2,2]]])
                    else: 
                        CC_= array([[CC[0,0],CC[0,1],CC[0,5]],[CC[1,0],CC[1,1],CC[1,5]],[CC[5,0],CC[5,1],CC[5,5]]])
                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                    if ipI in [3]: #[0,1,2,3]:
#                        print('XXX',Elem.Lch_,ipI,kap_,self.kapStrength,kap,D,'\n[',Eps__[0],Eps__[1],Eps__[2],Eps__[5],']', file=ff)
##                        print('XXX',Elem.Lch_,ipI,Elem.CrBwS,D,Hd,hdI,'\n[',Eps__[0],Eps__[1],Eps__[2],Eps__[5],']\n[',sig[0], sig[1], sig[2], sig[5],']', file=ff)
##                            print('XXX',(self.cc0+self.cc1*kap_/(2.0*J2S))*EpsD,'\n',self.cc2*kap_*nd_,'\n',self.cc3*kap_*array([1.,1.,1.,0.,0.,0.]),'\n',nd, file=ff)
##                            Material.Sig1_ = sig[0]
##                            Material.Sig2_ = sig[1]
##                            Material.Sig3_ = sig[2]
##                            Material.Eps1_ = Eps_[0]
##                            Material.CC_   = CC_[0,0]
##                            if Dps_[0]>ZeroD: Material.C1_ = (CC[0,0]*Dps_[0] + 2.*CC[0,1]*Dps_[1])/Dps_[0]
##                            else:             Material.C1_ = 0.
##                            Material.DEps_ = Dps_[0]
                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    return [sig[0],sig[1],sig[5]], CC_, [Eps__[0], Eps__[1], Eps__[2], Eps__[5], sig[0], sig[1], sig[2], sig[5],D,r] # r is random data
            elif Elem.dim == 4:
                CC_ = array([[CC[0,0],CC[0,1],CC[0,2],CC[0,5]],
                             [CC[1,0],CC[1,1],CC[1,2],CC[1,5]],
                             [CC[2,0],CC[2,1],CC[2,2],CC[2,5]],
                             [CC[5,0],CC[5,1],CC[5,2],CC[5,5]]])
#                if Elem.Label == 20 and ipI == 0: print('YYY', sig, '\n', CC)
#                print(Elem.Label, ipI, "%15.7e %15.7e %15.7e %15.7e"%(sig[0],sig[1],sig[2],sig[5]), file=ff)
                return [sig[0],sig[1],sig[2],sig[5]], CC_, [Eps__[0],Eps__[1],Eps__[2],Eps__[5],sig[0],sig[1],sig[2],sig[5],D,r]  # r is random data
            elif Elem.dim==3:
                if self.RType==1:                                           # gradient damage
                    ccc = 0.5*self.RegPar**2
                    fact = 0.5*self.Emod
                    CC_= array([[CC[0,0],CC[0,1],CC[0,2],CC[0,3],CC[0,4],CC[0,5],0,0,0],
                                [CC[1,0],CC[1,1],CC[1,2],CC[1,3],CC[1,4],CC[1,5],0,0,0],
                                [CC[2,0],CC[2,1],CC[2,2],CC[2,3],CC[2,4],CC[2,5],0,0,0],
                                [CC[3,0],CC[3,1],CC[3,2],CC[3,3],CC[3,4],CC[3,5],0,0,0],
                                [CC[4,0],CC[4,1],CC[4,2],CC[4,3],CC[4,4],CC[4,5],0,0,0],
                                [CC[5,0],CC[5,1],CC[5,2],CC[5,3],CC[5,4],CC[5,5],0,0,0],
                                [0,      0,      0,      0,      0,      0,      ccc*fact,0,0],
                                [0,      0,      0,      0,      0,      0,      0,ccc*fact,0],
                                [0,      0,      0,      0,      0,      0,      0,0,ccc*fact]])
                    CR = array([[ 0,             0,             0,             0,             0,             0,            -sig0[0]*hdI],
                                [ 0,             0,             0,             0,             0,             0,            -sig0[1]*hdI],
                                [ 0,             0,             0,             0,             0,             0,            -sig0[2]*hdI],
                                [ 0,             0,             0,             0,             0,             0,            -sig0[3]*hdI],
                                [ 0,             0,             0,             0,             0,             0,            -sig0[4]*hdI],
                                [ 0,             0,             0,             0,             0,             0,            -sig0[5]*hdI],
                                [-fact*nd[0]/Hd,-fact*nd[1]/Hd,-fact*nd[2]/Hd,-fact*nd[3]/Hd,-fact*nd[4]/Hd,-fact*nd[5]/Hd, fact]])
                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],ccc*Eps_[6]*fact,ccc*Eps_[7]*fact,ccc*Eps_[8]*fact], CC_, [0,0,0,0,0,0,(EpsR[6]-kap_)*fact], CR, [sig[0],sig[1],sig[2],Eps[0],Eps[1],Eps[2]]
                else:
                    return sig, CC, [Eps__[0],Eps__[1],Eps__[2],Eps__[3],Eps__[4],Eps__[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],D,r]
            elif Elem.dim==21:                                              # Continuum based shell
                CC_ = array([[CC[0,0],CC[0,1],CC[0,2],CC[0,3],CC[0,4],CC[0,5]],
                             [CC[1,0],CC[1,1],CC[1,2],CC[1,3],CC[1,4],CC[1,5]],
                             [0., 0., 0., 0., 0., 0.],
                             [CC[3,0],CC[3,1],CC[3,2],CC[3,3],CC[3,4],CC[3,5]],
                             [CC[4,0],CC[4,1],CC[4,2],CC[4,3],CC[4,4],CC[4,5]],
                             [CC[5,0],CC[5,1],CC[5,2],CC[5,3],CC[5,4],CC[5,5]]])
                if Elem.ShellRCFlag: return sig, CC_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps__[0],Eps__[1],Eps__[5], sig[0],sig[1],sig[5], D, svs]
#                if Elem.ShellRCFlag: return sig, CC_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps__[0],Eps__[1],Eps__[5], sig[0],sig[1],sig[5], D,kap_]
                else:                return sig, CC_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
            else: raise NameError ("ConFemMaterials::Isodam.Sig: not implemented for this element type")
        # C-Version
        else:
            DataOut = zeros((1),dtype=float)  # should have an eye on this whether if conforms to C-code
            if Elem.dim==1:
                sig  = zeros((2),dtype=float)
                sigR = zeros((2),dtype=float)
#                MatM_= zeros((4),dtype=float)
                MatM_= zeros((2,2),dtype=float)
                CR   = zeros((4),dtype=float)
                PlSt_ = False
            elif Elem.dim==2:
                sig  = zeros((3),dtype=float)  #zeros((5),dtype=float)
                sigR = zeros((1),dtype=float)  #zeros((4),dtype=float)
#                MatM_= zeros((9),dtype=float)  #zeros((25),dtype=float)
                MatM_ = zeros((3,3), dtype=float)
                CR   = zeros((1),dtype=float)  #zeros((16),dtype=float)
                PlSt_ = Elem.PlSt
            elif Elem.dim==3:
                sig  = zeros((6),dtype=float)  # zeros((9),dtype=float)
                sigR = zeros((1),dtype=float)  # zeros((7),dtype=float)
#                MatM_= zeros((36),dtype=float) # zeros((81),dtype=float)
                MatM_= zeros((6,6),dtype=float) # zeros((81),dtype=float)
                CR   = zeros((1),dtype=float)  # zeros((49),dtype=float)
                PlSt_ = False
            elif Elem.dim == 4:                                 # axisymmetric
                sig  = zeros((4), dtype=float)
                sigR = zeros((1), dtype=float)
#                MatM_= zeros((16), dtype=float)
                MatM_= zeros((4,4), dtype=float)
                CR = zeros((1), dtype=float)
                PlSt_ = False
            elif Elem.dim==21:
                sig  = zeros((6),dtype=float)  
                sigR = zeros((1),dtype=float)
#                MatM_= zeros((36),dtype=float)
                MatM_= zeros((6,6),dtype=float)
                CR   = zeros((1),dtype=float)
                PlSt_ = False
            if len(EpsR)==0: EpsR = zeros((1),dtype=float)
            rc = IsoDamC1( CalcType, Elem.dim, PlSt_, self.PrinStrains, Elem.Lch_, Elem.StateVar[ipI], Elem.StateVarN[ipI],\
                           Eps_, sig, MatM_, self.LiTy, cc0, cc1, cc2, cc3, self.RType, EpsR, self.kapStrength,\
                           Elem.CrBwS, self.gam2, self.kapUlt, self.edt, self.ed, self.gd, nu, Emod, Dps_, self.eta, self.RegPar, Elem.ScaleType,\
                           sigR, CR, Dt, self.svrTol, DataOut)
            D_   = Elem.StateVarN[ipI,0]
            kap_ = Elem.StateVarN[ipI,1]
#            print('YYY',Elem.Label,ipI,'\n',DataOut)
###############################                
#            if Elem.Label==6: print('XXX',ipI,D_,sig[0:3])
#            if Elem.Label==50 and ipI==0:
#                print 'ZZZ', Elem.Label, Eps_, sig, '__', sigR #  array([ [MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]])
###############################                
            if rc>110:
                raise NameError("ConFemMaterials::IsoDamage:sig:IsoDamC1 RC "+str(rc),Elem.Label,ipI, Eps_,D_)
            if Elem.dim==1:
                if self.RType==1:
#                    return [sig[0],sig[1]], array([[MatM_[0],MatM_[1]],[MatM_[2],MatM_[3]]]), [sigR[0],sigR[1]], array([[CR[0],CR[1]],[CR[2],CR[3]]]), [Eps_[0], sig[0], D_]
                    return sig, MatM_, [sigR[0],sigR[1]], array([[CR[0],CR[1]],[CR[2],CR[3]]]), [Eps_[0], sig[0], D_]
                else:
#                    return [sig[0],0],      array([[MatM_[0],MatM_[1]],[MatM_[2],MatM_[3]]]),                                                          [Eps_[0], sig[0], D_]
                    return sig, MatM_,                                                          [Eps_[0], sig[0], D_]
            elif Elem.dim==2:
                if self.RType==1:
                    pass
                else:
#                    return ([sig[0],sig[1],sig[2]]), array([ [MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]]), [Eps_[0],Eps_[1],0.,Eps_[2],sig[0],sig[1],0.,sig[2],D_,r]
                    return sig, MatM_, [Eps_[0],Eps_[1],0.,Eps_[2],sig[0],sig[1],0.,sig[2],D_,r]
            elif Elem.dim==3:
                if self.RType==1:
                    pass
                else:
#                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                           [Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],D_,r]
                    return sig, MatM_, [Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],D_,r]
            elif Elem.dim==4:
                if self.RType==1:
                    pass
                else:
#                    return ([sig[0],sig[1],sig[2],sig[3]]), \
#                            array([ [MatM_[0], MatM_[1], MatM_[2], MatM_[3] ],
#                                    [MatM_[4], MatM_[5], MatM_[6], MatM_[7] ],
#                                    [MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                                    [MatM_[12],MatM_[13],MatM_[14],MatM_[15]]]),\
#                    [Eps_[0],Eps_[1],Eps_[2],Eps_[3],sig[0],sig[1],sig[2],sig[3],D_,r]
#                    if Elem.Label==20 and ipI==0: print('XXX',sig,'\n',MatM_)
#                    print(Elem.Label, ipI, "%15.7e %15.7e %15.7e %15.7e" % (sig[0], sig[1], sig[2], sig[3]), file=ff)
                    return sig, MatM_, [Eps_[0],Eps_[1],Eps_[2],Eps_[3],sig[0],sig[1],sig[2],sig[3],D_,r]
            elif Elem.dim==21:
                if Elem.ShellRCFlag:
#                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps_[0],Eps_[1],Eps_[5], sig[0],sig[1],sig[5], D_, DataOut[0]]
                    return sig, MatM_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps_[0],Eps_[1],Eps_[5], sig[0],sig[1],sig[5], D_, DataOut[0]]
                else:
#                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
                    return sig, MatM_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):
            if Elem.StateVarN[j,1]>Elem.StateVar[j,1]:                      # damage measures
                Elem.StateVar[j,0] = Elem.StateVarN[j,0]
                Elem.StateVar[j,1] = Elem.StateVarN[j,1]
            Elem.StateVar[j,8] = Elem.StateVar[j,8] + Elem.StateVarN[j,8]   # "crack" energy
            Elem.StateVar[j,2] = Elem.StateVarN[j,2]                        # strain velocities - used in ViscExten3D
            Elem.StateVar[j,3] = Elem.StateVarN[j,3]
            Elem.StateVar[j,4] = Elem.StateVarN[j,4]
            Elem.StateVar[j,5] = Elem.StateVarN[j,5]
            Elem.StateVar[j,6] = Elem.StateVarN[j,6]
            Elem.StateVar[j,7] = Elem.StateVarN[j,7]

            Elem.StateVar[j,9] = Elem.StateVarN[j,9]                        # SDA crack quantities
            Elem.StateVar[j,10]= Elem.StateVarN[j,10]
            Elem.StateVar[j,11]= Elem.StateVarN[j,11]
            Elem.StateVar[j,12]= Elem.StateVarN[j,12]
            Elem.StateVar[j,13]= Elem.StateVarN[j,13]
            Elem.StateVar[j,14]= Elem.StateVarN[j,14]
            Elem.StateVar[j,15]= Elem.StateVarN[j,15]
            
#            Elem.StateVar[j,16:22] = 0.                                     # for strain increments of equilibrium iteration for algorithmic modulus
            
        return False

#class MicroPlane(Material):                                                 # microplane damage old Version
#    def __init__(self, PropMat, RotCrack,PrinStrain,ShearRetFac, S2ndCrack,S2ndCrackA, f6):
##                        (self, SymmetricVal, RTypeVal,   UpdateVal, Updat2Val, StateVarVal, NDataVal):
#        Material.__init__(self, False,        PropMat[8], False,     False,     None,        8,       "MicroPl")
##        self.Symmetric = False                                              # flag for symmetry of material matrices
##        self.RType = PropMat[4]                                             # type of regularization 0: without, 1: gradient 2: crack band
##        self.Update = False                                                 # has specific update procedure for update of state variables
##        self.Updat2 = False                                                 # no 2 stage update procedure
##        self.StateVar = 1                                                   # see below
##        self.NDataVal = 8
##        self.Type = "MicroPl"
#        self.PropMat= PropMat                                                # for wrappers, e.g. RCSHELL#
#
#        self.EE = PropMat[0]                                                # macroscopic Young's modulus
#        self.Emod = PropMat[0]                                              # used for SDA
#        self.nu = PropMat[1]                                                # macroscopic Poissons's ratio
#        self.type = PropMat[2]                                              # type of damage function
#        #
#        self.nInt = len(I21Points)
#        self.nState = 1                                                     # number of state variables per integration direction of unit sphere
#        self.iS       = self.nState*self.nInt+5                             # entry index for strain rate for, e.g. viscous extension
#        self.StateVar = self.nState*self.nInt+5+6+7    # -> 39              # number of state variables per integration point of element; 1st extra for dd_iso, 2nd for I1, 3rd for J2, 
#                                                                            # 4th currently not used, 5th for eps_33 in case of plane state, 6-11 for strain increment
#                                                                            # 12 largest principal value, 13-15 largest principal strain direction, 16-18 traction in largest principal strain direction 
#        ff = self.EE / ( ( 1. + self.nu ) * ( 1. - 2.*self.nu ) )
#        self.EMat = ff*array([[1.-self.nu,self.nu,self.nu,0,0,0],
#                              [self.nu,1.-self.nu,self.nu,0,0,0],
#                              [self.nu,self.nu,1.-self.nu,0,0,0],
#                              [0,0,0,(1.-2.*self.nu)/2.,0,0],
#                              [0,0,0,0,(1.-2.*self.nu)/2.,0],
#                              [0,0,0,0,0,(1.-2.*self.nu)/2.]])              # 3D isotropic linear elasticity
#        KK = self.EE/(3.*(1.-2.*self.nu))                                   # macroscopic bulk modulus
#        GG = self.EE/(2.*(1.+self.nu))                                      # macroscopic shear modulus
#        # V-D split
#        self.E_V= 3.*KK                                                     # moduli of microplane elasticity V-D split
#        self.E_D= 2.*GG
#        self.E_DM = array([[self.E_D,0,0],[0,self.E_D,0],[0,0,self.E_D]])
#        self.PlStressL = 0.01                                               # limit for plane stress iteration (dim = 2 (plane stress ), 21)
#        self.PlStressI = 10                                                 # max number of iterations
#        #
#        if self.type == 1:                                                  # V-D split Vree / Leukart damage function 
#            self.alphaV = PropMat[3] #0.9 #0.96
#            self.betaV  = PropMat[4] #3000. #300.
#            self.kap0V  = PropMat[5] #0.0001 #0.0005
#        elif self.type == 2:                                                # V-D split Vree / uhc damage function
#            self.alphaV = 0.
#            self.betaV  = 0.
#            self.kap0V  = 0.
#            self.fct= PropMat[3]                                            # tensile strength (PropMat[5,6] not used
#            self.eps_ct = 1.648721271*self.fct/self.EE                      # strain for uniaxial tensile strength to make e_0=0
#            self.gd = 2.
#            self.e0 = 0.                                                    # with eps_ct as above
##            self.e0 = 1.                                                    # to make it elastic for control purposes
#            self.ed = 2.331643981*PropMat[3]/self.EE                        # with eps_ct as above
#            self.gam2 = 12.*350.                                            # parameter for scaling of equivalent damage strain regularization
#            self.SpecCrEn = self.SpeCrEnergySample( 1.0, self.eps_ct, False,None,None, self.UniaxTensionUnit) # specific crack energy unscaled
#            self.fc = 0.7*self.fct*PropMat[4]                               # rough estimation of uniaxial compressive strength - required for pressure dependent bond
#            self.RType  = PropMat[7]                                        # indicator for regularization
#            self.RegPar = PropMat[8]                                        # crack energy for regularization (RType=2)
#            if self.RType==2: self.bw = self.RegPar/self.SpecCrEn 
#            else:             self.bw = 0.
##            Echo(f"Microplane: fct {self.fct:f}, fc {self.fc:f}, SpecCrEn {self.SpecCrEn:f}, bw {self.bw:f}", f6)
#            self.CrX, self.CrY = self.SpeCrEnergyData( 0.01, 100., self.CrBwN, self.eps_ct, self.UniaxTensionScaled) # arrays to determine scaling factor for given characteristic element length
#            self.CrX2,self.CrY2= self.SpeCrEnergyData( 0.01, 0.5,  self.CrBwN, self.eps_ct, self.UniaxTensionScaled2) # arrays to determine scaling factor for given characteristic element length
##            self.dDestr   = 1.-1.e-3                                        # for maximum damage allowed, subtractor should not be small in order to avoid problems in dividing by MMatS[2,2], see below  
#            self.dDestr   = 1.-1.e-4                                        # for maximum damage allowed, subtractor should not be small in order to avoid problems in dividing by MMatS[2,2], see below  
#            self.kapUlt   = exp(log(-log(1.-self.dDestr))/self.gd)*self.ed+self.e0
#        #
#        self.kV =         PropMat[4]                                        # (measure for) ratio of uniaxial compressive strength to tensile strength
#        self.kV0 = (self.kV-1.)/(2.*self.kV*(1.-2.*self.nu))
#        self.kV1 = self.kV0
#        self.kV2 = 3./(self.kV*(1.+self.nu)**2)
#        self.eta =        PropMat[9]                                        # for viscous regularization
#        self.etaCr = PropMat[10]                                            # artificial viscosity for SDA crack
#        self.svrTol = 0.01                                                  # minimum reference value for viscous stresses - for control computations only
#        self.ElCharLenBound = PropMat[11]                                   # related bound to distinguish between types of crack band regularization
#        self.Density =    PropMat[12]
#        self.RotCrack    = RotCrack                                         # rotating crack for CPS4 using this material      in case of SDA 
#        self.PrinStrains = PrinStrain                                       # SDA rankine criterion with principal strains / stresses 
#        self.S2ndCrack   = S2ndCrack                                        # allows for second crack
#        if self.S2ndCrack: self.S2ndCrackA = cos(S2ndCrackA*pi/180.)        # cos(deviation 2nd crack  from 1st crack in rad) should be smaller than this value, i.e. should be more or less perpendicular
#        if ShearRetFac > 0: self.ShearRetFac = ShearRetFac
#        else:               self.ShearRetFac = 1.0e-6
#        self.RankineFactor = 1.0
#    def UniaxTensionScaled(self, gam1, eps):
#        xx = eps - self.eps_ct
#        kap_ =gam1*xx + (1.-gam1)/(self.gam2) * (1.-exp(-self.gam2*xx)) + self.eps_ct
#        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps
#    def UniaxTensionScaled2(self, beta, eps):
#        kap = eps 
#        kap_=  (1-beta)*self.eps_ct * log( ( kap-beta*self.eps_ct )/( (1-beta)*self.eps_ct) ) + self.eps_ct
#        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps
#    def UniaxTensionUnit(self, dummy, eps):
#        kap_ = eps
#        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps
#    def C3(self, Emod, nu):
#        ff = Emod / ( ( 1. + nu ) * ( 1. - 2.*nu ) )
#        return ff*array([[1.-nu,nu,nu,0,0,0],[nu,1.-nu,nu,0,0,0],[nu,nu,1.-nu,0,0,0],[0,0,0,(1.-2.*nu)/2.,0,0],[0,0,0,0,(1.-2.*nu)/2.,0],[0,0,0,0,0,(1.-2.*nu)/2.]])
#    def DamFunc1(self, kapOld, eta):
#        kap = max(self.kap0V,kapOld)
#        if eta>kap: kap=eta
#        dd = 1.-self.kap0V/kap*(1.+self.alphaV*( exp( self.betaV*(self.kap0V-kap) )-1. ))
#        if dd>ZeroD: Dp = self.kap0V/kap**2*(1+self.alphaV*(exp(self.betaV*(self.kap0V-kap))-1.))+self.kap0V/kap*self.alphaV*self.betaV*exp(self.betaV*(self.kap0V-kap))
#        else:        Dp = 0.
#        return kap, dd, Dp
##    def DamFunc2(self, kapOld, eta):
##        kap = max(self.e0,kapOld)
##        if eta>kap: kap=eta
##        if kap<=self.e0: dd = 0.
##        else:            dd = 1.-exp(-((kap-self.e0)/self.ed)**self.gd)
##        if dd>ZeroD:     Dp = (1.-dd) * self.gd/self.ed**self.gd * (kap-self.e0)**(self.gd-1.) 
##        else:            Dp = 0.
##        return kap, dd, Dp
#    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
##
##        if Elem.Label==57 and ipI==0: Echo(f"AAA {Eps[0]:12.6e}, {Eps[1]:12.6e}, {Eps[2]:12.6e} ", ff)#
##
#        def DevStiffness():
#            DVD = zeros((6,6),dtype=float)
#            DVD[0,0] = m00*n02*tbt2+2.*m01*nn[0]*tbt*nn[1]*obt+2.*m02*nn[0]*obt*nn[2]*tbt+m11*n12*obt2+2.*m12*nn[1]*obt2*nn[2]+m22*n22*obt2
#            DVD[0,1] = m00*n02*tbt*obt+m01*nn[0]*tbt2*nn[1]+m02*nn[0]*obt*nn[2]*tbt+m01*nn[0]*obt2*nn[1]+m11*n12*obt*tbt+m12*nn[1]*obt2*nn[2]+m02*nn[0]*obt2*nn[2]+m12*nn[1]*obt*nn[2]*tbt+m22*n22*obt2
#            DVD[0,2] = m00*n02*tbt*obt+m01*nn[0]*tbt*nn[1]*obt+m02*nn[0]*tbt2*nn[2]+m01*nn[0]*obt2*nn[1]+m11*n12*obt2+m12*nn[1]*obt*nn[2]*tbt+m02*nn[0]*obt2*nn[2]+m12*nn[1]*obt2*nn[2]+m22*n22*obt*tbt 
#            DVD[0,3] = -0.50*m01*nn[0]*tbt*nn[2]-0.50*m02*nn[0]*tbt*nn[1]-0.50*m11*nn[2]*nn[1]*obt-0.50*m12*n12*obt-0.50*m12*n22*obt-0.50*m22*nn[2]*obt*nn[1] 
#            DVD[0,4] = -0.50*m00*nn[0]*tbt*nn[2]-0.50*m02*n02*tbt-0.50*m01*nn[2]*nn[1]*obt-0.50*m12*nn[0]*obt*nn[1]-0.50*m02*n22*obt-0.50*m22*nn[0]*obt*nn[2]
#            DVD[0,5] = -0.50*m00*nn[0]*tbt*nn[1]-0.50*m01*n12*obt-0.50*m02*nn[2]*nn[1]*obt-0.50*m01*n02*tbt-0.50*m11*nn[0]*obt*nn[1]-0.50*m12*nn[0]*obt*nn[2]
#            DVD[1,1] = m00*n02*obt2+2.*m01*nn[0]*tbt*nn[1]*obt+2.*m02*nn[0]*obt2*nn[2]+m11*n12*tbt2+2.*m12*nn[1]*obt*nn[2]*tbt+m22*n22*obt2
#            DVD[1,2] = m00*n02*obt2+m01*nn[0]*obt2*nn[1]+m02*nn[0]*obt*nn[2]*tbt+m01*nn[0]*tbt*nn[1]*obt+m11*n12*obt*tbt+m12*nn[1]*tbt2*nn[2]+m02*nn[0]*obt2*nn[2]+m12*nn[1]*obt2*nn[2]+m22*n22*obt*tbt
#            DVD[1,3] = -0.50*m01*nn[0]*obt*nn[2]-0.50*m02*nn[0]*obt*nn[1]-0.50*m11*nn[1]*tbt*nn[2]-0.50*m12*n12*tbt-0.50*m12*n22*obt-0.50*m22*nn[2]*obt*nn[1]
#            DVD[1,4] = -0.50*m00*nn[0]*obt*nn[2]-0.50*m02*n02*obt-0.50*m01*nn[2]*nn[1]*tbt-0.50*m12*nn[1]*tbt*nn[0]-0.50*m02*n22*obt-0.50*m22*nn[0]*obt*nn[2]
#            DVD[1,5] = -0.50*m00*nn[0]*obt*nn[1]-0.50*m01*n12*tbt-0.50*m02*nn[2]*nn[1]*obt-0.50*m01*n02*obt-0.50*m11*nn[1]*tbt*nn[0]-0.50*m12*nn[0]*obt*nn[2]
#            DVD[2,2] = m00*n02*obt2+2.*m01*nn[0]*obt2*nn[1]+2.*m02*nn[0]*obt*nn[2]*tbt+m11*n12*obt2+2.*m12*nn[1]*obt*nn[2]*tbt+m22*n22*tbt2
#            DVD[2,3] = -0.50*m01*nn[0]*obt*nn[2]-0.50*m02*nn[0]*obt*nn[1]-0.50*m11*nn[2]*nn[1]*obt-0.50*m12*n12*obt-0.50*m12*n22*tbt-0.50*m22*nn[2]*tbt*nn[1]
#            DVD[2,4] = -0.50*m00*nn[0]*obt*nn[2]-0.50*m02*n02*obt-0.50*m01*nn[2]*nn[1]*obt-0.50*m12*nn[0]*obt*nn[1]-0.50*m02*n22*tbt-0.50*m22*nn[0]*nn[2]*tbt
#            DVD[2,5] = -0.50*m00*nn[0]*obt*nn[1]-0.50*m01*n12*obt-0.50*m02*nn[1]*nn[2]*tbt-0.50*m01*n02*obt-0.50*m11*nn[0]*obt*nn[1]-0.50*m12*nn[0]*nn[2]*tbt
#            DVD[3,3] = 0.25*m11*n22+0.50*m12*nn[1]*nn[2]+0.25*m22*n12
#            DVD[3,4] = 0.25*m01*n22+0.25*m12*nn[0]*nn[2]+0.25*m02*nn[1]*nn[2]+0.25*m22*nn[0]*nn[1]
#            DVD[3,5] = 0.25*m01*nn[1]*nn[2]+0.25*m02*n12+0.25*m11*nn[0]*nn[2]+0.25*m12*nn[0]*nn[1]
#            DVD[4,4] = 0.25*m00*n22+0.50*m02*nn[0]*nn[2]+0.25*m22*n02
#            DVD[4,5] = 0.25*m00*nn[1]*nn[2]+0.25*m02*nn[0]*nn[1]+0.25*m01*nn[0]*nn[2]+0.25*m12*n02
#            DVD[5,5] = 0.25*m00*n12+0.50*m01*nn[0]*nn[1]+0.25*m11*n02
#            DVD[1,0] = DVD[0,1]
#            DVD[2,0] = DVD[0,2]
#            DVD[2,1] = DVD[1,2]
#            DVD[3,0] = DVD[0,3]
#            DVD[3,1] = DVD[1,3]
#            DVD[3,2] = DVD[2,3]
#            DVD[4,0] = DVD[0,4]
#            DVD[4,1] = DVD[1,4]
#            DVD[4,2] = DVD[2,4]
#            DVD[4,3] = DVD[3,4]
#            DVD[5,0] = DVD[0,5]
#            DVD[5,1] = DVD[1,5]
#            DVD[5,2] = DVD[2,5]
#            DVD[5,3] = DVD[3,5]
#            DVD[5,4] = DVD[4,5]
#            return DVD
#        if CalcType == 0: return [], [], []
#        if not ConFemMatCFlag: # False: # flag for not C-version
##        if True:
#            obt, obs, obn, tbt = 1./3., 1./6., 1./9., -2./3.
#            obt2 = obt**2
#            tbt2 = tbt**2
#            VV   = array([[obt,0,0],[0,obt,0],[0,0,obt]])                                       # projection tensor for volumetric part
#            ep2  = Elem.StateVarN[ipI,self.nState*self.nInt+4]                                  # value of last iteration taken, not from last converged step
#            if Elem.dim==2:
#                if Elem.PlSt: 
#                    epsT = array([ [Eps[0],0.5*Eps[2],0.] , [0.5*Eps[2],Eps[1],0.] , [0.,0.,ep2] ])
#                    Dps__ = array([ Dps[0], Dps[1], 0., 0., 0., Dps[2]])                        # plane stress --> strain Voigt notation --  used for viscous regularization only with diagonal stiffness
#                else:         
#                    epsT = array([ [Eps[0],0.5*Eps[2],0.] , [0.5*Eps[2],Eps[1],0.] , [0.,0.,0.] ])        #
#                    Dps__ = array( [Dps[0], Dps[1], 0., 0., 0., Dps[2] ])                       # plane strain --> strain Voigt notation
#            elif Elem.dim==3: 
#                epsT  = array([ [Eps[0],0.5*Eps[5],0.5*Eps[4]] , [0.5*Eps[5],Eps[1],0.5*Eps[3]] , [0.5*Eps[4],0.5*Eps[3],Eps[2] ]]) # strain tensor arrangement
#                Dps__ = array([ Dps[0], Dps[1], Dps[2], Dps[3], Dps[4], Dps[5]])                # triaxial strain increment Voigt notation
#            elif Elem.dim==21:                                                                  # continuum based shell
#                epsT  = array([ [Eps[0],0.5*Eps[5],0.5*Eps[4]] , [0.5*Eps[5],Eps[1],0.5*Eps[3]] , [0.5*Eps[4],0.5*Eps[3],ep2] ]) # strain tensor arrangement
#                Dps__ = array([ Dps[0], Dps[1], 0.,     Dps[3], Dps[4], Dps[5]])                # triaxial strain increment for continuum bases shell Voigt notation
#            else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented")
#            VVV = zeros((6,6),dtype=float)
#            VVV[0,0] = obn
#            VVV[0,1] = obn
#            VVV[0,2] = obn
#            VVV[1,1] = obn
#            VVV[1,2] = obn
#            VVV[2,2] = obn
#            VVV[1,0], VVV[2,0], VVV[2,1] = VVV[0,1], VVV[0,2], VVV[1,2]
##            epsVol = VV[0,0]*epsT[0,0]+VV[1,1]*epsT[1,1]+VV[2,2]*epsT[2,2]
#            ns = self.nState
#            PlStFlag = False
#            if (Elem.dim==2 and Elem.PlSt) or Elem.dim==21: PlStFlag = True
#            #
#            for ii in range(self.PlStressI):                                # for plane stress iteration
#                epsVol = VV[0,0]*epsT[0,0]+VV[1,1]*epsT[1,1]+VV[2,2]*epsT[2,2]  # moved from above due to update of epsT[2,2]
#                dd_iso, I1, J2, bufD = 0., 0., 0., []
#                MMatS, MMatT, DV1 = zeros((6,6),dtype=float), zeros((6,6),dtype=float), zeros((6,6),dtype=float)
#                sig = zeros((3,3),dtype=float)                              # stresses in tensor notation
##                if ipI==0: print('YYY',f"{epsT[0,0]:12.8f},{epsT[1,1]:12.8f},{epsT[2,2]:12.8f},{epsT[0,1]:12.8f}")
#                for i in range(self.nInt):                                  # microplane integration order
#                    kapOld = Elem.StateVar[ipI,ns*i]
#                    nn = array([I21Points[i,0],I21Points[i,1],I21Points[i,2]])
#                    # V-D-Split projection tensors
#                    D0 =     array([[nn[0]    -nn[0]*obt, 0.5*nn[1], 0.5*nn[2]],
#                                    [0.5*nn[1]         ,-nn[0]*obt , 0.       ],
#                                    [0.5*nn[2]         , 0.       , -nn[0]*obt]])
#                    D1 =     array([[-nn[1]*obt         , 0.5*nn[0], 0.       ],
#                                    [0.5*nn[0], nn[1]-nn[1]*obt    , 0.5*nn[2]],
#                                    [0.       , 0.5*nn[2]         , -nn[1]*obt]])
#                    D2 =     array([[-nn[2]*obt        ,       0.  , 0.5*nn[0]],
#                                    [0.       ,      -nn[2]*obt    , 0.5*nn[1]],
#                                    [0.5*nn[0], 0.5*nn[1],     nn[2]-nn[2]*obt]])
#                    # vector deviator strains
#                    epsDD = array([D0[0,0]*epsT[0,0]+D0[0,1]*epsT[0,1]+D0[0,2]*epsT[0,2]+\
#                                   D0[1,0]*epsT[1,0]+D0[1,1]*epsT[1,1]+D0[1,2]*epsT[1,2]+\
#                                   D0[2,0]*epsT[2,0]+D0[2,1]*epsT[2,1]+D0[2,2]*epsT[2,2],\
#                                   D1[0,0]*epsT[0,0]+D1[0,1]*epsT[0,1]+D1[0,2]*epsT[0,2]+\
#                                   D1[1,0]*epsT[1,0]+D1[1,1]*epsT[1,1]+D1[1,2]*epsT[1,2]+\
#                                   D1[2,0]*epsT[2,0]+D1[2,1]*epsT[2,1]+D1[2,2]*epsT[2,2],\
#                                   D2[0,0]*epsT[0,0]+D2[0,1]*epsT[0,1]+D2[0,2]*epsT[0,2]+\
#                                   D2[1,0]*epsT[1,0]+D2[1,1]*epsT[1,1]+D2[1,2]*epsT[1,2]+\
#                                   D2[2,0]*epsT[2,0]+D2[2,1]*epsT[2,1]+D2[2,2]*epsT[2,2]])
#                    # microplane strain invariants, state variable
#                    I1mp = 3.*epsVol
#                    J2mp = 3./2.*dot(epsDD,epsDD)
#                    eta_ = self.kV0*I1mp + sqrt( (self.kV1*I1mp)**2 + self.kV2*J2mp )
#                    # damage functions
#                    if   self.type==1: 
#                        kap, dd, Dp = self.DamFunc1( kapOld, eta_ )
#                    elif self.type==2: 
#                        if self.RType==2:                                                           # crack band regularization 
#                            if eta_>self.eps_ct:                                                    # limit strain exceeded
#                                beta= Elem.CrBwS
#                                if Elem.ScaleType==1:
#                                    xx  = eta_ - self.eps_ct
#                                    eta = beta*xx + (1.-beta)/(self.gam2) * (1.-exp(-self.gam2*xx)) + self.eps_ct
#                                    dkk = beta +    (1.-beta)             *     exp(-self.gam2*xx)
#                                else:
#                                    eta = (1-beta)*self.eps_ct * log( ( eta_-beta*self.eps_ct )/( (1-beta)*self.eps_ct) ) + self.eps_ct
#                                    dkk = (1-beta)*self.eps_ct/(eta_-beta*self.eps_ct)
#                            else:                                                                   # below limit strain
#                                eta = eta_
#                                dkk = 1.
#                        else:                                                                       # no regularization
#                            eta = eta_
#                            dkk = 1.
#                        if eta>self.kapUlt: eta = self.kapUlt                                       # to have some residual stiffness
#                        # damage function
#                        kap = max( self.e0, kapOld, eta)
#                        if kap<=self.e0: 
#                            dd = 0.
#                            Dp = 0.
#                        else:            
#                            dd = 1.-exp(-((kap-self.e0)/self.ed)**self.gd)
#                            Dp = (1.-dd) * self.gd/self.ed**self.gd * (kap-self.e0)**(self.gd-1.) 
#                        Dp = Dp*dkk                                                                 # regularization scaling of derivative dD/dkap tangential stiffness
#                    # stresses
#                    sigVV = (1.-dd)*    self.E_V*epsVol                                             # volumetric stress (scalar)
#                    sigDD = (1.-dd)*dot(self.E_DM,epsDD)                                            # deviatoric stress (vector)
#                    sig   = sig + 6.*I21Weights[i]*(sigVV*VV+sigDD[0]*D0+sigDD[1]*D1+sigDD[2]*D2)   # stress tensor integration # factor 6. calibrates the weighting factors
#                    # secant material stiffness
#                    n02 = nn[0]**2
#                    n12 = nn[1]**2
#                    n22 = nn[2]**2
#                    s00, s01, s02 = 1., 1., 1. 
#                    m00, m11, m22 = 1., 1., 1. 
#                    m01, m02, m12 = 0., 0., 0.
#                    DDD = DevStiffness()
#                    MMatS = MMatS + 6.*I21Weights[i]*( (1.-dd)*self.E_V*VVV + (1.-dd)*self.E_D*DDD ) # presumably symmetric, i.e. only upper right part builded
#                    # corrector for tangential material stiffness
#                    if kap>kapOld and dd>ZeroD:
#                        s00, s01, s02 = epsDD[0], epsDD[1], epsDD[2]
#                        m00, m11, m22 = s00**2,   s01**2,   s02**2
#                        m01, m02, m12 = s00*s01,  s00*s02,  s01*s02
#                        DDD = DevStiffness()
#                        DV1[0,0] = -obt*s00*nn[0]*tbt-obt2*s01*nn[1]-obt2*s02*nn[2] 
#                        DV1[0,1] = DV1[0,0] 
#                        DV1[0,2] = DV1[0,0]
#                        DV1[1,0] = -obt2*s00*nn[0]-obt*s01*nn[1]*tbt-obt2*s02*nn[2]
#                        DV1[1,1] = DV1[1,0]
#                        DV1[1,2] = DV1[1,0]
#                        DV1[2,0] = -obt2*s00*nn[0]-obt2*s01*nn[1]-obt*s02*nn[2]*tbt
#                        DV1[2,1] = DV1[2,0]
#                        DV1[2,2] = DV1[2,0]
#                        DV1[3,0] = obs*s01*nn[2]+obs*s02*nn[1]
#                        DV1[3,1] = DV1[3,0]
#                        DV1[3,2] = DV1[3,0]
#                        DV1[4,0] = obs*s00*nn[2]+obs*s02*nn[0]
#                        DV1[4,1] = DV1[4,0]
#                        DV1[4,2] = DV1[4,0]
#                        DV1[5,0] = obs*s00*nn[1]+obs*s01*nn[0]
#                        DV1[5,1] = DV1[5,0]
#                        DV1[5,2] = DV1[5,0]
#                        xx = sqrt( self.kV1**2 * I1mp**2 + self.kV2 * J2mp )
#                        thet =  9.*self.kV1**2/xx
#                        psi  =  1.5*self.kV2/xx
#                        alph = (3.*self.kV0 + thet*epsVol)*self.E_V*epsVol
#                        bet  = (3.*self.kV0 + thet*epsVol)*self.E_D
#                        gam  = psi*self.E_V*epsVol
#                        delt = psi*self.E_D
#                        MMatT = MMatT - 6.*I21Weights[i]*Dp*( alph * VVV + bet * DV1 + gam * transpose(DV1) + delt * DDD )
#                    #
##                    bufD += [dd]            
#                    dd_iso += 2.*I21Weights[i] * dd
#                    I1 += 2.*I21Weights[i] * I1mp
#                    J2 += 2.*I21Weights[i] * J2mp
#                    # microplane state variable updates
#                    Elem.StateVarN[ipI,ns*i] = kap
#                    
##                ep2 = -(MMatS[2,0]*epsT[0,0]+MMatS[2,1]*epsT[1,1]+MMatS[2,3]*epsT[1,2]*0.5+MMatS[2,4]*epsT[0,2]*0.5+MMatS[2,5]*epsT[0,1]) / MMatS[2,2] # lateral normal strain in case of zero lateral stress
#                ep2 = -( MMatS[2,0]*epsT[0,0] + MMatS[2,1]*epsT[1,1] + 2.0*MMatS[2,5]*epsT[0,1] ) / MMatS[2,2] # lateral normal strain in case of zero lateral stress
#                if abs(epsT[2,2])>ZeroD*1.e-3: xyz = abs((epsT[2,2]-ep2)/epsT[2,2])
#                else:                          xyz = 0.
#                if PlStFlag: epsT[2,2] = ep2
#                if not PlStFlag or xyz<self.PlStressL: break
##                if PlStFlag: epsT[2,2] = ep2
#    
#            if ii >= self.PlStressI-1: raise NameError("ConFemMaterials::Microplane.Sig: plane stress iteration", ii, sig[2,2])
#    #        I1_ = epsT[0,0]+epsT[1,1]+epsT[2,2]         # 1st invariant of strain tensor
#    #        epsTD = array([[Eps[0]-I1/3.,0.5*Eps[5],  0.5*Eps[4]],[0.5*Eps[5],   Eps[1]-I1/3.,0.5*Eps[3]],[0.5*Eps[4],   0.5*Eps[3],  Eps[2]-I1/3.]])
#    #        J2_ = 0.5*(epsTD[0,0]**2+epsTD[1,1]**2+epsTD[2,2]**2 + 2.*epsTD[0,1]**2 + 2.*epsTD[0,2]**2 + 2.*epsTD[1,2]**2)  # 2nd Invariant of strain deviator
#    #        print 'X', (I1-I1_)/I1, (J2-J2_)/J2
#    
#            Elem.StateVarN[ipI,self.nState*self.nInt] = dd_iso
#            Elem.StateVarN[ipI,self.nState*self.nInt+1] = I1
#            Elem.StateVarN[ipI,self.nState*self.nInt+2] = J2
#            Elem.StateVarN[ipI,self.nState*self.nInt+4] = ep2
#            MatM = MMatS + MMatT
#            #
#            #print('XXX',bufD)
#            #
#            if self.RType==3:                                               # before viscous contribution
#                if self.PrinStrains: laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([epsT[0,0],epsT[1,1],epsT[2,2],epsT[1,2],epsT[0,2],epsT[0,1]]) )
#                else:                laMax, eVec = self.Rankine( ConFemMatCFlag, Elem, array([ sig[0,0], sig[1,1], sig[2,2], sig[1,2], sig[0,2], sig[0,1]]) )
#            #
#            if self.eta>0.: 
#                zz, Veps = self.ViscExten3D( Dt, self.eta, Dps__, Elem, ipI, self.iS)           # for viscous regularization
#                sig[0,0] = sig[0,0] + self.eta*Veps[0]                                          # Voigt notation
#                sig[1,1] = sig[1,1] + self.eta*Veps[1]
#                sig[2,2] = sig[2,2] + self.eta*Veps[2]
#                sig[1,2] = sig[1,2] + self.eta*Veps[3]
#                sig[0,2] = sig[0,2] + self.eta*Veps[4]
#                sig[0,1] = sig[0,1] + self.eta*Veps[5]
#                for i in range(6): MatM[i,i] = MatM[i,i] + zz 
#            #
#            if self.RType==3: self.RankineUpd( Elem, ipI, array([sig[0,0], sig[1,1], sig[2,2], sig[1,2], sig[0,2], sig[0,1]]), eVec, laMax, self.iS) # iS = 26                             
#            #
#            if Elem.dim==2:                                                                     # 2D
#                if Elem.PlSt:                                                                   # plane stress
#                    cD = 1./MatM[2,2]
#                    MatM_ =array([[MatM[0,0]-MatM[0,2]*MatM[2,0]*cD,MatM[0,1]-MatM[0,2]*MatM[2,1]*cD,MatM[0,5]-MatM[0,2]*MatM[2,5]*cD],
#                                  [MatM[1,0]-MatM[1,2]*MatM[2,0]*cD,MatM[1,1]-MatM[1,2]*MatM[2,1]*cD,MatM[1,5]-MatM[1,2]*MatM[2,5]*cD], 
#                                  [MatM[5,0]-MatM[5,2]*MatM[2,0]*cD,MatM[5,1]-MatM[5,2]*MatM[2,1]*cD,MatM[5,5]-MatM[5,2]*MatM[2,5]*cD]])
#                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                    if ipI == 0:
#                        Material.Sig1_ = sig[0,0]
#                        Material.Sig2_ = sig[1,1]
#                        Material.Sig3_ = sig[2,2]
#                        Material.Eps1_ = Eps[0]
#                        Material.CC_   = MatM_[0,0]
#                        if Dps[0]>ZeroD: Material.C1_ = (MatM_[0,0]*Dps[0] + 2.*MatM_[0,1]*Dps[1])/Dps[0]
#                        else:            Material.C1_ = 0.
#                        Material.DEps_ = Dps[0]
#                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                else:                                                                           # plane strain
#                    MatM_ =array([[MatM[0,0],MatM[0,1],MatM[0,5]],[MatM[1,0],MatM[1,1],MatM[1,5]],[MatM[5,0],MatM[5,1],MatM[5,5]]])
##               return array([sig[0,0],sig[1,1],sig[0,1]]), MatM_, [Eps[0], Eps[1], Eps[2], sig[0,0], sig[1,1], sig[0,1]]
##
##                if Elem.Label==57 and ipI==0: print("BBB\n", sig[0,0],sig[1,1],sig[2,2],sig[1,2],'\n',MatM_)
##                
#                return array([sig[0,0],sig[1,1],sig[0,1]]), MatM_, [Eps[0], Eps[1], ep2, Eps[2], sig[0,0], sig[1,1], sig[2,2], sig[0,1]]
#            elif Elem.dim==3:                                                                   # 3D
#                return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]), MatM, [sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
#            elif Elem.dim==21:                                                                  # Continuum based shell
#                if abs(MatM[2,2])<ZeroD: 
#                    print('X\n', MMatS, '\n', MatM) 
#                    raise NameError("ConFemMaterials::Microplane.Sig: tangential Mat[2,2] to small")
#                cD = 1./MatM[2,2]
#                MatM_ = array([[MatM[0,0]-MatM[0,2]*MatM[2,0]*cD, MatM[0,1]-MatM[0,2]*MatM[2,1]*cD, 0., MatM[0,3]-MatM[0,2]*MatM[2,3]*cD, MatM[0,4]-MatM[0,2]*MatM[2,4]*cD, MatM[0,5]-MatM[0,2]*MatM[2,5]*cD],
#                               [MatM[1,0]-MatM[1,2]*MatM[2,0]*cD, MatM[1,1]-MatM[1,2]*MatM[2,1]*cD, 0., MatM[1,3]-MatM[1,2]*MatM[2,3]*cD, MatM[1,4]-MatM[1,2]*MatM[2,4]*cD, MatM[1,5]-MatM[1,2]*MatM[2,5]*cD],
#                               [0.,                               0.,                               0., 0.,                               0.,                               0.],
#                               [MatM[3,0]-MatM[3,2]*MatM[2,0]*cD, MatM[3,1]-MatM[3,2]*MatM[2,1]*cD, 0., MatM[3,3]-MatM[3,2]*MatM[2,3]*cD, MatM[3,4]-MatM[3,2]*MatM[2,4]*cD, MatM[3,5]-MatM[3,2]*MatM[2,5]*cD],
#                               [MatM[4,0]-MatM[4,2]*MatM[2,0]*cD, MatM[4,1]-MatM[4,2]*MatM[2,1]*cD, 0., MatM[4,3]-MatM[4,2]*MatM[2,3]*cD, MatM[4,4]-MatM[4,2]*MatM[2,4]*cD, MatM[4,5]-MatM[4,2]*MatM[2,5]*cD],
#                               [MatM[5,0]-MatM[5,2]*MatM[2,0]*cD, MatM[5,1]-MatM[5,2]*MatM[2,1]*cD, 0., MatM[5,3]-MatM[5,2]*MatM[2,3]*cD, MatM[5,4]-MatM[5,2]*MatM[2,4]*cD, MatM[5,5]-MatM[5,2]*MatM[2,5]*cD]])
#                if Elem.ShellRCFlag: return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]), MatM_, [sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1], Eps[0], Eps[1], Eps[2], sig[0,0], sig[1,1], sig[0,1], 0.,0.]
#                else:                return array([sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]), MatM_, [sig[0,0],sig[1,1],sig[2,2],sig[1,2],sig[0,2],sig[0,1]]
#            else: raise NameError ("ConFemMaterials::Microplane.Sig: not implemented for this element type")
#        # C-version
#        else:
#            DataOut  = zeros((1),dtype=float)
#            if Elem.dim==1:
#                sig  = zeros((2),dtype=float)
#                MatM_= zeros((4),dtype=float)
#                PlSt_ = False
#                raise NameError ("ConFemMaterials::Microplane.Sig: not implemented for this element type")
#            elif Elem.dim==2:
#                sig  = zeros((3),dtype=float)
#                MatM_= zeros((9),dtype=float)
#                PlSt_ = Elem.PlSt
#            elif Elem.dim==3:
#                sig  = zeros((6),dtype=float)  
#                MatM_= zeros((36),dtype=float)
#                PlSt_ = False
#            elif Elem.dim==21:
#                sig  = zeros((6),dtype=float)  
#                MatM_= zeros((36),dtype=float)
#                PlSt_ = False
#            rc = MicroPlaneC1( Elem.dim, PlSt_, self.PrinStrains, Elem.Lch_, Elem.StateVar[ipI], Elem.StateVarN[ipI],\
#                           Eps, sig, MatM_, self.type, self.E_V, self.E_D,\
#                           self.kV0, self.kV1, self.kV2, self.kap0V, self.alphaV, self.betaV, self.RType, Elem.CrBwS,self.gam2,self.kapUlt,Elem.ScaleType,\
#                           self.eps_ct, self.e0, self.ed, self.gd, self.nu, self.EE, Dps, self.eta,\
#                           Dt, self.svrTol, DataOut, self.nState, self.nInt, self.PlStressI, self.PlStressL, self.iS)
#            # not used: Elem.Lch_ but Elem.CrBwS for regularization
#            if rc> 110:  raise NameError("ConFemMaterials::MicroPlane:sig: RC "+str(rc))
#            if rc==110:  print(f"ConFemMaterials::MicroPlane:sig: no plain stress convergence El {Elem.Label:d}, IP {ipI:d}", file=ff)
#            D_   = Elem.StateVarN[ipI,self.nState*self.nInt]
#            if Elem.dim==1:
#                pass
##           elif Elem.dim==2:
##                
##                if ipI==0:
##                    print('YYY',f"{D_:12.9f}, {sig[0]:12.8f},{sig[1]:12.8f},{sig[2]:12.8f}")
##                    print('YYY',f"{D_:12.9f}, {sig[0]:12.8f},{sig[1]:12.8f},{sig[2]:12.8f}", file=ff)
##                    for jj in range(self.nInt): 
##                        ll_ = 7
##                        print(f"{int(DataOut[jj*ll_]):3d},{DataOut[jj*ll_+1]:16.9e},{DataOut[jj*ll_+2]:14.10f},{DataOut[jj*ll_+3]:14.10f},{DataOut[jj*ll_+4]:14.10f},{DataOut[jj*ll_+5]:14.10f},{DataOut[jj*ll_+6]:14.10f}",file=ff) 
##                    print('YYY',f"{sig[0]:14.10f},{sig[1]:14.10f},{sig[2]:14.10f}", file=ff)
##                
#                return ([sig[0],sig[1],sig[2]]), array([ [MatM_[0],MatM_[1],MatM_[2]],[MatM_[3],MatM_[4],MatM_[5]],[MatM_[6],MatM_[7],MatM_[8]]]), [Eps[0],Eps[1],0., Eps[2],sig[0],sig[1],0., sig[2]]
#            elif Elem.dim==3:
##                for i_ in range(len(DataOut)): print(i_,DataOut[i_], file=ff) 
##                if ipI==0: print(Elem.Label,ipI,D_,sig[0],sig[1], sig[2], file = ff)
#                return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                       [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                       [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                       [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                       [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                       [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                       [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
#            elif Elem.dim==21:
#                if Elem.ShellRCFlag:
#                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps[0],Eps[1],Eps[5], sig[0],sig[1],sig[5], D_]
#                else:
#                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
#                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
#                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
#                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
#                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
#                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
#                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
#                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
#    def UpdateStateVar(self, Elem, ff):
#        ns = self.nState
#        iS = self.iS                                                        # entry point for strain rate
#        for i in range(Elem.StateVar.shape[0]):
#            for j in range(self.nInt):                                      # loop over microplanes
#                kapOld = Elem.StateVar[i,ns*j]
#                kapNew = Elem.StateVarN[i,ns*j]
#                if kapNew>kapOld: Elem.StateVar[i,ns*j] = Elem.StateVarN[i,ns*j]
#            entry1 = ns*self.nInt
#            ddIsoOld = Elem.StateVar[i,entry1]
#            ddIsoNew = Elem.StateVarN[i,entry1]
#            if ddIsoNew>ddIsoOld: Elem.StateVar[i,entry1] = Elem.StateVarN[i,entry1] 
#            Elem.StateVar[i,entry1+1] = Elem.StateVarN[i,entry1+1]          # I1
#            Elem.StateVar[i,entry1+2] = Elem.StateVarN[i,entry1+2]          # J2
#            Elem.StateVar[i,entry1+3] = Elem.StateVarN[i,entry1+3]          # ep3
#            # 
#            Elem.StateVar[i,iS+0] = Elem.StateVarN[i,iS+0]                  # strain velocities - used in ViscExten3D
#            Elem.StateVar[i,iS+1] = Elem.StateVarN[i,iS+1]
#            Elem.StateVar[i,iS+2] = Elem.StateVarN[i,iS+2]
#            Elem.StateVar[i,iS+3] = Elem.StateVarN[i,iS+3]
#            Elem.StateVar[i,iS+4] = Elem.StateVarN[i,iS+4]
#            Elem.StateVar[i,iS+5] = Elem.StateVarN[i,iS+5]
        
class MicroPlaneDam(Material):                                                 # microplane damage actual version
    def __init__(self, PropMat, RotCrack,PrinStrain,ShearRetFac, S2ndCrack,S2ndCrackA, f6):
#                        (self, SymmetricVal, RTypeVal,   UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, False,        PropMat[8], False,     False,     None,        8,       'MICRODAMAGE')
#        self.Symmetric = False                                              # flag for symmetry of material matrices
#        self.RType = PropMat[4]                                             # type of regularization 0: without, 1: gradient 2: crack band
#        self.Update = False                                                 # has specific update procedure for update of state variables
#        self.Updat2 = False                                                 # no 2 stage update procedure
#        self.StateVar = 1                                                   # see below
#        self.NDataVal = 8
#        self.Type = "MicroPl"
        self.PropMat= PropMat                                                # for wrappers, e.g. RCSHELL
        self.EE = PropMat[0]                                                # macroscopic Young's modulus
#        self.Emod = PropMat[0]                                              # used for SDA
        self.nu = PropMat[1]                                                # macroscopic Poissons's ratio
#        self.type = PropMat[2]                                              # type of damage function
        #
        self.nInt = len(I21Points)
        self.nState = 1                                                     # number of state variables per integration direction of unit sphere
        self.iS       = self.nState*self.nInt+5                             # entry index for strain rate for, e.g. viscous extension
        self.StateVar = self.nState*self.nInt+5+6+7    # -> 39              # number of state variables per integration point of element; 1st extra for dd_iso, 2nd for I1, 3rd for J2, 
                                                                            # 4th currently not used, 5th for eps_33 in case of plane state, 6-11 for strain increment
                                                                            # 12 largest principal value, 13-15 largest principal strain direction, 16-18 traction in largest principal strain direction 
        KK = self.EE/(3.*(1.-2.*self.nu))                                   # macroscopic bulk modulus
        GG = self.EE/(2.*(1.+self.nu))                                      # macroscopic shear modulus
        # V-D split
        self.E_V= 3.*KK                                                     # moduli of microplane elasticity V-D split
        self.E_D= 2.*GG
        self.E_DM = array([[self.E_D,0,0],[0,self.E_D,0],[0,0,self.E_D]])
        self.PlStressL = 0.01                                               # limit for plane stress iteration (dim = 2 (plane stress ), 21)
        self.PlStressI = 10                                                 # max number of iterations
        self.kV =         PropMat[4]                                        # (measure for) ratio of uniaxial compressive strength to tensile strength
        self.kV0 = (self.kV-1.)/(2.*self.kV*(1.-2.*self.nu))
        self.kV1 = self.kV0
        self.kV2 = 3./(self.kV*(1.+self.nu)**2)
        self.fct= PropMat[3]                                                # tensile strength (PropMat[5,6] not used
        self.eps_ct = 1.648721271*self.fct/self.EE                          # strain for uniaxial tensile strength to make e_0=0
        self.gd = 2.
        self.e0 = 0.                                                        # with eps_ct as above
#            self.e0 = 1.                                        # to make it elastic for control purposes
        self.ed = 2.331643981*PropMat[3]/self.EE                            # with eps_ct as above
        self.gam2 = 12.*350.                                                # parameter for scaling of equivalent damage strain regularization
        self.SpecCrEn = self.SpeCrEnergySample( 1.0, self.eps_ct, False,None,None, self.UniaxTensionUnit) # specific crack energy unscaled
        self.fc = 0.7*self.fct*PropMat[4]                                   # rough estimation of uniaxial compressive strength - required for pressure dependent bond
        self.RType  = PropMat[7]                                            # indicator for regularization
        self.RegPar = PropMat[8]                                            # crack energy for regularization (RType=2)
        self.eta    = PropMat[9]                                            # for viscous regularization
        if self.RType==2: self.bw = self.RegPar/self.SpecCrEn
        else:             self.bw = 0.
        self.ElCharLenBound = 0.282 # PropMat[11]                                   # related bound to distinguish between types of crack band regularization
        self.Density =    PropMat[12]
        self.lamLow  = 0.003 # 0.01
        self.lamHigh = 100.
        self.lam2Low = 0.01
        self.lam2High= 0.5
        self.CrX, self.CrY = self.SpeCrEnergyData( "LargeEl", self.lamLow,self.lamHigh,  self.CrBwN, self.eps_ct, self.UniaxTensionScaled) # arrays to determine scaling factor for given characteristic element length
        self.CrX2,self.CrY2= self.SpeCrEnergyData( "SmallEl", self.lam2Low,self.lam2High,self.CrBwN, self.eps_ct, self.UniaxTensionScaled2) # arrays to determine scaling factor for given characteristic element length
        if self.CrX[-1]<self.CrX2[-1]:
            print('scaling factors arguments type 1,2 do not overlap ',self.CrX[0],self.CrX[-1],'__',self.CrX2[-1],self.CrX2[0])
        self.dDestr   = 1.-1.e-4                                            # for maximum damage allowed, subtractor should not be small in order to avoid problems in dividing by MMatS[2,2], see below
        self.kapUlt   = exp(log(-log(1.-self.dDestr))/self.gd)*self.ed+self.e0
    def UniaxTensionScaled(self, gam1, eps):
        xx = eps - self.eps_ct
        kap_ =gam1*xx + (1.-gam1)/(self.gam2) * (1.-exp(-self.gam2*xx)) + self.eps_ct
        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps
    def UniaxTensionScaled2(self, beta, eps):
        kap = eps 
        kap_=  (1-beta)*self.eps_ct * log( ( kap-beta*self.eps_ct )/( (1-beta)*self.eps_ct) ) + self.eps_ct
        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps
    def UniaxTensionUnit(self, dummy, eps):
        kap_ = eps
        return self.EE * exp(-pow(( kap_-self.e0)/self.ed ,self.gd)) * eps

    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps_, Eps_, dTmp, Temp, EpsR):
#
#        if Elem.Label==57 and ipI==0: Echo(f"XXX {Eps_[0]:12.6e}, {Eps_[1]:12.6e}, {Eps_[2]:12.6e} ", ff)
#        
        def VolStiffness():
            return array([[ob9, ob9, ob9, 0, 0, 0], 
                          [ob9, ob9, ob9, 0, 0, 0],
                          [ob9, ob9, ob9, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0]])
        def DevStiffness():
            return array([[ fb9*nn[0]**2+ob9*nn[1]**2+ob9*nn[2]**2, -tb9*nn[0]**2-tb9*nn[1]**2+ob9*nn[2]**2, -tb9*nn[0]**2+ob9*nn[1]**2-tb9*nn[2]**2, -ob3*nn[2]*nn[1],  ob6*nn[0]*nn[2],  ob6*nn[0]*nn[1]],
                          [-tb9*nn[0]**2-tb9*nn[1]**2+ob9*nn[2]**2,  ob9*nn[0]**2+fb9*nn[1]**2+ob9*nn[2]**2,  ob9*nn[0]**2-tb9*nn[1]**2-tb9*nn[2]**2,  ob6*nn[2]*nn[1], -ob3*nn[0]*nn[2],  ob6*nn[0]*nn[1]],
                          [-tb9*nn[0]**2+ob9*nn[1]**2-tb9*nn[2]**2,  ob9*nn[0]**2-tb9*nn[1]**2-tb9*nn[2]**2,  ob9*nn[0]**2+ob9*nn[1]**2+fb9*nn[2]**2,  ob6*nn[2]*nn[1],  ob6*nn[0]*nn[2], -ob3*nn[0]*nn[1]],
                          [-ob3*nn[2]*nn[1],  ob6*nn[2]*nn[1],  ob6*nn[2]*nn[1], ob4*nn[2]**2+ob4*nn[1]**2, ob4*nn[0]*nn[1],           ob4*nn[0]*nn[2]],
                          [ ob6*nn[0]*nn[2], -ob3*nn[0]*nn[2],  ob6*nn[0]*nn[2], ob4*nn[0]*nn[1],           ob4*nn[2]**2+ob4*nn[0]**2, ob4*nn[2]*nn[1]],
                          [ ob6*nn[0]*nn[1],  ob6*nn[0]*nn[1], -ob3*nn[0]*nn[1], ob4*nn[0]*nn[2],           ob4*nn[2]*nn[1],           ob4*nn[1]**2+ob4*nn[0]**2]])
        def CVStiffness():
            c0 =  tb3*epsDD[0]*nn[0]-ob3*epsDD[1]*nn[1]-ob3*epsDD[2]*nn[2]
            c1 = -ob3*epsDD[0]*nn[0]+tb3*epsDD[1]*nn[1]-ob3*epsDD[2]*nn[2]
            c2 = -ob3*epsDD[0]*nn[0]-ob3*epsDD[1]*nn[1]+tb3*epsDD[2]*nn[2]
            c3 =  ob2*epsDD[1]*nn[2]+ob2*epsDD[2]*nn[1]
            c4 =  ob2*epsDD[0]*nn[2]+ob2*epsDD[2]*nn[0]
            c5 =  ob2*epsDD[0]*nn[1]+ob2*epsDD[1]*nn[0]
            return array([[c0**2, c0*c1, c0*c2, c0*c3, c0*c4, c0*c5],
                          [c0*c1, c1**2, c1*c2, c1*c3, c1*c4, c1*c5],
                          [c0*c2, c1*c2, c2**2, c2*c3, c2*c4, c2*c5],
                          [c0*c3, c1*c3, c2*c3, c3**2, c3*c4, c3*c5],
                          [c0*c4, c1*c4, c2*c4, c3*c4, c4**2, c4*c5],
                          [c0*c5, c1*c5, c2*c5, c3*c5, c4*c5, c5**2]]),array ([[ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5], 
                                                                               [ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5], 
                                                                               [ob3*c0, ob3*c1, ob3*c2, ob3*c3, ob3*c4, ob3*c5], 
                                                                               [0, 0, 0, 0, 0, 0], 
                                                                               [0, 0, 0, 0, 0, 0], 
                                                                               [0, 0, 0, 0, 0, 0]])
        if CalcType == 0: return [], [], []
        if not ConFemMatCFlag: # False: # flag for not C-version
#        if True:
            ob3, ob6, ob9, tb3, ob2 = 1./3., 1./6., 1./9., 2./3., 0.5
            fb9, tb9, ob4 = 4./9., 2./9., 1./4.
            VVV = VolStiffness()
            ep2  = Elem.StateVarN[ipI,self.nState*self.nInt+4]              # value of last iteration taken, not from last converged step
            if Elem.dim==2:
                if Elem.PlSt:  Eps = array([ Eps_[0], Eps_[1], ep2, 0., 0., Eps_[2] ])
                else:          Eps = array([ Eps_[0], Eps_[1], 0., 0., 0., Eps_[2] ])
                Dps = array([ Dps_[0], Dps_[1], 0., 0., 0., Dps_[2]])                        # plane stress --> strain Voigt notation --  used for viscous regularization only with diagonal stiffness
            elif Elem.dim == 4:
                Eps = array([Eps_[0], Eps_[1], Eps_[2], 0., 0., Eps_[3]])
                Dps = array([Dps_[0], Dps_[1], Dps_[2], 0., 0., Dps_[3]])
            elif Elem.dim==3:
                Eps = array([Eps_[0], Eps_[1], Eps_[2], Eps_[3], Eps_[4], Eps_[5]])
                Dps = array([ Dps_[0], Dps_[1], Dps_[2], Dps_[3], Dps_[4], Dps_[5]]) # triaxial strain increment Voigt notation
            elif Elem.dim==21: 
                Eps = array([Eps_[0], Eps_[1], ep2,     Eps_[3], Eps_[4], Eps_[5]])  # continuum based shell
                Dps = array([ Dps_[0], Dps_[1], 0.,     Dps_[3], Dps_[4], Dps_[5]])  # triaxial strain increment for continuum bases shell Voigt notation
            else: raise NameError("ConFemMaterials::Microplane.Sig: not implemented")
            #
#            epsVol = ob3*( Eps[0]+Eps[1]+Eps[2] )
            ns = self.nState
            PlStFlag = False
            if (Elem.dim==2 and Elem.PlSt) or Elem.dim==21: PlStFlag = True
            #
            for ii in range(self.PlStressI):                                # for plane stress iteration
                epsVol = ob3*( Eps[0]+Eps[1]+Eps[2] )
                dd_iso, I1, J2 = 0., 0., 0.
                MMatS, MMatT = zeros((6,6),dtype=float), zeros((6,6),dtype=float)
                sig = zeros((6),dtype=float)                                # stresses in voigt notation
                for i in range(self.nInt):                                  # microplane integration order
                    kapOld = Elem.StateVar[ipI,ns*i]
                    nn = array([I21Points[i,0],I21Points[i,1],I21Points[i,2]])
                    # V-D-Split projection tensors
                    D0 =  array([ tb3*nn[0], -ob3*nn[0], -ob3*nn[0], 0,         ob2*nn[2], ob2*nn[1]])
                    D1 =  array([-ob3*nn[1],  tb3*nn[1], -ob3*nn[1], ob2*nn[2], 0,         ob2*nn[0]])   
                    D2 =  array([-ob3*nn[2], -ob3*nn[2],  tb3*nn[2], ob2*nn[1], ob2*nn[0], 0])
                    # vector deviator strains
                    epsDD = array([ dot(D0,Eps) , dot(D1,Eps), dot(D2,Eps)])
                    # microplane strain invariants, state variable
                    I1mp = 3.*epsVol
                    J2mp = 3./2.*dot(epsDD,epsDD)
                    eta_ = self.kV0*I1mp + sqrt( (self.kV1*I1mp)**2 + self.kV2*J2mp )
                    # equivalent strain
                    if self.RType==2:                                       # crack band regularization
                        if eta_>self.eps_ct:                                # limit strain exceeded
                            beta= Elem.CrBwS                                # --> ConFemMat::Material -- find scaling factors for given char length for both types
                            if Elem.ScaleType==1:
                                xx  = eta_ - self.eps_ct
                                eta = beta*xx + (1.-beta)/(self.gam2) * (1.-exp(-self.gam2*xx)) + self.eps_ct
                                dkk = beta +    (1.-beta)             *     exp(-self.gam2*xx)
                            else:
                                eta = (1-beta)*self.eps_ct * log( ( eta_-beta*self.eps_ct )/( (1-beta)*self.eps_ct) ) + self.eps_ct
                                dkk = (1-beta)*self.eps_ct/(eta_-beta*self.eps_ct)
                        else:                                               # below limit strain
                            eta = eta_
                            dkk = 1.
                    else:                                                   # no regularization
                        eta = eta_
                        dkk = 1.
                    if eta>self.kapUlt: eta = self.kapUlt                   # to have some residual stiffness
                    # damage function
                    kap = max( self.e0, kapOld, eta)
                    if kap<=self.e0: 
                        dd = 0.
                        Dp = 0.
                    else:            
                        dd = 1.-exp(-((kap-self.e0)/self.ed)**self.gd)
                        Dp = (1.-dd) * self.gd/self.ed**self.gd * (kap-self.e0)**(self.gd-1.) 
#                        Dp = (1.-dd) * self.gd/pow(self.ed,self.gd) * pow((kap-self.e0),(self.gd-1.))
                        Dp = Dp*dkk                                         # regularization scaling of derivative dD/dkap tangential stiffness
                    # stresses
                    sigVV = (1.-dd)*    self.E_V*epsVol                     # volumetric stress (scalar)
                    sigDD = (1.-dd)*dot(self.E_DM,epsDD)                    # deviatoric stress (vector)
                    sigV  = array([ ob3*sigVV, ob3*sigVV, ob3*sigVV, 0., 0., 0.])
                    sig   = sig + 6.*I21Weights[i]*( sigV + sigDD[0]*D0 + sigDD[1]*D1 + sigDD[2]*D2) # stress tensor integration # factor 6. calibrates the weighting factors
                    # secant stiffness
                    DDD = DevStiffness()
                    MMatS = MMatS + 6.*I21Weights[i]*( (1.-dd)*self.E_V*VVV + (1.-dd)*self.E_D*DDD ) # presumably symmetric, i.e. only upper right part builded
                    # corrector for tangential material stiffness
                    if kap>kapOld and dd>ZeroD:
                        Q = sqrt( (3.*self.kV1*epsVol)**2 + 1.5*self.kV2*(epsDD[0]*epsDD[0]+epsDD[1]*epsDD[1]+epsDD[2]*epsDD[2]))
                        Q = 1./Q
                        A = Dp * ( 3.*self.kV0 + Q*(3.*self.kV1)**2*epsVol )
                        B = Dp * ( 1.5*Q*self.kV2 )
                        sigV0 = self.E_V*epsVol
                        CC, VC = CVStiffness()
                        MMatT = MMatT - 6.*I21Weights[i] * ( sigV0*(A*VVV + B*VC) + self.E_D*(A*transpose(VC) + B*CC) )
                    #
                    dd_iso += 2.*I21Weights[i] * dd
                    I1 += 2.*I21Weights[i] * I1mp
                    J2 += 2.*I21Weights[i] * J2mp
                    # microplane state variable updates
                    Elem.StateVarN[ipI,ns*i] = kap
                    
                ep2 = -( MMatS[2,0]*Eps[0] + MMatS[2,1]*Eps[1] + MMatS[2,5]*Eps[5] ) / MMatS[2,2] # lateral normal strain in case of zero lateral stress
                if abs(Eps[2])>ZeroD*1.e-3: xyz = abs((Eps[2]-ep2)/Eps[2])
                else:                       xyz = 0.
                if PlStFlag: Eps[2] = ep2
                if not PlStFlag or xyz<self.PlStressL: break

            if ii >= self.PlStressI-1: raise NameError("ConFemMaterials::Microplane.Sig: plane stress iteration", ii, sig[2])
            #
            Elem.StateVarN[ipI,self.nState*self.nInt] = dd_iso
            Elem.StateVarN[ipI,self.nState*self.nInt+1] = I1
            Elem.StateVarN[ipI,self.nState*self.nInt+2] = J2
            Elem.StateVarN[ipI,self.nState*self.nInt+4] = ep2
            MatM = MMatS + MMatT
            #
            if self.eta>0.: 
                zz, Veps = self.ViscExten3D( Dt, self.eta, Dps, Elem, ipI, self.iS)           # for viscous regularization
                for i in range(6): 
                    sig[i] = sig[i] + self.eta*Veps[i]
                    MatM[i,i] = MatM[i,i] + zz 
            #
            if Elem.dim==2:                                                                     # 2D
                if Elem.PlSt:                                                                   # plane stress
                    cD = 1./MatM[2,2]
                    MatM_ =array([[MatM[0,0]-MatM[0,2]*MatM[2,0]*cD,MatM[0,1]-MatM[0,2]*MatM[2,1]*cD,MatM[0,5]-MatM[0,2]*MatM[2,5]*cD],
                                  [MatM[1,0]-MatM[1,2]*MatM[2,0]*cD,MatM[1,1]-MatM[1,2]*MatM[2,1]*cD,MatM[1,5]-MatM[1,2]*MatM[2,5]*cD], 
                                  [MatM[5,0]-MatM[5,2]*MatM[2,0]*cD,MatM[5,1]-MatM[5,2]*MatM[2,1]*cD,MatM[5,5]-MatM[5,2]*MatM[2,5]*cD]])
                else:                                                                           # plane strain
                    MatM_ =array([[MatM[0,0],MatM[0,1],MatM[0,5]],[MatM[1,0],MatM[1,1],MatM[1,5]],[MatM[5,0],MatM[5,1],MatM[5,5]]])
                return array([sig[0],sig[1],sig[5]]), MatM_, [Eps[0], Eps[1], ep2, Eps[5], sig[0], sig[1], sig[2], sig[5],dd_iso]
            elif Elem.dim == 4:
                MatM_ = array([[MatM[0, 0], MatM[0, 1], MatM[0, 2], MatM[0, 5]],
                               [MatM[1, 0], MatM[1, 1], MatM[1, 2], MatM[1, 5]],
                               [MatM[2, 0], MatM[2, 1], MatM[2, 2], MatM[2, 5]],
                               [MatM[5, 0], MatM[5, 1], MatM[5, 2], MatM[5, 5]]])
                return [sig[0], sig[1], sig[2], sig[5]], MatM_, [Eps[0], Eps[1], Eps[2], Eps[5], sig[0], sig[1], sig[2], sig[5], dd_iso]
            elif Elem.dim==3:                                                                   # 3D
                return array([sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]), MatM, [Eps[0],Eps[1],Eps[2],Eps[3],Eps[4],Eps[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],dd_iso]
            elif Elem.dim==21:                                                                  # Continuum based shell
                if abs(MatM[2,2])<ZeroD: 
                    print('X\n', MMatS, '\n', MatM) 
                    raise NameError("ConFemMaterials::Microplane.Sig: tangential Mat[2,2] to small")
                cD = 1./MatM[2,2]
                MatM_ = array([[MatM[0,0]-MatM[0,2]*MatM[2,0]*cD, MatM[0,1]-MatM[0,2]*MatM[2,1]*cD, 0., MatM[0,3]-MatM[0,2]*MatM[2,3]*cD, MatM[0,4]-MatM[0,2]*MatM[2,4]*cD, MatM[0,5]-MatM[0,2]*MatM[2,5]*cD],
                               [MatM[1,0]-MatM[1,2]*MatM[2,0]*cD, MatM[1,1]-MatM[1,2]*MatM[2,1]*cD, 0., MatM[1,3]-MatM[1,2]*MatM[2,3]*cD, MatM[1,4]-MatM[1,2]*MatM[2,4]*cD, MatM[1,5]-MatM[1,2]*MatM[2,5]*cD],
                               [0.,                               0.,                               0., 0.,                               0.,                               0.],
                               [MatM[3,0]-MatM[3,2]*MatM[2,0]*cD, MatM[3,1]-MatM[3,2]*MatM[2,1]*cD, 0., MatM[3,3]-MatM[3,2]*MatM[2,3]*cD, MatM[3,4]-MatM[3,2]*MatM[2,4]*cD, MatM[3,5]-MatM[3,2]*MatM[2,5]*cD],
                               [MatM[4,0]-MatM[4,2]*MatM[2,0]*cD, MatM[4,1]-MatM[4,2]*MatM[2,1]*cD, 0., MatM[4,3]-MatM[4,2]*MatM[2,3]*cD, MatM[4,4]-MatM[4,2]*MatM[2,4]*cD, MatM[4,5]-MatM[4,2]*MatM[2,5]*cD],
                               [MatM[5,0]-MatM[5,2]*MatM[2,0]*cD, MatM[5,1]-MatM[5,2]*MatM[2,1]*cD, 0., MatM[5,3]-MatM[5,2]*MatM[2,3]*cD, MatM[5,4]-MatM[5,2]*MatM[2,4]*cD, MatM[5,5]-MatM[5,2]*MatM[2,5]*cD]])
#            
#                if Elem.Label==16 and ipI==19:
#                    print('XXX\n',sig)
#            
                if Elem.ShellRCFlag: return array([sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]), MatM_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps[0], Eps[1], Eps[5], sig[0], sig[1], sig[5], 0.,0.]
                else:                return array([sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]), MatM_, [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]
            else: raise NameError ("ConFemMaterials::Microplane.Sig: not implemented for this element type")
        # C-version
        else:
            DataOut  = zeros((1),dtype=float) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if Elem.dim==1:
                sig  = zeros((2),dtype=float)
                MatM_= zeros((4),dtype=float)
                PlSt_ = False
                raise NameError ("ConFemMaterials::Microplane.Sig: not implemented for this element type")
            elif Elem.dim==2:
                sig  = zeros((3),dtype=float)
                MatM_= zeros((9),dtype=float)
                PlSt_ = Elem.PlSt
            elif Elem.dim==3:
                sig  = zeros((6),dtype=float)  
                MatM_= zeros((36),dtype=float)
                PlSt_ = False
            elif Elem.dim == 4:  # axisymmetric
                sig  = zeros((4), dtype=float)
                MatM_= zeros((16), dtype=float)
#                MatM_ = zeros((4, 4), dtype=float)
                PlSt_ = False
            elif Elem.dim==21:
                sig  = zeros((6),dtype=float)  
                MatM_= zeros((36),dtype=float)
                PlSt_ = False
            Dps = zeros((6),dtype=float)
            rc = MicroPlaneC2( Elem.dim, PlSt_, False, Elem.Lch_, Elem.StateVar[ipI], Elem.StateVarN[ipI],\
                           Eps_, sig, MatM_, 0, self.E_V, self.E_D,\
                           self.kV0, self.kV1, self.kV2, 0., 0., 0., self.RType, Elem.CrBwS,self.gam2,self.kapUlt,Elem.ScaleType,\
                           self.eps_ct, self.e0, self.ed, self.gd, self.nu, self.EE, Dps_, self.eta,\
                           Dt, 0., DataOut, self.nState, self.nInt, self.PlStressI, self.PlStressL, self.iS)
            # not used: Elem.Lch_ but Elem.CrBwS for regularization
            if rc> 110:  raise NameError("ConFemMaterials::MicroPlane:sig: RC "+str(rc))
            if rc==110:  print(f"ConFemMaterials::MicroPlane:sig: no plain stress convergence El {Elem.Label:d}, IP {ipI:d}", file=ff)
            D_   = Elem.StateVarN[ipI,self.nState*self.nInt]
            if Elem.dim==1:
                pass
            elif Elem.dim==2:
                return ([sig[0],sig[1],sig[2]]),\
                array([[MatM_[0],MatM_[1],MatM_[2]],
                       [MatM_[3],MatM_[4],MatM_[5]],
                       [MatM_[6],MatM_[7],MatM_[8]]]),\
                       [Eps_[0],Eps_[1],0., Eps_[2],sig[0],sig[1],0., sig[2],D_]
            elif Elem.dim==3:
                return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
                array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
                       [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
                       [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
                       [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
                       [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
                       [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
                       [Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5],sig[0],sig[1],sig[2],sig[3],sig[4],sig[5],D_]
            elif Elem.dim == 4:
                return [sig[0],sig[1],sig[2],sig[3]], \
                       array([[MatM_[0],  MatM_[1],  MatM_[2],  MatM_[3]],
                              [MatM_[4],  MatM_[5],  MatM_[6],  MatM_[7]],
                              [MatM_[8],  MatM_[9],  MatM_[10], MatM_[11]],
                              [MatM_[12], MatM_[13], MatM_[14], MatM_[15]]]), \
                       [ Eps_[0],Eps_[1],Eps_[2],Eps_[3], sig[0],sig[1],sig[2],sig[3], D_]
            elif Elem.dim==21:
                if Elem.ShellRCFlag:
                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps_[0],Eps_[1],Eps_[5], sig[0],sig[1],sig[5], D_]
                else:
                    return [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]],\
                    array([[MatM_[0], MatM_[1], MatM_[2], MatM_[3], MatM_[4], MatM_[5] ],
                           [MatM_[6], MatM_[7], MatM_[8], MatM_[9], MatM_[10],MatM_[11]],
                           [MatM_[12],MatM_[13],MatM_[14],MatM_[15],MatM_[16],MatM_[17]],
                           [MatM_[18],MatM_[19],MatM_[20],MatM_[21],MatM_[22],MatM_[23]],
                           [MatM_[24],MatM_[25],MatM_[26],MatM_[27],MatM_[28],MatM_[29]],
                           [MatM_[30],MatM_[31],MatM_[32],MatM_[33],MatM_[34],MatM_[35]]]),\
                           [sig[0],sig[1],sig[2],sig[3],sig[4],sig[5]]

    def UpdateStateVar(self, Elem, ff):
        ns = self.nState
        iS = self.iS                                                        # entry point for strain rate
        for i in range(Elem.StateVar.shape[0]):
            for j in range(self.nInt):                                      # loop over microplanes
                kapOld = Elem.StateVar[i,ns*j]
                kapNew = Elem.StateVarN[i,ns*j]
                if kapNew>kapOld: Elem.StateVar[i,ns*j] = Elem.StateVarN[i,ns*j]
            entry1 = ns*self.nInt
            ddIsoOld = Elem.StateVar[i,entry1]
            ddIsoNew = Elem.StateVarN[i,entry1]
            if ddIsoNew>ddIsoOld: Elem.StateVar[i,entry1] = Elem.StateVarN[i,entry1] 
            Elem.StateVar[i,entry1+1] = Elem.StateVarN[i,entry1+1]          # I1
            Elem.StateVar[i,entry1+2] = Elem.StateVarN[i,entry1+2]          # J2
            Elem.StateVar[i,entry1+3] = Elem.StateVarN[i,entry1+3]          # ep3
            # 
            Elem.StateVar[i,iS+0] = Elem.StateVarN[i,iS+0]                  # strain velocities - used in ViscExten3D
            Elem.StateVar[i,iS+1] = Elem.StateVarN[i,iS+1]
            Elem.StateVar[i,iS+2] = Elem.StateVarN[i,iS+2]
            Elem.StateVar[i,iS+3] = Elem.StateVarN[i,iS+3]
            Elem.StateVar[i,iS+4] = Elem.StateVarN[i,iS+4]
            Elem.StateVar[i,iS+5] = Elem.StateVarN[i,iS+5]

class Mises(Material):                              # elastoplastic Mises
    def __init__(self, PropMat, val):
#    def __init__(self,         SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,      False,     5,           9, 'MISES')
#        self.Symmetric = True                       # flag for symmetry of material matrices
#        self.Update = True                          # has specific update procedure for update of state variables
#        self.Updat2 = False                         # no 2 stage update procedure
#        self.StateVar = 5                           # number of state variables per integration point (may in the end be overruled by other types using mises, see RCBeam )
#                                                    # 0: current permanent strain upon unloading, 1: current yield stress, 2: current reference strain for smoothed uniaxial stress strain curve, 3: final stress of last time step used for smoothed version
#        self.NData = 8                              # number of data items

        self.Emod = PropMat[0]                      # Young's modulus
        self.nu = PropMat[1]                        # Poissons's ratio
        self.sigY = PropMat[2]                      # uniaxial yield stress
        self.sigU = PropMat[3]                      # strength
        self.epsU = PropMat[4]                      # limit strain
        self.alphaT = PropMat[5]                    # thermal expansion coefficient
        self.sfac = PropMat[6]                      # parameter to smooth the transition in the bilinear course
        if self.sfac>0.:
            denom = self.sfac*(-self.epsU*self.Emod+self.sigY)
            self.b0 = 0.25*(self.sfac**2*self.epsU*self.Emod-2*self.sfac*self.epsU*self.Emod+self.epsU*self.Emod-self.sigU+2*self.sfac*self.sigU-self.sfac**2*self.sigU)*self.sigY/denom
            self.b1 = 0.50*self.Emod*(-self.epsU*self.Emod+self.sigU-self.sfac*self.epsU*self.Emod-self.sfac*self.sigU+2*self.sigY*self.sfac)/denom
            self.b2 = 0.25*self.Emod**2*(self.epsU*self.Emod-self.sigU)/self.sigY/denom
        self.Density = PropMat[7]                   # specific mass
        self.fct = val[0]                           # concrete tensile strength for tension stiffening
        self.alpha = val[1]                         # tension stiffening parameter
        self.epsY = self.sigY/self.Emod             # uniaxial yield strain  ????
        self.Etan = (self.sigU-self.sigY)/(self.epsU-self.epsY)# tangential / hardening modulus
        self.H = self.Etan/(1-self.Etan/self.Emod)  # hardening modulus
        
        self.epsD  = self.sfac*self.epsY
        
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        sigy = max(self.sigY,Elem.StateVar[ipI][1]) # current state parameter - current uniaxial yield stress
        Elem.StateVarN[ipI][0] = 0                  # current values of state variables have to initialized again
        Elem.StateVarN[ipI][1] = 0
        if Elem.dim==2 or Elem.dim==3 or Elem.dim==21:           # plane stress/strain biaxial
            nu = self.nu                            # Poisson's ratio
            mu = self.Emod/(2.*(1.+nu))             # shear modulus
            C0 = self.Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),nu/(1-nu),0,0,0],
                                                           [nu/(1-nu),1,nu/(1-nu),0,0,0],
                                                           [nu/(1-nu),nu/(1-nu),1,0,0,0],
                                                           [0,0,0,(1-2*nu)/(2*(1-nu)),0,0],
                                                           [0,0,0,0,(1-2*nu)/(2*(1-nu)),0],
                                                           [0,0,0,0,0,(1-2*nu)/(2*(1-nu))]]) # triaxial isotropic elasticity 
            if Elem.dim==2:                         # plate plane stress / plane strain
                Sig = array( [Elem.DataP[ipI,4],Elem.DataP[ipI,5],Elem.DataP[ipI,6],0.,0.,Elem.DataP[ipI,7]] ) # stress of previous increment
                dEps = array([Dps[0],Dps[1],0.,0.,0.,Dps[2]]) # total strain increment
            elif Elem.dim==21:                      # cb shell -- dEps[2] --> 0, re-evaluated in the following
                Sig = array([Elem.DataP[ipI,0],Elem.DataP[ipI,1],Elem.DataP[ipI,2],Elem.DataP[ipI,3],Elem.DataP[ipI,4],Elem.DataP[ipI,5]] ) # stress of previous increment
#                Sig = array([Elem.DataP[ipI,6],Elem.DataP[ipI,7],Elem.DataP[ipI,8],Elem.DataP[ipI,9],Elem.DataP[ipI,10],Elem.DataP[ipI,11]] ) # stress of previous increment
                dEps = array([Dps[0],Dps[1],0.,Dps[3],Dps[4],Dps[5]]) # total strain increment
            elif Elem.dim==3:
#               Sig = array([Elem.DataP[ipI,0],Elem.DataP[ipI,1],Elem.DataP[ipI,2],Elem.DataP[ipI,3],Elem.DataP[ipI,4],Elem.DataP[ipI,5]] ) # stress of previous increment
                Sig = array([Elem.DataP[ipI,6],Elem.DataP[ipI,7],Elem.DataP[ipI,8],Elem.DataP[ipI,9],Elem.DataP[ipI,10],Elem.DataP[ipI,11]] ) # stress of previous increment
                dEps = array([Dps[0],Dps[1],Dps[2],Dps[3],Dps[4],Dps[5]]) # total strain increment
            Fn = sqrt(3.*(((Sig[0]-Sig[1])**2+(Sig[0]-Sig[2])**2+(Sig[1]-Sig[2])**2)/6.+Sig[3]**2+Sig[4]**2+Sig[5]**2))-sigy # distance to yield surface of previous step
            Eflag = False
            dEpp = zeros((6),dtype=float)           # initial value plastic strain increment
            dSig = zeros((6),dtype=float)           # initial value plastic strain increment
            dLam = 0.                               # initial value plastic multiplier increment
#            dsiy = 0.                               # initial value current yield stress increment
            ni = 20                                 # iteration limit
            for i in range(ni):
                if Elem.PlSt: dEps[2]=-( C0[2,0]*(dEps[0]-dEpp[0])
                                        +C0[2,1]*(dEps[1]-dEpp[1])
                                        +C0[2,3]*(dEps[3]-dEpp[3])
                                        +C0[2,4]*(dEps[4]-dEpp[4])
                                        +C0[2,5]*(dEps[5]-dEpp[5]))/C0[2,2]+dEpp[2] # lateral strain for plane stress
                dSig = dot(C0,dEps-dEpp)
                SigN = Sig + dSig
                J2 = ((SigN[0]-SigN[1])**2+(SigN[0]-SigN[2])**2+(SigN[1]-SigN[2])**2)/6.+SigN[3]**2+SigN[4]**2+SigN[5]**2 # 2nd stress deviator invariant
                if sqrt(3.*J2)-sigy<1.e-9:          # elastic loading, unloading or reloading
                    Eflag = True
                    break
                sm = (SigN[0]+SigN[1]+SigN[2])/3.   # 1st stress invariant / mean stress predictor stress
                rr = sqrt(3./(4.*J2))*array([SigN[0]-sm,SigN[1]-sm,SigN[2]-sm,SigN[3],SigN[4],SigN[5]]) # yield gradient predictor stress
                dL = (Fn + dot(rr,dSig)+rr[3]*dSig[3]+rr[4]*dSig[4]+rr[5]*dSig[5] - self.H*dLam)/(3.*mu+self.H) # with Voigt notation correction
                if dL<1.e-9: break
                dLam = dLam + dL                    # update plastic multiplier
                dEpp = dLam*array([rr[0],rr[1],rr[2],2.*rr[3],2.*rr[4],2.*rr[5]]) # plastic strain incremen with Voigt notation correction
            if i>=ni-1:  
                print(elI, ipI, i, ni, Sig, dEps, sigy, dLam, dEpp, dSig)  
                raise NameError ("ConFemMaterials::Mises.Sig: no convergence")
            Sig = SigN
            if Eflag:                               # elastic loading or unloading / reloading
                if Elem.dim==2:                     # plane stress / strain
                    if Elem.PlSt: MatM = self.Emod/(1-nu**2)*array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]]) # plane stress
                    else:         MatM = self.Emod*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),0],[nu/(1-nu),1,0],[0,0,(1-2*nu)/(2*(1-nu))]]) # plane strain
                elif Elem.dim==21:                  # continuum bases shell -- plane stress ???
                    MatM = self.Emod*array([[1./(1-nu**2), nu/(1-nu**2), 0., 0., 0., 0.],
                                            [nu/(1-nu**2), 1./(1-nu**2), 0., 0., 0., 0.],
                                            [0., 0., 0., 0., 0., 0.],
                                            [0., 0., 0., 1./(2+2*nu), 0., 0.],
                                            [0., 0., 0., 0., 1./(2+2*nu), 0],
                                            [0., 0., 0., 0., 0., 1./(2+2*nu)]]) 
                else: MatM = C0
            else:                                   # plastic loading
                Elem.StateVarN[ipI][1]= sqrt(3.*J2) # new equivalent yield limit 
                A = 1./(self.H + 3.*mu)             # --> 2.*mu*dot(rr,rr) --> dot(xx,rr) --> dot(C0,rr)
                CC = C0 -  A*4.*mu**2*outer(rr,rr)  # --> A*outer(xx,xx)            # tangential material stiffness
                cD = 1./CC[2,2]
                if Elem.dim==2:
                    if Elem.PlSt: MatM=array([[CC[0,0]-CC[0,2]*CC[2,0]*cD,CC[0,1]-CC[0,2]*CC[2,1]*cD,CC[0,5]-CC[0,2]*CC[2,5]*cD],
                                              [CC[1,0]-CC[1,2]*CC[2,0]*cD,CC[1,1]-CC[1,2]*CC[2,1]*cD,CC[1,5]-CC[1,2]*CC[2,5]*cD], 
                                              [CC[5,0]-CC[5,2]*CC[2,0]*cD,CC[5,1]-CC[5,2]*CC[2,1]*cD,CC[5,5]-CC[5,2]*CC[2,5]*cD]])
                    else:         MatM=array([[CC[0,0],CC[0,1],CC[0,5]],[CC[1,0],CC[1,1],CC[1,5]],[CC[5,0],CC[5,1],CC[5,5]]])
                elif Elem.dim==21:                               # shell -- plane stress
                    MatM=array([[CC[0,0]-CC[0,2]*CC[2,0]*cD,CC[0,1]-CC[0,2]*CC[2,1]*cD,0.,CC[0,3]-CC[0,2]*CC[2,3]*cD,CC[0,4]-CC[0,2]*CC[2,4]*cD,CC[0,5]-CC[0,2]*CC[2,5]*cD],
                                [CC[1,0]-CC[1,2]*CC[2,0]*cD,CC[1,1]-CC[1,2]*CC[2,1]*cD,0.,CC[1,3]-CC[1,2]*CC[2,3]*cD,CC[1,4]-CC[1,2]*CC[2,4]*cD,CC[1,5]-CC[1,2]*CC[2,5]*cD],
                                [0.,                        0.,                        0.,0.,                        0.,                        0.],
                                [CC[3,0]-CC[3,2]*CC[2,0]*cD,CC[3,1]-CC[3,2]*CC[2,1]*cD,0.,CC[3,3]-CC[3,2]*CC[2,3]*cD,CC[3,4]-CC[3,2]*CC[2,4]*cD,CC[3,5]-CC[3,2]*CC[2,5]*cD],
                                [CC[4,0]-CC[4,2]*CC[2,0]*cD,CC[4,1]-CC[4,2]*CC[2,1]*cD,0.,CC[4,3]-CC[4,2]*CC[2,3]*cD,CC[4,4]-CC[4,2]*CC[2,4]*cD,CC[4,5]-CC[4,2]*CC[2,5]*cD], 
                                [CC[5,0]-CC[5,2]*CC[2,0]*cD,CC[5,1]-CC[5,2]*CC[2,1]*cD,0.,CC[5,3]-CC[5,2]*CC[2,3]*cD,CC[5,4]-CC[5,2]*CC[2,4]*cD,CC[5,5]-CC[5,2]*CC[2,5]*cD]])
                elif Elem.dim==3:
                    MatM=array([[CC[0,0],CC[0,1],CC[0,2],CC[0,3],CC[0,4],CC[0,5]],
                                [CC[1,0],CC[1,1],CC[1,2],CC[1,3],CC[1,4],CC[1,5]],
                                [CC[2,0],CC[2,1],CC[2,2],CC[2,3],CC[2,4],CC[2,5]],
                                [CC[3,0],CC[3,1],CC[3,2],CC[3,3],CC[3,4],CC[3,5]],
                                [CC[4,0],CC[4,1],CC[4,2],CC[4,3],CC[4,4],CC[4,5]],
                                [CC[5,0],CC[5,1],CC[5,2],CC[5,3],CC[5,4],CC[5,5]]])
            if Elem.dim==2:
                sig = array([Sig[0],Sig[1],Sig[5]])
#                return sig, MatM, [Eps[0], Eps[1], Eps[2], Sig[0], Sig[1], Sig[5], Sig[2]] # data
                return sig, MatM, [Eps[0], Eps[1], 0., Eps[2], Sig[0], Sig[1], Sig[2], Sig[5],Elem.StateVarN[ipI][1]] # data
            elif Elem.dim==3:                                                           # 3D, shell
                sig = array([Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5]])
                return sig, MatM, [Eps[0],Eps[1],Eps[2],Eps[3],Eps[4],Eps[5], Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5], Elem.StateVarN[ipI][1]] # data
            elif Elem.dim==21:
                sig = array([Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5]])
                return sig, MatM, [ Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5]] # data to compute internal forces
        elif Elem.dim==1 or Elem.dim==10 or Elem.dim==11: 
            raise NameError("ConFemMat::Mises.Sig: dim mismatch") # uniaxial, Bernoulli beam, Timoshenko beam, MisesReMem
        else: raise NameError ("ConFemMat::Mises.Sig: not implemented for this element type")
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):    # loop over integration points same as for UpdateStateVar of RCBeam RList
            if Elem.StateVarN[j,1]>Elem.StateVar[j,1]:                 #Elem.StateVar[j] = Elem.StateVarN[j]
                Elem.StateVar[j,0] = Elem.StateVarN[j,0]
                Elem.StateVar[j,1] = Elem.StateVarN[j,1]
            Elem.StateVar[j,2] = Elem.StateVarN[j,2] # longitudinal strain
            if abs(Elem.StateVarN[j,4])>abs(Elem.StateVar[j,4]):
                Elem.StateVar[j,3] = Elem.StateVarN[j,3] # for smoothed uniaxial version of uniaxial mises for reinforcement
                Elem.StateVar[j,4] = Elem.StateVarN[j,4] # "
        return False

class MisesUniaxial( Mises ):                                       # elastoplastic Mises  for Elem.dim==1 or Elem.dim==10 or Elem.dim==11: # uniaxial, Bernoulli beam, Timoshenko beam, MisesReMem
    def __init__(self, PropMat, val):
        Mises.__init__(self, PropMat, val)
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        return self.SigU(ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR, -1)
    def SigU(self,       ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR, ipI2):
        #
        def sigEpsL(eps):                           # uniaxial stress strain curve smoothing between elastic and yielding branch within a range a ... b 
            Emod = self.Emod
            ET = self.Etan
            epsD = self.epsD
            epsY_= self.epsY
            sigY_= self.sigY
            gg  = 2.                                 # gg*epsD marks right side transition range
            g2  = gg**2
            g3  = gg**3
            ff = 1./(g3+3*gg+1+3*g2)
            if Elem.StateVar[ipI,3]==0: epsY, sigY = self.epsY, self.sigY  # initialization
            else:                       epsY, sigY = Elem.StateVar[ipI,3], Elem.StateVar[ipI,4] # actual values
            b0 =  ff*(sigY-2*g2*Emod*epsD+2*g2*epsD*ET+3*gg*sigY+Emod*epsY_*g3+3*g2*Emod*epsY_)
            b1 = -ff*(-4*gg*epsD*Emod-epsD*ET+ET*gg*epsD+g2*Emod*epsD-4*g2*epsD*ET+6*gg*Emod*epsY_-6*gg*sigY-g3*epsD*Emod)/epsD
            b2 = -ff*(2*g2*Emod*epsD-2*g2*epsD*ET+3*gg*Emod*epsY_-3*gg*sigY-2*gg*epsD*Emod+2*ET*gg*epsD-3*Emod*epsY_+3*sigY-2*epsD*ET+2*Emod*epsD)/epsD**2
            b3 =  ff*(2*Emod*epsY_-2*sigY+gg*epsD*Emod+epsD*ET-Emod*epsD-ET*gg*epsD)/epsD**3
            eps1 = eps+2*epsY_-epsY
            eps2 = eps        -epsY
            if   eps < epsY-2*epsY_-gg*epsD+1.e-9:
                sig = -sigY + ET*eps1
                Emo = ET
                Elem.StateVarN[ipI,3] =  eps + 2*epsY_ + gg*epsD
                Elem.StateVarN[ipI,4] = -sig - ET*gg*epsD
            elif eps < epsY-2*epsY_+epsD:
                sig =   b3*eps1**3 -  b2*eps1**2 + b1*eps1-b0
                Emo = 3*b3*eps1**2 -2*b2*eps1    + b1
            elif eps <   epsY-epsD:                                         # will presumably not work for cyclic loading with a change from compression to tension and vice versa
                sig = sigY_ + Emod*eps2 
                Emo = Emod
            elif eps < epsY+gg*epsD-1.e-9:
                sig =   b3*eps2**3 +  b2*eps2**2 + b1*eps2 + b0
                Emo = 3*b3*eps2**2 +2*b2*eps2    + b1
            else:
                sig = sigY + ET*eps2
                Emo = ET
                Elem.StateVarN[ipI,3] = eps - gg*epsD
                Elem.StateVarN[ipI,4] = sig - ET*gg*epsD
            return sig, Emo
        #
        if CalcType == 0: return [], [], []
        if ipI2==-1:                                                        # uniaxial constant across cross section
            sigy = max(self.sigY,Elem.StateVar[ipI][1])                     # current state parameter - current uniaxial yield stress
            Elem.StateVarN[ipI][0] = 0                                      # current values of state variables have to initialized again
            Elem.StateVarN[ipI][1] = 0
            epsP = Elem.StateVar[ipI][0]                                    # current permanent strain with zero stress
        else:
            sigy = max(self.sigY,Elem.StateVar[ipI][ipI2][1])                     # current state parameter - current uniaxial yield stress
            Elem.StateVarN[ipI][ipI2][0] = 0                                      # current values of state variables have to initialized again
            Elem.StateVarN[ipI][ipI2][1] = 0
            epsP = Elem.StateVar[ipI][ipI2][0]                                    # current permanent strain with zero stress
        eps = Eps[0] - self.alphaT*Temp                                     # stress inducing strain
        dps = Dps[0] - self.alphaT*dTmp                                     # stress inducing strain increment
        epc = 0                                                             # yielding correction strain in case of tension stiffening
        # 
        if ipI2==-1 and Elem.TensStiff:                                     # tension stiffening
            if Elem.Type=='SH4':                                            # Elem.dim was presumably changed temporarily
                if   Elem.nInt==2 and ipI>=16: jj = (ipI-16)//4             # floor division -> integer value -> index for reinforcement layer, 4 stands for number of integration points
                elif Elem.nInt==5 and ipI>=20: jj = (ipI-20)//5             # "
                rhoeff = Elem.Geom[jj+2,2]                                  # effective reinforcement ratio
                betat = Elem.Geom[jj+2,3]                                   # tension stiffening parameter betat
            else:
                rhoeff = Elem.Geom[elI+2,2]                                 # effective reinforcement ratio, elI presumably correct, see ConFemMaterials::RCBeam:sig - loop for reinforcement layers
                betat = Elem.Geom[elI+2,3]                                  # tension stiffening parameter betat
            if betat==0.: 
                eps_, Emod_, dsig = epsP, self.Emod, 0.
            else:
                if self.fct/rhoeff>sigy: raise NameError("ConFemMat::Mises.Sig: effective reinforcement ratio to less for minimum reinforcement")
                sigsr = self.fct/rhoeff
                eps_ = self.fct*(self.alpha-betat)/(self.Emod*rhoeff)
                Emod_ = self.alpha*sigsr/eps_
                dsig = betat*self.fct/rhoeff
                epc = dsig/self.Emod
        else:   
            eps_, Emod_, dsig = epsP, self.Emod, 0.                         # no tension stiffening
        #
        if self.sfac>0.:                                                    # elasto-plastic with smoothing of transition but WITHOUT elasto-plastic unloading / reloading
            sig_, Emod_ = sigEpsL(eps)
            MatM = array([[Emod_,0.],[0.,0.]])
            sig = array([sig_,0.])
            if abs(sig_) > sigy: Elem.StateVarN[ipI][1] = abs(sig_)         # Elem.StateVarN[ipI][1] for output only
            else:Elem.StateVarN[ipI][1] = sigy
        #
        elif eps<=(epsP-self.epsY):                                         # plastic compressive
            MatM = array([[self.Etan,0.],[0.,0.]])                          # tangential material stiffness
            sig  = array([-sigy + self.Etan*dps,0.])                        # stress
            if ipI2==-1:                                                    # uniaxial constant across cross section
                Elem.StateVarN[ipI][0]= eps+self.epsY                       # update state variables
                Elem.StateVarN[ipI][1]=-sig[0]
            else:                                                           
                Elem.StateVarN[ipI][ipI2][0]= eps+self.epsY                 # update state variables
                Elem.StateVarN[ipI][ipI2][1]=-sig[0]
        elif eps<(epsP):                                                    # elastic compressive
            MatM = array([[self.Emod,0.],[0.,0.]])
            sig  = array([self.Emod*(eps-epsP),0.])
        elif eps<(epsP+self.epsY-epc-ZeroD):                                # elastic tension with some tolerance
            if eps<eps_:                                                    # to consider tension stiffening 
                MatM = array([[Emod_,0.],[0.,0.]])
                sig  = array([Emod_*(eps-epsP),0.])
            else:
                MatM = array([[self.Emod,0.],[0.,0.]])
                sig  = array([self.Emod*(eps-epsP) + dsig,0.])
        else:                                                               # plastic tensile
            MatM = array([[self.Etan,0.],[0.,0.]])
            sig  = array([ sigy + self.Etan*dps,0.])
            if ipI2==-1:                                                    # uniaxial constant across cross section
                Elem.StateVarN[ipI][0]= eps-self.epsY                       # permanent strain
                Elem.StateVarN[ipI][1]= sig[0]                              # yield stress
            else:
                Elem.StateVarN[ipI][ipI2][0]= eps-self.epsY                 # permanent strain
                Elem.StateVarN[ipI][ipI2][1]= sig[0]                        # yield stress
#        return sig, MatM, [ eps, sig[0], epsP] 
        if ipI2 == -1: return sig, MatM, [ eps, sig[0], max(sigy,Elem.StateVar[ipI][1]) ] # ipI2 = -1 by default, see calling routine 
        else:          return sig, MatM, [ eps, sig[0], epsP]

class MisesBeam2D_(MisesUniaxial):                                       # elastoplastic Mises  for  Elem.dim==10 Bernoulli beam
    def __init__(self, PropMat, val):
        Mises.__init__(self, PropMat, val)
        self.nI = 50
#        self.epsLim = PropMat[2]/PropMat[4]                # uhc ???
    def Sig( self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        r = 0.5*Elem.Geom[1,2] # Diam/2
        z1, z2 = Elem.zLow, Elem.zUpp
        nI = self.nI
        dz = (z2-z1)/(nI-1)                                                 # cross sectional height / number of integration points
        nn = zeros((nI), dtype=double)
        mm = zeros((nI), dtype=double)
        nde= zeros((nI), dtype=double)                                      # material stiffness
        ndk= zeros((nI), dtype=double)
        mdk= zeros((nI), dtype=double)
        z = z1                                                              # cross section / strain coordinate, initial value
        for i in range(nI):                                                 # numerical integration 
            bb = 2*sqrt(r**2 - z**2)                                        # width of cross section 
            eps = Eps[0] - z*Eps[1]                                         # strain 
            dep = Dps[0] - z*Dps[1]                                         # strain increment
            sig, dsig, _ = self.SigU( ff, CalcType, Dt, elI, ipI, Elem, [dep,0], [eps,0], 0,0, EpsR, i)
            sig_   = bb*sig[0]
            dsig_  = bb*dsig[0][0]
            nn[i]  =    sig_                                                # normal force
            mm[i]  = -z*sig_                                                # moment
            nde[i] =       dsig_                                            # normal force grad eps
            ndk[i] = -z   *dsig_                                            # normal force grad kappa
            mdk[i] =  z**2*dsig_                                            # moment grad kappa
            z = z + dz
            if abs(z-r) < ZeroD : z = r
            if i ==    0: sig_lo = sig[0]
            if i == nI-1: sig_up = sig[0]
        NN  = integrate.trapz( nn, x=None, dx=dz)                           # numerical integration normal force
        MM  = integrate.trapz( mm, x=None, dx=dz)                           # numerical integration moment (Script -> last term in Eq. (3.18))
        NDE = integrate.trapz( nde,x=None, dx=dz)                           # numerical normal force grad eps (Script -> last term in Eq. (3.20)1)
        NDK = integrate.trapz( ndk,x=None, dx=dz)                           # numerical normal force grad kappa (Script -> last term in Eq. (3.20)2)
        MDK = integrate.trapz( mdk,x=None, dx=dz)                           # numerical moment grad kappa (Script -> last term in Eq. (3.21))
        MDE = NDK
        MatM= array( [[NDE,NDK],[MDE,MDK]] )                                # material tangential stiffness (Script Eq. (3.21))
        Sig = array( [NN,MM] )                                              # stress vector
        return Sig, MatM, [Eps[0], Eps[1], Sig[0], Sig[1],sig_lo,sig_up]
    def UpdateStateVar(self, Elem, ff):
        for j in range(len(Elem.StateVar)):                                 # loop over all element integration points
            for i in range(self.nI):                         # loop over integration points same as for UpdateStateVar of RCBeam RList
                if Elem.StateVarN[j][i,1]>Elem.StateVar[j][i,1]:                  # yield stress condition
                    Elem.StateVar[j][i,0] = Elem.StateVarN[j][i,0]
                    Elem.StateVar[j][i,1] = Elem.StateVarN[j][i,1]
                for k in range(len(Elem.DataPi[j][i])): Elem.DataPi[j][i][k] = Elem.Datai[j][i][k]
#                for k in range(len(Elem.DataP[j][i])): Elem.DataP[j][i][k] = Elem.Data[j][i][k]
        return False


class MisesBeam2D(Material):                                       # elastoplastic Mises  for  Elem.dim==10 Bernoulli beam
    def __init__(self, PropMat, val):
        #    def __init__(self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,      False,     15,           10, 'MISESBEAM')
                                            # 0  lower yield stress abs       -- non-reversible state variable        0
                                            # 1  lower plastic strain signed  -- reversible state variable            1
                                            # 2  lower stress abs upon switching from plastic to elastic              2
                                            # 3  lower plastic strain upon switching from plastic to elastic (signed) 3
                                            # 4  largest plastic coordinate (signed) reached                              4
                                            # 5  index for current lower elastic (0) or plastic range (1)             5

                                            # 6  upper yield stress abs                                               6
                                            # 7  upper plastic strain signed                                          7
                                            # 8  upper switch stress                                                  8
                                            # 9  upper plastic strain upon switching from plastic to elastic (signed) 9
                                            # 10 smallet plastic coordinate (signed) reached                           10
                                            # 11 index for current upper elastic or plastic range                     11

                                            # 12 normal force                                                         12
                                            # 13 bending moment                                                       13
                                            # 14 shear force
                                            # respect also IniBeam
        self.PropMat = PropMat
        self.Emod = PropMat[0]
        self.nu   = PropMat[1]
        self.fY   = PropMat[2]
        self.epsY = self.fY/self.Emod
        self.fU   = PropMat[3]
        self.epsU = PropMat[4]
        self.EmodT= (self.fU-self.fY)/(self.epsU-self.epsY)
        if self.EmodT < 0.01:
            raise NameError("ConFemMat::MisesBeam2D_.__init__: some hardening > 0.01 required",self.EmodT)
        self.Density = PropMat[7]                                           # specific mass
    def Sig( self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR): # Eps is current strain, Dps is strain increment from previous step strain to current strain?
                                                                            # explicit strain intergraion - error prone to strain incement size
        TolA = 1.e-6
        def CrossSecRect( bb, zL, zY1, zY2, zU):
            AAlow = bb * (zY1 - zL)
            AAcen = bb * (zY2 - zY1)
            AAupp = bb * (zU - zY2)
            SSlow = 0.5 * bb * (zY1 ** 2 - zL ** 2)
            SScen = 0.5 * bb * (zY2 ** 2 - zY1 ** 2)
            SSupp = 0.5 * bb * (zU ** 2 - zY2 ** 2)
            JJlow = 1. / 3. * bb * (zY1 ** 3 - zL ** 3)
            JJcen = 1. / 3. * bb * (zY2 ** 3 - zY1 ** 3)
            JJupp = 1. / 3. * bb * (zU ** 3 - zY2 ** 3)
            return AAlow, AAcen, AAupp, SSlow, SScen, SSupp, JJlow, JJcen, JJupp
        # see CircleCrossSecValues.mws
        def CircArea( R, z_):
            if 1.-abs(z_)<TolA: ar = sign(z_)*0.5*pi
            else:               ar = arctan(z_/sqrt(1-z_**2))
            return R**2 * ( z_*sqrt(1-z_**2) + ar )
        def Circ1stMom( R, z_):
            return 2./3.* R**3 * (-1+z_**2) * sqrt(1-z_**2)
        def Circ2ndMom( R, z_):
            if 1.-abs(z_)<TolA: ar = sign(z_)*0.5*pi
            else:               ar = arctan(z_/sqrt(1-z_**2))
            return R**4 * ( -0.25*z_*sqrt(1-z_**2) + 0.5*z_**3*sqrt(1-z_**2) + 0.25*ar )
        def CrossSecCircle( rr, zL, zY1, zY2, zU):
            zL  =  zL/rr
            zY1 = zY1/rr
            zY2 = zY2/rr
            zU  = zU/rr
            AAlow =   CircArea( rr, zY1) -   CircArea( rr, zL)
            AAcen =   CircArea( rr, zY2) -   CircArea( rr, zY1)
            AAupp =   CircArea( rr, zU)  -   CircArea( rr, zY2)
            SSlow = Circ1stMom( rr, zY1) - Circ1stMom( rr, zL)
            SScen = Circ1stMom( rr, zY2) - Circ1stMom( rr, zY1)
            SSupp = Circ1stMom( rr, zU)  - Circ1stMom( rr, zY2)
            JJlow = Circ2ndMom( rr, zY1) - Circ2ndMom( rr, zL)
            JJcen = Circ2ndMom( rr, zY2) - Circ2ndMom( rr, zY1)
            JJupp = Circ2ndMom( rr, zU)  - Circ2ndMom( rr, zY2)
            return AAlow, AAcen, AAupp, SSlow, SScen, SSupp, JJlow, JJcen, JJupp
        # update section coordinates and state variables
        def LimPlas( z,z_,i, eps,dps, Eps):             # z current edge coordinaten, z_ opposite edged corrdinate, eps,dps edge strain; epsXps plastic strain upon switching from plastic to elastic (signed)
            Tol = 1.0e-12
            sigy = Elem.StateVar[ipI, i]                                    # lower/upper yield stress
            epsp = Elem.StateVar[ipI, i+1]                                  # lower/upper plastic strain
            sigs = Elem.StateVar[ipI, i+2]                                  # switch stress
            epsSp= Elem.StateVar[ipI, i+3]                                  # plastic strain belonging to switch
            zY_  = Elem.StateVar[ipI, i+4]                                  # latest coordinate of plastic zone
            delE = sigy/Emod + Tol                                          # elastic strain range -- sigy is updated below which introduces an error
            if (epsp-delE <= eps) and (eps <= epsp+delE):                   # elastic range
                zY = z                                                      # elastic coordinate with edge
                Elem.StateVarN[ipI, i+5]=0.                                 # again in elastic range
            else:                                                           # plastic range
                Elem.StateVarN[  ipI,i+5]=1.                                # going into plastic range in the following
                if Elem.StateVar[ipI,i+5]<1.: sigs = sigy                   # switch from elastic into plastic -- stress
                if abs(Eps[1]) > ZeroD:                                     # some bending there
                    epsS = eps - sign(dps)*(sigy-sigs)/EmodT                # switch strain
                    zY__ = ((Eps[0]) - (epsS-epsSp)) / (Eps[1])
                    dzY = 1./Eps[1] * Dps[0] -(Eps[0]-epsS+epsSp)/(Eps[1]**2) * Dps[1]
                    dzY2 = dzY - 1./(Eps[1]**2) *Dps[0]*Dps[1] + 2.*(Eps[0]-epsS+epsSp)/(Eps[1]**3) * Dps[1]**2 # 2nd order taylor expansion
                    zY = zY_ + dzY2
                    if z<0.:
#                        zY= max( z, zY)                                # plastic zone shouild not fall below lower edge
                        if zY<z:  zY = z                                    # current edge
                        if zY>z_: zY = z_                                   # opposite to current edge
                    else:
#                        zY= min( z, zY)
                        if zY>z:  zY = z                                    # current edge
                        if zY<z_: zY = z_                                   # opposite to current edge
                else:
                    zY = 0.                                                 # plastic coordinate down to reference axis
                epsp = epsp + dps*(1.-EmodT/Emod)
                sigy = sigy + abs(EmodT * dps)
            return zY, epsp, sigy, sigs, delE                               # coordinate for plastic area, plastic strain, largest stress reached, switch stress?, elastic strain range
        # end of def
        if CalcType == 0: return [], [], []
        Emod = self.Emod
        EmodT= self.EmodT
        if Elem.CrossSecType == "RECT":
            bb = Elem.Geom[1, 1]
            hh = Elem.Geom[1, 2]
            zL   = -0.5*hh
            zU   =  0.5*hh
        elif Elem.CrossSecType=="CIRCLE":
            rr = 0.5*Elem.Geom[1,2]                                         # Diam/2
            zL = -rr
            zU =  rr
        else:
            raise NameError("ConFemMat::MisesBeam2D.Sig: cross section type not yet implemented",Elem.CrossSecType )
        epsL = Eps[0] - zL*Eps[1]                                           # current lower strain
        epsU = Eps[0] - zU*Eps[1]                                           # current upper strain
        DpsL = Dps[0] - zL*Dps[1]                                           # current lower strain increment
        DpsU = Dps[0] - zU*Dps[1]                                           # current upper strain increment
        #
        zY1, epsLp, sigLy,sigLs, delEL = LimPlas( zL,zU,0, epsL,DpsL, Eps)# returns coordinate for plastic area, plastic strain, largest stress reached, switch stress?, elastic strain range
        zY2, epsUp, sigUy,sigUs, delEU = LimPlas( zU,zL,6, epsU,DpsU, Eps)
        # cross sectional values
        if Elem.CrossSecType=="RECT":
            AAlow, AAcen, AAupp, SSlow, SScen, SSupp, JJlow, JJcen, JJupp = CrossSecRect( bb, zL, zY1, zY2, zU)
        elif Elem.CrossSecType=="CIRCLE":
            AAlow, AAcen, AAupp, SSlow, SScen, SSupp, JJlow, JJcen, JJupp = CrossSecCircle(rr, zL, zY1, zY2, zU)
        # tangential material stiffness, internal forces
        if Elem.dim in [10, 12]:                                            # bernoulli beam, axisymmetric bernoulli beam
            MatM = array([[  EmodT*AAlow+Emod*AAcen+EmodT*AAupp,                 -(EmodT*SSlow+Emod*SScen+EmodT*SSupp)],
                          [-(EmodT*SSlow+Emod*SScen+EmodT*SSupp),  self.matbStiff*(EmodT*JJlow+Emod*JJcen+EmodT*JJupp)]])
        elif Elem.dim in [11, 13]:                                          # timoshenko beam, axisymmetric timoshenko beam
            alpha = 0.8
            AA = Elem.Geom[1, 3]                                            # cross sectional area
            GG = AA*alpha*0.5*Emod/(1+self.nu)
            MatM = array([[  EmodT*AAlow+Emod*AAcen+EmodT*AAupp,                 -(EmodT*SSlow+Emod*SScen+EmodT*SSupp), 0],
                          [-(EmodT*SSlow+Emod*SScen+EmodT*SSupp),  self.matbStiff*(EmodT*JJlow+Emod*JJcen+EmodT*JJupp),  0],
                          [0, 0, GG]])
        else:
            raise NameError("ConFemMat::MisesBeam2D.Sig: element type not allowed",Elem.Type )
        DelSig = dot(MatM,Dps)
        if Elem.dim in [10, 12]: Sig = array([Elem.StateVar[ipI,12]+DelSig[0], Elem.StateVar[ipI,13]+DelSig[1]])
        else:                    Sig = array([Elem.StateVar[ipI,12]+DelSig[0], Elem.StateVar[ipI,13]+DelSig[1], Elem.StateVar[ipI,14]+DelSig[2]])
        # state variables
        Elem.StateVarN[ipI, 0] = sigLy                                      # lower yield stress
        Elem.StateVarN[ipI, 1] = epsLp                                      # lower plastic strain
        Elem.StateVarN[ipI, 2] = sigLs                                      # lower switch stress
        # 3 is updated in UpdateStateVar
        Elem.StateVarN[ipI, 4] = zY1
        Elem.StateVarN[ipI, 6] = sigUy                                      # upper yield stress
        Elem.StateVarN[ipI, 7] = epsUp                                      # upper plastic strain
        Elem.StateVarN[ipI, 8] = sigUs                                      # upper switch stress
        # 9 is updated in UpdateStateVar
        Elem.StateVarN[ipI,10] = zY2
        Elem.StateVarN[ipI,12] = Sig[0]                                     # normal force
        Elem.StateVarN[ipI,13] = Sig[1]                                     # bending moment
        if Elem.dim in [11, 13]: Elem.StateVarN[ipI,14] = Sig[2]
        #
        if Elem.dim in [10, 12]: return Sig, MatM, [ Eps[0],Eps[1], Sig[0],Sig[1], sigLy,sigUy, zY1,zY2]
        else:                    return Sig, MatM, [ Eps[0],Eps[1],Eps[2], Sig[0],Sig[1],Sig[2], sigLy,sigUy, zY1,zY2 ]

    def UpdateStateVar(self, Elem, ff):
        for ipI in range(Elem.StateVar.shape[0]):                           # loop over integration points

            if Elem.StateVarN[ipI,0] > Elem.StateVar[ipI,0]:                # lower hardening with increasing yield stress (unsigned)
                Elem.StateVar[ipI,0] = Elem.StateVarN[ipI,0]                # update of lower yield stress
                Elem.StateVar[ipI,1] = Elem.StateVarN[ipI,1]                # update of lower plastic strain
            Elem.StateVar[ipI,2] = Elem.StateVarN[ipI,2]                    # lower switch stress
            if Elem.StateVarN[ipI,4] > Elem.StateVar[ipI,4]:
                Elem.StateVar[ipI,4] = Elem.StateVarN[ipI,4]                #
            if (Elem.StateVar[ipI,5]>0.) and (Elem.StateVarN[ipI,5]<1.):    # lower from plastic to elastic
                Elem.StateVar[ipI,3] = Elem.StateVarN[ipI, 1]               # refers to lower plastic strain
            Elem.StateVar[ipI,5] = Elem.StateVarN[ipI,5]                    # index for current lower elastic or plastic range

            if Elem.StateVarN[ipI,6] > Elem.StateVar[ipI,6]:                # upper hardening with increasing yield stress (unsigned)
                Elem.StateVar[ipI,6] = Elem.StateVarN[ipI,6]
                Elem.StateVar[ipI,7] = Elem.StateVarN[ipI,7]
            Elem.StateVar[ipI,8] = Elem.StateVarN[ipI,8]                    # upper switch stress
            if Elem.StateVarN[ipI,10] < Elem.StateVar[ipI,10]:
                Elem.StateVar[ipI,10] = Elem.StateVarN[ipI,10]          #
            if (Elem.StateVar[ipI,11]>0.) and (Elem.StateVarN[ipI,11]<1.):  # upper from plastic to elastic
                Elem.StateVar[ipI, 9] = Elem.StateVarN[ipI, 7]              # refers to upper plastic strain
            Elem.StateVar[ipI,11] = Elem.StateVarN[ipI,11]                  # index for current upper elastic or plastic range

            Elem.StateVar[ipI,12] = Elem.StateVarN[ipI,12]                  # normal force
            if False: #Elem.StateVar[ipI,13]*Elem.StateVarN[ipI,13]<0:
                raise NameError("ConFemMat::MisesBeam2D_.Sig: moment sign change not allowed",Elem.StateVar[ipI,13],Elem.StateVarN[ipI,13])
            else:
                Elem.StateVar[ipI,13] = Elem.StateVarN[ipI,13]                  # bending moment
            Elem.StateVar[ipI, 14] = Elem.StateVarN[ipI, 14]                # shear force
        return False

class Lubliner(Material):                              # elastoplastic Mises
    def __init__(self, PropMat):
#    def __init__(self,         SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, False,        None,     True,      False,     17,           8, "Lubliner")
#        self.Symmetric = True                       # flag for symmetry of material matrices
#        self.Update = True                          # has specific update procedure for update of state variables
#        self.Updat2 = False                         # no 2 stage update procedure
#        self.StateVar = 17                           # number of state variables per integration point
                                                    # 16: PFlag of latest step: False -> 0, True -> 1, strictly there should be one for compression and one for tension? not quite shure about this
#        self.NData = 8                             # number of data items

        self.E_0  = PropMat[0]                      # initial Young's modulus
        self.Emod = PropMat[0]                      # for compatibility - initial Young's modulus
        self.nu   = PropMat[1]                      # Poissons's ratio
        self.f_c0 = PropMat[2]                      # elastic uniaxial compressive yield (?) (unsigned)
        self.f_cm = PropMat[3]                      # uniaxial compressive strength (unsigned)
        self.f_t0 = PropMat[4]                      # uniaxial tensile strength
        alpha_    = PropMat[5]                      # ratio of biaxial compressive strength to uniaxial compressive strength
        K_c       = PropMat[6]                      # ratio of strength on tensile meridian compared to compressive meridian
        self.G_F  = PropMat[7]                      # tensile cracking energy
        self.G_ch = PropMat[8]                      # compressive cracking energy
        self.alpha_p = PropMat[9]*pi/180.           # angle of dilatancy (input in degrees, transformed to rad)
        self.ecc  = PropMat[10]                     # eccentricity of plastic potential surface
        self.CalDamC = PropMat[11]                  # assumed damage for maximum uniaxial compression stress
        self.CalDamT = PropMat[12]                  # assumed damage for half maximum uniaxial tensile stress
#        self.CalDamC = 0.4                          # assumed damage for maximum uniaxial compression stress
#        self.CalDamT = 0.5                          # assumed damage for half maximum uniaxial tensile stress
        self.Density = PropMat[13]                  # specific mass
        self.ff = self.E_0*(1.-self.nu)/((1.+self.nu)*(1.-2*self.nu))
        # uniaxial
        xx = self.f_cm/self.f_c0
        self.a_c  = 2.*xx - 1.+2.*sqrt(xx**2-xx)
        self.dbbC = log(1.-self.CalDamC) / log((1.+self.a_c)/(2.*self.a_c)) # Eq. (79)
        self.a_t  = 1.0                                                     # xx=1 in case of tension
        self.dbbT = log(1.-self.CalDamT) / log((1.+self.a_t-sqrt(1.+self.a_t**2))/(2.*self.a_t)) # Eq. (84)
        # multiaxial
        self.alpha = (alpha_-1.)/(2.*alpha_-1.)
        self.gamma = 3.*(1.-K_c)/(2.*K_c-1.)
        # elastic contributions
        nu = self.nu                                # Poisson's ratio
        self.C0 = self.E_0*(1-nu)/((1+nu)*(1-2*nu))*array([[1,nu/(1-nu),nu/(1-nu),0,0,0],
                                                       [nu/(1-nu),1,nu/(1-nu),0,0,0],
                                                       [nu/(1-nu),nu/(1-nu),1,0,0,0],
                                                       [0,0,0,(1-2*nu)/(2*(1-nu)),0,0],
                                                       [0,0,0,0,(1-2*nu)/(2*(1-nu)),0],
                                                       [0,0,0,0,0,(1-2*nu)/(2*(1-nu))]]) # triaxial isotropic elasticity
        self.ZeroTol = 1.e-3
        self.KapTol  = 1.e-4
#        self.KapTol  = 1.e-6                                                # parameter for plastic strain, kappa iteration with relatively high sensitivity
    def fC(self, kappa, ipI,Label):                                         # used by Lubliner::CalcKapC
        phi = 1. + self.a_c*(2.+self.a_c)*kappa
        fC = self.f_c0 *( (1.+self.a_c)*sqrt(phi) - phi )/self.a_c
        dfC= self.f_c0*( 0.5*(1.+self.a_c)*(2.*self.a_c+self.a_c**2)/sqrt(phi)-self.a_c*(2.+self.a_c) )/self.a_c
#        if kappa>=1.: raise NameError(f"ConFemMaterials::Lubliner.fC: unvalid compressive damage, {Label:d}, {ipI:d}__{phi:f},{fC:f},{dfC:.4e}")
        return fC, dfC
    def fT(self, kappa, ipI,Label):                                         # used by Lubliner::CalcKapC
        phi = 1. + self.a_t*(2.+self.a_t)*kappa
        fT = self.f_t0 *( (1.+self.a_t)*sqrt(phi) - phi )/self.a_t
        dfT= self.f_t0*( 0.5*(1.+self.a_t)*(2.*self.a_t+self.a_t**2)/sqrt(phi)-self.a_t*(2.+self.a_t) )/self.a_t
#        if kappa>=1.: raise NameError(f"ConFemMaterials::Lubliner.fC: unvalid tensile damage, {Label:d}, {ipI:d}__{phi:f},{fT:f},{dfT:.4e}")
        return fT, dfT
    def FF(self, I1, J2, sig_max, fC, fT):                                  # Yield function, used by Lubliner::GP, Lubliner::sig
        if   sig_max<-self.ZeroTol: 
            return (sqrt(3.*J2) + self.alpha*I1 + self.gamma*sig_max)/(1.-self.alpha)
        elif sig_max< self.ZeroTol: 
            return (sqrt(3.*J2) + self.alpha*I1 )                    /(1.-self.alpha)
        else:
            beta_ =  fC/fT*(1.-self.alpha) - (1.+self.alpha)
            return (sqrt(3.*J2) + self.alpha*I1 + beta_*sig_max)/(1.-self.alpha)
    def dFF(self, J2, sig_dev, nd_, sig_max, fC, dfC, fT, dfT, kapC, kapT, DC, DT):  # gradient of yield function, used by Lubliner::GP
        x = 0.5*sqrt(3./J2)
        y = 1./(1.-self.alpha)
        a_c = self.a_c
        a_t = self.a_t
        dbbC = self.dbbC
        dbbT = self.dbbT
        f_c0 = self.f_c0
        phiC  = 1.+a_c*(2.+a_c)*kapC
        RphiC = sqrt(phiC)
        xx = 1. + a_c - RphiC
        yy = 0.5*f_c0*(2.+a_c)/(1.-DT)*pow(xx/a_c,-dbbC)/(xx*RphiC) * ( (2.-dbbC)*phiC + (dbbC+dbbC*a_c-3.-3.*a_c)*RphiC + (1.+a_c)**2 ) 
        phiT  = 1.+a_t*(2.+a_t)*kapT
        RphiT = sqrt(phiT)
        xx = 1. + a_t - RphiT 
        zz = 0.5*fC/(1.-DC) * dbbT    *pow(xx/a_t,-dbbT)/(xx*RphiT) * a_t*(2.+a_t)
        dfC_eff = array([yy,zz])                                                    # gradient of effective compression "cohesion" with respect to internal state variables kappa_c, kappa_t
        if   sig_max<-self.ZeroTol:
#            z =  self.alpha +self.gamma
            return array([ y*(x*sig_dev[0]+self.alpha +self.gamma*nd_[0]),
                           y*(x*sig_dev[1]+self.alpha +self.gamma*nd_[1]),
                           y*(x*sig_dev[2]+self.alpha +self.gamma*nd_[2]),
                           y*(x*sig_dev[3]            +self.gamma*nd_[3]),
                           y*(x*sig_dev[4]            +self.gamma*nd_[4]),
                           y*(x*sig_dev[5]            +self.gamma*nd_[5])]), array([ -dfC_eff[0], -dfC_eff[1]])  
        elif sig_max< self.ZeroTol: 
            return array([ y*(x*sig_dev[0]+self.alpha),
                           y*(x*sig_dev[1]+self.alpha),
                           y*(x*sig_dev[2]+self.alpha),
                           y*(x*sig_dev[3]),
                           y*(x*sig_dev[4]),
                           y*(x*sig_dev[5])]), array([ -dfC_eff[0], -dfC_eff[1]])  
        else:                       
            beta_ = fC/fT*(1.-self.alpha) - (1.+self.alpha)
            return array([ y*(x*sig_dev[0]+self.alpha +beta_*nd_[0]),
                           y*(x*sig_dev[1]+self.alpha +beta_*nd_[1]),
                           y*(x*sig_dev[2]+self.alpha +beta_*nd_[2]),
                           y*(x*sig_dev[3]            +beta_*nd_[3]),
                           y*(x*sig_dev[4]            +beta_*nd_[4]),
                           y*(x*sig_dev[5]            +beta_*nd_[5])]), array([ dfC*sig_max/fT - dfC_eff[0], -sig_max*dfT*fC/(fT**2) - dfC_eff[1]]) # ???
    def dGG(self, J2, sig_dev):                                 # gradient of flow potential, used by Lubliner::GP
        x = 1.5/sqrt( (self.ecc*self.f_t0*tan(self.alpha_p))**2 + 3.*J2)
        y = tan(self.alpha_p)/3.
        return array([x*sig_dev[0] + y,
                      x*sig_dev[1] + y,
                      x*sig_dev[2] + y,
                      x*sig_dev[3],
                      x*sig_dev[4],
                      x*sig_dev[5]])
    def GP(self, sig, fC, dfC, fT, dfT, kapC, kapT, DC, DT, Flag):                  # whole plastic stuff for given stress state
        I1  = sig[0] + sig[1] + sig[2]
        ee, ev  = eigh(array([ [sig[0], sig[5], sig[4]] , [sig[5], sig[1], sig[3]] , [sig[4], sig[3], sig[2]] ])) # to find largest principal value
        Psig_max = ee[2]                                                            # ascending order, largest (signed!) on last position 
        J2  = ((ee[0]-ee[1])**2+(ee[0]-ee[2])**2+(ee[1]-ee[2])**2)/6.
        FF  = self.FF( I1, J2, Psig_max, fC, fT)       # Yield function
        if not Flag: return FF
        pp  = I1/3.
#        if abs(ee[2]-ee[1])<self.ZeroTol:                       # same two largest principal stresses
#            ev[0,2] = ev[0,1]+ev[0,2]                           # second index indicates eigenvector
#            ev[1,2] = ev[1,1]+ev[1,2]
#            ev[2,2] = ev[2,1]+ev[2,2]
#            eL = sqrt( ev[0,2]**2 + ev[1,2]**2 + ev[2,2]**2 )
#            ev[0,2] = ev[0,2]/eL                           # second index indicates eigenvector
#            ev[1,2] = ev[1,2]/eL
#            ev[2,2] = ev[2,2]/eL
        sig_dev  = array([sig[0]-pp, sig[1]-pp, sig[2]-pp, sig[3], sig[4], sig[5]])
        nd_      = array([ev[0,2]*ev[0,2],ev[1,2]*ev[1,2],ev[2,2]*ev[2,2],ev[1,2]*ev[2,2],ev[0,2]*ev[2,2],ev[0,2]*ev[1,2]]) # gradient of largest principal stress corresponding to Voigt notation, largest eigenvalue in last position (ascending order)
        dff, dkk = self.dFF( J2, sig_dev, nd_, Psig_max, fC, dfC, fT, dfT, kapC, kapT, DC, DT)      # gradient of yield function with respect to stress
        dgg      = self.dGG( J2, sig_dev)                                           # gradient of flow potential
        return FF, dff, dkk, dgg
    def CalcKapC(self, g_ch, g_F, r, dEpsP, kapC_old, kapT_old, ipI,Label):
        ni, i_, j_  = 10, 0, 0
        # compression
        kapC_, DkC = kapC_old, 0.
        if dEpsP[0]<0.:
            for i_ in range(ni):
                fC_, dfC_ = self.fC(kapC_, ipI,Label)
                CC = (1.-r)/g_ch * dEpsP[0]
                RRC= DkC + CC*fC_
                JJ = 1. + CC*dfC_
                dkC = -RRC/JJ
                DkC = DkC + dkC
                kapC_ = kapC_old + DkC 
                if kapC_>=1.: raise NameError(f"ConFemMaterials::Lubliner.CalcKapC: invalid compressive damage, {Label:d}, {ipI:d}__{i_},{kapC_old:f},{DkC:f},{dEpsP[0]:.4e}")
                if abs(RRC)<self.KapTol: break
        else:
            fC_, dfC_ = self.fC(kapC_, ipI,Label)
        # tension
        kapT_, DkT = kapT_old, 0.
        if dEpsP[2]>0:
            for j_ in range(ni):
                fT_, dfT_ = self.fT(kapT_, ipI,Label)
                CC = r/g_F*dEpsP[2]
                RRT= DkT - CC*fT_
                JJ = 1. - CC*dfT_
                dkT = - RRT/JJ
                DkT = DkT + dkT
                kapT_ = kapT_old + DkT
                if kapT_>=1.: raise NameError(f"ConFemMaterials::Lubliner.CalcKapT: invalid tensile damage, {Label:d}, {ipI:d}__{j_:d},{kapT_old:f},{DkT:f},{dEpsP[2]:.4e}")
                if abs(RRT)<self.KapTol: break
        else:
            fT_, dfT_ = self.fT(kapT_, ipI,Label)
        #
        if (i_>=ni-1) or (j_>=ni-1): raise NameError ("ConFemMaterials::Lubliner.CalcKapC: no convergence A")
        return kapC_, kapT_, fC_, dfC_, fT_, dfT_
    def Dam(self, kapC, kapT):                                              # used by Lubliner::sig
        a_c = self.a_c
        phiC= 1. + a_c*(2.+a_c)*kapC
        RpC= sqrt(phiC)
        z  = 1. + a_c - RpC
        x  = z / a_c
        y  = pow(x,self.dbbC)
        DC = 1. - y 
        GC = self.dbbC * y/z * (2.*a_c+a_c**2)/(2.*RpC) 
        #
        a_t = self.a_t
        phiT= 1. + a_t*(2.+a_t)*kapT
        RpT= sqrt(phiT)
        z  = 1. + a_t - RpT
        x  = z / a_t
        y = pow(x,self.dbbT)
        DT = 1. - y 
        GT = self.dbbT * y/z * (2.*a_t+a_t**2)/(2.*RpT)
        #
        return DC, DT, GC, GT               

    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps_, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        kapT_old = Elem.StateVar[ipI][0]
        kapT_    = kapT_old
        kapC_old = Elem.StateVar[ipI][1]
        kapC_    = kapC_old
        fc_old   = Elem.StateVar[ipI][2] + self.f_c0
        ft_old   = Elem.StateVar[ipI][3] + self.f_t0
        DC       = Elem.StateVar[ipI][4]                                            # scalar damage
        DT       = Elem.StateVar[ipI][5]
        DD       = 1.-(1.-DC)*(1.-DT)                                               # update damage
        Lch      = Elem.Lch
        g_ch     = self.G_ch/Lch                                                    # given dissipated energy per surface (!) modified to yield energy per volume
        g_F      = self.G_F/Lch
        EpsP_old = array([Elem.StateVar[ipI][6],Elem.StateVar[ipI][7],Elem.StateVar[ipI][8],Elem.StateVar[ipI][9],Elem.StateVar[ipI][10],Elem.StateVar[ipI][11]]) # old plastic strain in global system
        LdEpsP   = zeros((6), dtype=float)
        EpsP     = zeros((6), dtype=float)
        Eps2     = Elem.StateVar[ipI][12]                                           # lateral strain in case of plane stress
        #
        PlStFlag = False
        if (Elem.dim==2 and Elem.PlSt) or Elem.dim==21: 
            PlStFlag = True
            Eps2     = Elem.StateVar[ipI][12]                                       # lateral strain in case of plane stress
        if Elem.dim==2:                                                             # plate plane stress / plane strain
            Eps  = array([Eps_[0],Eps_[1],0.,0.,0.,Eps_[2]])                        # total strain increment
            dEps = array([Dps[0], Dps[1], 0.,0.,0.,Dps[2]])                         # total strain increment
            if PlStFlag: 
                Eps[2]=-( self.C0[2,0]*(Eps[0]-EpsP_old[0])
                         +self.C0[2,1]*(Eps[1]-EpsP_old[1])
                         +self.C0[2,3]*(Eps[3]-EpsP_old[3])
                         +self.C0[2,4]*(Eps[4]-EpsP_old[4])
                         +self.C0[2,5]*(Eps[5]-EpsP_old[5]))/self.C0[2,2]+EpsP_old[2] # lateral strain for plane stress -- (1-D) falls out
        elif Elem.dim==3:
            Eps  = array([Eps_[0],Eps_[1],Eps_[2],Eps_[3],Eps_[4],Eps_[5]])         # total strain
            dEps = array([Dps[0], Dps[1], Dps[2], Dps[3], Dps[4], Dps[5]])          # total strain increment
        else: raise NameError ("CaeFemMaterials::Lubliner.Sig: not implemented for this element type")
        #
        SigEff = dot(self.C0,(Eps-EpsP_old))                                        # trial stress
        fc_eff = fc_old/(1.-DD)                                                     # current effective uniaxial compressive strength
        FF     = self.GP( SigEff, fc_old, None, ft_old, None, None, None, DC, DT, False) # yield function
        PFlag  = False                                                              # initial setting for flag for plastic range
        Elem.StateVarN[ipI][16] = 0                                                 # holds latest PFlag for update control -> initialize
        if norm(Dps)>0.001*ZeroD and FF>=fc_eff:
            PFlag = True                                                            # plasticity / damage occurs in current iteration
            dEpsP = zeros((6), dtype=float)                                        # to sum up plastic strain sub-increments used for next plane stress iteration 
            np_ = 30 # 30
            for j_ in range(np_):                                                   # lateral strain iteration for plane stress
                kapT_old = Elem.StateVar[ipI][0]                                    # this should be updated regarding sub-incrementing and has to be reinitialized
                kapC_old = Elem.StateVar[ipI][1]                                    # this should be updated regarding sub-incrementing and has to be reinitialized
                EpsP_old = array([Elem.StateVar[ipI][6],Elem.StateVar[ipI][7],Elem.StateVar[ipI][8],Elem.StateVar[ipI][9],Elem.StateVar[ipI][10],Elem.StateVar[ipI][11]]) # re-initialize old plastic strain in global system
                if PlStFlag:
                    dEps[2]=-( self.C0[2,0]*(dEps[0]-dEpsP[0])
                              +self.C0[2,1]*(dEps[1]-dEpsP[1])
                              +self.C0[2,3]*(dEps[3]-dEpsP[3])
                              +self.C0[2,4]*(dEps[4]-dEpsP[4])
                              +self.C0[2,5]*(dEps[5]-dEpsP[5]))/self.C0[2,2]+dEpsP[2] # lateral strain for plane stress -- (1-D) falls out
                    Eps[2]=-( self.C0[2,0]*(Eps[0]-EpsP_old[0]-dEpsP[0])
                             +self.C0[2,1]*(Eps[1]-EpsP_old[1]-dEpsP[1])
                             +self.C0[2,3]*(Eps[3]-EpsP_old[3]-dEpsP[3])
                             +self.C0[2,4]*(Eps[4]-EpsP_old[4]-dEpsP[4])
                             +self.C0[2,5]*(Eps[5]-EpsP_old[5]-dEpsP[5]))/self.C0[2,2]+EpsP_old[2]+dEpsP[2] # lateral strain for plane stress -- (1-D) falls out
                dEpsPr, gg_ = eigh(array([ [dEps[0], dEps[5], dEps[4]] , [dEps[5], dEps[1], dEps[3]] , [dEps[4], dEps[3], dEps[2]] ]))  # principal values and eigenvectors of strain increment tensor, smallest (signed) 1st!
                gg_[0,2] = gg_[1,0]*gg_[2,1]-gg_[2,0]*gg_[1,1]                      # to have a right handed system
                gg_[1,2] = gg_[2,0]*gg_[0,1]-gg_[0,0]*gg_[2,1]
                gg_[2,2] = gg_[0,0]*gg_[1,1]-gg_[1,0]*gg_[0,1]
                gg = array([[gg_[0,0],gg_[1,0],gg_[2,0]],[gg_[0,1],gg_[1,1],gg_[2,1]],[gg_[0,2],gg_[1,2],gg_[2,2]]])
                TTS = array([[ gg[0,0]**2, gg[0,1]**2, gg[0,2]**2, 2.*gg[0,1]*gg[0,2], 2.*gg[0,0]*gg[0,2], 2.*gg[0,0]*gg[0,1]],
                             [ gg[1,0]**2, gg[1,1]**2, gg[1,2]**2, 2.*gg[1,1]*gg[1,2], 2.*gg[1,0]*gg[1,2], 2.*gg[1,0]*gg[1,1]],
                             [ gg[2,0]**2, gg[2,1]**2, gg[2,2]**2, 2.*gg[2,1]*gg[2,2], 2.*gg[2,0]*gg[2,2], 2.*gg[2,0]*gg[2,1]],
                             [ gg[1,0]*gg[2,0], gg[1,1]*gg[2,1], gg[1,2]*gg[2,2], gg[1,2]*gg[2,1]+gg[1,1]*gg[2,2], gg[1,2]*gg[2,0]+gg[1,0]*gg[2,2], gg[1,1]*gg[2,0]+gg[1,0]*gg[2,1]],
                             [ gg[0,0]*gg[2,0], gg[0,1]*gg[2,1], gg[0,2]*gg[2,2], gg[0,2]*gg[2,1]+gg[0,1]*gg[2,2], gg[0,2]*gg[2,0]+gg[0,0]*gg[2,2], gg[0,1]*gg[2,0]+gg[0,0]*gg[2,1]],
                             [ gg[0,0]*gg[1,0], gg[0,1]*gg[1,1], gg[0,2]*gg[1,2], gg[0,2]*gg[1,1]+gg[0,1]*gg[1,2], gg[0,2]*gg[1,0]+gg[0,0]*gg[1,2], gg[0,1]*gg[1,0]+gg[0,0]*gg[1,1]]])
                TTD = array([[ gg[0,0]**2, gg[0,1]**2, gg[0,2]**2,    gg[0,1]*gg[0,2],    gg[0,0]*gg[0,2],    gg[0,0]*gg[0,1]],
                             [ gg[1,0]**2, gg[1,1]**2, gg[1,2]**2,    gg[1,1]*gg[1,2],    gg[1,0]*gg[1,2],    gg[1,0]*gg[1,1]],
                             [ gg[2,0]**2, gg[2,1]**2, gg[2,2]**2,    gg[2,1]*gg[2,2],    gg[2,0]*gg[2,2],    gg[2,0]*gg[2,1]],
                             [ 2.*gg[1,0]*gg[2,0], 2.*gg[1,1]*gg[2,1], 2.*gg[1,2]*gg[2,2], gg[1,2]*gg[2,1]+gg[1,1]*gg[2,2], gg[1,2]*gg[2,0]+gg[1,0]*gg[2,2], gg[1,1]*gg[2,0]+gg[1,0]*gg[2,1]],
                             [ 2.*gg[0,0]*gg[2,0], 2.*gg[0,1]*gg[2,1], 2.*gg[0,2]*gg[2,2], gg[0,2]*gg[2,1]+gg[0,1]*gg[2,2], gg[0,2]*gg[2,0]+gg[0,0]*gg[2,2], gg[0,1]*gg[2,0]+gg[0,0]*gg[2,1]],
                             [ 2.*gg[0,0]*gg[1,0], 2.*gg[0,1]*gg[1,1], 2.*gg[0,2]*gg[1,2], gg[0,2]*gg[1,1]+gg[0,1]*gg[1,2], gg[0,2]*gg[1,0]+gg[0,0]*gg[1,2], gg[0,1]*gg[1,0]+gg[0,0]*gg[1,1]]])
                
                # loop over strain sub-increment size
                kList = [1,2,4,8,12,16]     # > 1 currently does not work for uniax2d reloading compression transition elastic plastic
#                kList = [12]
                for nk_ in kList:
                    kf = 1./nk_
                    LdEps  = kf*array([dEpsPr[0],dEpsPr[1],dEpsPr[2],0.,0.,0.]) # extend principal strains to length 6 for local strain
                    LdEpsP[:] = LdEps[:]
                    dEpsP[:] = 0.                                          # to sum up plastic strain sub-increments used for next plane stress iteration 
                    # loop over strain sub-increments
                    for k_ in range(nk_):
                        XXX = []
                        NotReady = False
                        ni_ = 10 # 10 #15 #10
                        # loop over inner iteration 
                        for i_ in range(ni_):                               # iteration loop for implicit computations of kap, plastic strain and stress
                            LSigEff   = dot(TTS,SigEff)                     # transform stresses to actual local system
                            rr = (max(LSigEff[0],0.)+max(LSigEff[1],0.)+max(LSigEff[2],0.)) / (abs(LSigEff[0])+abs(LSigEff[1])+abs(LSigEff[2])) # stress weight factor
                            if 1.-abs(rr)<1.e-2: rr=1.0                     # to avoid spurious effects from intermediate iterated stresses
                            LdEpsP[0] = (1.-rr)*LdEpsP[0]                   # compression - presumably to avoid spurious lateral strain influence
                            LdEpsP[2] =     rr *LdEpsP[2]                   # tension - presumably to avoid spurious lateral strain influence
                            kapC, kapT, fC, dfC, fT, dfT = self.CalcKapC(g_ch, g_F, rr, LdEpsP, kapC_old, kapT_old, ipI,Elem.Label)
                            DC, DT, GC, GT = self.Dam(kapC, kapT)
                            # compute plastic strain increment
                            _, dff, dkk, dgg = self.GP( LSigEff, fC, dfC, fT, dfT, kapC, kapT, DC, DT, True)   # gradients should be evaluated in local system
                            GF = array([[dgg[0]*dff[0],dgg[0]*dff[1],dgg[0]*dff[2],dgg[0]*dff[3],dgg[0]*dff[4],dgg[0]*dff[5]],
                                        [dgg[1]*dff[0],dgg[1]*dff[1],dgg[1]*dff[2],dgg[1]*dff[3],dgg[1]*dff[4],dgg[1]*dff[5]],
                                        [dgg[2]*dff[0],dgg[2]*dff[1],dgg[2]*dff[2],dgg[2]*dff[3],dgg[2]*dff[4],dgg[2]*dff[5]],
                                        [dgg[3]*dff[0],dgg[3]*dff[1],dgg[3]*dff[2],dgg[3]*dff[3],dgg[3]*dff[4],dgg[3]*dff[5]],
                                        [dgg[4]*dff[0],dgg[4]*dff[1],dgg[4]*dff[2],dgg[4]*dff[3],dgg[4]*dff[4],dgg[4]*dff[5]],
                                        [dgg[5]*dff[0],dgg[5]*dff[1],dgg[5]*dff[2],dgg[5]*dff[3],dgg[5]*dff[4],dgg[5]*dff[5]]])
                            zz = dot(self.C0,dgg)
                            yy = dff[0]*zz[0]+dff[1]*zz[1]+dff[2]*zz[2]+dff[3]*zz[3]+dff[4]*zz[4]+dff[5]*zz[5]
                            h0 = -fC*(1.-rr)/g_ch
                            h2 =  fT*    rr /g_F
                            zz_= -dkk[0]*h0*dgg[0] -dkk[1]*h2*dgg[2] + yy
                            EP = 1./zz_ * dot(GF,self.C0)                   # plastic operator
                            LdEpsP[:] = dot(EP,LdEps)[:]                    # local plastic strain increment from local strain increment (constant within this loop) and plastic operator
                            dEpsP_  = dot(transpose(TTS),LdEpsP)            # transform local plastic strain increment back to global system
                            EpsP[:]= (EpsP_old+dEpsP_)[:]
                            SigEff = dot(self.C0,(Eps - EpsP))              # compute effective stress
                            
#                            if Elem.Label in [421] and ipI in [0] and nk_>0: # and k_>0:
#                                print(f'j{j_:d},k{k_:d}{nk_:d},i{i_:d}, LEpsP   {", ".join(zzz for zzz in [f"{x:11.4e}" for x in LEpsP[0:6]])}, res {", ".join(zzz for zzz in [f"{x:11.4e}" for x in XXX[-5:-1]])}', file=ff)
#                                print(f'j{j_:d},k{k_:d}{nk_:d},i{i_:d}, LdEpsP {", ".join(zzz for zzz in [f"{x:11.4e}" for x in LdEpsP[0:6]])}', file=ff)
#                                print(f'j{j_:d},k{k_:d}{nk_:d},i{i_:d},rr{rr:4f}, LEpsP  {",".join(z for z in [f"{x:11.4e}" for x in LEpsP[0:6]])}, resx {",".join(z for z in [f"{x:11.4e}" for x in XXX[-5:-1]])},{DC:f},{DT:f}')
        ##                        print(f'sigL {", ".join(zzz for zzz in [f"{x:8.4f}" for x in LSigEff_[0:6]])}', file=ff)
        ##                        print(f'ff {", ".join(zzz for zzz in [f"{x:8.4f}" for x in dff[0:6]])}', file=ff)
        ##                        print(f'gg {", ".join(zzz for zzz in [f"{x:8.4f}" for x in dgg[0:6]])}', file=ff)
        ##                        print(f'PE0 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[0,0:6]])}', file=ff)
        ##                        print(f'PE1 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[1,0:6]])}', file=ff)
        ##                        print(f'PE2 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[2,0:6]])}', file=ff)
        ##                        print(f'PE3 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[3,0:6]])}', file=ff)
        ##                        print(f'PE4 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[4,0:6]])}', file=ff)
        ##                        print(f'PE5 {", ".join(zzz for zzz in [f"{x:8.4f}" for x in EP[5,0:6]])}', file=ff)
        #                        print(f'{Elem.Label:d},{ipI:d},{j_:d},{i_:d}', file=ff)
        #                        print(f'prinstrains {", ".join(zzz for zzz in [f"{x:11.4e}" for x in LdEps[0:6]])}', file=ff)
        ##                        print(f'strain {", ".join(zzz for zzz in [f"{x:11.4e}" for x in Eps[0:6]])}', file=ff)
        ##                        print(f'strainP{", ".join(zzz for zzz in [f"{x:11.4e}" for x in EpsP[0:6]])}', file=ff)
        ##                        print(f'sigEff {", ".join(zzz for zzz in [f"{x:8.4f}" for x in SigEff[0:6]])}', file=ff)
                            
                            if sqrt( (kapC_-kapC)**2 + (kapT_-kapT)**2 ) < self.KapTol*sqrt( kapC**2 + kapT**2 ): break # convergence control
                            dkapC, dkapT = kapC_-kapC, kapT_-kapT           # for log purposes only
                            kapC_, kapT_ = kapC, kapT 
                            # end loop inner iteration
                        if i_>=ni_-1:                                       # no convergence in inner loop
                            if k_>=kList[-1]-1:
                                print(f'-- ConFemMat::Lubliner:Sig: no convergence B, {Elem.Label:d},{ipI:d},j{j_:d},k{k_:d},i{i_:d},{kapC:11.4e},{dkapC:11.4e},{kapT:11.4e},{dkapT:11.4e}, res {", ".join(zzz for zzz in [f"{x:11.4e}" for x in XXX[-5:-1]])}')
                                print(f'-- ConFemMat::Lubliner:Sig: no convergence B, {Elem.Label:d},{ipI:d},j{j_:d},k{k_:d},i{i_:d},{kapC:11.4e},{dkapC:11.4e},{kapT:11.4e},{dkapT:11.4e}, res {", ".join(zzz for zzz in [f"{x:11.4e}" for x in XXX[-5:-1]])}', file=ff)
                            NotReady = True
                            break                                           # divergent break intermediate loop
                        else:                                               # next k_
                            kapC_old = kapC
                            kapT_old = kapT
                            dEpsP[:]  = (dEpsP+dEpsP_)[:]                   # for next plane stress iteration
                            EpsP_old[:]= (EpsP_old+dEpsP_)[:]               # update of old plastic stress with sub-increment
                        # convergent end intermediate loop over LdEps-increments
           
                    if NotReady:                                            # next nk_
                        pass                                                # further restores for next iteration of outer loop
                        continue                                            # iteration with next LdEps-increment due to previous divergence
                    else: 
                        break                                               # convergent break variable LdEps-increments loops
                    # end loop variable LdEps-increments
                 
                # general updates
                DD = 1.-(1.-DC)*(1.-DT)
                DM = array([ (1.-DT)*GC*h0 , 0. , (1.-DC)*GT*h2, 0., 0., 0. ]) #
                LSigEff  = dot(TTS,SigEff)
                X = array([[LSigEff[0]*DM[0],LSigEff[0]*DM[1],LSigEff[0]*DM[2],LSigEff[0]*DM[3],LSigEff[0]*DM[4],LSigEff[0]*DM[5]],
                           [LSigEff[1]*DM[0],LSigEff[1]*DM[1],LSigEff[1]*DM[2],LSigEff[1]*DM[3],LSigEff[1]*DM[4],LSigEff[1]*DM[5]],
                           [LSigEff[2]*DM[0],LSigEff[2]*DM[1],LSigEff[2]*DM[2],LSigEff[2]*DM[3],LSigEff[2]*DM[4],LSigEff[2]*DM[5]],
                           [LSigEff[3]*DM[0],LSigEff[3]*DM[1],LSigEff[3]*DM[2],LSigEff[3]*DM[3],LSigEff[3]*DM[4],LSigEff[3]*DM[5]],
                           [LSigEff[4]*DM[0],LSigEff[4]*DM[1],LSigEff[4]*DM[2],LSigEff[4]*DM[3],LSigEff[4]*DM[4],LSigEff[4]*DM[5]],
                           [LSigEff[5]*DM[0],LSigEff[5]*DM[1],LSigEff[5]*DM[2],LSigEff[5]*DM[3],LSigEff[5]*DM[4],LSigEff[5]*DM[5]]])
                Y = dot(X,EP)
                if not PlStFlag: break
#                if abs(SigEff[2])<1.e-2: break
                if abs(SigEff[2])<1.e-3: break
            if j_>=np_-1: 
                print(f'-- CaeFemMaterials::Lubliner:MatC: no convergence C {Elem.Label:d} {ipI:d} - {j_:d} {k_:d} {i_:d} - {", ".join(z for z in [f"{x:8.4f}" for x in SigEff[0:6]])}')
                print(f'-- CaeFemMaterials::Lubliner:MatC: no convergence C {Elem.Label:d} {ipI:d} - {j_:d} {k_:d} {i_:d} - {", ".join(z for z in [f"{x:8.4f}" for x in SigEff[0:6]])}', file=ff)
            
        Elem.StateVarN[ipI][12] = Eps[2]
        if not PFlag and DD>ZeroD and norm(Dps)<=1.0e-3*ZeroD: # to reconstruct equilibrium control of last loading step in case of elasticity - no plastic increment
            if PlStFlag: Eps[2]=Eps2
            SigEff = dot( self.C0,(Eps - EpsP_old) )
        if PFlag: 
            CP  = (1.-DD)*(self.C0 - dot(self.C0,EP)) - Y                   # tangential material stiffness in local system
            MatM = dot(transpose(TTD),dot(CP,TTD))                          # transform local stiffness into global system
            Elem.StateVarN[ipI][0]  = kapT
            Elem.StateVarN[ipI][1]  = kapC
            Elem.StateVarN[ipI][2]  = fC - self.f_c0
            Elem.StateVarN[ipI][3]  = fT - self.f_t0
            Elem.StateVarN[ipI][4]  = DC
            Elem.StateVarN[ipI][5]  = DT
            Elem.StateVarN[ipI][6]  = EpsP[0]                               # update of plastic strains
            Elem.StateVarN[ipI][7]  = EpsP[1] 
            Elem.StateVarN[ipI][8]  = EpsP[2] 
            Elem.StateVarN[ipI][9]  = EpsP[3] 
            Elem.StateVarN[ipI][10] = EpsP[4] 
            Elem.StateVarN[ipI][11] = EpsP[5] 
#            Elem.StateVarN[ipI][13] = LdEpsPPr[0]                           # principal plastic strain increment
#            Elem.StateVarN[ipI][14] = LdEpsPPr[1]                           # principal plastic strain increment
#            Elem.StateVarN[ipI][15] = LdEpsPPr[2]                           # principal plastic strain increment
            Elem.StateVarN[ipI][16] = 1                                     # -> PFlag for update procedure
        else:    
            MatM = (1.-DD)*self.C0
        Sig = (1.-DD)*SigEff
        #
#        if Elem.Label in [421] and ipI in [0]:
#            if PFlag: print(f'AAA, {Elem.Label:d},{ipI:d},{Elem.StateVar[ipI][12]:12.5e},{Eps2:12.5e},{Elem.StateVarN[ipI][12]:12.5e}x{", ".join(zzz for zzz in [f"{x:12.5e}" for x in Eps[0:6]])}',file=ff)
#            else:     print(f'BBB, {Elem.Label:d},{ipI:d},{Elem.StateVar[ipI][12]:12.5e},{Eps2:12.5e},{Elem.StateVarN[ipI][12]:12.5e}x{", ".join(zzz for zzz in [f"{x:12.5e}" for x in Eps[0:6]])}',file=ff) 
#            if PFlag: print(f'AAA, {Elem.Label:d},{ipI:d},{Eps2:12.5e},{", ".join(zzz for zzz in [f"{x:11.7f}" for x in Sig[0:6]])},x\
#{", ".join(zzz for zzz in [f"{x:12.5e}" for x in Eps[0:6]])},x\
#{", ".join(zzz for zzz in [f"{x:12.5e}" for x in EpsP[0:6]])}', file=ff)
#            else:     print(f'BBB, {Elem.Label:d},{ipI:d},{Eps2:12.5e},{", ".join(zzz for zzz in [f"{x:11.7f}" for x in Sig[0:6]])},x\
#{", ".join(zzz for zzz in [f"{x:12.5e}" for x in Eps[0:6]])},x\
#{", ".join(zzz for zzz in [f"{x:12.5e}" for x in EpsP_old[0:6]])}', file=ff)
#        
        #
        if Elem.dim==2:                                                             # plane stress / strain
            if Elem.PlSt: 
                cD = 1./MatM[2,2]
                MatM_=array([[MatM[0,0]-MatM[0,2]*MatM[2,0]*cD, MatM[0,1]-MatM[0,2]*MatM[2,1]*cD, MatM[0,5]-MatM[0,2]*MatM[2,5]*cD],
                             [MatM[1,0]-MatM[1,2]*MatM[2,0]*cD, MatM[1,1]-MatM[1,2]*MatM[2,1]*cD, MatM[1,5]-MatM[1,2]*MatM[2,5]*cD], 
                             [MatM[5,0]-MatM[5,2]*MatM[2,0]*cD, MatM[5,1]-MatM[5,2]*MatM[2,1]*cD, MatM[5,5]-MatM[5,2]*MatM[2,5]*cD]])
            else:         
                MatM_=array([[MatM[0,0],MatM[0,1],MatM[0,5]],[MatM[1,0],MatM[1,1],MatM[1,5]],[MatM[5,0],MatM[5,1],MatM[5,5]]])
            return array([Sig[0],Sig[1],Sig[5]]), MatM_, [Eps[0],Eps[1],Eps[2],Eps[5], Sig[0],Sig[1],Sig[2],Sig[5]] # data
        else: 
            sig = array([Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5]])
            return sig, MatM, [Sig[0],Sig[1],Sig[2],Sig[3],Sig[4],Sig[5], 0., 0.]   # data
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):    # loop over integration points 
            if Elem.StateVarN[j][16] > 0:                                   # -> PFlag for update procedure
                if Elem.StateVarN[j][0] > Elem.StateVar[j][0]:
#                    
#                    if j==0: print('YYY')
#                    
                    Elem.StateVar[j][0] = Elem.StateVarN[j][0]              # kapT
                    Elem.StateVar[j][3] = Elem.StateVarN[j][3]              # fT
                    Elem.StateVar[j][5] = Elem.StateVarN[j][5]              # DT
                    Elem.StateVar[j][6] = Elem.StateVarN[j][6]              # plastic strain
                    Elem.StateVar[j][7] = Elem.StateVarN[j][7]
                    Elem.StateVar[j][8] = Elem.StateVarN[j][8]
                    Elem.StateVar[j][9] = Elem.StateVarN[j][9]
                    Elem.StateVar[j][10]= Elem.StateVarN[j][10]
                    Elem.StateVar[j][11]= Elem.StateVarN[j][11]
#                    Elem.StateVar[j][13]= Elem.StateVarN[j][13]             # latest plastic principal strain increment
#                    Elem.StateVar[j][14]= Elem.StateVarN[j][14]
#                    Elem.StateVar[j][15]= Elem.StateVarN[j][15]
                if Elem.StateVarN[j][1] > Elem.StateVar[j][1]:
#                    
#                    if j==0: print('ZZZ',Elem.StateVar[j][1],Elem.StateVarN[j][1])
#                    
                    Elem.StateVar[j][1] = Elem.StateVarN[j][1]              # kapC
                    Elem.StateVar[j][2] = Elem.StateVarN[j][2]              # fC
                    Elem.StateVar[j][4] = Elem.StateVarN[j][4]              # DC
                    Elem.StateVar[j][6] = Elem.StateVarN[j][6]              # plastic strain
                    Elem.StateVar[j][7] = Elem.StateVarN[j][7]
                    Elem.StateVar[j][8] = Elem.StateVarN[j][8]
                    Elem.StateVar[j][9] = Elem.StateVarN[j][9]
                    Elem.StateVar[j][10]= Elem.StateVarN[j][10]
                    Elem.StateVar[j][11]= Elem.StateVarN[j][11]
#                    Elem.StateVar[j][13]= Elem.StateVarN[j][13]             # latest plastic principal strain increment
#                    Elem.StateVar[j][14]= Elem.StateVarN[j][14]
#                    Elem.StateVar[j][15]= Elem.StateVarN[j][15]
            Elem.StateVar[j][12]    = Elem.StateVarN[j][12]                 # lateral strain for plane stress
        return False

class RCBeam(Material):                                             # Reinforced concrete beam cross section
    def __init__(self, PropMat):
                                                                    # concrete data in PropMat[0], reinforcement data in PropMat[1]
        if PropMat[0][4]<1.e-6: raise NameError ("ConFemMaterials::RCBEAM.__init__: non zero tensile strength has to be defined for tension stiffening")
        self.PropMat = PropMat                                      # [0] concrete (see below, [0][4] tensile strength-tension stiff only, [0][5] cross section integration), [1][] reinforcement
        self.NData = 8                                              # number of data items in Elem.Data / DataP
        self.NullD = 1.*10**(-6)
        self.Reinf = MisesUniaxial( PropMat[1], [PropMat[0][4],1.3] )# initialization material for reinforcement + tension stiffening parameters (tensile strength, factor for stabilized cracking)
        self.ReinfTC = None
        self.StateVar = 5                                           # number of state variables per integration point per concrete/reinforcement layer
                                                                    # 0: current permanent strain upon unloading, 1: current yield stress, 2: current strain (just obsolete)
                                                                    # 2: current reference strain for uniaxial stress strain curve, 3: final stress of last time step, see also mises for reinforcement
                                                                    # StateVar of Mises should be the same
        self.Update = True                                          # has specific update procedure for update of state variables
        self.Updat2 = False                                         # no two-stage update procedure
        self.alphaT = PropMat[0][6]                                 # thermal expansion coefficient
        self.phi    = PropMat[0][7]                                 # creep number
        self.zeta   = PropMat[0][8]                                 # viscosity parameter
        self.psi = (1+self.phi)*self.zeta                           # auxiliary value
        self.Density = PropMat[0][9]                                # specific mass
        if self.PropMat[0][5]>0 : self.Symmetric = True             # flag for symmetry of material matrices, PropMat[0][5] indicates number of integration over cross section, 
                                                                    # value 0 (no numerical integration) might be creeping which probably is not symmetric
        if len(PropMat[0])==11:
            self.fct = PropMat[0][10]                               # introduce concrete tensile strength for internal forces by back door 
            self.epsRef = self.fct/PropMat[0][0]                    # reference strain for tensile stress-strain relation --> MatUniaaxCoTens
            self.epsFctCoeff = 1.5 # 2.            
            self.epsLimCoeff = 3. # 3. # 5. # 6.        
            self.epsLim = self.epsLimCoeff*self.epsRef              # strain to zero
            self.epsFct = self.epsFctCoeff*self.epsRef              # strain to tensile strength
        else: 
            self.epsLim = 0.
        self.Type = 'RCBEAM'
    
    def MatUniaxCo( self, PropM, eps):                              # uniaxial concrete (see DIN1045-1 Sec. 9.1.5)
        Ec0     = PropM[0]                                          # Concrete Young's modulus
        fc      = PropM[1]                                          # compressive strength
        epsc1   = PropM[2]                                          # strain at compressive strength
        epsc1u  = PropM[3]                                          # ultimate compressive strain
        Case = 0
        if 0. < eps <= self.epsLim:
            sig, dsig = self.MatUniaxCoTens( eps, Ec0)
#            sig = Ec0*eps
#            dsig = Ec0
        elif epsc1u<=eps and eps<=0.:
            if Case == 0:
                eta = eps/epsc1                                     #DIN1045-1 Eq. (63)
                k = -Ec0*epsc1/fc                                   #DIN1045-1 Eq. (64)
                sig = -fc * (k*eta-eta**2) / (1.+(k-2)*eta)         #DIN1045-1 Eq. (62)
                dsig = -fc * ( (k-2*eta)/(1+(k-2)*eta) - (k*eta-eta**2)/(1+(k-2)*eta)**2*(k-2) ) / epsc1
            elif Case == 1:                                         # spring spline
                s1 = -epsc1
                F1 = fc
                s2 = -epsc1u
                F2 = 1.0*fc
                s0 = 0.
                d0 = Ec0
                sig, dsig = SpringSpline( s1, F1, s2, F2, s0, d0, eps)
            elif Case == 2:                                         # parabola square
                fcd = fc
                eps_c2 = epsc1
                eps_c3 = epsc1u
                if   eps_c2<=eps:
                    xyz = 1.-(eps)/eps_c2
                    sig =  -fcd*(1.-xyz*xyz)
                    dsig=  -2*fcd*(1+eps/eps_c2)/eps_c2
                elif eps_c3<=eps:
                    sig = -fcd
                    dsig = 0.
            else:
                raise NameError("no uniaxial concrete law") 
        else:
            sig = 0
            dsig = 0
        return sig, dsig
    def MatUniaxCoTens(self, s, E0):
        s1 = self.epsFctCoeff*self.epsRef 
        F1 = self.fct
        s2 = self.epsLimCoeff*self.epsRef 
        F2 = 0.
        s0 = 0. 
        d0 = E0
        ssig, dsig = SpringSpline( s1, F1, s2, F2, s0, d0, s)
        return ssig, dsig
    def IntegrateConc( self, Elem, Eps, Tr, Tg, z1, z2, nI):
        dz = (z2-z1)/nI                                         # z1 local lower, z2 local upper; cross sectional height / number of integration points
        nn = zeros((nI+1), dtype=double)
        mm = zeros((nI+1), dtype=double)
        nde= zeros((nI+1), dtype=double)
        ndk= zeros((nI+1), dtype=double)
        mdk= zeros((nI+1), dtype=double)
        bb = Elem.Geom[1,1]                                     # width of cross section, may eventually be overridden
        z = z1                                                  # cross section / strain coordinate, initial value
        for i in range(nI+1):                                   # numerical integration of concrete contribution
            if Elem.CrossSecType=="POLYLINE": bb = Elem.WidthOfPolyline(z) # width of cross section defined by polyline
            ept  = self.alphaT*(Tr - z*Tg)                      # temperature strain in concrete fiber
            eps  = Eps[0]-z*Eps[1]
            sig, dsig =self.MatUniaxCo( self.PropMat[0], eps-ept) # current fiber stress
            sig_ = bb*sig
            dsig_= bb*dsig
            nn[i]  =    sig_                                    # normal force
            mm[i]  = -z*sig_                                    # moment
            nde[i] =       dsig_                                # normal force grad eps
            ndk[i] = -z   *dsig_                                # normal force grad kappa
            mdk[i] =  z**2*dsig_                                # moment grad kappa
            z = z + dz
        NN  = integrate.trapz( nn, x=None, dx=dz)           # numerical integration normal force
        MM  = integrate.trapz( mm, x=None, dx=dz)           # numerical integration moment
        NDE = integrate.trapz( nde, x=None, dx=dz)        # numerical normal force grad eps
        NDK1= integrate.trapz( ndk, x=None, dx=dz)        # numerical normal force grad kappa
        MDK = integrate.trapz( mdk, x=None, dx=dz)        # numerical moment grad kappa
        return NN, MM, NDE, NDK1, MDK
#    @profile
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        self.PropMat[0][1]=RandomField_Routines().get_property(Elem.Set,Elem.Label,self.PropMat[0][1])[0]
        if CalcType == 0: return [], [], []
        epsc = Eps[0]                                               # strain on compressive side -- initializations
        epss = Eps[0]                                               # maximum reinforcement strain
        NN   = 0.                                                   # resulting normal force
        MM   = 0.                                                   # resulting moment
        NDE  = 0.                                                   # normal force gradient strain
        NDK  = 0.                                                   # normal force gradient curvature
        MDK  = 0.
        NDK1, NDK2 = 0., 0.                                                   # moment gradient curvature
        if Elem.CrossSecType=="POLYLINE": zc1, zc3, hh = Elem.Geom[1,1], Elem.Geom[1,2], Elem.Geom[1,2]-Elem.Geom[1,1]
        else:                     zc1, zc3, hh = -Elem.Geom[1,2]/2, Elem.Geom[1,2]/2, Elem.Geom[1,2] # coordinate of bottom fibre, top fibre, height of cross section    
        Tr   = 0.5*(Temp[0]+Temp[1])                                # temperature of reference / middle axis
        Tg   = (Temp[0]-Temp[1])/hh                                 # temperature gradient
        dTr  = 0.5*(dTmp[0]+dTmp[1])                                # temperature increment of reference / middle axis
        dTg  = (dTmp[0]-dTmp[1])/hh                                 # temperature increment gradient
        # reinforcement contributions
        nR =Elem.Geom.shape[0]-2                                    # number of reinforcement layers
        for i in range(nR):                                         # reinforcement contribution
            ipL = (1+nR)*ipI+(1+i)                                  # every IP has state variables for concrete and nR reinforcement layers
            As = Elem.Geom[i+2,0]                                   # reinforcement cross section
            ys = Elem.Geom[i+2,1]                                   # coordinates of reinforcement
            eps = Eps[0] - ys*Eps[1]                                # strain in reinforcement
            dep = Dps[0] - ys*Dps[1]                                # strain increment in reinforcement
            tmp = Tr - ys*Tg                                        # temperature in reinforcement
            dtp = Tr - ys*Tg                                        # temperature increment in reinforcement
            if Elem.RType[i]=='TC':
                if len(Eps)==4: eps = (Eps[0]-Eps[2]) - ys*(Eps[1]-Eps[3])   # for material tester
                Sig, Dsig, _ = self.ReinfTC.Sig( ff, CalcType, Dt, i, ipL, Elem, [dep, 0], [eps, 0], dtp, tmp, None)
            else:                   
                Sig, Dsig, _ =   self.Reinf.Sig( ff, CalcType, Dt, i, ipL, Elem, [dep, 0], [eps, 0], dtp, tmp, None) # -> MisesUniaxial
            NN = NN   + Sig[0] *As                                  # reinforcement contribution to normal force
            MM = MM   - Sig[0] *As*ys                               # reinforcement contribution to moment
            NDE = NDE + Dsig[0,0]*As                                # reinf. contrib. to normal force grad eps
            NDK = NDK - Dsig[0,0]*As*ys                             # reinf. contrib. to normal force grad kappa
            MDK = MDK + Dsig[0,0]*As*ys**2                          # reinf. contrib. to moment grad kappa
            Elem.StateVarN[ipL,2] = eps                             # 
            if eps>epss: epss=eps                                   # maximum reinforcement strain

        BZ = zeros((2,2))                                           # for linear compressive concrete without integration
        if fabs(Eps[1])>self.NullD:                                 # contributing concrete part
            z0  = (Eps[0])/Eps[1]                                   # zero line
            z0_ = (Eps[0]-self.epsLim)/Eps[1]                       # 'zero' line - with contribution of concrete tensile strength
            if Eps[1]>0:                                            # UPward curvature
                epsc = Eps[0]-zc3*Eps[1]                            # minimum compressive strain
                if z0_ > zc3:                                       # predominant tension
                    zc2 = zc3
                    zc1 = zc3                                       # no concrete zone
                elif z0_ >zc1:                                      # predominant bending
                    zc2 = min(z0, zc3)                              # concrete zero line coordinate, should be z1_ >= zc1
                    zc1 = min(z0_,zc3)                              # concrete zone lower coordinate
                    assert zc1 <= zc2
                    BZ=array([[1/Eps[1],-Eps[0]/(Eps[1]**2)],[0,0]])
            else:                                                   # DOWNward curvature
                epsc = Eps[0]-zc1*Eps[1]                            # minimum compressive strain
                if z0_ < zc1:                                       # predominant tension
                    zc3 = zc1                                       # no concrete zone
                    zc2 = zc1
                elif z0_ < zc3:                                     # predominant bending 
                    zc3 = max(z0_,zc1)                              # concrete zone upper coordinate
                    zc2 = max(z0, zc1)                              # concrete zero line coordinate, should be z2_ <= zc3
                    assert zc2 <= zc3
                    BZ=array([[0,0],[1/Eps[1],-Eps[0]/(Eps[1]**2)]])
            
        if epsc<self.PropMat[0][3]:
            Echo(f"el {Elem.Label:d}, el-ind {elI} ip {ipI:d} ConFemMaterials::RCBeam.Sig: allowed concrete strain {self.PropMat[0][3]:f} exceeded with {epsc:f}", ff)
            print((Elem.Label,elI,ipI,'__',self.PropMat[0][3],epsc,"ConFemMaterials::RCBeam.Sig: allowed concrete stress exceeded"))

        Emod= self.PropMat[0][0]
        nI  = int(self.PropMat[0][5])                               # number of integration points (for integration over cross section)
        ipL = (1+nR)*ipI                                            # every IP has state variables for concrete and nR reinforcement layers
        Elem.StateVarN[ipL,2] = epsc                                # currently used by reference strain for uniaxial smoothed stress strain curve
        # concrete contribution
        if nI>0:
            xx = zc3 - zc1                                          # concrete zone height
            if xx > ZeroD:                                          # nonlinear behavior in compressive range
                NN_, MM_, NDE_, NDK1_, MDK_ = self.IntegrateConc( Elem, Eps, Tr, Tg, zc1, zc3, nI )
                NN  += NN_
                MM  += MM_
                NDE += NDE_
                MDK += MDK_
                NDK1 = NDK + NDK1_
                NDK2 = NDK1
        else:                                                       # linear behavior in compressive zone, features creep
            psiI= 1/(1+Dt*self.psi)
            eps1 = Eps[0]-zc1*Eps[1] - self.alphaT*(Tr-zc1*Tg)        # effective strain in lower concrete fiber
            eps2 = Eps[0]-zc3*Eps[1] - self.alphaT*(Tr-zc3*Tg)        # effective strain in upper concrete fiber
            dep1 = Dps[0]-zc1*Dps[1] - self.alphaT*(dTr-zc1*dTg)      # effective strain increment in lower concrete fiber
            dep2 = Dps[0]-zc3*Dps[1] - self.alphaT*(dTr-zc3*dTg)      # effective strain increment in upper concrete fiber
            if abs(eps1)>self.NullD: sic1 = psiI*( Elem.StateVar[ipL,0] + Emod*dep1 + Dt*self.zeta*Emod*(eps1) )# implicit Euler integration !!! zeta inverse to zeta in book
            else: sic1 = 0.
            if abs(eps2)>self.NullD: sic2 = psiI*( Elem.StateVar[ipL,1] + Emod*dep2 + Dt*self.zeta*Emod*(eps2) )# implicit Euler integration !!!
            else: sic2 = 0.
            if Elem.CrossSecType=="POLYLINE":
                AA, SS, JJ, hh = Elem.CrossSecGeom(zc1,zc3)           # area, statical moment, inertial moment
                AS = 1./(zc3-zc1)*array([[AA*zc3-SS,-AA*zc1+SS],[-SS*zc3+JJ,SS*zc1-JJ]]) # A_sigma in book
                NN = NN + AS[0,0]*sic1 + AS[0,1]*sic2
                MM = MM + AS[1,0]*sic1 + AS[1,1]*sic2
                AZ = Elem.CrossSecGeomA(zc1,zc3,sic1,sic2)
            else:
                bb = Elem.Geom[1,1]                                 # width of cross section
                AS = array([[bb*(zc3-zc1)/2,bb*(zc3-zc1)/2],[-bb*(zc3+2*zc1)*(zc3-zc1)/6,-bb*(2*zc3+zc1)*(zc3-zc1)/6]]) # A_sigma in book
                NN =  NN + bb*(zc3-zc1)*(sic1+sic2)/2
                MM =  MM - bb*( (zc3+2*zc1)*(zc3-zc1)*sic1 + (2*zc3+zc1)*(zc3-zc1)*sic2 )/6
                AZ = array([[-bb*(sic1+sic2)/2,bb*(sic1+sic2)/2],[bb*( (4*zc1-zc3)*sic1 + (2*zc1+zc3)*sic2 )/6, -bb*( (zc1+2*zc3)*sic1 + (-zc1+4*zc3)*sic2 )/6]])
            BE = array([[1,-zc1],[1,-zc3]])
            SE = Emod*dot(AS,BE)                                    # material tangential stiffness 
            SZ = dot(AZ,BZ)                                         # geometrical tangential stiffness
            NDE = NDE + SE[0,0] + SZ[0,0]
            NDK1= NDK + SE[0,1] + SZ[0,1]
            NDK2= NDK + SE[1,0] + SZ[1,0]
            MDK = MDK + SE[1,1] + SZ[1,1]
            Elem.StateVarN[ipL,0] = sic1                            # update state variables
            Elem.StateVarN[ipL,1] = sic2                            #
        if Elem.dim==10:                                            # Bernoulli beam
            MatM = array( [[NDE,NDK1],[NDK2,MDK]] )                 # material tangential stiffness (Script Eq. (3.21))
            sig = array( [NN,MM] )                                  # stress vector
            return sig, MatM, [Eps[0],Eps[1],sig[0],sig[1],1000*epss,1000*epsc]
        elif Elem.dim==11:                                          # Timoshenko beam
            if Elem.CrossSecType=="POLYLINE": bb = Elem.WidthOfPolyline(0.5*(zc3-zc1)) # width of cross section defined by polyline
            else:                             bb = Elem.Geom[1,1]           # width of cross section, may eventually be overridden
            MatM = array( [[NDE,NDK1,0],[NDK2,MDK,0],[0,0,bb*hh*0.9*Emod/4]] )
            sig = array( [NN,MM,MatM[2,2]*Eps[2]] )                 # stress vector
            return sig, MatM, [Eps[0],Eps[1],Eps[2],sig[0],sig[1],sig[2],1000*epss,1000*epsc]
    def UpdateStateVar(self, Elem, ff):
        nR =Elem.Geom.shape[0]-2                                    # number of reinforcement layers
        RList = [ (1+nR)*j+(1+i) for j in range(Elem.nIntL) for i in range(nR)] # list of integration point indices for reinforcement
        for j in range(Elem.StateVar.shape[0]):                     # loop over all concrete and reinforcement integration point contributions
            if j in RList:                                          # reinforcement integration contribution
                if Elem.StateVarN[j,1]>Elem.StateVar[j,1]:          # plastic loading - unsigned highest stress increased - in case of yielding
                    Elem.StateVar[j,0] = Elem.StateVarN[j,0]        # new permanent strain
                    Elem.StateVar[j,1] = Elem.StateVarN[j,1]        # new highest stress (unsigned) - in case of yielding
                Elem.StateVar[j,2] = Elem.StateVarN[j,2]            # strain of fiber
                Elem.StateVar[j,3] = Elem.StateVarN[j,3]            # for smoothed uniaxial version of uniaxial mises for reinforcement -> activation of sigEpsL
                Elem.StateVar[j,4] = Elem.StateVarN[j,4]            # "
            else: 
                Elem.StateVar[j] = Elem.StateVarN[j]                # concrete integration contribution
        return False

class TextileR(Material):
    def __init__(self, PropMat):
#    def __init__(self,         SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,      False,     1,           2, "TexR")
#        self.Symmetric = True                       # flag for symmetry of material matrices
#        self.Update = True                          # has specific update procedure for update of state variables
#        self.Updat2 = False                         # no 2 stage update procedure
#        self.StateVar = 1                           # number of state variables per integration point
#        self.NData = 2                              # number of data items
        self.PropMat = PropMat
#        self.StateVar = None
#        self.NData = 2                              # number of data items
        self.alphaT = 0                             # thermal expansion coefficient
        self.Density = 0                            # specific mass
        self.Emod = PropMat[0]                      # Young's modulus
        self.eps_tund = PropMat[1]                  # PropMat[1] is strain 
        self.sig_tund = self.eps_tund*self.Emod     
        self.sig_tu = PropMat[2]                    # failure stress
        self.eps_tu = PropMat[3]                    # failure strain
        if self.eps_tu==self.eps_tund:
            self.linearity=1
            self.Etan = self.Emod 
        else:
            self.linearity=2
            self.Etan = (self.sig_tu-self.sig_tund)/(self.eps_tu-self.eps_tund) 
        
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
#        self.sig_tu=RandomField_Routines().get_property(Elem.Set,Elem.Label,self.sig_tu)[-1] #sig_tu is overwritten by the random value
#                                                                                            #chose of the last value!!!
#        if self.eps_tu==self.eps_tund:
#            self.linearity=1
#            self.Etan = self.Emod 
#        else:
#            self.linearity=2
#            self.Etan = (self.sig_tu-self.sig_tund)/(self.eps_tu-self.eps_tund) 
   
        if Eps[0]<=self.eps_tund:
            MatM = array([[self.Emod,0],[0,0]])          # material stiffness uniaxial
            sig = [self.Emod*Eps[0],0.]
        elif Eps[0]<=self.eps_tu and self.linearity==2:
            MatM = array([[self.Etan,0],[0,0]]) 
            sig  = [(Eps[0]-self.eps_tund)*self.Etan+self.sig_tund,0]
        else: 
            MatM = array([[self.Etan,0],[0,0]])          # material stiffness uniaxial
            sig  = [(Eps[0]-self.eps_tund)*self.Etan+self.sig_tund,0] 
        Elem.StateVarN[ipI][0] = sig[0] 
        return sig, MatM, [Eps[0], sig[0], 0.]          # last value to make this consistent with WraRCShell

class WraTCBeam(RCBeam):                                # Textile Reinforced concrete beam
    def __init__(self, PropMat):
        RCBeam.__init__(self, [PropMat[0], PropMat[1]])
        self.ReinfTC = TextileR( PropMat[2] )           # used as additional option by RCBeam

class WraRCShell(Material):                                                 # Wrapper for reinforced shell
    def __init__(self, PropMat, MatPointer, ReMatType, ff):
        if isinstance( MatPointer, ElasticLT):
            self.Conc = ElasticLT(PropMat[0])
            self.ConcS = 'ELASTICLT' 
            fct = PropMat[0][2]                                             # for tension stiffening
        elif isinstance( MatPointer, IsoDamage):
            self.Conc = IsoDamage(PropMat[0], False, False, 0., False, 0.)
            self.ConcS = 'ISODAMAGE' 
            fct = self.Conc.fct                                             # for tension stiffening
        elif isinstance( MatPointer, Elastic):
            self.Conc = Elastic(PropMat[0])
            self.ConcS = 'ELASTIC' 
            self.Conc.NData = 8                                             # to make this consistent with ElasticLT and IsoDamage
            fct = 0.
        elif isinstance( MatPointer, MicroPlaneDam):
            self.Conc = MicroPlaneDam(PropMat[0], False, False, 0., False, 0., ff)
            self.ConcS = "MicroPl"  
            self.Conc.NData = 8                                             # to make this consistent with ElasticLT and IsoDamage
            fct = 0.
        elif isinstance( MatPointer, MicroPlaneDam):
            self.Conc = MicroPlaneDam(PropMat[0], False, False, 0., False, 0., ff)
            self.ConcS = 'MICRODAMAGE' 
            self.Conc.NData = 8                                             # to make this consistent with ElasticLT and IsoDamage
            fct = 0.
        else: raise NameError ("ConFemMaterials::WraRCShell.__init__: material type not regarded")
        self.StateVar  = self.Conc.StateVar                                 # number of state variables per integration point
        self.NData     = self.Conc.NData                                    # number of data items, subject to change in ConFemInOut::ReadInputFile for input of SH4.
        self.Update    = self.Conc.Update
        self.Updat2    = self.Conc.Updat2
        self.RType     = self.Conc.RType                                    # regularization type for concrete
        self.ReinfType = None                                               # reinforcement type
        self.Density   = PropMat[0][-1]                                     # density of concrete component
        if self.RType==2:
            self.bw    = self.Conc.bw
            self.CrX   = self.Conc.CrX
            self.CrY   = self.Conc.CrY
            self.CrX2  = self.Conc.CrX2
            self.CrY2  = self.Conc.CrY2
            self.CrBwN = self.Conc.CrBwN                                    # crack band regularization: number of intervals for scaling factor interpolation
            self.ElCharLenBound = self.Conc.ElCharLenBound
        if ReMatType == "TR":                                               # leads to TextileR, but input syntax slightly different with corresponding input syntax of TCBEAM
            E_0, sig_1, E_1, sig_u = PropMat[1][0], PropMat[1][2], PropMat[1][3], PropMat[1][4]
            eps_1 = sig_1/E_0
            eps_u = eps_1 + (sig_u-sig_1)/E_1  
            self.Reinf  = TextileR( [ E_0 , eps_1 , sig_u , eps_u ] )       # to make this compatible with init of TextileR
            self.ReinfType = 'TR'
        else:                 
            self.Reinf  = MisesUniaxial( PropMat[1], [fct,1.3])             # initialization material for reinforcement + tension stiffening parameters
            self.ReinfType = 'MISES'
        self.Type = 'RCSHELL'
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if ipI<Elem.nIntLi:                                                 # initial value of number of bulk integration points before extension for reinforcement 
            SigL, MatL, dataL = self.Conc.Sig( ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, []) # 4/5 integration points over cross section height
        else:                                                               # reinforcement contributions -- every reinforcement sheet should have its own IpI
            if CalcType == 0: return [], [], []
            As = SampleWeightRCShell[Elem.Set,Elem.IntT,Elem.nInt-1,ipI]
            if As>ZeroD:
                Elem.dim=1                                                  # for uniaxial reinforcement constitutive law
                Dps_ = array([Dps[0],Dps[1],Dps[5]])
                Eps_ = array([Eps[0],Eps[1],Eps[5]])
                if   Elem.nInt==2: j = (ipI-16)//4                          # floor division -> integer value -> index for reinforcement layer, 4 stands for number of integration points in base area
                elif Elem.nInt>=5: j = (ipI-20)//4                          # " corresponding
                phi = Elem.Geom[2+j,4]*pi/180.                              # angle to direction of reinforcement
                Trans=array([[cos(phi)**2,sin(phi)**2, cos(phi)*sin(phi)],
                             [sin(phi)**2,cos(phi)**2,-cos(phi)*sin(phi)],
                             [-2*cos(phi)*sin(phi),2*cos(phi)*sin(phi),cos(phi)**2-sin(phi)**2]])# transformation matrix for strains from DIRECTOR to directed
                MStiff=array([[cos(phi)**4,cos(phi)**2*sin(phi)**2,cos(phi)**3*sin(phi)],
                              [cos(phi)**2*sin(phi)**2,sin(phi)**4,cos(phi)*sin(phi)**3],
                              [cos(phi)**3*sin(phi),cos(phi)*sin(phi)**3,cos(phi)**2*sin(phi)**2]])
                dpsL = dot(Trans,Dps_)                                      # directed strain increment
                epsL = dot(Trans,Eps_)                                      # directed strain
                SigL, MatL, dataL = self.Reinf.Sig( ff, CalcType, Dt, elI, ipI, Elem, dpsL, epsL, dTmp, Temp, None)
                MatM = MatL[0,0]*MStiff                                     # anisotropic tangential material stiffness in DIRECTOR system
                sig = dot(Trans.transpose(),[SigL[0],0.,0.])                # rotate back to DIRECTOR system ???
                Elem.dim=21                                                 # undo
                SigL = array([sig[0],sig[1],0.,0.,0.,sig[2]])               # 3 comp --> 6 comp 
                MatL = array([[MatM[0,0],MatM[0,1],0.,0.,0.,MatM[0,2]],     # "
                              [MatM[1,0],MatM[1,1],0.,0.,0.,MatM[1,2]],
                              [0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.],
                              [0.,0.,0.,0.,0.,0.],
                              [MatM[2,0],MatM[2,1],0.,0.,0.,MatM[2,2]]])
            else:
                MatL = zeros((6,6), dtype=float)
                SigL = zeros((6), dtype=float)
                dataL = zeros((3), dtype=float)
            dataL = [SigL[0],SigL[1],SigL[2],SigL[3],SigL[4],SigL[5],dataL[0],dataL[1],dataL[2],0.,0.,0.,0.,0.]# trailing zeros to make this consistent with data format of ElasticLT 2d
        # different dataL for bulk and reinforcement !
        return SigL, MatL, dataL
    def UpdateStateVar(self, Elem, ff):
        Flag = self.Conc.UpdateStateVar(Elem, ff)
        if Elem.Type=='SH4': 
            if   Elem.nInt==2: nn = list(range(16,Elem.StateVar.shape[0])) # for reinforcement type mises/RC only
            elif Elem.nInt==5: nn = list(range(20,Elem.StateVar.shape[0])) # "
            if self.ReinfType == 'MISES':
                for j in nn:                          # loop over integration points for reinforcement only
                    if Elem.StateVarN[j,1]>Elem.StateVar[j,1]: Elem.StateVar[j] = Elem.StateVarN[j]
            if self.ReinfType == "TR":
                for j in nn:
                    if abs(Elem.StateVarN[j,0])>abs(Elem.StateVar[j,0]): Elem.StateVar[j] = Elem.StateVarN[j]
                    if abs(Elem.StateVarN[j,0])>self.Reinf.sig_tu: print('WraRCShell::TR: strength exceeded', Elem.Label, j, Elem.StateVar[j,0], file=ff)  
        return Flag
    def UpdateStat2Var(self, Elem, ff, Flag, LoFl):
        self.Conc.UpdateStat2Var(Elem, ff, Flag, LoFl)
    
class WraMisesReMem(Material):                      # anisotropic reinforcement membrane with elastoplastic Mises
    def __init__(self, PropMat, val):
        self.PropMat = PropMat
        self.Density = PropMat[7]                   # specific mass
        self.NData = 8                              # number of data items
#        self.Reinf = Mises( PropMat, [0.0,1.3] )    # initialization material for reinforcement + tension stiffening parameters
        self.Reinf = MisesUniaxial( PropMat, [0.0,1.3] )    # initialization material for reinforcement + tension stiffening parameters
        self.StateVar = self.Reinf.StateVar
        self.Update = True                          # has specific update procedure for update of state variables
        self.Updat2 = False                         # no 2 stage update procedure
        phi = PropMat[1]*pi/180.                    # angle of orientation of reinforcement
        self.Trans=array([[cos(phi)**2,sin(phi)**2, cos(phi)*sin(phi)],
                          [sin(phi)**2,cos(phi)**2,-cos(phi)*sin(phi)],
                          [-2*cos(phi)*sin(phi),2*cos(phi)*sin(phi),cos(phi)**2-sin(phi)**2]])# transformation matrix for strains from global to local
        self.MStiff=array([[cos(phi)**4,cos(phi)**2*sin(phi)**2,cos(phi)**3*sin(phi)],
                          [cos(phi)**2*sin(phi)**2,sin(phi)**4,cos(phi)*sin(phi)**3],
                          [cos(phi)**3*sin(phi),cos(phi)*sin(phi)**3,cos(phi)**2*sin(phi)**2]])
        self.Type = 'RESHEET'
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
#        Emod = self.PropMat[0]
#        nu = 0.2                                    # Poissons ratio fixed here !
        if Elem.dim==2:                             # plane strain biaxial
            dpsL = dot(self.Trans,Dps)              # local strain increment
            epsL = dot(self.Trans,Eps)              # local strain
            Elem.dim=1
            SigL, MatL, _ = self.Reinf.Sig( ff, CalcType, Dt, elI, ipI, Elem, dpsL, epsL, dTmp, Temp, EpsR)
            Elem.dim=2
            MatM = MatL[0,0]*self.MStiff            # anisotropic local tangential material stiffness
            sig = dot(self.Trans.transpose(),[SigL[0],0.,0.]) # global stress times thickness
#            return sig, MatM, [Eps[0], Eps[1], Eps[2], sig[0], sig[1], sig[2]] # !!!!!!!
            return sig, MatM, [Eps[0], Eps[1], 0., Eps[2], sig[0], sig[1], 0., sig[2]] # !!!!!!!
        else: raise NameError ("ConFemMaterials::WraMisesReMem.Sig: misesReMem does not match element type")
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):# loop over all integration and integration sub points
            if abs(Elem.StateVarN[j,4])>abs(Elem.StateVar[j,4]): Elem.StateVar[j] = Elem.StateVarN[j]
        return False

class NLSlab(Material):                                         # Simple nonlinear Kirchhoff slab
    def __init__(self, PropMat):
        self.posKX  = PropMat[0]                                # pos moment initial bending stiffness x
        self.posKY  = PropMat[1]                                # pos moment initial bending stiffness y
        self.Mx_y   = PropMat[2]                                # pos moment initial yield moment x
        self.My_y   = PropMat[3]                                # pos moment initial yield moment y
        self.posKTX = PropMat[4]                                # pos moment hardening stiffness x
        self.posKTY = PropMat[5]                                # pos moment hardening stiffness y
        self.Mx_t   = PropMat[6]                                # pos moment failure moment x
        self.My_t   = PropMat[7]                                # pos moment failure moment y
        self.alpha  = PropMat[8]                                # factor for twisting stiffness
        if len(PropMat) > 9 :
            self.negKX   = PropMat[9]                           # neg moment initial bending stiffness x
            self.negKY   = PropMat[10]                          # neg moment initial bending stiffness y
            self.negMx_y = PropMat[11]                          # neg moment initial yield moment x - signed
            self.negMy_y = PropMat[12]                          # neg moment initial yield moment y - signed
            self.negKTX  = PropMat[13]                          # neg moment hardening stiffness x
            self.negKTY  = PropMat[14]                          # neg moment hardening stiffness y
            self.negMx_t = PropMat[15]                          # neg moment failure moment x - signed
            self.negMy_t = PropMat[16]                          # neg moment failure moment y - signed
        else:
            self.negKX   =  self.posKX
            self.negKY   =  self.posKY
            self.negMx_y = -self.Mx_y
            self.negMy_y = -self.My_y
            self.negKTX  =  self.posKTX
            self.negKTY  =  self.posKTY
            self.negMx_t = -self.Mx_t
            self.negMy_t = -self.My_t
            
        self.StateVar = 9                                       # pos. moment: 0 plastic curv x, 1 current moment x, 2 yield mom x, 3-5 same y,
                                                                # 6 actual twisting moment
                                                                # neg. moment: 7 yield mom x, 8 same y
        self.NData = 6                                          # number of data items, extended by 2 in SB3 initialization
        self.Update = True                                      # has specific update procedure for update of state variables
        self.Updat2 = False                                     # no 2 stage update procedure
        self.Type = 'NLSLAB'
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        if CalcType == 0: return [], [], []
        if Elem.dim==20:                                        # Kirchhoff slab
            # previous time step values yield moments
            if Elem.StateVar[ipI][2]==0: posMxx_y = self.Mx_y
            else:                        posMxx_y = Elem.StateVar[ipI][2]
            if Elem.StateVar[ipI][5]==0: posMyy_y = self.My_y
            else:                        posMyy_y = Elem.StateVar[ipI][5] 
            if Elem.StateVar[ipI][7]==0: negMxx_y = self.negMx_y
            else:                        negMxx_y = Elem.StateVar[ipI][7]
            if Elem.StateVar[ipI][8]==0: negMyy_y = self.negMy_y
            else:                        negMyy_y = Elem.StateVar[ipI][8] 
            # previous time step values
            EpsP = array([Elem.StateVar[ipI][0], Elem.StateVar[ipI][3],0.]) # permanent curvature at the end of of previous time increment
#            if abs(Eps[0])>ZeroD and abs(Eps[1])>ZeroD: EpsP[2] = Eps[2]*0.5*(EpsP[0]/Eps[0]+EpsP[1]/Eps[1]) # strong assumption
            sigOld = array([Elem.StateVar[ipI][1],Elem.StateVar[ipI][4],Elem.StateVar[ipI][6]]) # moments at the end of of previous time increment
            # actual tangential stiffness
            MatM = zeros((3,3), dtype=float)
            # x-direction
            if   Eps[0]<=(EpsP[0]+negMxx_y/self.negKX+ZeroD): MatM[0,0] = self.negKTX # negative plastic range 
            elif Eps[0]<= EpsP[0]:                            MatM[0,0] = self.negKX  # negative elastic range
            elif Eps[0]<=(EpsP[0]+posMxx_y/self.posKX-ZeroD): MatM[0,0] = self.posKX  # positive elastic range
            else:                                             MatM[0,0] = self.posKTX # positive plastic range               
            # y-direction
            if   Eps[1]<=(EpsP[1]+negMyy_y/self.negKY+ZeroD): MatM[1,1] = self.negKTY # negative plastic range
            elif Eps[1]<= EpsP[1]:                            MatM[1,1] = self.negKY  # negative elastic range
            elif Eps[1]<=(EpsP[1]+posMyy_y/self.posKY-ZeroD): MatM[1,1] = self.posKY  # positive elastic range 
            else:                                             MatM[1,1] = self.posKTY # positive plastic range               
            MatM[2,2] = self.alpha*sqrt(MatM[0,0]*MatM[1,1])    # twisting stiffness
            #
            sig = sigOld + dot(MatM,Dps)
            #
            # update state variables
            Elem.StateVarN[ipI][1],Elem.StateVarN[ipI][4],Elem.StateVarN[ipI][6] = sig[0],sig[1],sig[2] # new stress as state variable
            if   sig[0] > posMxx_y:
                Elem.StateVarN[ipI][0] = Eps[0]-sig[0]/self.posKX   # new permanent curvature 
                Elem.StateVarN[ipI][2] = sig[0]                     # new yield moment
            elif sig[0] < negMxx_y:
                Elem.StateVarN[ipI][0] = Eps[0]-sig[0]/self.negKX   # new permanent curvature 
                Elem.StateVarN[ipI][7] = sig[0]                     # new yield moment
            if   sig[1] > posMyy_y:
                Elem.StateVarN[ipI][3] = Eps[1]-sig[1]/self.posKY   # new permanent curvature 
                Elem.StateVarN[ipI][5] = sig[1]                     # new yield moment
            elif sig[1] < negMyy_y:
                Elem.StateVarN[ipI][3] = Eps[1]-sig[1]/self.negKY   # new permanent curvature 
                Elem.StateVarN[ipI][8] = sig[1]                     # new yield moment
            # control for failure
            if sig[0] > self.Mx_t:    raise NameError("ConFemMat::NlSlab:sig failure pos x El", Elem.Label,ipI,"moment",sig[0],"curv",Eps[0])
            if sig[0] < self.negMx_t: raise NameError("ConFemMat::NlSlab:sig failure neg x El", Elem.Label,ipI,"moment",sig[0],"curv",Eps[0])
            if sig[1] > self.My_t:    raise NameError("ConFemMat::NlSlab:sig failure pos y El", Elem.Label,ipI,"moment",sig[1],"curv",Eps[1])
            if sig[1] < self.negMy_t: raise NameError("ConFemMat::NlSlab:sig failure neg y El", Elem.Label,ipI,"moment",sig[1],"curv",Eps[1])
            #    
            return sig, MatM, [Eps[0],Eps[1],Eps[2],sig[0],sig[1],sig[2]]
#            return sig, MatM, [Eps[0],Eps[1],Eps[2],sig[0],sig[1],sig[2],EpsP[0],EpsP[1],EpsP[2]]
        else: raise NameError ("ConFemMaterials::NLSlab.Sig: not implemented for this element type")
    def UpdateStateVar(self, Elem, ff):
        for j in range(Elem.StateVar.shape[0]):                 # loop over all integration and integration sub points
            Elem.StateVar[j] = Elem.StateVarN[j]
        return False

class Spring(Material):                                             # Nonlinear spring for bond
    def __init__(self, PropMat):
        self.PropMat = PropMat
        self.StateVar = None
        self.NData = 2                                              # number of data items
        self.Type = 'SPRING'
    def Sig(self, ff, CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
        s1 = self.PropMat[0]
        F1 = self.PropMat[1]
        s2 = self.PropMat[2]
        F2 = self.PropMat[3]
        s0 = self.PropMat[4]
        d0 = self.PropMat[5]
        #
        ssig, dsig = SpringSpline( s1, F1, s2, F2, s0, d0, Eps[0])
        #
        sig =  array([ssig,0], dtype=double)
        MatM = array([[dsig,0],[0,0]], dtype=double)                # material tangential stiffness
        return sig, MatM, [Eps[0], sig[0]]
        
class Bond(Material):
    def __init__(self, PropMat, PropMat2, TraPre ):
#            def __init__(self, SymmetricVal, RTypeVal, UpdateVal, Updat2Val, StateVarVal, NDataVal):
        Material.__init__(self, True,         None,     True,     False,     5,           6, 'BOND')                 # should be symmetric as only diagonal of MatM is  used
#        self.Symmetric = True                                              # flag for symmetry of material matrices
#        self.RType                                                         # type of regularization 0: SDA, 1: gradient 2: crack band, 3: SDA
#        self.Update = False                                                # has specific update procedure for update of state variables
#        self.StateVar = 4                                                  # [0] smallest slip ever reached; [1] corresponding bond stress value; [2] largest slip ever reached; [3] corr. bond stress
                                                                            # [4] bond status: 1 loading positive, 2 unloading positive, -1 loading negative, -2 unloading negative  
#        self.NData = 6
        self.StiffLong = PropMat[0]
        self.StiffLate = PropMat[1]
        self.PropMat   = PropMat2
        self.TraPre    = TraPre                                             # transverse pressure (True or False) - Ahmad
    def Sig(self, ff, CalcType, Dt, elI, ipI,   Elem, ds,   sl,   dtmpI, tempI, EpsR):
        if CalcType == 0: return [], [], []
        if self.TraPre:
            sigper, fc, fct = Elem.transverse_stress()                      # sigper = mean stress in the concrete orthogonal to the bar axis 
            if sigper > 0:
                mult = max(0.,1.0 - 0.3*sigper/fct)                         # according to modelcode 2010 6.1.1.3.2
            else:
                mult = 1.0 - tanh(0.2*sigper/(0.1*fc))
        else: mult = 1.
        # loading
        if sl[0]>Elem.StateVar[ipI,0]-1.0e-12 or sl[0]<Elem.StateVar[ipI,2]+1.0e-12:
            Elem.StateVarN[ipI,4] = 0
            s1 = self.PropMat[0]                                            # PropMat[0] is 2nd value in input line, 1st currently ignored
            F1 = self.PropMat[1]
            s2 = self.PropMat[2]
            F2 = self.PropMat[3]
            s0 = self.PropMat[4]                                            # up to this slip initial stiffness is constant
            d0 = self.PropMat[5]                                            # initial stiffness / Young's modulus
            ssig, dsig = SpringSpline( s1, F1, s2, F2, s0, d0, sl[0])
        # unloading positive
        elif sl[0]>0:
#            dsig = Elem.StateVar[ipI,1]/Elem.StateVar[ipI,0]
            dsig = Elem.StateVarN[ipI,1]/Elem.StateVarN[ipI,0]
            ssig = dsig*sl[0]
            Elem.StateVarN[ipI,4] =  1 
        # unloading negative
        else:
#            dsig = Elem.StateVar[ipI,3]/Elem.StateVar[ipI,2]
            dsig = Elem.StateVarN[ipI,3]/Elem.StateVarN[ipI,2]
            ssig = dsig*sl[0] 
            Elem.StateVarN[ipI,4] = -1 
        #
        ssig = mult*ssig                                                    # bond stress multiplied by transverse pressure factor - Ahmad
        dsig = mult*dsig                                                    # tang stiffness multiplied - uhc
        if   sl[0] > Elem.StateVar[ipI,0]:                                  # [ipI,0] largest slip (signed) ever reached
            Elem.StateVarN[ipI,0] = sl[0]
            Elem.StateVarN[ipI,1] = ssig                                    # corr. bond stress
        elif sl[0] < Elem.StateVar[ipI,2]:                                  # [ipI,2] smallest slip ever reached
            Elem.StateVarN[ipI,2] = sl[0] 
            Elem.StateVarN[ipI,3] = ssig                                    # corr. bond stress
        #
        if Elem.Type in   ["Bond2D2","Bond2D3","BondAX2",'BondAX3']:
            MatM = array([[dsig, 0.],[ 0., self.StiffLate      ]])
            tt   = array( [ssig,           self.StiffLate*sl[1]])
            return tt, MatM, [sl[0], sl[1], tt[0], tt[1]]                   # returns slip, bond stress
        elif Elem.Type in ["Bond3D2","Bond3D3"]:
            MatM = array([[dsig, 0., 0.],[ 0., self.StiffLate, 0. ],[ 0., 0., self.StiffLate] ])
            tt   = array( [ssig, self.StiffLate*sl[1], self.StiffLate*sl[2]])
            return tt, MatM, [sl[0], sl[1], sl[2], tt[0], tt[1], tt[2]]     # # returns slip, bond stress
        else: raise NameError("ConFemMaterials::Bond.MatC: unknown element material dimension "+Elem.Type)
    def UpdateStateVar(self, Elem, ff):
#        for j in range(Elem.StateVar.shape[0]):
#            if Elem.StateVarN[j,4] != Elem.StateVar[j,4]: 
#                print("bond status changed ",Elem.Label,j,Elem.StateVar[j,4],"->",Elem.StateVarN[j,4])
#                print("bond status changed ",Elem.Label,j,Elem.StateVar[j,4],"->",Elem.StateVarN[j,4], file=ff)
        for j in range(Elem.StateVar.shape[0]):
            Elem.StateVar[j,:] = Elem.StateVarN[j,:]


class MatTester():
    def __init__(self ):
        self.ZeroD = 1.e-9 
    def Run(self, Name, MName, ElSet, end=0.005):                                                # Smallest float for division by Zero
        from ConFemInOut import DataInput
        from ConFemBasics import EvNodeC2D, EvNodeEl
#        from ConFemElem import TxDxI
        import scipy.spatial as spatial
        f1=open( Name+".in.txt", 'r')
        f6=open( Name+".MatTester.txt", 'w')
#        NodeList, ElList, MatList, _, NoLabToInd, SecDict = ReadInputFile(f1, f6, False)
        NodeList, ElList, MatList, SecDict, NoLabToInd, _, _, _ = DataInput( f1, f6, False)
        NoIndToCMInd = [i for i in range(len(NodeList))]                    # maps identical before Cuthill McKee
        f1.close()
        EvNodeEl( NodeList, ElList, NoLabToInd)                             # determine element indices belonging to each node and more - update required for Cuthill McKee
        ContinuumNodes, ConNoToNoLi = EvNodeC2D( NodeList )                 # build list of node coordinates for 2D/3D continuum nodes only and build list to map them to NodeList indices
        if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes )   # for search purposes, e.g. for EFG or aggregates or embedded truss elements
        else:                     CoorTree = None
        for el in ElList:
            if el.Type in ["T2D2I","T2D3I","T3D2I","T3D3I","B23I","B23EI","TAX2I","BAX21EI",'BAX23EI']:
                el.CreateBond( ElList,NodeList,SecDict, CoorTree,ConNoToNoLi, MatList, NoLabToInd,NoIndToCMInd)    # determine continuum elements connected via bond
        Echo(f"ConFemMat: {Name:s}, Mat name: {MName:s}", f6)
        MatLi = MatList[MName.upper()]
#        end, N = 0.005, 400
        N = 400
        Found = False
        for Elem in ElList:

            print('XXX',Elem.Set.upper(),ElSet.upper(),Elem.MatN.upper())

            if Elem.Set.upper() == ElSet.upper():
                if Elem.MatN.upper()==MName.upper():
                    Found = True
                    Echo(f"Found specified element set: {ElSet:s}", f6)
                    Echo(f"Found specified material: {MName:s} of type: {MatLi.Type:s}", f6)
                    Elem.Ini2( NodeList,NoIndToCMInd, MatList, SecDict )
                    self.RunMatT( MatLi, Elem, 0, True, f6, end, N )
                    break
        if not Found: raise NameError("ConMat: Specified material or element set not found",ElSet,MName)
        Echo("finished", f6)
        f6.close()
        plt.show()

    def Cplot(self, x, y, z, label, col, xLegend, yLegend ): 
        FonSizAx='large'     # fontsize axis
        LiWiCrv=2            # line width curve in pt
        PP = plt.figure()
        P0 = PP.add_subplot(111,title=label)
        P0.title.set_fontsize('x-large')
        if col !=None:  
            for i in range(len(x)): P0.plot(x[i], y[i], col, linewidth=LiWiCrv)
        else: 
            for i in range(len(x)): P0.plot(x[i], y[i], linewidth=LiWiCrv)
        P0.tick_params(axis='x', labelsize=FonSizAx)
        P0.tick_params(axis='y', labelsize=FonSizAx) #'x-large')
        P0.grid()
        if z != None:
            P1 = P0.twinx()
            for i in range(len(z)): P1.plot(x[i], z[i])
        if len(xLegend)>0: P0.set_xlabel( xLegend, fontsize='large' )      # medium, large, x-large
        if len(yLegend)>0: P0.set_ylabel( yLegend, fontsize='large' )
        PP.autofmt_xdate()

    def AmpVal(self, Time, Data):
        nD = len(Data)
        if Time>=Data[nD-1][0]:                         # time exceeds last value
            DelT = Data[nD-1][0]-Data[nD-2][0]
            if DelT<self.ZeroD: raise NameError ("ConFemMaterials::MatTester.AmpVal: something wrong with amplitude 1")
            return (Time-Data[nD-2][0])/DelT * (Data[nD-1][1]-Data[nD-2][1]) + Data[nD-2][1]
        for i in range(len(Data)-1):
            if Data[i][0]<=Time and Time<=Data[i+1][0]: # time within interval
                DelT = Data[i+1][0]-Data[i][0]
                if DelT<self.ZeroD: raise NameError ("ConFemMaterials::MatTester.AmpVal: something wrong with amplitude 2")
                return (Time-Data[i][0])/DelT * (Data[i+1][1]-Data[i][1]) + Data[i][1]
    def AmpVal2(self, Time):                        # exponential function with decreasing derivative
        if Time<=1.:
            return Time
        else:      
            a, b  = 5.687010672, .2133556056           # see maple preprocessor
            return (a-1.)*(1-exp(-b*(Time-1.)))+1.
    def RunMatT(self, MatLi, Elem, PreL, PlotF, ff, Scale, N):
        import ConFemMat
        # RCBeam
        if isinstance( MatLi, ConFemMat.RCBeam) or isinstance( MatLi, ConFemMat.WraTCBeam) or\
           isinstance( MatLi, RCBeam) or isinstance( MatLi, WraTCBeam):
            T = TestRCBeam( MatLi, Elem, ff) 
            T.BendingTest( [MatLi.PropMat[0][3], 1.5*MatLi.epsLim], [[ 0. , -1., -2.]], PreL)    
        # Mises beam
        if isinstance(MatLi, ConFemMat.MisesBeam2D_) or isinstance(MatLi, ConFemMat.MisesBeam2D):
            T = TestMisesBeam( MatLi, Elem, ff)
            T.BendingTest()

        # ElasticLT
        elif isinstance( MatLi, ConFemMat.ElasticLT):
            print('ElasticLT', MatLi.Emod, MatLi.nu, MatLi.fct, MatLi.Gf, MatLi.epsct) 
            N = 1000 # 200
            de = MatLi.epsct
            deps= [de/N,0,0]
            Eps,  Sig,  Emo,  WW,  Eps_,  Sig_,  Emo_,  WW_, Tim, Tim_ = [], [], [], [], [], [], [], [], [], []
            Eps1, Sig1, Emo1, WW1, Eps1_, Sig1_, Emo1_, WW1_           = [], [], [], [], [], [], [], []
#           Amp = [[0.,0.],[9.,9.]]
#            Amp = [[0.,0.],[4.,4.],[8.5,-0.5],[9.,0.]]
            Amp = [[0.,0.],[2.,2.],[4.5,-0.5],[5.,0.]]
            Amp1 = [[0.,0.],[23.,0.],[30.,7.],[32.,-0.5],[33.,0.]]
            State_, State1_  = 0, 0
            for i in range(int(23*N)):                                       # 5, 9, 23, 35 -- this rules length of loading path 
                fac, fac1 = self.AmpVal(float(i)/N, Amp), self.AmpVal(float(i)/N, Amp1)
                eps = [fac*de,fac1*de,0]
                Dt = 0.001
#                                     Sig  ff,CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR):
                sigI, C, Data = MatLi.Sig( ff, 0, Dt,   0,0,      Elem, deps,eps, 0, 0,       None)
                sigI, C, Data = MatLi.Sig( ff, 1, Dt,   0,0,      Elem, deps,eps, 0, 0,       None)# stress, stress incr, tangential stiffness (Script Eq. (1.20) -> linear-elastic material;
                C_ = array([[C[0,0],C[0,1]],[C[1,0],C[1,1]]])
                la,v = eigh(C_) 
                Flag_ = MatLi.UpdateStateVar(Elem, ff)
                MatLi.UpdateStat2Var(        Elem, ff, True, True)
                for k in range(len(Data)): Elem.DataP[0][k] = Data[k]

                State  = int(Elem.StateVarN[0,0])
                State1 = int(Elem.StateVarN[0,6])
                if State != State_:                                         # new state - new color 
                    Eps.append(Eps_); Sig.append(Sig_); Emo.append(Emo_); WW.append(WW_); Tim.append(Tim_) 
                    Eps_, Sig_, Emo_, WW_, Tim_ = [], [], [], [], []
                    State_ = State
                if State1!=State1_:                                         # new state - new color
                    Eps1.append(Eps1_); Sig1.append(Sig1_); Emo1.append(Emo1_); WW1.append(WW1_)
                    Eps1_, Sig1_, Emo1_, WW1_ = [], [], [], []
                    State1_ = State1
                    
                Sig_.append( sigI[0]); Eps_.append( eps[0]); Emo_.append( C[0,0]); WW_.append( Data[8]); Tim_.append(float(i)/N)
                Sig1_.append(sigI[1]); Eps1_.append(eps[1]); Emo1_.append(C[1,1]); WW1_.append(Data[9])
                Echo(f"{i:4d} State {State:d}, {fac:6.3f}, {eps[0]:8.3e}, {eps[1]:8.3e}, s {sigI[0]:8.3e}, {C[0,0]:8.3e}, EV {la[0]:8.3e}, {la[1]:8.3e}", ff) 
            
            Eps.append(Eps_); Sig.append(Sig_); Emo.append(Emo_); WW.append(WW_); Tim.append(Tim_)
            Eps1.append(Eps1_); Sig1.append(Sig1_); Emo1.append(Emo1_); WW1.append(WW1_)
            if PlotF: 
                self.Cplot(Eps, Sig, None,'eps0-sig0',None,"","")
                self.Cplot(Eps1,Sig1,None,'eps1-sig1',None,"","")
                self.Cplot(Eps, Emo, None,'eps-Emod', None,"","")
                self.Cplot(WW,  Sig, None,'w-sig',    None,"","")
                self.Cplot(Eps,  WW, None,'eps-w',    None,"","")
                self.Cplot(Tim, Eps, None,'time-eps', None,"","")
#                plt.show()
            print('Finished ElasticLT')
            
        elif isinstance( MatLi, ConFemMat.Spring):
            T = TestSpring( MatLi, Elem, ff)
            end = 1.2*MatLi.PropMat[2]
            T.UniaxTest( end=end )
        elif isinstance( MatLi, ConFemMat.Bond):
            T = TestBond( MatLi, Elem, Scale, ff)
            T.UniaxTest()
        elif isinstance( MatLi, ConFemMat.IsoDamage) or isinstance( MatLi, ConFemMat.MicroPlaneDam) or isinstance(MatLi, ConFemMat.ElasticSDA):
            T = TestContinuum( MatLi, Elem, ff )
            #
#            epsT, sigT = T.UniaxTest( PlotF, N=10000, RangeEpsN= 2 )        # stability tension
#            epsC, sigC = T.UniaxTest( PlotF, N=10000, RangeEpsN=15, Tension=False ) # stability compression
            #
            epsT, sigT = T.UniaxTest( PlotF, N=200 )
#            epsC, sigC = T.UniaxTest( PlotF, N=200, RangeEpsN=15, Tension=False )
#            pp = plt.figure().add_subplot();      pp.plot(epsT, sigT);            pp.plot(epsC, sigC);             pp.grid()   # quick and dirty plot of both
        else:
            print('Unknown material type',MatLi.Type)
        return 0

class TestMisesBeam( MatTester):
    def __init__(self, MatLi, Elem, ff):
        self.MatLi   = MatLi
        self.Elem    = Elem
        self.LogFile = ff
    def BendingTest(self):
        eps, kap = 0., 0.
        epsUlt = self.MatLi.epsU
        rr = 0.5*self.Elem.Geom[1,2]                                        # Diam/2
        kk_, ee_, N_, M_, EJ_, sig = [0.], [0.], [0.], [0.], [], [0, 0]
        eps_, kap_ = 1.0e-5, 1.0e-4                                         # nominal values
        x = 1.2                                                             # factor for strain reference axis
        #
#        si, epT, brc, Eps  = [ 1., -1.,  1.], [ 0.25, -0.15, 0.1 ], ["ge", "le", "ge"], [[1.0e-5, 0.],[1.0e-5,0.],[1.0e-5,0.]] # normal force tension
#        si, epT, brc, Eps = [-1.,  1., -1.], [-0.25,  0.15, -0.1], ["le", "ge", "le"], [[1.0e-5,0.],[1.0e-5,0.],[1.0e-5,0.]] # normal force compression
        #
        si, epT, brc, Eps = [ 1., -1.,  1.], [0.12, 0.065, 0.22], ["ge", "le", "ge"], [[0., 1.0e-4], [0., 1.0e-4],[0., 1.0e-4]]  #      pos bending with un- and reloading
        si, epT, brc, Eps = [-1.,  1., -1.], [-0.12,-0.065, -0.22], ["le", "ge", "le"], [[0.,1.0e-4],[0.,1.0e-4],[0.,1.0e-4]]            # neg bending with un- and reloading
        #
        si, epT, brc, Eps =  [  1., -1.,  1.], [  0.17, 0.165, 0.30],["ge", "le", "ge"], [[x*eps_,kap_],[x*eps_,kap_],[x*eps_,kap_]]            # pos bending with normal force
#        si, epT, brc, Eps =  [ -1.,  1., -1.], [   -0.17, -0.15, -0.30], ["le", "ge", "le"], [[x*eps_,kap_],[x*eps_,kap_],[x*eps_,kap_]]            # pos bending with normal force
        #
        for si_, ept_, brc_, Eps_ in zip( si, epT, brc, Eps):
            print(si_,ept_,brc_,Eps_)
            while True:
                eps = eps + si_*Eps_[0]
                kap = kap + si_*Eps_[1]
                epsL = eps + rr*kap
                epsU = eps - rr*kap
    #            sig, MatM, dataS = self.MatLi.Sig(ff,  CalcType, Dt, elI, ipI, Elem, Dps, Eps, dTmp, Temp, EpsR)
                sig, MatM, data = self.MatLi.Sig(None, 1, 0, 0, 0, self.Elem, [si_*Eps_[0], si_*Eps_[1]], [eps, kap], [0., 0.], [0., 0.], None)
                Echo(f"eps {eps:9.6f} kap {kap:9.6f} N {sig[0]:11.8f} M {sig[1]:11.8f} EA {MatM[0,0]:11.5e} EJ {MatM[1,1]:11.5e}\
 - {data[0]:9.6f} {data[1]:9.6f} {data[2]:9.6f} {data[3]:9.6f} {data[4]:9.6f} {data[5]:9.6f} {data[6]:3.3f} {data[7]:3.3f}", self.LogFile)
                kk_ += [kap]
                ee_ += [eps]
                N_  += [sig[0]]
                M_  += [sig[1]]
                EJ_ += [MatM[1, 1]]
                self.MatLi.UpdateStateVar(self.Elem, None)
                if brc_=="ge":
                    if epsL >= ept_*epsUlt: break
                else:
                    if epsL <= ept_*epsUlt: break
        self.Cplot([kk_],[M_], None, 'moment-curvature', None,"x","y")
        self.Cplot([ee_],[N_], None, 'normal force-strain', None, "x", "y")

class TestRCBeam( MatTester):
    def __init__(self, MatLi, Elem, ff):
        self.MatLi   = MatLi
        self.Elem    = Elem
        self.LogFile = ff
    def BendingTest(self, Data0, Data1, PreL, PlotF=True, n=200):
        # stress strain curve concrete
        e1, e2 = Data0[0], Data0[1] #self.MatLi.PropMat[0][3], 1.5*self.MatLi.epsLim
        de = (e2-e1)/n
        ep, si, ds = [], [], []
        for i in range(n+1):
            eps = e1 + i*de
            sig, dsig = self.MatLi.MatUniaxCo( self.MatLi.PropMat[0], eps)
            ep += [1000.*eps]
            si += [sig]
            ds += [dsig]
        if PlotF: self.Cplot( [ep], [si], [ds], 'eps-sig','r-',"curvature","moment")
        # moment curvature
        Unload = False #True # False        
        kk, N, M, EJ = [], [], [], []
        kapD_= 0.0002
#            Deps0 = array([0.000001, 0,      0, 0])
#            Deps1 = array([0,        0.00001, 0, 0])
        for NorT in Data1[0]: 
            kapD = kapD_
            kk_, N_, M_, EJ_, sig  = [], [], [], [], [0,0]
            kap, epsL, epsLold, kap_, epsL_ = -kapD, 0, 0, 0, 0
            self.Elem.StateVar[:,:] = 0.
            Flag = True
            while Flag:
#                    Mold= sig[1]
                kap = kap + kapD
                epsLold = epsL
                NIter = 20
                for j in range(NIter):
                    eps = array([epsL, kap, epsL_, kap_])
                    #                         Sig( self.LogFile, CalcType, Dt, elI, ipI, self.Elem, Dps,               Eps,       dTmp,    Temp,    EpsR):
                    sig, MatM, dataS  = self.MatLi.Sig(None,1,        0,  0,   0,   self.Elem,[epsL-epsLold,kapD],eps,       [0.,0.], [0.,0.], None)
                    # sense of the following not quite clear uhc 200301 -- looks like numerical differentiation
#                    sig0,MatM0,dataS0 = self.MatLi.Sig(None,1,        0,  0,   0,   self.Elem,[epsL-epsLold,kapD],eps-Deps0, [0.,0.], [0.,0.], None)
#                    sig1,MatM1,dataS1 = self.MatLi.Sig(None,1,        0,  0,   0,   self.Elem,[epsL-epsLold,kapD],eps-Deps1, [0.,0.], [0.,0.], None)
#                    MM00 = (sig[0]-sig0[0])/Deps0[0]
#                    MM01 = (sig[0]-sig1[0])/Deps1[1]
#                    MM10 = (sig[1]-sig0[1])/Deps0[0]
#                    MM11 = (sig[1]-sig1[1])/Deps1[1]
                    #
                    Res = sig[0] - NorT
                    if fabs(Res)<1.e-5: break
                    epsL = epsL - (1/MatM[0,0]) * Res
                if j == NIter: Echo(f"convergence interrupted with kappa {kap:f}, iter {j:d}", self.LogFile)
#                    Flag_ = self.MatLi.UpdateStateVar(self.Elem, None)                        # update for sigy in MisesUniaxial, Flag_ is dummy
                if abs(sig[1])<PreL: epsL_, kap_ = epsL, kap
                # termination control - StateVarN update in Sig        
                if self.Elem.StateVarN[0,2]<0.95*self.MatLi.PropMat[0][3]: Flag,TVal,TType = False,self.Elem.StateVarN[0,2],"concrete failure"   # control of concrete - keep sign in mind 0.98
                for i in range(self.Elem.Geom.shape[0]-2):                           # control of reinforcement strains over number of reinforcement layers
                    if self.Elem.RType[i]=='RC' and abs(self.Elem.StateVarN[i+1,2])>0.99*self.MatLi.PropMat[1][4]: Flag,TVal,TType = False,self.Elem.StateVarN[i+1,2],"rebar failure" # ultimate rebar strain
                    if self.Elem.RType[i]=='TC' and abs(self.Elem.StateVarN[i+1,2])>0.99*self.MatLi.ReinfTC.PropMat[3]: Flag,TVal,TType = False,self.Elem.StateVarN[i+1,2],"TRCr failure"
                if not Flag:
                    if not Unload: 
                        Echo(f"{TType:s} {TVal:f}", self.LogFile) 
                        break                                        # terminate for next normal force 
                    else:          
                        Echo(f"unloading", self.LogFile) 
                        kapD, Flag = -kapD_, True                               # unloading
                if (kapD<0 and sig[1]<0.05):                                    # terminate unloading for next normal force
                    kapD, Flag = kapD_, False 
                    break #continue
                # numerical differentiation control
#                    dMom1 = sig[1]-Mold
#                    dMom2 = MatM[1,0]*(epsL-epsLold)+MatM[1,1]*kapD
#                    dMomE = ((dMom1-dMom2)/dMom1)
#                    print(j, PreL, eps, sig, '__',   MatM[0,0], MM00, '_', MatM[0,1], MM01,'_', MatM[1,0], MM10, '_', MatM[1,1], MM11, '__', dMom1, dMom2, dMomE)
#                    if abs(dMomE)>0.10: print(dMomE)
                if abs(eps[1])>ZeroD:
                    Echo(f"{j:d}, kap {eps[1]:f}, N {sig[0]:10.6f}, M {sig[1]:f}, eps_s {dataS[4]:8.3f}, eps_c {dataS[5]:8.3f}, EJ_s {sig[1]/eps[1]:f}, EJ_t {MatM[1,1]:f}", self.LogFile)
                kk_ += [kap]
                N_  += [sig[0]]
                M_  += [sig[1]]
                EJ_ += [MatM[1,1]]
            kk += [kk_]
            N  += [N_]
            M  += [M_]
            EJ += [EJ_]
        if PlotF: 
            self.Cplot(kk,N, None, 'normal force-curvature', None,"x1","y1")
#            self.Cplot(kk,M,EJ,'moment-curvature', None,"x","y")
            self.Cplot(kk,M, None, 'moment-curvature', None,"x","y")
        for m in self.MatLi.PropMat: Echo([f"Par {v:f}" for v in m], self.LogFile)
        if self.MatLi.ReinfTC!=None: 
            for m in self.MatLi.ReinfTCt: Echo([f"Par {v:f}" for v in m], self.LogFile)
    
class TestBond( MatTester):
    def __init__(self, MatLi, Elem, Scale, ff):
        self.MatLi   = MatLi
        self.Elem    = Elem
        self.Scale = Scale
        self.LogFile = ff
    def UniaxTest(self, PlotF=True, N=200):
        end = self.Scale
        ep, si, ds = [], [], []
        de = end/N
        for i in range(N+1):      
            eps = i*de
            sig, MatM, _ = self.MatLi.Sig( self.LogFile, 1, 0, 0, 0,  self.Elem, 0,   [eps,0],   None, None, None)
            ep += [eps]
            si += [sig[0]]
            ds += [MatM[0,0]]
            Echo(f"{i:4d}, e {eps:f}, s {sig[0]:f}, ds {MatM[0,0]:e}", self.LogFile)
            self.Elem.StateVar[0,:] = self.Elem.StateVarN[0,:]
        Echo([f"Par {v:f}" for v in self.MatLi.PropMat], self.LogFile)
        if PlotF: 
            self.Cplot( [ep], [si], [ds], 'Bond-Slip','darkred','slip [m]', 'bond stress [MPa]')
#                self.Cplot( [xi], [si], None, 'Bond-Slip Steel',None,'slip [m]', 'bond stress [MPa]')
#                self.Cplot( [xi], [si], None, 'Bond-Slip Carbon','darkred','slip [m]', 'bond stress [MPa]')
class TestSpring( MatTester):
    def __init__(self, MatLi, Elem, ff):
        self.MatLi   = MatLi
        self.Elem    = Elem
        self.LogFile = ff
    def UniaxTest(self, PlotF=True, N=200, end=0.002 ):
        ep, si, ds = [], [], []
        de = end/N
        for i in range(N+1):                           # stress strain curve
            eps = i*de
            sig, dsig, _ = self.MatLi.Sig( self.LogFile, None, None, None, None, None, [de,0], [eps,0], None, None, None)
            ep += [1000.*eps]
            si += [sig[0]]
            ds += [dsig[0,0]]
            Echo(f"{i:4d}, displ {eps:f}, force {sig[0]:f}, df {dsig[0,0]:e}", self.LogFile)
        Echo([f"Par {v:f}" for v in self.MatLi.PropMat], self.LogFile)
        if PlotF: 
            self.Cplot( [ep], [si], [ds], 'eps-sig','r-', "_", "_")

class TestContinuum( MatTester):
    def __init__(self, MatLi, Elem, ff):
        self.MatLi   = MatLi
        self.Elem    = Elem
        self.LogFile = ff
        self.buffer  = None
    # from  https://stackoverflow.com/questions/39474254/cardanos-formula-not-working-with-numpy
    def find_cubic_roots(self, a,b,c,d, bp = False):
        from numpy import isreal, real 
        from scipy.special import cbrt
        a,b,c,d = a+0j, b+0j, c+0j, d+0j
        all_ = (a != pi)
        Q = (3*a*c - b**2)/ (9*a**2)
        R = (9*a*b*c - 27*a**2*d - 2*b**3) / (54 * a**3)
        D = Q**3 + R**2
        S = 0 #NEW CALCULATION FOR S STARTS HERE
        if isreal(R + sqrt(D)): S = cbrt(real(R + sqrt(D)))
        else:                   S = (R + sqrt(D))**(1/3)
        T = 0 #NEW CALCULATION FOR T STARTS HERE
        if isreal(R - sqrt(D)): T = cbrt(real(R - sqrt(D)))
        else:                   T = (R - sqrt(D))**(1/3)
        result = zeros(tuple(list(a.shape) + [3])) + 0j
        result[all_,0] = - b / (3*a) + (S+T)
        result[all_,1] = - b / (3*a)  - (S+T) / 2 + 0.5j * sqrt(3) * (S - T)
        result[all_,2] = - b / (3*a)  - (S+T) / 2 -  0.5j * sqrt(3) * (S - T)
        #if bp:
            #pdb.set_trace()
        return result    
    def MatStability(self, C, Emod, Tens):
        from numpy import real, isreal
        def d( x ):                                                         # value of determinant
            return f*( a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0 )
        def dd( x ):                                                        # value of 1st derivative of determinant
            return f*( 4.*a4*x**3 + 3.*a3*x**2 + 2.*a2*x + a1 )
        f = 1./sqrt(Emod)
        C_ = f*C                                                         # to normalize large numbers
        a4 = (2*C_[2,2]*C_[1,1]-2*C_[2,1]*C_[1,2]) 
        a3 = (2*C_[0,2]*C_[1,1]+C_[2,0]*C_[1,1]-2*C_[0,1]*C_[1,2]-C_[2,1]*C_[1,0])
        a2 = (C_[0,0]*C_[1,1]+2*C_[0,2]*C_[2,1]-2*C_[0,1]*C_[2,2]-C_[0,1]*C_[1,0]+2*C_[2,0]*C_[1,2]-2*C_[2,2]*C_[1,0])
        a1 = (C_[0,0]*C_[2,1]+2*C_[0,0]*C_[1,2]-2*C_[0,2]*C_[1,0]-C_[0,1]*C_[2,0])
        a0 = (2*C_[0,0]*C_[2,2]-2*C_[0,2]*C_[2,0])
        rr = []
        cr = self.find_cubic_roots(4.*a4, 3.*a3, 2.*a2, a1)
        ro = isreal(cr)                                                     # returns bool 
        for i, b in enumerate(ro):                                          # solution given as tan of orientation angle
            if b: rr += [ real(cr[i]) ] 
        qq = [ d(i)      for i in rr ]
        dq = [ dd(i)     for i in rr ]
        rp = [ arctan(i) for i in rr ]                                      # orientation angle from tan
        Q, ax, epsJ = [], [], []
        for i, q in enumerate(qq):
            if q < 0.:
                n = array([cos(rp[i]),sin(rp[i])])
                # from MatStability2D.mws
                QQ = array([[ n[0]*n[0]*C_[0,0]+2*n[0]*C_[0,2]*n[1]+n[1]*C_[2,0]*n[0]+2*n[1]*n[1]*C_[2,2], 2*n[0]*n[0]*C_[0,2]+n[0]*C_[0,1]*n[1]+2*n[1]*C_[2,2]*n[0]+n[1]*n[1]*C_[2,1] ],\
                            [ n[0]*n[0]*C_[2,0]+2*n[1]*C_[2,2]*n[0]+n[1]*C_[1,0]*n[0]+2*n[1]*n[1]*C_[1,2], 2*n[0]*n[0]*C_[2,2]+n[0]*C_[2,1]*n[1]+2*n[1]*C_[1,2]*n[0]+n[1]*n[1]*C_[1,1] ]])
                self.StabilityFlag = True 
                mL= sqrt(1.+(QQ[0,1]/QQ[0,0])**2)                           # length of vector m 
                m = array([1.,-QQ[0,1]/QQ[0,0]])/mL
#                m = array([-1.,QQ[0,1]/QQ[0,0]])/mL
                epsJ = array([ [n[0]*m[0],0.5*(n[0]*m[1]+n[1]*m[0])] , [0.5*(m[0]*n[1]+n[0]*m[1]),n[1]*m[1]] ]) # normalized strain jump
                if not Tens: epsJ = -epsJ
                la,v = eigh( epsJ)
                Q += [[ 180.*self.buffer[0][i]/pi,self.buffer[1][i],self.buffer[2][i],self.buffer[3],self.buffer[4],self.buffer[5], self.buffer[6],self.buffer[7],self.buffer[8],self.buffer[9] ],\
                      [ 180.*rp[i]/pi, dq[i], qq[i], epsJ, la, v, C, QQ, n, m ]]
        self.buffer = [rp, dq, qq, array([[0,0],[0,0]]), array([0,0]) , array([[0,0],[0,0]]), C, [], [], [] ] 
        return Q
    def UniaxTest(self, PlotF, N=100, RangeEpsN=5, Tension=True):
        Sig1, Sig2, Sig3, Eps1, DSig, DSig_, CC = [0], [0], [0], [0], [0], [0], [0]
        EpsEnd = RangeEpsN * self.MatLi.eps_ct
        if Tension: de, eps1 =  EpsEnd/N, 0
        else:       de, eps1 = -EpsEnd/N, 0    
        nu = self.MatLi.nu
        self.StabilityFlag = False
        for i in range(N):
            eps1 = eps1+de 
            eps  = array([eps1,-nu*eps1,0,0,0,0])
            deps = array([de  ,-nu*de,  0,0,0,0])
            sigI, C, _ = self.MatLi.Sig( self.LogFile, 1, de,0,0, self.Elem, deps, eps, [], [], [])
            dSig = C[0,0] - nu*( C[0,1] ) #+ C[0,1])
            Echo(f"{i:4d}, {eps[0]:10.3e}, {eps[1]:10.3e}, {sigI[0]:10.3e}, {sigI[1]:10.3e}, _ , {C[0,0]:10.3e}, {C[0,1]:10.3e}, {dSig:10.3e}", self.LogFile)
            if not self.StabilityFlag: 
                Q = self.MatStability( C, self.MatLi.Emod, Tension)
                for j in Q:
                    Echo(f"Mat Instability: orientation {j[0]:f}, dq {j[1]:f}, q {j[2]:f}, epsJump {j[3][0,0]:f},{j[3][0,1]:f},{j[3][1,0]:f},{j[3][1,1]:f}, prin eps 1 {j[4][0]:f}_{j[5][0,0]:f},{j[5][0,1]:f}, 2 {j[4][1]:f}_{j[5][1,0]:f},{j[5][1,1]:f}", self.LogFile)
                    C, Q, n, m = j[6], j[7], j[8], j[9]
                    Echo(f'C 0 {", ".join(zzz for zzz in [f"{x:8.2f}" for x in C[0,0:3]])}', self.LogFile)
                    Echo(f'C 1 {", ".join(zzz for zzz in [f"{x:8.2f}" for x in C[1,0:3]])}', self.LogFile)
                    Echo(f'C 2 {", ".join(zzz for zzz in [f"{x:8.2f}" for x in C[2,0:3]])}', self.LogFile)
                    if len(Q)>0:
                        Echo(f'Q 0 {", ".join(zzz for zzz in [f"{x:8.3f}" for x in Q[0,0:2]])}', self.LogFile)
                        Echo(f'Q 1 {", ".join(zzz for zzz in [f"{x:8.3f}" for x in Q[1,0:2]])}', self.LogFile)
                        Echo(f'n   {", ".join(zzz for zzz in [f"{x:8.3f}" for x in n[0:2]])}', self.LogFile)
                        Echo(f'm   {", ".join(zzz for zzz in [f"{x:8.3f}" for x in m[0:2]])}', self.LogFile)
            Sig1 += [sigI[0]]
            Sig2 += [sigI[1]]
            Sig3 += [sigI[2]]
            Eps1 += [eps[0]]
            CC   += [C[0,0]]
            DSig_+= [dSig]
            if i>1: DSig += [0.5*(Sig1[i]-Sig1[i-2])/de]
        DSig += [0.5*(Sig1[N]-Sig1[N-2])/de, 0.]
        if PlotF:
            ZeLi = zeros((N+1),dtype=float)                                 # list with zeros
#            self.Cplot([Eps1],[Sig1],None,'Sigma-Epsilon',None,'epsilon [-]','sigma [Mpa]')
            self.Cplot([Eps1,Eps1,Eps1],[Sig1,Sig2,Sig3],None,'sig-eps',None,'','')
#            self.Cplot([Eps1,Eps1,Eps1],[Sig1,Sig2,Sig3],[DSig,CC,ZeLi],'sig-eps',None,'','')
            self.Cplot([Eps1,Eps1,Eps1],[Sig1,Sig2,Sig3],[DSig,DSig_,ZeLi],'sig-eps',None,'','')
#                self.Cplot([Eps1],[DSig],None,'dsig-eps',None)
        return Eps1, Sig1
    
if __name__ == "__main__":
    Name, MatName, ElSet = "../DataExamples/E03/E3-04_", 'MAT5', 'bond'
    Name, MatName, ElSet = "../DataExamples/E08/E8-04", 'BOND1', 'bar1'     # elset relates to truss / beam
#    Name, MatName, ElSet = "../DataExamples/E04/E4-01", 'mat1', 'prop1'         #beam
    
#    Name, MatName, ElSet = "../_DataTmp/phi-45coarse", 'bond2', 'E2'
#    Name, MatName, ElSet ="../DataExamples/E06/E6-02", 'ISOD', 'EL1'        # input dataset must also be adapted
#    Name, MatName, ElSet ="../_DataSpecimen/One2D/WillamsTest2D", 'ISOD', 'EL1'
#    Name, MatName, ElSet ="C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/LSP/Coarse/ConFem/Job-CU", "ISOD",  "Part-1-1__PickedSet5"
#    Name, MatName, ElSet ="C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/OneElement2D/Large/ConFem/Job-1", "IsoDam",  "_PickedSet4"
#    Name, MatName, ElSet ="../DataExamples/E06/Uniax2D", 'ELLT', 'EL1'
#    Name, MatName, ElSet ="../DataExamples/E06/Uniax2D", 'MIPL', 'EL1'
#    Name, MatName, ElSet ="C:/Users/uhc/Documents/Work/FoilPap/C3/L6/Abschlussbericht/ConFem/1-one_Steel_Bar_Fig-41_42/Deep_beam", 'ISOD', 'CONC'
#    Name, MatName, ElSet ="C:/Users/uhc/Documents/Work/FoilPap/C3/L6/Abschlussbericht/ConFem/1-one_Steel_Bar_Fig-41_42/Deep_beam", 'BOND1', 'bar1'
#    Name, MatName, ElSet ="C:/Users/uhc/Documents/Work/FoilPap/C3/L6/Abschlussbericht/ConFem/3-one_Carbon_Bar_Fig-45_46/Deep_beam", 'BOND1', 'bar1'
#    Name, MatName, ElSet = "                     ../DataExamples/E08-03/E6-03", 'MAT2', 'EL1'                                # choose material type !
    Name, MatName, ElSet ="../_DataBond/PulloutAxiSym", 'bond', 'STEELBAR'
    Name, MatName, ElSet = "../_DataTmp/AxiSlab", 'matbond', 'REINF'
    Name, MatName, ElSet = "C:/Users/uhc/Documents/Work/FoilPap/2024/Note_FlatSlab/ConFem/PSlabN_2/PSlabN_2", 'matbond',  'UReinfR'
 #                           C:\Users\uhc\Documents\Work\FoilPap\2024\Note_FlatSlab\ConFem\PSlabN_2
#    Name, MatName, ElSet = "../_DataTmp/BeamMises", 'mat12', 'prop1'
#    Name, MatName, ElSet = "../_DataTmp/BeamMises", 'mat2', 'prop1'
#    Name, MatName, ElSet = "../_DataTmp/AxiSlab", 'matbond', 'REINF'
#    Name, MatName, ElSet = "C:/Users/uhc/Documents/Work/Edu/VorlesungenUebungen/CE_RC_2022TUM/ConFem/UniaxTensionBar2/UniaxTenBar_12R", 'MATB', 'ELB',

    X = MatTester()
    X.Run( Name, MatName, ElSet, end=2.0)
