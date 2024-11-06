import matplotlib.pyplot as plt
from numpy import array, fabs,zeros
from scipy import integrate
ZeroD = 1.e-9                                               # to control division by zero

def Echo(*args):
    print(args[0])
    if len(args)>1: print(args[0], file=args[1])

def Cplot( x, y, z, label, col, xLegend, yLegend):
    FonSizAx = 'large'                                  # fontsize axis
    LiWiCrv = 2                                         # line width curve in pt
    PP = plt.figure()
    P0 = PP.add_subplot(111, title=label)
    P0.title.set_fontsize('x-large')
    if col != None:
        for i in range(len(x)): P0.plot(x[i], y[i], col, linewidth=LiWiCrv)
    else:
        for i in range(len(x)): P0.plot(x[i], y[i], linewidth=LiWiCrv)
    P0.tick_params(axis='x', labelsize=FonSizAx)
    P0.tick_params(axis='y', labelsize=FonSizAx)        # 'x-large')
    P0.grid()
    if z != None:
        P1 = P0.twinx()
        for i in range(len(z)): P1.plot(x[i], z[i])
    if len(xLegend) > 0: P0.set_xlabel(xLegend, fontsize='large')  # medium, large, x-large
    if len(yLegend) > 0: P0.set_ylabel(yLegend, fontsize='large')
    PP.autofmt_xdate()

class PolyMat():                                        # stress-strain from polygon
    def __init__(self, Label, eMin, eMax, PropE, PropS):
        self.Label = Label
        # following two values must be defined with values
        self.epsCU = eMin                               # limit for compressive strain (signed)
        self.epsTU = eMax                               # limit for tensile strain  (signed)
        if len(PropE) != len(PropS): raise NameError("PolyMat 1")
        self.epsL  = PropE
        self.sigL  = PropS
    def SigEps(self, eps):
        for i, e in enumerate(self.epsL):
            if i>0:
                Emod = (self.sigL[i]-self.sigL[i-1])/(self.epsL[i]-self.epsL[i-1])
                sig = self.sigL[i-1] + Emod * (eps-self.epsL[i-1])
                if eps<=e: break
        return sig, Emod

class ConcLaw1():                                           # uniaxial concrete (see DIN1045-1 Sec. 9.1.5)
    def __init__(self, Label, PropM):
        self.Label = Label
        # following two values must be defined with values
        self.epsCU = PropM[0]                               # limit for compressive strain (signed)
        self.epsTU = PropM[1]                               # limit for tensile strain  (signed)
        # rest is free according to particular uniaxial stress-strain law
        self.Ec0     = PropM[2]                             # Concrete Young's modulus
        self.fc      = PropM[3]                             # compressive strength
        self.epsc1   = PropM[4]                             # strain at compressive strength
        self.epsDelt = 1.e-6                                # delta for numerical integration
    # this header is mandatory as is
    def SigEps(self, eps):
        def Sig( eps):
#            eta = eps / self.epsc1
            return -self.fc * (k*eta - eta**2) / (1. + (k-2) * eta)
        k = -self.Ec0 * self.epsc1 / self.fc
        # numerical differentiation
#        if eps<=self.epsCU:
#            dsig = (Sig( eps+self.epsDelt) - Sig( eps))/self.epsDelt
#        elif eps>=self.epsTU:
#            dsig = (Sig( eps) - Sig( eps-self.epsDelt))/self.epsDelt
#        else:
#            dsig = 0.5*(Sig( eps+self.epsDelt) - Sig( eps-self.epsDelt))/self.epsDelt
        eta = eps / self.epsc1
        dsig = -self.fc * ((k - 2*eta) / (1 + (k-2)*eta) - (k*eta - eta**2) / (1+(k-2)*eta)**2*(k-2)) / self.epsc1
        return Sig( eps), dsig

class ReinfLaw():
    def __init__(self, Label, PropM):
        self.Label = Label
        # following two values must be defined with values
        self.epsCU = PropM[0]                               # limit for compressive strain (signed)
        self.epsTU = PropM[1]                               # limit for tensile strain  (signed)
        # rest is free according to particular uniaxial stress-strain law
        self.Emod  = PropM[2]                               # Young's modulus
        self.fy    = PropM[3]                               # yield limit
        self.ft    = PropM[4]                               # tensile strength
        self.epsY  = self.fy / self.Emod
        self.Etan  = (self.ft-self.fy)/(self.epsTU-self.epsY)
    def SigEps(self, eps):
        if   eps> self.epsTU:
#            raise NameError("Failure: reinforcement tensile failure",self.epsTU,eps)
            return 0., 0.                                   # catched by IntReinf
        elif eps> self.epsY:
            sig = self.fy + self.Etan * (eps-self.epsY)
            return sig, self.Etan
        elif eps>-self.epsY:
            return self.Emod*eps, self.Emod
        else:
            sig = -self.fy + self.Etan * (eps+self.epsY)
            return sig, self.Etan

class CrossSection():
    def __init__(self, PropC):
        self.height = PropC[0]
        self.width  = PropC[1]

class ReinfCrossSection():
    def __init__(self, PropRC):
        self.PropCross = PropRC

class MatTest():
    def __init__(self, LogFile):
        self.LogFile = LogFile
        self.NIter = 20
        self.IterTarg = 1.0e-5                              # target for strain iteration

    def UniaxTest(self, Mat, PlotF=True, n=200):
        e1, e2 = Mat.epsCU, Mat.epsTU
        de = (e2-e1)/n                                      # negative by definition
        ep, si, ds = [], [], []
        for i in range(n+1):
            eps = e1 + i*de
            sig, dsig = Mat.SigEps( eps)
            ep += [1000.*eps]
            si += [sig]
            ds += [dsig]                                    # tangential material stiffness
        if PlotF:
#            Cplot( [ep], [si], [ds], Mat.Label+' sig-eps','r-',"eps %o","sigma")
            Cplot( [ep], [si], None, Mat.Label+' sig-eps',None,"eps %o","sigma")

    def IntCrossSec(self, epsRef, kap, z1, z2, Mat, CrossSec, nI=50):
        dz = (z2-z1)/nI                                     # z1 local lower, z2 local upper; cross sectional height / number of integration points
        nn = zeros((nI+1), dtype=float)
        mm = zeros((nI+1), dtype=float)
        nde= zeros((nI+1), dtype=float)
        ndk= zeros((nI+1), dtype=float)
        mdk= zeros((nI+1), dtype=float)
        bb = CrossSec.width                                 # width of cross section
        z = z1                                              # cross section / strain coordinate, initial value
        for i in range(nI+1):                               # numerical integration of concrete contribution
            eps  = epsRef - z*kap                           # fiber strain
            sig, dsig = Mat.SigEps( eps)                    # uniaxial stress strain definition
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

    def IntReinf(self, epsRef, kap, MatR, CrossSecR):
        NR,MR, NDE,NDK,MDK, Fail = 0.,0., 0.,0.,0., []
        for i, r in enumerate(CrossSecR.PropCross):
            As = r[0]
            zs = r[1]
            eps = epsRef - zs*kap                           # fiber strain
            if eps > MatR.epsTU: Fail +=["tensile failure reinforcement layer "+"%i"%i+" eps %f"%eps ]
            sig, dsig = MatR.SigEps(eps)                    # uniaxial stress strain definition
            NR +=  As*sig
            MR += -As*sig*zs
            NDE = NDE + dsig*As                             # reinf. contrib. to normal force grad eps
            NDK = NDK - dsig*As*zs                          # reinf. contrib. to normal force grad kappa
            MDK = MDK + dsig*As*zs**2                       # reinf. contrib. to moment grad kappa
        return NR,MR, NDE,NDK,MDK, Fail

    def BendingTest(self, ff, Mat,CrossSec, MatR, CrossSecR, NList, kapStart=0.,kapD=0.0002, PlotF=True):
        z1 = -0.5*CrossSec.height                           # z-coordinate lower edge
        z2 =  0.5*CrossSec.height                           # z-coordinate upper edge
        kL, NL, ML, EJL = [], [], [], []
        for N in NList:                                     # loop for prescribed normal forces
            ic = 0                                          # increment counter
            RFail = []                                      # gets string in case of reinforcement failure
            kap = kapStart-kapD                             # lower starter for curvature
            epsR = 0.                                       # iteration starter for longitudinal strain of reference axis
            kk_, NN_, MM_ = [], [], []                      # lists for respective normal force
            while True:                                     # loop for curvature increments
                kap += kapD
                if fabs(kap) > ZeroD:                       # division by kap in the following
                    ic += 1
                    for j in range(self.NIter):             # iterate for longitudinal strain
                        zCu = (epsR - Mat.epsCU) / kap      # position of compressive strain limitation
                        zTu = (epsR - Mat.epsTU) / kap      # position of tensile strain limitation
                        if zTu < z1: zTu = z1               # tensile limit line below cross section - upward bending compression
                        if zTu > z2: zTu = z2               # tensile limit line above cross section - downward bending compression
                        if zCu < z1: zCu = z1               # compressive limit below cross ssection - downward bending tension
                        if zCu > z2: zCu = z2               # compressive limit above cross section  - upward bending tension
                        # concrete / bulk contribution
                        NN, MM, NDE, NDK1, MDK = self.IntCrossSec( epsR, kap, zTu, zCu, Mat, CrossSec )
                        # reinforcement contribution
                        if CrossSecR != None:
                            NR, MR, NDER, NDKR, MDKR, RFail = self.IntReinf( epsR, kap, MatR, CrossSecR)
                        else:
                            NR, MR, NDER = 0., 0., 0.
                        # iteration control
                        Res = NN+NR - N
                        if fabs(Res) < self.IterTarg:
                            Echo(f"{ic:d}, {j:d}, eps {epsR:f}, kap {kap:f}, N {NN+NR:10.6f}, M {MM+MR:f}, zt {zTu:f}, zc {zCu:f}", ff)
                            kk_ += [kap]
                            NN_ += [NN+NR]
                            MM_ += [MM+MR]
                            break                           # break longitudinal strain iteration
                        epsR = epsR - (1./(NDE+NDER)) * Res # corrector for longitudinal strain
                    # termination criteria
                    if j==self.NIter-1:
                        Echo(f"{ic:d}, kap {kap:f}, iteration not finished zt {zTu:f} zc {zCu:f}", ff)
                        break
                    if (zCu>0. and zCu<z2) or (zCu<0. and zCu>z1):  # break while-loop due to concrete compressive failure upward / downward bending
                        Echo(f"compression failure zt {zTu:f} zc {zCu:f}", ff)
                        break
                    if len(RFail) > 0:                      # reinforcement failure
                        Echo(f"{RFail}", ff)
                        break
#                    if ic>1: break                         # for test purposes
                else:                                       # pass zero kap
                    continue
                #
            kL += [kk_]
            NL += [NN_]
            ML += [MM_]
        if PlotF:
            Cplot( kL, NL, None, 'N-kappa', None,"kappa","N")
            Cplot( kL, ML, None, 'M-kappa', None,"kappa","M")

    def Run(self, Name, MName, ElSet, NormalFList):
        if Name==None:
            f6 = open(self.LogFile, 'w')
            # concrete material parameters: strain range for stresses - starts with lower value,end with upper value (both signed);
            # initial Young's modulus; compressive strength (unsigned); strain for compressive strength (signed)
#            MatC = ConcLaw1("Conc", [ -0.002, 0.0004, 30000., 30., -0.0015 ])
            MatC = ConcLaw1("Conc", [-0.0035, 0.0, 35000., 48., -0.0023])
            # reinforcement material parameters: strain range for stresses - starts with lower value,end with upper value (both signed);
            # elastic Young's modulus; yield limit; ultimate strength
            MatR = ReinfLaw("Steel", [-0.025, 0.025, 200000., 550., 575.])
            # polygonal uniaxial stress-strain behavior
            MatA = PolyMat("CrazyMat", -0.03, 0.03, [-0.02, -0.015, -0.01, 0.01, 0.015, 0.02], [-550., -100., -500., 500., 100., 550.])
            # concrete cross section parameters
            CrC = CrossSection([0.4, 0.2])
            # reinforcement cross section parameters
#            CrR  = ReinfCrossSection([[0.1256637062e-2, -0.15],[0.6283185310e-3,  0.15]])
            CrR = ReinfCrossSection([[0.1256637062e-2, -0.15], [0.1256637062e-2, 0.15]])
        else:
            from ConFemMat import RCBeam, WraTCBeam
            from ConFemInOut import DataInput
            from ConFemBasics import EvNodeC2D, EvNodeEl
            import scipy.spatial as spatial
            f1 = open(Name + ".in.txt", 'r')
            f6 = open(Name + ".MatTester.txt", 'w')
            #        NodeList, ElList, MatList, _, NoLabToInd, SecDict = ReadInputFile(f1, f6, False)
            NodeList, ElList, MatList, SecDict, NoLabToInd, _, _, _ = DataInput(f1, f6, False)
            NoIndToCMInd = [i for i in range(len(NodeList))]  # maps identical before Cuthill McKee
            f1.close()
            EvNodeEl(NodeList, ElList,NoLabToInd)  # determine element indices belonging to each node and more - update required for Cuthill McKee
            ContinuumNodes, ConNoToNoLi = EvNodeC2D( NodeList)  # build list of node coordinates for 2D/3D continuum nodes only and build list to map them to NodeList indices
            if len(ContinuumNodes) > 0:
                CoorTree = spatial.cKDTree( ContinuumNodes)  # for search purposes, e.g. for EFG or aggregates or embedded truss elements
            else:
                CoorTree = None
            for el in ElList:
                if el.Type in ["T2D2I", "T2D3I", "T3D2I", "T3D3I", "B23I", "B23EI"]:
                    el.CreateBond(ElList, NodeList, CoorTree, ConNoToNoLi, MatList, NoLabToInd,NoIndToCMInd)  # determine continuum elements connected via bond
            Echo(f"ConFemMat: {Name:s}", f6)
            MatLi = MatList[MName.upper()]
            Found = False
            for Elem in ElList:
                if Elem.Set.upper() == ElSet.upper() and Elem.MatN.upper() == MName.upper():
                    Found = True
                    Echo(f"Found specified element set: {ElSet:s}\nFound specified material: {MName:s} of type: {MatLi.Type:s}", f6)
                    Elem.Ini2(NodeList, NoIndToCMInd, MatList, SecDict)
                    if isinstance(MatLi, RCBeam) or isinstance(MatLi, WraTCBeam):
                        M_ = MatLi.PropMat[0]
                        MatC = ConcLaw1("Conc", [M_[3], 0., M_[0], M_[1], M_[2]])
                        M_ = MatLi.PropMat[1]
                        MatR = ReinfLaw("Steel", [-M_[4], M_[4], M_[0], M_[2], M_[3]])
                        if Elem.CrossSecType=="POLYLINE":
                            raise NameError("ConStressStrain: polygonal cross section not yet implemented ")
#                                zc1, zc3, hh = Elem.Geom[1, 1], Elem.Geom[1, 2], Elem.Geom[1, 2] - Elem.Geom[1, 1]
                        else:
                            zc1, zc3, hh = -Elem.Geom[1, 2] / 2., Elem.Geom[1, 2] / 2., Elem.Geom[1, 2]  # coordinate of bottom fibre, top fibre, height of cross section
                            bb = Elem.Geom[1, 1]
                            CrC = CrossSection([ hh, bb])
                            CrR_= []
                            nR = Elem.Geom.shape[0] - 2  # number of reinforcement layers
                            for i in range(nR):  # reinforcement contribution
#                                    ipL = (1 + nR) * ipI + (1 + i)  # every IP has state variables for concrete and nR reinforcement layers
                                As = Elem.Geom[i + 2, 0]  # reinforcement cross section
                                ys = Elem.Geom[i + 2, 1]  # coordinates of reinforcement
                                CrR_ += [[As,ys]]
                                CrR = ReinfCrossSection(CrR_)
                    break
            if not Found: raise NameError("ConMat: Specified material or element set not found in data set", MName, ElSet)
        self.UniaxTest(MatC)
        self.BendingTest(f6, MatC, CrC, MatR, CrR, NormalFList)
        #    X.BendingTest( MatA,CrossSection([0.2, 0.2]), MatR,None, [0.], kapStart=-0.15)
        Echo("finished", f6)
        f6.close()
        plt.show()

if __name__ == "__main__":
    Name, MatName, ElSet = "../DataExamples/E04/E4-01", 'mat1', 'prop1'
    NormalFList = [0., -1., -2.]                                            # normal forces

#    X = MatTest( 'log.txt' )
#    X.Run(None, MatName, ElSet, NormalFList)

    X = MatTest( Name + ".MatTester.txt")
    X.Run( Name, MatName, ElSet, NormalFList)


