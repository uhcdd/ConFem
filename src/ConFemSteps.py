# ConFemSteps -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from numpy import transpose, sqrt, zeros, array, double, dot, cos, ones

from ConFemBasics import FindIndexByLabel, SamplePoints, SampleWeight
from ConFemMat import MisesUniaxial
from ConFemBasics import Echo

class Boundary(object):
    def __init__(self, NodeLabel, Dof, Val, NoList, AmpLbl, AddVal):
        self.Dof = Dof
        self.Val = Val
        self.NodeLabel = NodeLabel
        self.Amplitude = AmpLbl
        self.AddVal    = AddVal
        self.ValOffset = 0.
class CLoad(object):
    def __init__(self, NodeLabel, Dof, Val, NoList, AmpLbl):
        self.Dof = Dof
        self.Val = Val
        self.NodeLabel = NodeLabel
        self.Amplitude = AmpLbl
class Temperature(object):
    def __init__(self, NodeLabel, Val, NoList, AmpLbl):
        self.Val = Val
        self.NodeLabel = NodeLabel
        self.Amplitude = AmpLbl
class PreStress(object):                                            # this is for general and geometry data -- initialization and use in later steps
    def __init__(self, Type, PrePar, GeomL, AmpLbl):
        self.For = PrePar[0]                                        # prestressing force
        self.Ap = PrePar[1]                                         # cross section of tendon
        self.GeomL = GeomL                                          # elements and geometry
        self.Amplitude = AmpLbl
        self.PreLength = 1.                                         # prestressing length
        self.PreLenZ = 1.                                           # length of prestressing after 1st step                      
        self.Reinf = MisesUniaxial( [PrePar[2],.2,PrePar[3],PrePar[4],PrePar[5],PrePar[6],PrePar[7],PrePar[8]], [0.,1.] )# initialization material for prestressing
        self.StateVar = zeros((1,3), dtype=float)                   # state variables for mises material + tendon stress during tensioning 
        self.StateVarN= zeros((1,3), dtype=float)                   # updated state variables
        self.TensStiff = False
        self.Geom = array([[0.],[1.]])
        self.dim = 1                                                # for compatibility reasons
        self.nint = 3
        if   Type == "POST-BONDED":   self.BondLess = False
        elif Type == "POST-UNBONDED": self.BondLess = True
        else: raise NameError ("Unknown type of prestressing")
        if not self.BondLess: self.PrStrain = zeros((len(GeomL),self.nint),dtype=float) # storage for strain of prestressed fiber - after tensioning 
class PreStreSt(object):                                            # prestressing data related to step
    def __init__(self, Name, AmpLbl):
        self.Name = Name                                            # prestressing force
        self.Amplitude = AmpLbl
class DLoad(object):
    def __init__(self, ElLabel, Dof, Val, AmpLbl):
        self.ElLabel = ElLabel
        self.Dof = Dof
        self.Val = Val
        self.Amplitude = AmpLbl
#class ElFile(object):
#    def __init__(self, OutTime):
#        self.OutTime = float(OutTime)                                      # time interval output
#class NoFile(object):
#    def __init__(self, OutTime):
#        self.OutTime = float(OutTime)                                      # time interval output
class Step(object):
    current = 0                                                     # initialization of current step
    PreSDict = {}                                                   # dictionary for all prestressings
    def __init__(self):
        self.IterTol  = 1.e-3                                       # default values
        self.IterNum  = 10
        self.TimeStep = 0.                                          # should be initialized within confem
        self.TimeTarg = 0.
        self.SolType  ="NR"                                         # NR: Newton Raphson, BFGS: BFGS
        self.BoundList = []                                         # boundary conditions
        self.CLoadList = []                                         # nodal loads
        self.DLoadList = []                                         # distributed loads
        self.TempList  = []                                         # nodal temperatures
        self.PrestList = []                                         # prestressing data with reference to step - defined in ConFemInOut:DataInput, used by Step.NodalPrestress
        self.ElFilList = []
        self.NoFilList = []
        self.ReFilList = [1.]                                       # in case no "restart file" defined in input
        self.AmpDict   = {}                                         # dictionary for Amplitudes
        self.AmpDict['Default'] = [[0.,0.],[1.,1.]]
        self.ZeroD     = 1.e-9                                      # smallest float for division by Zero
        self.Dyn       = False                                      # flag for implicit dynamics
        self.Eigenmode2TS = False                             # largest natural period to timstep
        self.Eigenmode2TSrel = 0.1                                  # default ratio time step to largest natural period if
        self.Buckl     = False                                      # flag for buckling analysis
        self.Eigenmodes = False                                     # flag for computation of dynamic eigenmodes
        self.EigenmodesN = 5                                        # !!!!! default for number of comuted eigemodes in case of Eigenmodes=True - also used for buckling
        self.EigenmodeIter = 100                                    # !!!!! max number of iterations for eigenmode iteration - also used for buckling
        self.Eigenvalues = zeros((10), dtype=double)                # array for eigenvalues -- hard coded length
        self.NMbeta    = 0.25                                       # Newmark parameter beta
        self.NMgamma   = 0.5                                        # Newmark parameter gamma
        self.Damp      = False                                      # flag for artificial / Rayleigh damping
        self.RaAlph    = 0.0                                        # default Rayleigh damping parameter with stiffness
        self.RaBeta    = 0.0                                        # default Rayleigh damping parameter with mass
        self.EigenVal2Beta = False                                  # derive beta (stiffness) from eigenmode
        self.ZetaBeta  = 0.0                                        # damping beta (stiffness) related to critical damping
        self.EigenVal2Alph = False                                  # derive alpha (mass) from eigenmode
        self.ZetaAlph  = 0.0                                        # damping alpah (mass) related to critical damping
        self.ArcLen    = False                                      # flag for arc length control
        self.ArcLenV   = 0                                          # arc length parameter
        self.NLGeom    = False                                      # flag for large deformations
        self.varTimeSteps= False                                    # flag for variable TimestepSizes
        self.ScaleMass = 1.0                                        # for explicit to manipulate critical time step
    def AmpVal(self, Time, Label):                                  # determine value for amplitude
        Data = self.AmpDict[Label]
        nD = len(Data)
        if Time>=Data[nD-1][0]:                                     # time exceeds last value
            DelT = Data[nD-1][0]-Data[nD-2][0]
            if DelT<self.ZeroD: raise NameError ("something wrong with amplitude 1",Label,Time,Data)
            return (Time-Data[nD-2][0])/DelT * (Data[nD-1][1]-Data[nD-2][1]) + Data[nD-2][1]
        elif Time<0:
            return 0.
#            raise NameError("Time less zero", Time)
        for i in range(len(Data)-1):
            if Data[i][0]<=Time and Time<=Data[i+1][0]:             # time within interval
                DelT = Data[i+1][0]-Data[i][0]
                if DelT<self.ZeroD: raise NameError ("something wrong with amplitude 2",Time,Data)
                return (Time-Data[i][0])/DelT * (Data[i+1][1]-Data[i][1]) + Data[i][1]
        raise NameError ("something wrong with amplitude 3",Time,Data)
    # !!! needs to be tested again
    def BoundOffset(self, NodeList,NoLabToNoInd,NoIndToCMInd, VecU):        # add offset for prescribed displacement from current displacement for OPT=ADD
        for i in self.BoundList:                                    # loop over all boundary conditions of step
            if i.AddVal:
#                nI = FindIndexByLabel( NodeList, i.NodeLabel)       # node index of bc
                nI = NoIndToCMInd[ NoLabToNoInd[i.NodeLabel] ]      # node index of bc
                Node = NodeList[nI]
                for j, jj in enumerate(Node.DofT):
                    if jj==i.Dof: break                             # dof type equals prescribed dof type 
                k = Node.GlobDofStart + j                           # constrained dof global index
                i.ValOffset = VecU[k]
        return 0
    def BoundCond(self, N, Time, TimeS, TimeTar, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, Kmat,KVecU,KVecL,Skyline,SDiag, CalcType, SymS, ff):# introduce loads and boundary conditions into system
        # P is actual load, P0 is target load for step
        for i in self.BoundList:                                    # loop over all boundary conditions of step
            val = self.AmpVal(Time, i.Amplitude)*i.Val + i.ValOffset  # prescribed value for actual time
            valT = self.AmpVal(TimeTar, i.Amplitude)*i.Val + i.ValOffset # prescribed value for final target
            valS = self.AmpVal(TimeS, i.Amplitude)*i.Val + i.ValOffset # prescribed value for step target
            valT = valT - valS
            nI = FindIndexByLabel( NodeList, i.NodeLabel)           # node index of bc
            if nI < 0: raise NameError ("ConFemSteps.Step.BoundCond: node",i.NodeLabel,"does not exist")
            no = NodeList[nI]
            for j, jj in enumerate(no.DofT):                                    # loop over all dofs of node
                if jj==i.Dof:                                       # dof type equals prescribed dof type 
                    break
            k = no.GlobDofStart + j                                 # constrained dof global index
            if CalcType==2:                                         # only for NR, not for MNR, BFGS
                if Kmat!=None:
                    for j in range(N):                              # loop over rows and columns simul
                        VecI[j] = VecI[j] - VecU[k]*Kmat[j,k]       # internal forces due to prescribed displacement
                        VecP[j] = VecP[j] - val*    Kmat[j,k]       # external forces to enforce prescribed displacement
                        VecP0[j] = VecP0[j] - valT* Kmat[j,k]       # nominal external forces corresponding to prescribed displacements
                        Kmat[j,k] = 0.                              # modification stiffness matrix
                        Kmat[k,j] = 0.
                    Kmat[k,k] = 1.                                  # modification stiffness matrix
                elif len(KVecU)>0:
                    jj = SDiag[k]                                   # same for stiffness vectors
                    # modify values connected with system matrix column
                    for j in range(Skyline[k]):
                        VecI[k-j]  = VecI[k-j]  - VecU[k]*KVecU[jj+j] # internal forces due to prescribed displacement
                        VecP[k-j]  = VecP[k-j]  - val    *KVecU[jj+j] # external forces to enforce prescribed displacement
                        VecP0[k-j] = VecP0[k-j] - valT   *KVecU[jj+j] # external nominal forces corresponding to prescribed displacements
                    for j in range(Skyline[k]): KVecU[jj+j] = 0.    # upper right part of system matrix 
                    if not SymS:
                        for j in range(Skyline[k]): KVecL[jj+j] = 0.# lower left part of system matrix
                    # modify values connected with system matrix row -- "diff" within loops unfortunately is not constant due to banded vector storage
                    if not SymS:                                    # unsymmetric system
                        for j in range(k,N):
                            diff = j-k
                            if diff<Skyline[j]:
                                valK = KVecL[SDiag[j]+diff]
                                VecI[j] = VecI[j] - VecU[k]*valK
                                VecP[j] = VecP[j] - val*    valK
                                VecP0[j] = VecP0[j] - valT* valK
                        for j in range(k,N):
                            diff = j-k
                            if diff<Skyline[j]:
                                KVecL[SDiag[j]+diff] = 0.           # stiffness off diagonal
                    else:                                           # use KVecU instead KVecL in case of symmetric system
                        for j in range(k,N):
                            diff = j-k
                            if diff<Skyline[j]:
                                valK = KVecU[SDiag[j]+diff]
                                VecI[j] = VecI[j] - VecU[k]*valK 
                                VecP[j] = VecP[j] - val*    valK 
                                VecP0[j] = VecP0[j] - valT* valK 
                    for j in range(k,N):
                        diff = j-k
                        if diff<Skyline[j]: KVecU[SDiag[j]+diff] = 0. # stiffness off diagonal
                    #
                    KVecU[jj] = 1.                                  # stiffness diagonal
            # k is index of constrained dof
            VecI[k] = VecU[k]                                       # modification internal forces with actual displacement
            VecP[k] = val                                           # modification actual load vector with prescribed displacement -- difference to VecI yields prescribed displacment increment
            VecP0[k] = valT                                         # modification step target load vector
            BCIn[k] = 0                                             # index 0 for prescribed displacement instead default 1
            BCIi[k] = 1                                             # index 1 for prescribed displacement instead default 0
        return 0
    def BoundCondCheck(self, NodeList, ff):                                 # echo for mismatch of boundary dof types
        for b in self.BoundList:
            Found = False
            i = FindIndexByLabel( NodeList, b.NodeLabel)                    # node index of bc
            if i < 0: raise NameError ("ConFemSteps.Step.BoundCond: node",b.NodeLabel,"does not exist")
            no = NodeList[i]
            for jj in no.DofT:                                              # loop over all dofs of node
                if jj==b.Dof:                                               # dof type equals prescribed dof type 
                    Found = True
                    break
            if not Found: 
                Echo(f"ConFemSteps.Step.BoundCond: dof type {b.Dof:d} missing correspondence for node {b.NodeLabel:d}, skipped this", ff) #
#                raise NameError("exit") 

    def NodalLoads(self, N, Time, TimeTar, NodeList,NoLabToNoInd,NoIndToCMInd, VecP, VecP0):  # introduce concentrated loads into system
        VecP[:]  = 0.
        VecP0[:] = 0.
        for i, CL in enumerate(self.CLoadList):                     # loop over all concentrated loads of ste
            val  = self.AmpVal(Time, CL.Amplitude)*CL.Val           # prescribed value
            valT = self.AmpVal(TimeTar, CL.Amplitude)*CL.Val        # prescribed value
#            nI = FindIndexByLabel( NodeList, CL.NodeLabel)          # node index of concentrated load
            nI = NoIndToCMInd[ NoLabToNoInd[CL.NodeLabel] ]                 # node index of bc
            for j, j0 in enumerate(NodeList[nI].DofT):
                if j0 == CL.Dof: break                              # dof type equals prescribed dof type
            k = NodeList[nI].GlobDofStart + j                       # loaded dof global index
            VecP[k]  = VecP[k] + val                                # global load vector
            VecP0[k] = VecP0[k] + valT                              # global load vector
        return 0
    def ElementLoads(self, Time, TimeTar, ElList, VecP, VecP0):# introduce distributed loads into system
        for i, DL in enumerate(self.DLoadList):                     # loop over all distributed/element loads of step
            Label= DL.ElLabel
            Dof  = DL.Dof
            Val  = self.AmpVal(Time, DL.Amplitude)*DL.Val           # prescribed value
            ValT = self.AmpVal(TimeTar, DL.Amplitude)*DL.Val        # prescribed value
            for j, Elem in enumerate(ElList):
                if Elem.Set==Label or str(Elem.Label)==Label:
                    nfie = Elem.nFie                                # number of loading degrees of freedom
                    if Elem.Rot: TraM = Elem.Trans[0:nfie,0:nfie]
                    ev  = zeros((nfie), dtype=double)               # array for load value
                    ev0 = zeros((nfie), dtype=double)
                    if Elem.Type=='SB3':
                        ev[1]  = Val
                        ev0[1] = ValT
                    else:
                        ev[Dof-1]  = Val                            # Dof: input index for load dof 
                        ev0[Dof-1] = ValT
                    if Elem.Rot:
                        ev = dot(TraM,ev)
                        ev0 = dot(TraM,ev0)
                    pv  = zeros(( Elem.DofE ),dtype=double)         # element nodal loads
                    pv0 = zeros(( Elem.DofE ),dtype=double)         # element nodal loads
                    InT  = Elem.IntT                                # integration type for element
                    nint = Elem.nInt                                # integration order
                    nintL= Elem.nIntLi                              # total number of integration points
                    for k in range(nintL):                          # integration point loop
                        r = SamplePoints[InT,nint-1,k][0]
                        s = SamplePoints[InT,nint-1,k][1]
                        t = SamplePoints[InT,nint-1,k][2]
                        if Elem.Type=='SH4': t=0

                        a = Elem.JacoD(r,s,t)
                        b = Elem.Geom[0,0]
                        c = SampleWeight[InT,nint-1,k]

                        f = Elem.JacoD(r,s,t)*Elem.Geom[0,0]*SampleWeight[InT,nint-1,k]# weighting factor, Elem.Geom[0,0] is a measure for length, area ...
                        N = Elem.FormN(r,s,t)                       # shape function
                        pv  = pv  + f*dot( N.T, ev)
                        pv0 = pv0 + f*dot( N.T, ev0)
                    if Elem.Rot:
                        pv  = dot(Elem.Trans.transpose(),pv)        # transform local to global forces (Script Eq. (3.68))
                        pv0 = dot(Elem.Trans.transpose(),pv0)       # transform local to global forces (Script Eq. (3.68))
                    ndof = 0
                    for k in range(Elem.nNod):                      # assemble
                        for l in range(Elem.DofN[k]):
                            ii = Elem.DofI[k,l]
                            VecP[ii] = VecP[ii] + pv[ndof+l] 
                            VecP0[ii] = VecP0[ii] + pv0[ndof+l] 
                        ndof = ndof + Elem.DofN[k]
        return 0
    def NodalPrestress(self, N, Time, ElList, VecP, VecU, NLg):  # introduce nodal prestress
        for i, PL in enumerate(self.PrestList):                     # loop over all prestressings of step
            PreStreD = self.PreSDict[PL.Name]
            if PreStreD.BondLess: pF = PreStreD.PreLength/PreStreD.PreLenZ # relative change of prestressing length
            else:                 pF = 1
            Val =  self.AmpVal(Time, PL.Amplitude) * PreStreD.For   # prescribed value
            Emod= 2.e5                                              # Young's modulus
            if self.current == 0: PreStreD.StateVarN[0,2] = Val/(PreStreD.Ap*Emod) # prestressing strain during tensioning
            GeomL = PreStreD.GeomL                                  # other prestressing geometry data
            LP = 0.                                                 # initialization arc length of prestressing
            for j, GL in enumerate(GeomL):                          # loop over all elements of prestressing
                Elem = ElList[GL[0]]                                # GL[0] -> element index
                uvec = zeros((Elem.DofE), dtype=double)             # element displacements
                ndof0 = 0
                for k in range(Elem.nNod):                          # loop over all element nodes to determine element displacements
                    for k1 in range(Elem.DofN[k]):                  # loop over all node dofs
                        uvec[ndof0+k1] = VecU[Elem.DofI[k,k1]]      # element displacements from global displacements
                    ndof0 = ndof0 + Elem.DofN[k]                    # update entry index for element displacements
                nfie = Elem.nFie                                    # number of loading degrees of freedom
                if Elem.Rot: TraM = Elem.Trans[0:nfie,0:nfie]
                pv =    zeros(( Elem.DofE ),dtype=double)           # element nodal loads
                InT  = Elem.IntT                                    # integration type for element
                nint = Elem.nInt                                    # integration order
                nintL= Elem.nIntLi                                  # total number of integration points
                for k in range(nintL):                              # integration point loop
                    r = SamplePoints[InT,nint-1,k][0]               # local coordinate
                    f =  Elem.Geom[0,0]*SampleWeight[InT,nint-1,k]  # weighting factor
                    B, jd, TM = Elem.FormB(r, 0, 0, NLg)            # shape function derivative
                    epsI = dot( B, uvec)                            # integration point generalized strains
                    P = Elem.FormP(r, 0)                            # form function for prestressing 
                    DatP1 = dot(P[0,  :],GL[1])                     # interpolation of tendon force
                    DatP2 = dot(P[1:3,:],GL[1])                     # interpolation of tendon geometry
                    DatP3 = dot(P[2,  :],uvec)                      # interpolation of inclination from deformation - used for arc deformed length
                    # prestressing with bond
                    if not PreStreD.BondLess:
                        epsP = epsI[0] - DatP2[0]*epsI[1]           # strain of beam fiber in height of tendon
                        if self.current==0: 
                            PreStreD.PrStrain[j,k] = epsP
                        else:
                            Sig, Dsig, dataL = PreStreD.Reinf.Sig(None, 1,0.,j,0,PreStreD, [0.,0], [epsP,0],0.,0.,None)#CalcType, Dt, i, ipL, Elem, [dep, 0], [eps, 0], dtp, tmp
                            pF = ( epsP - PreStreD.PrStrain[j,k] + PreStreD.StateVarN[0,2] )/PreStreD.StateVarN[0,2] # local change of prestressing length
                    # pF for bondless set above
                    sigP = pF*Val*DatP1*cos(DatP2[1]) * array([1,-DatP2[0]])# internal forces due to prestressing
                    pv = pv + f*dot( B.T, sigP)                     # element nodal forces due to prestressing
                    LP = LP + f*sqrt((1+epsI[0])**2+(DatP2[1]+DatP3)**2)# arc length of deformed prestressing
                if Elem.Rot: pv = dot(Elem.Trans.transpose(),pv)    # transform local to global forces
                ndof = 0
                for k in range(Elem.nNod):                          # assemble
                    for l in range(Elem.DofN[k]):
                        ii = Elem.DofI[k,l]
                        VecP[ii] = VecP[ii] - pv[ndof+l] 
                    ndof = ndof + Elem.DofN[k]
            # end of element loop
            PreStreD.PreLength = LP                                 # current tendon length        
            if self.current==0: PreStreD.PreLenZ=LP                 # tendon length upon application of prestressing in 1st step
        return 0

    def NodalTemp(self, N, Time, NodeList,NoLabToNoInd,NoIndToCMInd, VecT): # nodal temperatures
        VecT[:] = 0.
        for i, TL in enumerate(self.TempList):
#            nI = FindIndexByLabel( NodeList, TL.NodeLabel)          # node index of temperature
            nI = NoIndToCMInd[ NoLabToNoInd[TL.NodeLabel] ]                 # node index of bc
            k  = NodeList[nI].GlobDofStart                          # global index
            tL = TL.Val                                             # list with temperature values
            tA = self.AmpVal(Time, TL.Amplitude)                    # amplitude value
            for j in range(len(tL)): 
                VecT[k+j] = VecT[k+j] + tA*tL[j]                    # global temperature vector
        return len(self.TempList)                                   # to indicate whether there is temperature in system

    def MaskForArclen(self, NodeList, NoIndToCMInd, NoLabToNoInd, N): # mask for displacment arc length 
        if len(self.ArcLenNodes)>0:
            Mask = zeros((N), dtype=int)
            for no in self.ArcLenNodes:
                No  = NodeList[ NoIndToCMInd[NoLabToNoInd[no]] ]
                ind = No.GlobDofStart
                for k, k0 in enumerate(No.DofT):
                    if k0 != 7:                                     # not dof type gradient damage 
                        Mask[ind+k] = 1
        else:
            Mask = ones((N),dtype=int)
            for No in NodeList:
                ind = No.GlobDofStart
                for k, k0 in enumerate(No.DofT):
                    if k0 == 7:                                     # dof type gradient damage 
                        Mask[ind+k] = 0
        return Mask

