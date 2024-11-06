# ConFem -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import process_time
from scipy import sparse
from scipy.sparse.linalg import aslinearoperator
from os import path as pth
import os
import scipy.spatial as spatial
import pickle

from ConFemBasics import *
import ConFemMat
from ConFemElem import *
from ConFemSteps import *
from ConFemInOut import *
from ConFemData import *
try:
    import LinAlg2
    LinAlgFlag = True
except ImportError:
    LinAlgFlag = False

def WriteXData(TimeX, Time, StepFinishedFlag):
    if Time >= TimeX or StepFinishedFlag:
        return True
    else:
        return False

from ConFemBasics import SamplePoints, _Truss2DAll_, _TriNodes_,_Bond2D_, _Length2D_, _FourNodes2D_,_Shells3D_,_BeamsAll_,__Length1D__

class ConExpliFem:
    def __init__(self):
        pass
    def TimeStep(self, SecDic, MatList, ElList, ScaleMass,RaAlph,RaBeta, NLGeom, ff):
        Echo(f"Critical time steps -- mass scaling {ScaleMass:.2e}", ff)
        SecDT =  {}
        for S in SecDic:
            Material = MatList[SecDic[S].Mat]
            Material.Density = Material.Density * ScaleMass
            rho = Material.Density                                          # apply before mass matrix computation
            if rho < ZeroD:
                Echo(f"Material {SecDic[S].Mat:s} with zero mass",ff)
                continue                                                    # for bond elements -- there are enough sections to contribute to mass, hopefully
#                raise NameError("ConExpliFem::TimeStep: 10")
            if isinstance(Material, ConFemMat.Elastic) \
                    or isinstance(Material, ConFemMat.IsoDamage) \
                    or isinstance(Material, ConFemMat.WraMisesReMem) \
                    or isinstance(Material, ConFemMat.ElasticR1D) \
                    or isinstance(Material, ConFemMat.ElasticC1D) \
                    or isinstance(Material, ConFemMat.MicroPlaneDam) \
                    or isinstance(Material, ConFemMat.MisesBeam2D):
                Emod = Material.PropMat[0]                                  # Emod
            else:
                raise NameError("ConExpliFem::TimeStep: 20", Material.Type)
            cc = sqrt(Emod / rho)
            dTMin = 1.0e6
            lamMax = 0.
            for e in SecDic[S].Elems:
                el = ElList[e]
                # following is estimation
                Lc = el.Lch_
                dT = Lc / cc
                if dT < dTMin: dTMin = dT
                # element matrices
                nint = el.nInt  # integration order
                InT  = el.IntT  # integration type of element
                mmat = np.zeros((el.DofEini, el.DofEini), dtype=np.double)  # element mass
                mvec = np.zeros((el.DofEini),             dtype=np.double)  # element mass into vector
                kmat = np.zeros((el.DofEini, el.DofEini), dtype=np.double)  # element stiffness
                CalcType = 11
                for j in range(el.nIntL):                                   # build element matrices with integration loop
#                    r = SamplePoints[InT, nint - 1, j][0]
#                    s = SamplePoints[InT, nint - 1, j][1]
#                    t = SamplePoints[InT, nint - 1, j][2]
#                    f = el.Geom[0, 0] * SampleWeight[InT, nint - 1, j]      # weighting factor

                    r, s, t, f = RSTf( CalcType,el, InT,nint, j)

                    if False: #ConFemElemCFlag and el.Type in ['SH4']:
                        B = np.zeros((6, 20), dtype=float)
                        Data_ = np.zeros((1), dtype=float)
                        TM = np.zeros((6, 6), dtype=float)
                        _ = SH4FormBC(r, s, t, B, Data_, el.XX, el.a, el.Vn, el.EdgeDir, el.gg[0], el.gg[1], TM)
                        det = Data_[0]
                    elif False: #ConFemElemCFlag and el.Type in ['CPS4','CPE4']:
                        B = np.zeros((3, 8), dtype=float)
                        det_ = np.zeros((1), dtype=float)
                        Data_ = np.zeros((1), dtype=float)
                        _ = CPS4FormBC(el.X0, el.Y0, el.X1, el.Y1, el.X2, el.Y2, el.X3, el.Y3, r, s, t, B, det_, Data_)
                        det = det_[0]
                    else:
                        B, det, TM = el.FormB(r, s, t, NLGeom)              # shape function derivative, alternatively Geom[1][0] for det
                    N  = el.FormN(r, s, t)  # shape function
                    mm = MatList[SecDic[S].Mat].Mass(el)                    # element mass in integration point
#                    factor = ScaleMass*det*f*el.Geom[1][0]
                    f_ =           det*f*el.Geom[1][0]
                    neps = B.shape[0]
                    dpsI, epsI = np.zeros((neps), dtype=float), np.zeros((neps), dtype=float)
                    if el.Type in _BeamsAll_:
                        tempI, dtmpI = np.zeros((2), dtype=float), np.zeros((2), dtype=float)
                    else:
                        tempI, dtmpI = 0., 0.
                    _, C, _ = Material.Sig(ff, 2, 1.0, None, j, el, dpsI, epsI, dtmpI, tempI, []) # CalcType --> 2
                    if False: #ConFemElemCFlag and el.Type not in ["B23E","B23EI","BAX23E","BAX23EI"]:
                        X = np.zeros((1, 1), dtype=float)  # currently dummy
                        Data_ = np.zeros((1), dtype=float)
                        rc = BTxCxB(mmat, f_, N, mm, Data_, X)
                        if rc != 0: raise NameError("ConFemExpliFem::TimeStep:NTxmxN", rc)
                        rc = BTxCxB(kmat, f_, B, C,  Data_, X)
                        if rc != 0: raise NameError("ConFemExpliFem::TimeStep:BTxCxB", rc)
                    else:
                        kmat = kmat + f_ * np.dot(np.transpose(B), np.dot(C, B))  # element stiffness
                        IPMassMatrix( el, mm, j,r,s,t,f, mmat, True)        # element mass matrix from ConFemBasics

                for i in range(el.DofEini):                                 # contract mass matrix to vector / diagonal matrix
                    for j in mmat[i]:
                        mvec[i] += ScaleMass * j
                # determine largest eigenvalue of element
                # https://www.geeksforgeeks.org/power-method-determine-largest-eigenvalue-and-eigenvector-in-python/
                x = np.random.rand(kmat.shape[1]).T
                tol, max_iter, lam_prev = 1e-6, 100, 0
                if min(np.absolute(mvec)) < 1.0e-12:  raise NameError("ConFemExpliFem::TimeStep: division by", min(np.absolute(mvec)))
                for i in range(max_iter):                                   # iterate for largest ev of element
                    y  = kmat @ x
                    y  = np.divide( y, mvec)
                    yN = np.linalg.norm(y)
                    x  = y / yN
                    lam = (x.T @ kmat @ x) / (x.T @ np.multiply( mvec, x))  # lam_ is eigenvalue of K * u = lam_ M u
                    if np.abs(lam - lam_prev) < tol: break
                    lam_prev = lam
                if lam > lamMax: lamMax = lam                               # check for magnitude within set
            # end of element loop
            omMax = sqrt(lamMax)                                            # max eigen angular velocity of elements within set
            xi = 0.5*( RaAlph/omMax ) # + RaBeta*omMax)                     # damping ratio
            dtMin = 2./omMax * (sqrt(1.+xi*xi)-xi)                          # critical time step for set
            y = lamMax * dtMin*dtMin / (4.-2*RaAlph*dtMin)                  # scaling factor for mass -- should be 1 here
            Echo(f"{S:s}, {Material.Type:s}, {Emod:.2e}, {rho:.2e}, {cc:.2e}, {Lc:.2e}, {dTMin:.4e} -- {dtMin:.6e} -- {RaAlph:.2e},{xi:.2e} - {y:f}  ",ff)
            SecDT[S] = [ lamMax, None, RaAlph, dtMin]                       # dictionary
        #
        dTMinAll = []
        for S in SecDT:
            S_ = SecDT[S]
            dTMinAll += [S_[3]]
        if ScaleMass>1.:
            dT = max(dTMinAll)
            for S in SecDT:
                S_ = SecDT[S]
                y = S_[0] * dT*dT / (4.-2.*S_[2]*dT)                        # scaling factor for mass
                Material = MatList[SecDic[S].Mat]
                Material.Density = Material.Density * y
                Echo(f"Superposed section mass scaling {S:s}: {y:.2f}", ff)
        else:
            dT = min(dTMinAll)
        Echo(f"derived time step {dT:.4e}", ff)
        return dT

#    @profile
    def Run(self, Name, LogName, PloF, LinAlgFlag, Restart, ResType, EigenvalPar, ElPlotTimes, StressStrainOut, VTK):
        StressStrainOutNames = SelOutIni( Name, ".elemout.", StressStrainOut) # initialize files for single element output
        if Restart:
            fd = open(Name+'.pkl', 'rb')
            # has to be in sync with pickle.dumo 
            NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);ElsticLTFlag=pickle.load(fd);\
                VecU0=pickle.load(fd);VecU1=pickle.load(fd);VecI1=pickle.load(fd);VecP1=pickle.load(fd);VecPN=pickle.load(fd);VecP0=pickle.load(fd);VecBm=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
                VecA=pickle.load(fd);VecV0=pickle.load(fd);VecAm=pickle.load(fd);VecVm=pickle.load(fd);VecD=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
                TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);StepCounter=pickle.load(fd);                     Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);\
                NoLabToNoInd=pickle.load(fd);NoIndToCMInd=pickle.load(fd);ContinuumNodes=pickle.load(fd);ConNoToNoLi=pickle.load(fd);SecDic=pickle.load(fd);LinAlgFlag=pickle.load(fd);ResultTypes=pickle.load(fd);Header=pickle.load(fd);
            fd.close()
            f6=open( Name+".protocol.txt", 'a')                             #
            Echo(f"\nConFem restart {Name:s}", f6)
            f1=open( Name+".in.txt", 'r')
            MatList, StepList = DataInput(f1, f6, Restart)              # read input file 
            f1.close()
            f2=open( Name+".elemout.txt", 'a')                              #
            f3=open( Name+".nodeout.txt", 'a')                              #
            f5=open( Name+".timeout.txt", 'a')                              #
            NameElemOut_ = Name+".elemout_.txt"
            NameNodeOut_ = Name+".nodeout_.txt"
            f7, MaxType =None, None
            ElemDataAll = {}                                                # skip this in case it is retrieved from pickle !!!
        else:
            DirName, FilName = DirFilFromName(Name)
            try: f1=open( Name+".in.txt", 'r')
            except: raise NameError(Name,"not found")
            f6=open( Name+".protocol.txt", 'w')
            Echo(f"ConFem  {Name:s}", f6)
            NodeList, ElList, MatList, SecDic, NoLabToNoInd, StepList, ResultTypes, Header = DataInput( f1, f6, Restart)
            f1.close()
            NoIndToCMInd = [i for i in range(len(NodeList))]                # maps identical before Cuthill McKee
            EvNodeEl( NodeList, ElList, NoLabToNoInd)                       # determine element indices belonging to each node and more - update required for Cuthill McKee
            ContinuumNodes, ConNoToNoLi = EvNodeC2D( NodeList )             # build list of node coordinates for 2D/3D continuum nodes only and build list to map them to NodeList indices
            if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes )   # for search purposes, e.g. for EFG or aggregates or embedded truss elements
            else:                     CoorTree = None
            for el in ElList:
                if el.Type in ["T2D2I","T2D3I","T3D2I","T3D3I","TAX2I","TAX3I","B23I","B23EI","BAX23I","BAX23EI","BAX21I","BAX21EI"]:
                    el.CreateBond( ElList, NodeList, SecDic, CoorTree,ConNoToNoLi, MatList, NoLabToNoInd, NoIndToCMInd)    # determine continuum elements connected via bond
            EvNodeEl( NodeList, ElList, NoLabToNoInd)                       # determine element indices belonging to each node after element creation - required for Cuthill McKee - assign type to continuum nodes

            if LinAlgFlag:
#                NoIndToCMInd = WithoutBandwidthOpti( NodeList )  # NoIndToCMInd#
#                NoIndToCMInd = CuthillMckee(NodeList, ElList)       # to reduce the skyline
                NoIndToCMInd = Sloan(NodeList, ElList)                      # to reduce the skyline
                NodeList.sort(key=lambda t: t.CMIndex)
            N, Skyline, SDiag, SLen, ElemDataAll = AssignGlobalDof( NodeList, ElList, MatList, SecDic, NoIndToCMInd) # assign degrees of freedom (dof) to nodes and elements -> see above
            f2=open( Name+".elemout.txt", 'w')
            NameElemOut_ = Name+".elemout_.txt"
            NameNodeOut_ = Name+".nodeout_.txt"
            if pth.exists(NameElemOut_): os.remove(NameElemOut_)
            if pth.exists(NameNodeOut_): os.remove(NameNodeOut_)
            f3=open( Name+".nodeout.txt", 'w')
            WrNodes, LineS = None, None
            MaxType = [] 
            f5, f7 = None, None
            if pth.isfile(Name+".opt.txt"):                                 # read options file if there is any
                f4=open( Name+".opt.txt", 'r')
                WrNodes, LineS, _, MaxType, _ = ReadOptionsFile(f4, NodeList,NoLabToNoInd,NoIndToCMInd)
                f4.close()
                f5=open( Name+".timeout.txt", 'w')
                if len(MaxType)>0: f7=open( Name+".elemmax.txt", 'w')
            for i in list(MatList.values()):                                # check, whether particular material types are in system
                ElsticLTFlag = ( isinstance(i,ConFemMat.ElasticLT) or isinstance( i.Conc, ConFemMat.ElasticLT)) and i.Used
                if ElsticLTFlag: break

            # Initializations
            VecU0 = zeros((N), dtype=double)                                # displacement vector time t
            VecU1 = zeros((N), dtype=double)                                # displacement vector time t+1
            VecVm = zeros((N), dtype=double)                                # velocity vector time t-1
            VecV0 = zeros((N), dtype=double)                                # velocity vector time t
            VecA  = zeros((N), dtype=double)                                # acceleration vector time t
#            VecI1 = zeros((N), dtype=double)                                # internal nodal forces vector time t+1
            VecI0 = zeros((N), dtype=double)                                # internal nodal forces vector time t
            VecP1 = zeros((N), dtype=double)                                # load vector time t+1
            VecP0 = zeros((N), dtype=double)                                # load vector time t
            VecPN = zeros((N), dtype=double)                                # nominal load vector for step time target  -- for compatibility only
            VecT  = zeros((N), dtype=double)                                # current temperatures vector - for compatibility only
            VecS  = zeros((N), dtype=double)                                # temperatures vector of previous time step  - for compatibility only
            BCIn  = ones((N),  dtype=int)                                   # indices for dofs with prescribed displacements
            BCIi  = zeros((N), dtype=int)                                   # indices for dofs with prescribed displacements --> 1, --> 0 otherwise
#            VecB  = zeros((N),dtype=double)                                 # reaction forces
            VecBm = zeros((N),dtype=double)                                 # previous time step t-1 reaction forces of - presumably not used, for compatibility
            VecAm = zeros((N),dtype=double)                                 # previous time step t-1 accelerations - presumably not used, for compatibility
            Time  = 0.
            TimeOld = 0.
            StepCounter = 0                                                        # step counter
            SymSys = IsSymSys( MatList )                                    # if there is at least one un-symmetric material the whole system is un-symmetric

        Eint  = 0
        Eext  = 0
        VecD  = zeros((N), dtype=double)
        VecV0_= zeros((N), dtype=double)
        VecVm_= zeros((N), dtype=double)
        VecV1 = zeros((N), dtype=double)
        VecR  = zeros((N), dtype=double)
        VecM  = zeros((N), dtype=double)                                    # lumped masses
        VecMI = zeros((N), dtype=double)                                    # lumped masses inverted
        VecD1 = zeros((N), dtype=double)                                    # nodal damping forces vector
        VecD0 = zeros((N), dtype=double)                                    # nodal damping forces vector
        VecI1 = zeros((N), dtype=double)  # internal nodal forces vector time t+1
        VecI0 = zeros((N), dtype=double)  # internal nodal forces vector time t
        VecB = zeros((N), dtype=double)  # reaction forces
        Node.ActiveNodes = 0
        for no in NodeList: 
            if len(no.NodeEl)>0: 
                Node.ActiveNodes +=1
                no.used = True
        EchoSys( f6, SymSys, LinAlgFlag, ConFemMatCFlag, ConFemElemCFlag, N, ElList, Node.ActiveNodes, SLen, MatList, Element, ElemDataAll) # ElemDataAll will currently not work with restart !!!

        # Calculations
        stime = process_time()
        while StepCounter<len(StepList):
            StLi = StepList[StepCounter]
            Echo(f"{StepCounter:d} step starts, solution type {StLi.SolType:s}", f6)
            dt = self.TimeStep(SecDic, MatList, ElList, StLi.ScaleMass,StLi.RaAlph,StLi.RaBeta, StLi.NLGeom,f6)
            StLi.BoundCondCheck( NodeList, f6)
            TimeTarg = StLi.TimeTarg                                        # time target for step
            if StepCounter>0: TimeS = StepList[StepCounter-1].TimeTarg      # time target of previous step
            else:      TimeS = 0.
            TimeEl = TimeS                                                  # base for elemout times
            TimeNo = TimeS                                                  # base for nodoout times
            TimeRe = TimeS
            StLi.current = StepCounter                                      # used in ConFemSteps
            MatM = sparse.lil_matrix((N, N))                                # sparse mass matrix initialization
            sysMass = IntForces( MatList, ElList, Time-TimeOld, VecU1,VecU0,VecU0,VecS,VecT,VecI1, MatM,None,None,None,None,       11, f6, StLi.NLGeom, SymSys,False,0) # mass matrix
            if StLi.Damp:
                MatD = StLi.RaAlph*MatM                                     # damping matrix
            for ii in range(N):
                for jj in range(N):                                         # mass vector
                    VecM[ii]=VecM[ii]+MatM[ii,jj]
            Echo(f"System mass {sysMass:f}", f6)
            TotalMass = 0
            for ii in range(N):
                TotalMass = TotalMass + VecM[ii]
                if VecM[ii] > ZeroD:
                    VecMI[ii] = 1. / VecM[ii]
                else:
                    raise NameError("Zero mass detected",ii,VecM[ii])
            if Time>TimeTarg-1.e-6: StepFinishedFlag = True    # step time target
            else:                   StepFinishedFlag = False
            StLi.BoundOffset( NodeList,NoLabToNoInd,NoIndToCMInd, VecU0)                               # add offset for prescribed displacement from current displacement for OPT=ADD
            
            IncCounter = 0
            while not StepFinishedFlag:                                     # loop over time increments
                IncCounter += 1
                TimeOld = Time
                if Time>=TimeEl and len(StLi.ElFilList)>0: TimeEl=Time + StLi.ElFilList[-1]# set time for element output
                if Time>=TimeNo and len(StLi.NoFilList)>0: TimeNo=Time + StLi.NoFilList[-1]# set time for element output
                if Time>=TimeRe:                           TimeRe=Time + StLi.ReFilList[-1]     # update time time for restart output
                Time = Time + dt                                            # new time for dynamic calculation

                VecV0_[:]= VecV0[:] + 0.5 * dt * VecA[:]                    # 1st update velocity
                VecD[:]  = dt * VecV0_[:]                                   # displacement increment vector
                VecU1[:] = VecU0[:] + VecD[:]
                VecI0[:] = VecI1[:]                                         # store current internal forces for energy before update
                VecP0[:] = VecP1[:]                                         # store current external forces for energy before update
                StLi.NodalLoads(N,Time, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecP1, VecPN)  # introduce concentrated loads into system -- Time = 0 initially
                StLi.ElementLoads(Time, TimeTarg, ElList, VecP1, VecPN)     # introduce distributed loads into system
                StLi.BoundCond(N, Time, 0., TimeTarg, NodeList, VecU1, VecI1, VecP1, VecPN, BCIn, BCIi, None,[],None,None,None, 2, None,f6 )  # presumably only VecP1, BCIn needed in the follwoing
                for jj in range(N):                                         # prescribe new displacements in case of boundary dofs
                    if BCIn[jj] == 0:
                        VecU1[jj] = VecP1[jj]
                if ElsticLTFlag:
                    IntForces( MatList, ElList, Time-TimeOld, VecU0,VecU1,VecD, VecS,VecT, VecI1, None,None,None,None,None, 0,f6,StLi.NLGeom, None, None, 0)  # check system state for certain materials
                IntForces(     MatList, ElList, Time-TimeOld, VecU0,VecU1,VecD, VecS,VecT, VecI1, None,None,None,None,None, 1,f6,StLi.NLGeom, None, None, 0)  # internal nodal forces / stiffness matrix
                VecB[:] = VecI1[:]                  # has to be cleared ?
                VecR[:] = VecP1[:] - VecI1[:]                               # residual vector
                if StLi.Damp:
                    VecD0[:] = VecD1[:]
                    VecVm_[:]= 0.5 * (VecV0[:] + VecVm[:])                          # velocity time t-1/2
                    VecD1 = aslinearoperator(MatD).matvec(VecVm_)           #                    VecD1 = aslinearoperator(MatD).matvec(VecV0)
                    VecR[:] = VecR[:] - VecD1[:]                                     # residual vector with damping
                for jj in range(N):
                    if BCIn[jj] == 0: VecR[jj] = 0                          # clear residuum for boundary conditions
                VecA = VecMI * VecR                                         # acceleration
                VecV1[:] = VecV0_[:] + 0.5 * dt * VecA[:]                            # 2nd update velocity

                Ekin = 0.5 * dot(VecV1, (VecM * VecV1))
                Eint = Eint + 0.5 * dot(VecI1 + VecI0, VecD)
                Eint = Eint + 0.5 * dot(VecD1 + VecD0, VecD)                # for damping
                Eext = Eext + 0.5 * dot(VecP0 + VecP1, VecD)
                Etot = Ekin + Eint - Eext
                Eref = max(Ekin, Eint, Eext)
                if Eref > ZeroD: Etot = Etot / Eref
                Resi = norm(VecR)                                           # residual norm
                Echo(f"{StepCounter:d} {IncCounter:3d} {TimeTarg:.3f} {Time:.5f} {Resi:e} Eint {Eint:e} Ekin {Ekin:e} Eext {Eext:e} Etot {Etot:e}", f6)
                VecVm[:] = VecV0[:]                                         # store current velocities as previous for next time step
                VecV0[:] = VecV1[:]                                         # store ahead velocities as current for next time step
                VecU0[:] = VecU1[:]                                         # store ahead displacements as current for next time step
                FinishEquilibIteration( MatList, ElList, NodeList,NoIndToCMInd, f6, StLi.NLGeom, True)# update state variables etc.
                SelOutWrite( "elem", StressStrainOut, StressStrainOutNames, Time, ElList, None, None, None, StepCounter, IncCounter)   # output for single elements over all time steps
                
                for el in ElList:
                    if el.Type in ["CPE4S","CPE3S","CPS4S","CPS3S"] and el.Active and el.NoCracks==1: 
                        el.CreateSDA2ndCrack( MatList, NodeList, NoIndToCMInd, f6, VecU0,VecD, dt)  # must precede CreateSDARankine and not in a common ElemList loop
                for el in ElList:
                    if el.RegType==3 and el.Type in ["CPS4","CPS3", "CPE4","CPE3","C3D8","CPS6","CPE6"] and el.Active:
                        nea = len(ElList)  # actual number of all elements (not only active)
                        nea = el.CreateSDARankine( nea, MatList, NodeList, ElList, SecDic, NoLabToNoInd, NoIndToCMInd, f6, VecU0,VecD, dt )
                
                if Time > TimeTarg-1.e-6:                                   # step time target
                    StepFinishedFlag = True                                       # time target reached, finalize computation for step
#                def WriteXData(TimeX):
#                    if Time+1.e-6>=TimeX or StepFinishedFlag: return True
#                    else:                                     return False
                if WriteXData( TimeEl,Time,StepFinishedFlag):
                    DataOutStress(NameElemOut_, ElList, NodeList, NoIndToCMInd, Time, "a", f6)
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    DataOut(NameNodeOut_, NodeList, VecU0, VecB, VecR, Time, "a" )
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                    WriteElemData( f2, f7, Time, ElList, NodeList,NoIndToCMInd, MatList, MaxType, ResultTypes)# write element data
                    f2.flush()
                    Echo(f"Element Data written {Time:f} {TimeEl:f}", f6)
                if WriteXData( TimeNo,Time,StepFinishedFlag):
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    WriteNodalData( f3, Time, NodeList, VecU0, VecB)         # write nodal data
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                    f3.flush()
                    Echo(f"Nodal Data written {Time:f} {TimeNo:f}", f6)
                if WriteXData( TimeRe,Time,StepFinishedFlag):
                    fd = open(Name+'.pkl', 'wb')                            # Serialize data and store for restart
                    pickle.dump(NodeList,fd);pickle.dump(ElList,fd);pickle.dump(MatList,fd);pickle.dump(StepList,fd);pickle.dump(N,fd);pickle.dump(WrNodes,fd);pickle.dump(LineS,fd);pickle.dump(ElsticLTFlag,fd);\
                        pickle.dump(VecU0,fd);pickle.dump(VecU1,fd);pickle.dump(VecI1,fd);pickle.dump(VecP1,fd);pickle.dump(VecPN,fd);pickle.dump(VecP0,fd);pickle.dump(VecBm,fd);pickle.dump(VecT,fd);pickle.dump(VecS,fd);pickle.dump(VecA,fd);pickle.dump(VecV0,fd);pickle.dump(VecAm,fd);pickle.dump(VecVm,fd);pickle.dump(VecD,fd);\
                        pickle.dump(BCIn,fd);pickle.dump(BCIi,fd);pickle.dump(Time,fd);pickle.dump(TimeOld,fd);pickle.dump(TimeEl,fd);pickle.dump(TimeNo,fd);pickle.dump(TimeS,fd);pickle.dump(StepCounter,fd);                     pickle.dump(Skyline,fd);pickle.dump(SDiag,fd);\
                    pickle.dump(SLen,fd);pickle.dump(SymSys,fd);pickle.dump(NoLabToNoInd,fd);pickle.dump(NoIndToCMInd,fd);pickle.dump(ContinuumNodes,fd);pickle.dump(ConNoToNoLi,fd);pickle.dump(SecDic,fd);pickle.dump(LinAlgFlag,fd);pickle.dump(ResultTypes,fd);pickle.dump(Header,fd);
                    fd.close()
                    Echo(f"Restart data written {Time:f}", f6)
                if f5!=None:
                    WriteNodes( f5, WrNodes, Time, VecU0, VecB, VecP1, StepCounter, IncCounter)
                    f5.flush()
            # end of time increment loop for current step
            StepCounter += 1                                                       # next step
        # end of step loop 

        Echo(f"total comp time {process_time()-stime:.0f} seconds", f6)
        f2.close()
        f3.close()
        if f5!=None: f5.close()
        f6.close()
        if f7!=None: f7.close()
#        fX.close()
        RC = FinishAllStuff(PloF, DirName, FilName, Name, ResType, VTK)
        return RC

if __name__ == "__main__":
    LogName="../LogFiles"                                                   # to log temporary data
#    numpy.seterr(all='raise')
    Name, Plot, Restart, HashSource, Eig, ElPlotTimes, StressStrainOut, VTK = DefData()
#    Name=str(raw_input('Filename without extension: '))                    # input data name
    ConFem_ = ConExpliFem()
    RC=ConFem_.Run(Name, LogName, Plot, LinAlgFlag, Restart, HashSource, Eig, ElPlotTimes, StressStrainOut, VTK)
    print(RC)
