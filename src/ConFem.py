# ConFem -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import process_time
from numpy import ones
#from scipy.linalg import norm
from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import splu
from scipy.sparse.linalg import aslinearoperator
from os import path as pth
from collections import deque
import os
#from scipy.spatial import KDTree
import scipy.spatial as spatial
#import sys 
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

def CheckInz(ElemList, Label):
    for el in ElemList:
        if el.Label==Label:
            print("CheckInz ",el.Label,el.Type,el.Inzi)

class ConFem:
    def __init__(self):
        pass
#    @profile
    def Run(self, Name, LogName, PloF, LinAlgFlag, Restart, ResType, EigenvalPar, ElPlotTimes, StressStrainOut, VTK):
        StressStrainOutNames = SelOutIni( Name, ".elemout.", StressStrainOut) # initialize files for single element output
        DirName, FilName = DirFilFromName(Name)
        if Restart:
            fd = open(Name+'.pkl', 'rb')
            # has to be in sync with pickle.dumo 
            NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);ElasticLTFlag=pickle.load(fd);\
                VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
                VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
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
#            DirName, FilName = DirFilFromName( Name )
            try: f1=open( Name+".in.txt", 'r')
            except: raise NameError(Name,"not found")
            f6=open( Name+".protocol.txt", 'w')
            Echo(f"ConFem  {Name:s}", f6)
            NodeList, ElList, MatList, SecDic, NoLabToNoInd, StepList, ResultTypes, Header = DataInput( f1, f6, Restart)
            f1.close()
            NoIndToCMInd = [i for i in range(len(NodeList))]                # maps identical before Cuthill McKee
            EvNodeEl( NodeList, ElList, NoLabToNoInd)                       # determine element indices belonging to each node and more - update required for Cuthill McKee
            ContinuumNodes, ConNoToNoLi = EvNodeC2D( NodeList )             # build list of [[nodal coordinates]] for 2D/3D continuum nodes only and list to map them to NodeList indices
            if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes )   # for search purposes, e.g. for EFG or aggregates or embedded truss elements
            else:                     CoorTree = None
            for el in ElList:
                if el.Type in ["T2D2I","T2D3I","T3D2I","T3D3I","TAX2I","TAX3I","B23I","B23EI","BAX23I","BAX23EI","BAX21I","BAX21EI"]:
                    el.CreateBond( ElList,NodeList,SecDic, CoorTree,ConNoToNoLi, MatList, NoLabToNoInd,NoIndToCMInd)    # determine continuum elements connected via bond
            EvNodeEl( NodeList, ElList, NoLabToNoInd)                       # determine element indices belonging to each node after element creation - required for Cuthill McKee - assign type to continuum nodes

            if LinAlgFlag:
#                NoIndToCMInd = WithoutBandwidthOpti( NodeList )  # NoIndToCMInd#
#                NoIndToCMInd = CuthillMckee(NodeList, ElList)       # to reduce the skyline
                NoIndToCMInd = Sloan(NodeList, ElList)                      # to reduce the skyline - initial node index (not label) to index according to bandwidth optimization
                NodeList.sort(key=lambda t: t.CMIndex)                      # CMIndex holds new index
            N, Skyline, SDiag, SLen, ElemDataAll = AssignGlobalDof( NodeList, ElList, MatList, SecDic, NoIndToCMInd) # assign degrees of freedom (dof) to nodes and elements -> see above
            f2=open( Name+".elemout.txt", 'w')
            NameElemOut_ = Name+".elemout_.txt"
            NameNodeOut_ = Name+".nodeout_.txt"
            if pth.exists(NameElemOut_): os.remove(NameElemOut_)
            if pth.exists(NameNodeOut_): os.remove(NameNodeOut_)
            f3=open( Name+".nodeout.txt", 'w')
            WrNodes, LineS, MaxEquiIter_ = None, None, None
            MaxType = [] 
            f5, f7 = None, None
            if pth.isfile(Name+".opt.txt"):                                 # read options file if there is any
                f4=open( Name+".opt.txt", 'r')
                WrNodes, LineS, _, MaxType, MaxEquiIter_ = ReadOptionsFile(f4, NodeList,NoLabToNoInd,NoIndToCMInd)
                f4.close()
                f5=open( Name+".timeout.txt", 'w')
                if len(MaxType)>0: f7=open( Name+".elemmax.txt", 'w')
            for i in list(MatList.values()):                                # check, whether particular material types are in system
                ElasticLTFlag = ( isinstance(i,ConFemMat.ElasticLT) or isinstance( i.Conc, ConFemMat.ElasticLT)) and i.Used
                if ElasticLTFlag: break
#            fX = open( Name+".stuff.txt", 'w')
            # Initializations
            VecU = zeros((N),dtype=double)                                  # current displacement vector
            VecC = zeros((N),dtype=double)                                  # displacement vector of previous time increment after equilibrium iteration
            VecY = zeros((N),dtype=double)                                  # displacement increment of previous time increment after equilibrium iteration
            VevU = zeros((N),dtype=double)                                  # current velocities
            VeaU = zeros((N),dtype=double)                                  # current accelerations
            VevC = zeros((N),dtype=double)                                  # previous time step velocities
            VeaC = zeros((N),dtype=double)                                  # previous time step accelerations
            VecI = zeros((N),dtype=double)                                  # internal nodal forces vector
            VecP = zeros((N),dtype=double)                                  # load vector
            VecP0= zeros((N),dtype=double)                                  # nominal load vector
            VecP0old= zeros((N),dtype=double)                               # nominal load vector of previous calculation step
            VecBold= zeros((N),dtype=double)                                # reaction forces of previous step
            VecT = zeros((N),dtype=double)                                  # current temperatures vector
            VecS = zeros((N),dtype=double)                                  # temperatures vector of previous time step
            BCIn = ones((N),dtype=int)                                      # indices for dofs with prescribed displacements --> 0, --> 1 otherwise
            BCIi = zeros((N),dtype=int)                                     # indices for dofs with prescribed displacements --> 1, --> 0 otherwise
            Time = 0.
            TimeOld = 0.
            SymSys = IsSymSys( MatList )                                    # if there is at least one un-symmetric material the whole system is un-symmetric
            StepCounter = 0                                                 # step counter
        VecP0i= zeros((N),dtype=double)                                  # intermediate storage for nominal load vector of step
        VecUP = zeros((N),dtype=double)                                     # displacement vector of latest equilibrium iteration -- for line search
        VecP0_= zeros((N),dtype=double)                                     # 
        VecB = zeros((N), dtype=double)                                     # reaction forces
        Node.ActiveNodes = 0
        for no in NodeList: 
            if len(no.NodeEl)>0: 
                Node.ActiveNodes +=1
                no.used = True
        EchoSys( f6, SymSys, LinAlgFlag, ConFemMatCFlag, ConFemElemCFlag, N, ElList, Node.ActiveNodes, SLen, MatList, Element, ElemDataAll) # ElemDataAll will currently not work with restart !!! 

        # Calculations
        stime = process_time()
        while StepCounter<len(StepList):
            ndeq = 10
            if f5 != None: timeoutQueues  = [ deque(maxlen=ndeq), deque(maxlen=ndeq), deque(maxlen=ndeq) ]
            equiiterQueues                = [ deque(maxlen=ndeq), deque(maxlen=ndeq), deque(maxlen=ndeq) ]
            #
            StLi = StepList[StepCounter]
            Echo(f"{StepCounter:d} step starts, solution type {StLi.SolType:s}", f6)
            StLi.BoundCondCheck( NodeList, f6)                              # whether bc dof has counterpart in dof of node (->node.DofT)
            if StLi.ArcLen: Mask = StLi.MaskForArclen( NodeList, NoIndToCMInd,NoLabToNoInd, N)
            if StLi.varTimeSteps:
                TimeVarL = len(StLi.TimeTargVar)
                TimeTarg = StLi.TimeTargVar[-1]
                if StepCounter>0: TimeS = StepList[StepCounter-1].TimeTargVar[-1]         # time target of previous step
                else:   TimeS = 0.
            else:
                TimeTarg = StLi.TimeTarg                                    # time target for step
                if StepCounter>0: TimeS = StepList[StepCounter-1].TimeTarg  # time target of previous step
                else:   TimeS = 0.
            TS = TimeTarg-TimeS
            if not StLi.Buckl:
                if TS<ZeroD: raise NameError("ConFem: TimeTarg <= Times",StepCounter,TimeS,TimeTarg," -- desired?")
                else: DTS = 1./TS
            if len(StLi.ElFilList)>0: TimeEl = TimeS                        # TimeEl updated below with StLi.ElFilList[-1]
            else:                     TimeEl = TimeTarg
            if len(StLi.NoFilList)>0: TimeNo = TimeS
            else:                     TimeNo = TimeTarg                     #
            if len(StLi.ReFilList)>0: TimeRe = TimeS + StLi.ReFilList[-1]
            else:                     TimeRe = TimeTarg
            StLi.current = StepCounter                                      # used in ConFemSteps
            if StLi.Dyn:                                                    # preparation of dynamic calculation
                dt = StLi.TimeStep                                          # time step for dynamic calculation
                if LinAlgFlag:
                    if StLi.Damp:
                        DVecU = zeros(SLen, dtype=float)                    # Initialize upper right part of damping vector
                        if not SymSys: DVecL = zeros(SLen, dtype=float)     # Initialize lower left part of damping vector
                        else:          DVecL = None
                    MVecU = zeros(SLen, dtype=float)                        # Initialize upper right part of mass vector
                    if not SymSys:     MVecL = zeros(SLen, dtype=float)     # Initialize lower left part of mass vector
                    else:              MVecL = None
                    sysMass = IntForces( MatList, ElList, Time-TimeOld, VecC,VecU,VecU,VecS,VecT, VecI, None,None,MVecU,MVecL,SDiag, 10, f6, StLi.NLGeom, SymSys,False,0) # mass matrix
                else:
                    MatM = sparse.lil_matrix((N, N))                        # sparse mass matrix initialization
                    sysMass = IntForces( N, MatList, ElList, Time-TimeOld, VecC,VecU,VecU,VecS,VecT, VecI, MatM,None,None,None,None,       10, f6, StLi.NLGeom, SymSys,False,0) # mass matrix
                a_ = StLi.NMgamma/StLi.NMbeta                               # Newmark auxiliary value
                a0 = 1./(StLi.NMbeta*dt**2)                                 # Newmark auxiliary value
                a1 = a_/dt                                                  # Newmark auxiliary value for Rayleigh damping
                Echo(f"System mass {sysMass:f}", f6)
            elif StLi.Buckl: 
                MatM = sparse.lil_matrix((N, N))                            # used for geometric stiffness matrix
                LinAlgFlag = False
            else: 
                MatM = None
            if Time>TimeTarg-1.e-6 and not StLi.Buckl: StepFinishedFlag = True    # step time target
            else:                                      StepFinishedFlag = False
            if LineS!=None and LineS[0]>1:                                  # parameters for line search iteration
                LinS = LineS[0]
                LinT = LineS[1]
            else: 
                LinS=1
            StLi.BoundOffset( NodeList,NoLabToNoInd,NoIndToCMInd, VecU)     # add offset for prescribed displacement from current displacement for OPT=ADD
            
            IncCounter = 0                                                  # counter for increments within step
            if MaxEquiIter_ != None:
                EquiFailedMax = MaxEquiIter_
            else:
                EquiFailedMax = 3
            EquiFailedCounter = 0                                           # to count failing equilibrium iterations - initialization
            while not StepFinishedFlag:                                     # loop over time steps
                CalcType = 2                                                # 0: check system, 1: internal forces only, 2: internal forces and tangential stiffness matrix
                IncCounter += 1
                if Time+1.e-6>=TimeEl: TimeEl=Time + StLi.ElFilList[-1]     # set time for element output
                if Time+1.e-6>=TimeNo: TimeNo=Time + StLi.NoFilList[-1]     # set time for nodeal output
                A_BFGS, B_BFGS, rho_BFGS = [], [], []                       # A_BFGS, B_BFGS: BFGS auxiliary list of vectors, rho_BFGS: BFGS auxiliary list of scalars
                VecD = zeros((N),dtype=double)                              # displacement increment vector
                VecDI= zeros((N),dtype=double)                              # 
                VecR = zeros((N),dtype=double)                              # residual nodal forces vector
                VecRP= zeros((N),dtype=double)                              # residual nodal forces vector
                if StLi.Dyn:                                                # implicit dynamic calculation
                    Time = Time + dt                                        # new time for dynamic calculation
                    VecX = VecC + dt*VevC + 0.5*dt**2*(1-2*StLi.NMbeta)*VeaC # Newmark auxiliary vector
                    VecX_= -a1*VecC +(1-a_)*VevC + dt*(1-0.5*a_)*VeaC       # Newmark auxiliary vector
                else:
                    dt = 0                                                  # initial value for time step in case of quasistatic computation -- actual value later determined
                tt = 0                                                      # initialize time increment for 1st equilibrim iteration
                TimeTargetActiveFlag = False                                # flag to proceed with time / loading - initialized to check for initial equilibrium
                En0, En1 = 1., 1.                                           # initialization initial energy

                # equilibrium iteration loop
                for j in range(StLi.IterNum):
                    S = [0.]                                                # line search auxiliary value
                    G = [dot(VecRP,VecD)]                                   # line search auxiliary value
                    sds = 1.0                                               # line search scalar
                    VecUP[:] = VecU[:]                                      # remember displacement of last load increment
                    k_ = -1                                                 # auxiliary counter for line search
                    # line search iteration -- generally broken for j > 1 ???
                    for k in range(LinS):                                   
                        VecU = VecUP + sds*VecD                             # update displacements
                        if StLi.ArcLen: dU = norm( (VecU-VecC) * Mask )     # dU is for output only
                        else:           dU = norm(VecU-VecC)                # displacement increment compared to last step
                        if StLi.Dyn:                                        # auxiliary values in case of dynamic calculation
                            VevU = a1*VecU + VecX_                          # actual velocities
                            VeaU = a0*(VecU-VecX)                           # actual accelerations
                        nTmp = StLi.NodalTemp( N, Time, NodeList,NoLabToNoInd,NoIndToCMInd, VecT)            # introduce nodal temperatures into system with computation of nodal temperatures -> in SimFemSteps.py

                        # internal nodal forces and system stiffness (CalcType==2)
                        if ElasticLTFlag and j==0:                          # check system state for certain materials -- might lead to 'system change'
                            IntForces( MatList,ElList,Time-TimeOld, VecC,VecU,VecD,VecS,VecT,VecI, None,None,None,None,None, 0,f6,StLi.NLGeom,SymSys,False, j)
                        if LinAlgFlag:
                            if CalcType==2:                                 # 0: check system, 1: internal forces only, 2: internal forces and tangential stiffness matrix
                                KVecU = zeros(SLen, dtype=float)            # Initialize upper right part of stiffness vector
                                if not SymSys: KVecL = zeros(SLen, dtype=float) # Initialize lower left part of stiffness vector
                                else:               KVecL = None
                            # ! Time-TimeOld initially different for dynamics and quasi statics
                            IntForces( MatList, ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, None,None,KVecU,KVecL,SDiag, CalcType,f6,StLi.NLGeom,SymSys,False, j,nTmp=nTmp)# internal nodal forces / stiffness matrix
                        else:
                            if CalcType==2:
                                MatK = sparse.lil_matrix((N, N))            # sparse stiffness matrix initialization
                            IntForces( MatList, ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, MatK,MatM,None,None,None,    CalcType,f6,StLi.NLGeom,SymSys,StLi.Buckl, j,nTmp=nTmp)# internal nodal forces / stiffness matrix

                        # eigenmode analysis
                        if StLi.Buckl:
                            if not SymSys: raise NameError("ConFem::Run: Symmetric system required for eigenforms")
                            Eigen( EigenvalPar[0], EigenvalPar[1], N, NodeList, ElList, StLi, MatK, MatM, NoLabToNoInd,NoIndToCMInd, PloF, f6)
                            StepFinishedFlag = True
                            break                                           # 1st break to come out for buckling analysis
                        
                        # external loading -- before application of displacement boundary conditions
                        StLi.NodalLoads( N, Time, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecP, VecP0)# introduce concentrated loads into system
                        StLi.ElementLoads( Time, TimeTarg, ElList, VecP, VecP0)# introduce distributed loads into system
                        StLi.NodalPrestress( N, Time, ElList, VecP, VecU, StLi.NLGeom)# introduce prestressing
                        StLi.NodalPrestress( N, TimeTarg, ElList, VecP0, VecU, StLi.NLGeom)# nominal prestressing (maybe this works not optimal with P0-approach)
                        VecP0i[:] = VecP0[:]                                # remember current total load vector for use after step completion
                        VecP0[:]  = VecP0[:]-VecP0old[:]                    # only nominal load change compared to last step is relevant
                        VecB[:]   = VecI[:]-VecP[:]                         # boundary forces from internal / external forces difference -- before displacement boundary conditions
                        
                        # Newmark dynamics
                        if StLi.Dyn:
                            if LinAlgFlag: 
                                if CalcType==2:                             # internal forces and tangential stiffness matrix
                                    if StLi.Damp:
                                        DVecU                = StLi.RaBeta*KVecU + StLi.RaAlph*MVecU
                                        if not SymSys: DVecL = StLi.RaBeta*KVecL + StLi.RaAlph*MVecL
                                        KVecU                = KVecU + a1*DVecU
                                        if not SymSys: KVecL = KVecL + a1*DVecL
                                    KVecU                    = KVecU + a0 * MVecU
                                    if not SymSys:     KVecL = KVecL + a0 * MVecL
                                if SymSys:     LinAlg2.sim1_mmul(VecY, MVecU,        VeaU, SDiag, Skyline, N)
                                else:          LinAlg2.sim0_mmul(VecY, MVecU, MVecL, VeaU, SDiag, Skyline, N)
                                VecP = VecP - VecY
                                if StLi.Damp:
                                    if SymSys: LinAlg2.sim1_mmul(VecY, DVecU,        VevU, SDiag, Skyline, N)
                                    else:      LinAlg2.sim0_mmul(VecY, DVecU, DVecL, VevU, SDiag, Skyline, N)
                                    VecP = VecP - VecY
                            else:
                                if StLi.Damp:
                                    MatD = StLi.RaBeta*MatK + StLi.RaAlph*MatM  
                                    MatK = MatK + a1*MatD
                                    VecP = VecP - aslinearoperator(MatD).matvec(VevU)
                                MatK = MatK + a0*MatM                       # modify system matrix
                                MatK = MatK.tolil()                         # became CSR
                                VecP = VecP - aslinearoperator(MatM).matvec(VeaU) # modify load vector
                                
                        # displacement boundary conditions and residual vector
                        if LinAlgFlag:
                            StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, None,KVecU,KVecL,Skyline,SDiag, CalcType, SymSys, f6)# introduce boundary conditions
                        else:
                            StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, MatK,[],None,None,None, CalcType, SymSys, f6)# introduce boundary conditions
                        if CalcType==2: VecP0_[:]= VecP0[:]                 # remember nominal load in case stiffness matrix is not updated upcoming -- stiffness updating relevant for prescribed displacements
                        else:           VecP0[:] = VecP0_[:]                # use previous value of nominal load in case stiffness matrix was not updated
                        VecR[:] = VecP[:] - VecI[:]                         # residual vector - entries should be zero for prescribed displacement for j>0

                        # line search scaling, if prescribed
                        if j>1 and LinS>1:
                            G_ = dot(VecR,VecD)
                            G0 = dot(VecRP,VecD)
                            if G_*G[k_]<0.:
                                k_ = k_+1
                                S += [sds]
                                G += [G_]
                            if abs(G_)>LinT*abs(G0):
                                sds = S[k_] - G[k_]*(sds-S[k_])/(G_-G[k_])
                                if sds<0.0:   sds = 1.0
                                elif sds<0.2: sds = 0.2
                                elif sds>5.0: sds = 5.0
                                print('line search',k_,'G/G0',G_/G0,'step',sds)
                                print('line search',k_,'G/G0',G_/G0,'step',sds, file=f6)
                            else: break                                     # end of line search iteration
                        else: break                                         # no line search iteration
                    # end of line search loop
                    if StLi.Buckl: break                                    # 2nd break to come out for buckling analysis

                    # equilibrium control
                    Resi = norm(VecR)                                       # residual norm
                    if StLi.SolType=='BFGS': VecRD = VecRP - VecR + tt*DTS*VecP0 # BFGS auxiliary vector
                    VecRP[:] = VecR[:]                                      # residual nodal forces vector of previous iteration
                    if j>0:
                        En1 = dot(VecD, (VecR + BCIi * (VecB-VecBold)))     # residual energy, BCIi masks dofs NOT prescribed
                        if j==1: En0 = En1                                  # for convergence control via energy criterion
                        VecBold[:] = VecB[:]
                    if abs(En0)>ZeroD: Resi_ = abs(En1/En0)                 # energy based convergence indicator
                    else:              Resi_ = 1.
                    Echo(f"{StepCounter:d} {IncCounter:3d} {TimeTarg:.3f} {Time:.5f} iter {j:d} {Resi:e} {Resi_:e}    {dU:e}", f6)
                    f6.flush()
                    
                    # termination control
                    if Resi<StLi.IterTol:                                   # convergence criterion reached
                        if StLi.Dyn:
                            break                                           # continue after equilibrium iteration loop - independent from j
                        if j==0:
                            TimeTargetActiveFlag = True                     # initial system is in static equilibrium, proceed with time / load increment to iterate for equilibrium
                        else:                                               # successfully iterated for equilibrium
                            break                                           # continue after equilibrium iteration loop

                    # determine new solution - equation solver
                    if j==0 or StLi.SolType=='NR':                          # Newton Raphson
                        if LinAlgFlag:
                            if SymSys: LinAlg2.sim1_lu(      KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_lu(      KVecU, KVecL, SDiag, Skyline, N)
                            VecD[:]  = VecR[:]
                            if SymSys: LinAlg2.sim1_so(VecD, KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_so(VecD, KVecU, KVecL, SDiag, Skyline, N)
                            VecDI[:] = VecP0[:]
                            if SymSys: LinAlg2.sim1_so(VecDI,KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_so(VecDI,KVecU, KVecL, SDiag, Skyline, N)
                        else:
                            K_LU = splu(MatK.tocsc(),permc_spec=3) #triangulization of stiffness matrix
                            VecD = K_LU.solve( VecR )                       # solution of K*u=R -> displacement increment
                            VecDI= K_LU.solve( VecP0 )                      # solution contribution for arc length control
                        if StLi.SolType != 'NR': CalcType = 1               # stiffness matrix not built anymore 
                    elif StLi.SolType=='BFGS':                              # BFGS according to Matthies & Strang 1979 (S. 1617) for indefinite matrices
                        VecRD = VecRD * BCIn                                # mask prescribed dofs
                        VecD = VecD * BCIn                                  # mask prescribed dofs
                        A_BFGS += [VecRD]  
                        B_BFGS += [sds*VecD]
                        rho_BFGS += [1./dot(transpose(VecRD),sds*VecD)]
                        alpha_BFGS = []
                        for k in range(j-1,-1,-1):
                            alpha = rho_BFGS[k]*dot(B_BFGS[k],VecR)
                            VecR  = VecR - alpha * A_BFGS[k]
                            alpha_BFGS += [alpha]
                        alpha_BFGS.reverse()
                        if LinAlgFlag:
                            for jj in range(N): VecD[jj] = VecR[jj]
                            if SymSys: LinAlg2.sim1_so(VecD, KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_so(VecD, KVecU, KVecL, SDiag, Skyline, N)
                        else:
                            VecD = K_LU.solve(VecR)
                        for k in range(j):
                            beta_BFGS = rho_BFGS[k]*dot(A_BFGS[k],VecD)
                            VecD = VecD + (alpha_BFGS[k] - beta_BFGS) * B_BFGS[k]# end BFGS
                    else:
                        if LinAlgFlag:
                            if SymSys: LinAlg2.sim1_so(VecR, KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_so(VecR, KVecU, KVecL, SDiag, Skyline, N)
                            VecD[:] = VecR[:]
                            if SymSys: LinAlg2.sim1_so(VecP0,KVecU,        SDiag, Skyline, N)
                            else:      LinAlg2.sim0_so(VecP0,KVecU, KVecL, SDiag, Skyline, N)
                            VecDI[:] = VecP0[:]
                        else:
                            VecD = K_LU.solve( VecR )                       # solution of K*u=R -> displacement increment, modified Newton-Raphson
                            VecDI= K_LU.solve( VecP0 )                      # solution contribution for arc length control
                    # end equation solver

                    tt = 0                                                  # default time corrector
                    # determine time increment dt in case of quasi-statics -- has already been done for dynamics
                    if not StLi.Dyn and TimeTargetActiveFlag:               # TimeTargetActiveFlag set to False before equilibrium iteration loop -- set to True if equilbrium in j=0 iteration
                        if StLi.ArcLen: 
                            tt = TS*ArcLength(StLi.ArcLenV, VecDI, VecD, VecU-VecC, VecY, Mask) # time increment with arc length control
                        elif j==0:                                          # time increment determined for the following quilibrium j-iterations
                            if StLi.varTimeSteps:
                                for targ_i in range(TimeVarL): # 
                                    if Time < StLi.TimeTargVar[targ_i]: 
                                        tt = StLi.TimeStepVar[targ_i]
                                        break                               # breaks current targ_i for-loop
                            else:
                                tt=StLi.TimeStep                            # time increment prescribed 
                        dt = dt + tt                                        # update time increment corresponding to loading increment - dt initially 0 for j==0
                        Time = TimeOld + dt                                 # update time

                    VecD = VecD + tt*DTS*VecDI                              # new displacement increment with two contributions; tt=0 in case of dynamics and for quasi-statics j>0 & fixed time step
                    if j >= (StLi.IterNum - 2): LogResiduals(LogName, IncCounter, j, NodeList, VecRP, VecD)
                # end of equilibrium iteration loop -- either by reaching convergence criteria or by interation counter limit
                equiiterQueues[0].append( IncCounter )                      # load increment counter
                equiiterQueues[1].append( Time )                            # actual time
                equiiterQueues[2].append( j )                               # last iteration

                # premature termination cases -- if so
                if not StLi.Dyn and j==StLi.IterNum-1:
                    if EquiFailedCounter==EquiFailedMax-1:
                        Echo(f"{Name:s}\n premature termination - {EquiFailedCounter+1:d} iteration attempts - time {Time:.4f} of target {TimeTarg:.4f}", f6)
                        break                                               # continues after time increment loop
                    EquiFailedCounter += 1
                else:
                    EquiFailedCounter = 0
                if f5!=None:
                    q1 = list(timeoutQueues[1])                             # reaction forces or so
                    if len(q1) > 3:
                        if abs(q1[-1]) > 1.0e+01 * 0.5*(abs(q1[-2])+abs(q1[-3])):
                            Echo(f"{Name:s}\n premature termination - unstable reaction {q1[-1]:.4f} - time {Time:.4f} of target {TimeTarg:.4f}",f6)
                            break                                               # continues after time increment loop
                # book keeping for iterated time step
                TimeOld = Time
                VecY = VecU-VecC                                            # store displacement increment in time increment
                VecC[:] = VecU[:]                                           # store current displacements as previous for next time step
                VecS[:] = VecT[:]                                           # store current temperatures as previous for next time step
                VevC[:] = VevU[:]                                           # store current velocities as previous for next time step
                VeaC[:] = VeaU[:]                                           # store current accelerations as previous for next time step
                FinishEquilibIteration( MatList, ElList, NodeList,NoIndToCMInd, f6, StLi.NLGeom, TimeTargetActiveFlag)# update state variables etc.
                SelOutWrite( "elem", StressStrainOut, StressStrainOutNames, Time, ElList, None, None, None, StepCounter, IncCounter)   # output for single elements over all time steps
                
                nea = len(ElList)                                           # actual number of all elements (not only active)
                for el in ElList:
                    if el.Type in ["CPE4S","CPE3S","CPS4S","CPS3S"] and el.Active and el.NoCracks==1: 
                        el.CreateSDA2ndCrack( MatList, NodeList, NoIndToCMInd, f6, VecU,VecY, dt)  # must precede CreateSDARankine and not in a common ElemList loop
                for el in ElList:
                    if el.RegType==3 and el.Type in ["CPS4","CPS3", "CPE4","CPE3","C3D8","CPS6","CPE6"] and el.Active:
                        nea = el.CreateSDARankine( nea, MatList, NodeList, ElList, SecDic, NoLabToNoInd, NoIndToCMInd, f6, VecU,VecY, dt )
                
                if Time > TimeTarg-1.e-6:                                   # step time target
                    StepFinishedFlag = True                                 # time target reached, finalize computation for step
                if StLi.ArcLen: 
                    try:
                        pass
                        if norm(VecB)/norm(VecU) < 0.001: StepFinishedFlag = True                         # basically to stop snap back
                    except:
                        pass
                def WriteXData(TimeX):
                    EqF = not StLi.Dyn and j == StLi.IterNum - 1            # equilibrium iteration failed
                    if not EqF and not StLi.Buckl and (Time+1.e-6>=TimeX or StepFinishedFlag): return True
                    else:                                                                      return False
                if WriteXData(TimeEl):
                    DataOutStress(NameElemOut_, ElList, NodeList, NoIndToCMInd, Time, "a", f6)
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    DataOut(NameNodeOut_, NodeList, VecU, VecB, VecR, Time, "a" )
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                    WriteElemData( f2, f7, Time, ElList, NodeList,NoIndToCMInd, MatList, MaxType, ResultTypes)# write element data
                    f2.flush()
                    Echo(f"Element data written {Time:f}", f6)
                if WriteXData(TimeNo):
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    WriteNodalData( f3, Time, NodeList, VecU, VecB)         # write nodal data
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                    f3.flush()
                    Echo(f"Nodal data written {Time:f}", f6)
                if (Time+1.e-6>=TimeRe and j<StLi.IterNum-1) or StepFinishedFlag:
                    if j<(StLi.IterNum-1):
                        fd = open(Name+'.pkl', 'wb')                            # Serialize data and store for restart
                        pickle.dump(NodeList,fd);pickle.dump(ElList,fd);pickle.dump(MatList,fd);pickle.dump(StepList,fd);pickle.dump(N,fd);pickle.dump(WrNodes,fd);pickle.dump(LineS,fd);pickle.dump(ElasticLTFlag,fd);\
                            pickle.dump(VecU,fd);pickle.dump(VecC,fd);pickle.dump(VecI,fd);pickle.dump(VecP,fd);pickle.dump(VecP0,fd);pickle.dump(VecP0old,fd);pickle.dump(VecBold,fd);pickle.dump(VecT,fd);pickle.dump(VecS,fd);pickle.dump(VeaU,fd);pickle.dump(VevU,fd);pickle.dump(VeaC,fd);pickle.dump(VevC,fd);pickle.dump(VecY,fd);\
                            pickle.dump(BCIn,fd);pickle.dump(BCIi,fd);pickle.dump(Time,fd);pickle.dump(TimeOld,fd);pickle.dump(TimeEl,fd);pickle.dump(TimeNo,fd);pickle.dump(TimeS,fd);pickle.dump(StepCounter,fd);                     pickle.dump(Skyline,fd);pickle.dump(SDiag,fd);\
                        pickle.dump(SLen,fd);pickle.dump(SymSys,fd);pickle.dump(NoLabToNoInd,fd);pickle.dump(NoIndToCMInd,fd);pickle.dump(ContinuumNodes,fd);pickle.dump(ConNoToNoLi,fd);pickle.dump(SecDic,fd);pickle.dump(LinAlgFlag,fd);pickle.dump(ResultTypes,fd);pickle.dump(Header,fd);
                        fd.close()
                        Echo(f"Restart data written {Time:f}", f6)
                        if len(StLi.ReFilList)>0: TimeRe += StLi.ReFilList[-1]  # update time for restart output
                    else :
                        Echo(f"No restart data writing {Time:f}, residual {Time:e}", f6)
                if f5!=None:
                    WriteNodes( f5, WrNodes, Time, VecU, VecB, VecP, StepCounter, IncCounter, timeoutQueues)
                    f5.flush()
            # end of time increment loop for current step
                    
            VecP0old[:] = VecP0i[:]                                         # remember nominal load of this step
            StepCounter += 1                                                # next step
            #
            Echo("last iteration data", f6)
            for u, p, s in zip(list(equiiterQueues[0]), list(equiiterQueues[1]), list(equiiterQueues[2])):
                Echo(f" {u:4d}, {p:.4f}, {s:4d}", f6)
            if f5 != None:
                Echo("last timeout data", f6)
                for u, p, s in zip(list( timeoutQueues[0]), list(timeoutQueues[1]), list(timeoutQueues[2])):
                    Echo(f" {s:4d}, {u:.6f}, {p:.2f}", f6)
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
    Name, Plot, Restart, HashSource, Eig, ElPlotTimes, StressStrainOut, VTK = DefData() # VTK currently not used
#    Name=str(raw_input('Filename without extension: '))                    # input data name
    ConFem_ = ConFem()
    RC=ConFem_.Run(Name, LogName, Plot, LinAlgFlag, Restart, HashSource, Eig, ElPlotTimes, StressStrainOut, VTK)
    print(RC)