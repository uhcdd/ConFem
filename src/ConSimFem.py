# ConSimFem -- 2022-09-17
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import *
#from numpy import *
from numpy import ones, copy
from scipy.linalg import norm
from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import splu
from scipy.sparse.linalg import aslinearoperator
#import matplotlib.pyplot as plt
from os import path as pth
import pickle
#import sys
from collections import deque

from ConFemBasics import *
from ConFemMat import *
import ConFemElem
import ConFemSteps
from ConFemInOut import *
try:
    import LinAlg2
    LinAlgFlag = True
except ImportError:
    LinAlgFlag = False 

class ConSimFem:
    def __init__(self):
        self.SymSys = True                                                          # flag for symmetric system
#    @profile
    def Run(self, Name, PloF, LinAlgFlag, ResType, ElPlotTimes=['1.0000']):
        DirName, FilName = DirFilFromName(Name)
        f1=open( Name+".in.txt", 'r')
        f6=open( Name+".protocol.txt", 'w')
        Echo(f"SimFem  {Name:s}", f6)
        f5, f7 = None, None
        NameElemOut_ = Name + ".elemout.txt"
        NameNodeOut_ = Name + ".nodeout.txt"
        if pth.exists(NameElemOut_): os.remove(NameElemOut_)
        if pth.exists(NameNodeOut_): os.remove(NameNodeOut_)
        #        NodeList, ElList, MatList, StepList, NoLabToNoInd, SecDic = ReadInputFile(f1, f6, False)  # read input file and create node, element, material and step lists -> in SimFemInOut.py
        NodeList, ElList, MatList, SecDic, NoLabToNoInd, StepList, ResultTypes, Header = DataInput( f1, f6, False)
        f1.close()
        NoIndToCMInd = [i for i in range(len(NodeList))]                # maps identical before Cuthill McKee
        EvNodeEl( NodeList, ElList, NoLabToNoInd)                                     # determine element indices belonging to each node - required for Cuthill McKee
        if LinAlgFlag:
#                NoIndToCMInd = WithoutBandwidthOpti( NodeList )  # NoIndToCMInd
#                NoIndToCMInd = CuthillMckee(NodeList, ElList)       # to reduce the skyline
            NoIndToCMInd = Sloan(NodeList, ElList)                      # to reduce the skyline
            NodeList.sort(key=lambda t: t.CMIndex)
        N, Skyline, SDiag, SLen, ElemDataAll = AssignGlobalDof( NodeList, ElList, MatList, SecDic, NoIndToCMInd) # assign degrees of freedom (dof) to nodes and elements -> see above
        WrNodes, LineS = None, None
        MaxType = [] 
        if pth.isfile(Name+".opt.txt"):                                             # read options file if there is any
            f4=open( Name+".opt.txt", 'r')
            WrNodes, ReDes, MaxType, _,_,_,_,_ = ReadOptionsFile(f4, NodeList, NoLabToNoInd,NoIndToCMInd)
            f4.close()
            f5=open( Name+".timeout.txt", 'w')

        # Initializations
        VecU = zeros((N),dtype=double)                                              # current displacement vector
        VecC = zeros((N),dtype=double)                                              # displacement vector of previous time step
        VecI = zeros((N),dtype=double)                                              # internal nodal forces vector
        VecP = zeros((N),dtype=double)                                              # load vector
        VecP0= zeros((N),dtype=double)                                              # nominal load vector
        VecT = zeros((N),dtype=double)                                              # current temperatures vector
        VecS = zeros((N),dtype=double)                                              # temperatures vector of previous time step
        BCIn = ones((N),dtype=int)                                                  # indices for dofs with prescribed displacements --> 0, --> 1 otherwise
        BCIi = zeros((N),dtype=int)                                                 # indices for dofs with prescribed displacements --> 1, --> 0 otherwise 
        VecB = zeros((N),dtype=double)                                              # nodal boundary reactions
        Time = 0.
        TimeS = 0.                                                                  # Target time of previous step
        TimeOld = 0.
        TimeEl = 0.
        TimeNo = 0.
        Step = 0
    
        SymSys = IsSymSys( MatList )                                                # if there is at least one un-symmetric material the whole system is un-symmetric
        Node.ActiveNodes = 0
        for no in NodeList: 
            if len(no.NodeEl)>0: Node.ActiveNodes +=1
#        Echo(f"symmetric system {SymSys}, used LinAlg2 {LinAlgFlag}, used ConFemMatC {ConFemMatCFlag}, used CaeFemElemC {ConFemElemCFlag}", f6)
#        Echo(f"ndofs {N:d}, elems {len(ElList):d}, active nodes {Node.ActiveNodes:d}, abs.size {1.*SLen/1024:.2f} MItems, rel.size {SLen/(N**2):.4f} (abs.size/(ndofs**2)", f6)
        EchoSys( f6, SymSys, LinAlgFlag, ConFemMatCFlag, ConFemElemCFlag, N, ElList, Node.ActiveNodes, SLen, MatList, Element, ElemDataAll)
        # Calculations
        i = 0
        stime = process_time()
        for i in range(len(StepList)):                                              # loop over steps
            if f5 != None:
                ndeq = 20
                queues = [ deque(maxlen=ndeq), deque(maxlen=ndeq), deque(maxlen=ndeq) ]
            StLi = StepList[i]
            Echo(f"{i:d} step starts, solution type {StLi.SolType:s}", f6)
            TimeTarg = StLi.TimeTarg                                                # time target for step
            if i>0: TimeS = StepList[i-1].TimeTarg                                  # time target of previous step
            StLi.current = i
            Stop = False                                                            # flag to stop computation
            counter = 0
            StLi.BoundOffset( NodeList,NoLabToNoInd,NoIndToCMInd, VecU)                                       # add offset for prescribed displacement from current displacement for OPT=ADD
            if len(StLi.ElFilList)>0: TimeEl = TimeS
            else:                     TimeEl = TimeTarg
            if len(StLi.NoFilList)>0: TimeNo = TimeS
            else:                     TimeNo = TimeTarg

            maxWriteNodes = zeros((3), dtype=float)                                 # see ConFemInOut::WriteNodes
            while not Stop:                                                         # loop over time steps
                counter = counter+1
                if Time+1.e-6>=TimeEl: TimeEl=Time + StLi.ElFilList[-1]             # set time for element output
                if Time+1.e-6>=TimeNo: TimeNo=Time + StLi.NoFilList[-1]             # set time for nodeal output
                VecD = zeros((N),dtype=double)                                      # displacement increment vector
                Time = Time + StLi.TimeStep
                for j in range(StLi.IterNum):                                       # equilibrium iteration loop
                    VecU = VecU + VecD                                              # update displacements
                    if LinAlgFlag:
                        KVecU = zeros(SLen, dtype=float)                            # Initialize upper right part of stiffness vector
                        KVecL = zeros(SLen, dtype=float)                            # Initialize lower left part of stiffness vector
                        IntForces( MatList,ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, None,None,KVecU,KVecL,SDiag, 2,f6,StLi.NLGeom,SymSys,False, j)# internal nodal forces / stiffness matrix
                    else:
                        MatK = sparse.lil_matrix((N, N))                            # sparse stiffness matrix initialization
                        IntForces( MatList,ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, MatK,None,None,None,None,    2,f6,StLi.NLGeom,SymSys,False, j)# internal nodal forces / stiffness matrix
                    StLi.NodalLoads( N, Time, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecP, VecP0)      # introduce concentrated loads into system -> in SimFemSteps.py
                    StLi.ElementLoads( Time, TimeTarg, ElList, VecP, VecP0)# introduce distributed loads into system -> in SimFemSteps.py
                    VecB[:] = VecI[:]-VecP[:] #copy(VecI-VecP)                                          # boundary reactions
                    if LinAlgFlag:
                        StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, None,KVecU,KVecL,Skyline,SDiag, 2, SymSys, f6)# introduce boundary conditions
                    else:
                        StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, MatK,[],None,None,None, 2, SymSys, f6) # introduce boundary conditions
                    VecR = VecP - VecI                                              # residual vector
#                        for jj in xrange(len(VecI)): sys.stdout.write('%5i%16.8f%16.8f%16.8f\n'%(jj,VecI[jj],VecP[jj],VecR[jj]))
                    Resi = norm(VecR)                                               # residual norm
                    Echo(f"{i:d} {counter:3d} {TimeTarg:.3f} {Time:.5f} iter {j:d} {Resi:e}", f6)
                    if j>0 and Resi<StLi.IterTol: break                             # convergence criterium reached, continue after for-loop
                    if LinAlgFlag:
                        if SymSys:
                            LinAlg2.sim1_lu(      KVecU, SDiag, Skyline, N)
                            LinAlg2.sim1_so(VecR, KVecU, SDiag, Skyline, N)
                        else:
                            LinAlg2.sim0_lu(      KVecU, KVecL, SDiag, Skyline, N)
                            LinAlg2.sim0_so(VecR, KVecU, KVecL, SDiag, Skyline, N)
                        VecD[:] = VecR[:]
                    else:
#                        K_LU = linsolve.splu(MatK.tocsc(),permc_spec=3)             #triangulization of stiffness matrix
                        K_LU = splu(MatK.tocsc(),permc_spec=3)             #triangulization of stiffness matrix
                        VecD = K_LU.solve( VecR )                                   # solution of K*u=R -> displacement increment

                if Time > StLi.TimeTarg-1.e-6: Stop = True                          # time target reached, finalize computation, eventually with one equilibrium iteration
                TimeOld = Time
                VecC[:] = VecU[:]                                                   # store current displacements as previous for next time step
                FinishEquilibIteration( MatList, ElList, NodeList,NoIndToCMInd, f6, StLi.NLGeom, True)     # update state variables etc.
                if Time+1.e-6>=TimeEl: 
                    WriteElemData(NameElemOut_, ElList, NodeList, NoIndToCMInd, Time, "a", f6)
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    WriteNodeData(NameNodeOut_, NodeList, VecU, VecB, VecR, Time, "a")
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                    fd = open(Name+'.pkl', 'wb')                                    # Serialize data and store for restart
                    Flag, VecP0old, VecBold, VeaU, VevU, VeaC, VevC, VecY, FlElasticLT, ContinuumNodes, ConNoToNoLi = None, None, None, None, None, None, None, None, None, [], None # for compatibility with confem
                    pickle.dump(NodeList,fd);pickle.dump(ElList,fd);pickle.dump(MatList,fd);pickle.dump(StepList,fd);pickle.dump(N,fd);pickle.dump(WrNodes,fd);pickle.dump(LineS,fd);pickle.dump(FlElasticLT,fd);\
                        pickle.dump(VecU,fd);pickle.dump(VecC,fd);pickle.dump(VecI,fd);pickle.dump(VecP,fd);pickle.dump(VecP0,fd);pickle.dump(VecP0old,fd);pickle.dump(VecBold,fd);pickle.dump(VecT,fd);pickle.dump(VecS,fd);pickle.dump(VeaU,fd);pickle.dump(VevU,fd);pickle.dump(VeaC,fd);pickle.dump(VevC,fd);pickle.dump(VecY,fd);\
                        pickle.dump(BCIn,fd);pickle.dump(BCIi,fd);pickle.dump(Time,fd);pickle.dump(TimeOld,fd);pickle.dump(TimeEl,fd);pickle.dump(TimeNo,fd);pickle.dump(TimeS,fd);pickle.dump(Step,fd);                     pickle.dump(Skyline,fd);pickle.dump(SDiag,fd);\
                    pickle.dump(SLen,fd);pickle.dump(SymSys,fd);pickle.dump(NoLabToNoInd,fd);pickle.dump(NoIndToCMInd,fd);pickle.dump(ContinuumNodes,fd);pickle.dump(ConNoToNoLi,fd);pickle.dump(SecDic,fd);pickle.dump(LinAlgFlag,fd);pickle.dump(ResultTypes,fd);pickle.dump(Header,fd);
                    fd.close()

                if Time+1.e-6>=TimeNo:
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                    WriteNodeData(NameNodeOut_, NodeList, VecU, VecB, VecR, Time, "a")
                    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                if f5!=None:
                    WriteNodes( f5, WrNodes, Time, VecU,VecB,VecP, i, counter, queues,ndeq,maxWriteNodes)
                    f5.flush()
            Step += 1                                                       # next step
        # end of step loop 
        Echo(f"total comp time {process_time()-stime:.0f} seconds", f6)
        if f5!=None: f5.close()
        f6.close()
        RC = FinishAllStuff(PloF, DirName, FilName, Name, ResType, False)
#        StiffAlloc(MatK)
        return RC
 
if __name__ == "__main__":
#    numpy.seterr(all='raise')
    Plot = True                                   
    Name,Plot ="../DataExamples/E08/E8-01", True                                                # linear elastic plate
    Name="../DataExamples/E09/E9-01"                                               # input data name - linear elastic slab
    Name="../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load"                  #  "Shell"  "bridge_el05m"
    Name, Plot = "../_DataBeams/E3-02_B23", True                                                                           # TestSuiteTrusses
#    Name=str(raw_input('Filename without extension: '))                                 # input data name
    ConSimFem_ = ConSimFem()
    RC = ConSimFem_.Run(Name, Plot, LinAlgFlag, "elemout") # flag for plot
    print(RC)

