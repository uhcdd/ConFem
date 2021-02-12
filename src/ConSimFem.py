# ConSimFem -- 2014-01-13
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import *
#from numpy import *
from numpy import ones, copy
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import aslinearoperator
#import matplotlib.pyplot as plt
from os import path as pth
import pickle
#import sys

import ConFemBasics
import imp
imp.reload(ConFemBasics)
from ConFemBasics import *
import ConFemMat
imp.reload(ConFemMat)
from ConFemMat import *
import ConFemElem
imp.reload(ConFemElem)
from ConFemElem import *
import ConFemSteps
imp.reload(ConFemSteps)
from ConFemSteps import *
import ConFemInOut
imp.reload(ConFemInOut)
from ConFemInOut import *
try:
    import LinAlg2
    imp.reload(LinAlg2)
    from LinAlg2 import *
    Linalg_possible=True
except ImportError:
    Linalg_possible=False

class ConSimFem:
    def __init__(self):
        self.SymSys = True                                                          # flag for symmetric system
#    @profile
    def Run(self, Name, PloF, LinAlgFlag, ResType, ElPlotTimes=['1.0000']):
        f1=open( Name+".in.txt", 'r')
        f2=open( Name+".elemout.txt", 'w')
        f3=open( Name+".nodeout.txt", 'w')
        f6=open( Name+".protocol.txt", 'w')
        Echo(f"SimFem  {Name:s}", f6)
        f5, f7 = None, None
        NodeList, ElList, MatList, StepList, NoLabToNoInd, SecDic = ReadInputFile(f1, f6, False)  # read input file and create node, element, material and step lists -> in SimFemInOut.py
        f1.close()
        NoIndToCMInd = [i for i in range(len(NodeList))]                # maps identical before Cuthill McKee
        EvNodeEl( NodeList, ElList, NoLabToNoInd)                                     # determine element indices belonging to each node - required for Cuthill McKee
        if LinAlgFlag:
#                NoIndToCMInd = WithoutBandwidthOpti( NodeList )  # NoIndToCMInd
#                NoIndToCMInd = CuthillMckee(NodeList, ElList)       # to reduce the skyline
            NoIndToCMInd = Sloan(NodeList, ElList)                      # to reduce the skyline
            NodeList.sort(key=lambda t: t.CMIndex)
        N, Mask, Skyline, SDiag, SLen = AssignGlobalDof( NodeList, ElList, MatList, NoIndToCMInd) # assign degrees of freedom (dof) to nodes and elements -> see above
        WrNodes, LineS = None, None
        MaxType = [] 
        if pth.isfile(Name+".opt.txt"):                                             # read options file if there is any
            f4=open( Name+".opt.txt", 'r')
            WrNodes, LineS, ReDes, MaxType = ReadOptionsFile(f4, NodeList, NoLabToNoInd,NoIndToCMInd)
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
    
        SymSys = IsSymSys( MatList )                                                # if there is at least one un-symmetric material the whole system is un-symmetric
        Echo(f"symmetric system {SymSys}, used LinAlg2 {LinAlgFlag}, used ConFemMatC {ConFemMatCFlag}, used CaeFemElemC {ConFemElemCFlag}", f6)
#        Echo(f"ndofs {N:d}, elems {len(ElList):d}, nodes {len(NodeList):d}, rel size {SLen/(N**2):.4f}, abs size {1.*SLen/1024:.4f} MWords", f6)
        Node.ActiveNodes = 0
        for no in NodeList: 
            if len(no.NodeEl)>0: Node.ActiveNodes +=1
        Echo(f"ndofs {N:d}, elems {len(ElList):d}, active nodes {Node.ActiveNodes:d}, rel size {SLen/(N**2):.4f}, abs size {1.*SLen/1024:.4f} MWords", f6)
        # Calculations
        i = 0
        stime = process_time()
        for i in range(len(StepList)):                                              # loop over steps
            StLi = StepList[i]
            Echo(f"{i:d} step starts, solution type {StLi.SolType:s}", f6)
            TimeTarg = StLi.TimeTarg                                                # time target for step
            if i>0: TimeS = StepList[i-1].TimeTarg                                  # time target of previous step
            StLi.current = i
            Stop = False                                                            # flag to stop computation
            counter = 0
            StLi.BoundOffset( NodeList,NoLabToNoInd,NoIndToCMInd, VecU)                                       # add offset for prescribed displacement from current displacement for OPT=ADD

            while not Stop:                                                         # loop over time steps
                counter = counter+1
                if len(StLi.ElFilList)>0 and Time+1.e-6>=TimeEl: TimeEl=Time + StLi.ElFilList[0].OutTime# set time for element output
                if len(StLi.NoFilList)>0 and Time+1.e-6>=TimeNo: TimeNo=Time + StLi.NoFilList[0].OutTime# set time for element output
                VecD = zeros((N),dtype=double)                                      # displacement increment vector
                VecR = zeros((N),dtype=double)                                      # residual nodal forces vector
                Time = Time + StLi.TimeStep
                for j in range(StLi.IterNum):                                       # equilibrium iteration loop
                    VecU = VecU + VecD                                              # update displacements
                    if LinAlgFlag:
                        KVecU = zeros(SLen, dtype=float)                            # Initialize upper right part of stiffness vector
                        KVecL = zeros(SLen, dtype=float)                            # Initialize lower left part of stiffness vector
                        IntForces( N, MatList,ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, None,None,KVecU,KVecL,Skyline,SDiag, 2,f6,StLi.NLGeom,SymSys,False, j)# internal nodal forces / stiffness matrix
                    else:
                        MatK = sparse.lil_matrix((N, N))                            # sparse stiffness matrix initialization
                        IntForces( N, MatList,ElList, Time-TimeOld, VecC,VecU,VecD,VecS,VecT, VecI, MatK,None,None,None,None,None,       2,f6,StLi.NLGeom,SymSys,False, j)# internal nodal forces / stiffness matrix
                    StLi.NodalLoads( N, Time, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecP, VecP0)      # introduce concentrated loads into system -> in SimFemSteps.py
                    StLi.ElementLoads( Time, TimeTarg, ElList, VecP, VecP0)# introduce distributed loads into system -> in SimFemSteps.py
                    VecB[:] = VecI[:]-VecP[:] #copy(VecI-VecP)                                          # boundary reactions
                    if LinAlgFlag:
                        StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecU, VecI, VecP, VecP0, BCIn, BCIi, None,KVecU,KVecL,Skyline,SDiag, 2, SymSys)# introduce boundary conditions
                    else:
                        StLi.BoundCond( N, Time, TimeS, TimeTarg, NodeList,NoLabToNoInd,NoIndToCMInd, VecU, VecI, VecP, VecP0, BCIn, BCIi, MatK,[],None,None,None, 2, SymSys)# introduce boundary conditions
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
                        K_LU = linsolve.splu(MatK.tocsc(),permc_spec=3)             #triangulization of stiffness matrix
                        VecD = K_LU.solve( VecR )                                   # solution of K*u=R -> displacement increment

                if Time > StLi.TimeTarg-1.e-6: Stop = True                          # time target reached, finalize computation, eventually with one equilibrium iteration
                TimeOld = Time
                VecC[:] = VecU[:]                                                   # store current displacements as previous for next time step
                FinishEquilibIteration( MatList, ElList, NodeList,NoIndToCMInd, f6, StLi.NLGeom, True)     # update state variables etc.
                if Time+1.e-6>=TimeEl: 
                    WriteElemData( f2, f7, Time, ElList, NodeList,NoLabToNoInd,NoIndToCMInd, MatList, MaxType)# write element data
                    fd = open(Name+'.pkl', 'wb')                                    # Serialize data and store for restart
                    Flag, VecP0old, VecBold, VeaU, VevU, VeaC, VevC, VecY = None, None, None, None, None, None, None, None # for compatibility with confem
                    pickle.dump(NodeList,fd);pickle.dump(ElList,fd);pickle.dump(MatList,fd);pickle.dump(StepList,fd);pickle.dump(N,fd);pickle.dump(WrNodes,fd);pickle.dump(LineS,fd);pickle.dump(Flag,fd);\
                        pickle.dump(VecU,fd);pickle.dump(VecC,fd);pickle.dump(VecI,fd);pickle.dump(VecP,fd);pickle.dump(VecP0,fd);pickle.dump(VecP0old,fd);pickle.dump(VecBold,fd);pickle.dump(VecT,fd);pickle.dump(VecS,fd);pickle.dump(VeaU,fd);pickle.dump(VevU,fd);pickle.dump(VeaC,fd);pickle.dump(VevC,fd);pickle.dump(VecY,fd);\
                        pickle.dump(BCIn,fd);pickle.dump(BCIi,fd);pickle.dump(Time,fd);pickle.dump(TimeOld,fd);pickle.dump(TimeEl,fd);pickle.dump(TimeNo,fd);pickle.dump(TimeS,fd);pickle.dump(i,fd);pickle.dump(Mask,fd);pickle.dump(Skyline,fd);pickle.dump(SDiag,fd);
                    pickle.dump(SLen,fd);pickle.dump(SymSys,fd);pickle.dump(NoLabToNoInd,fd);pickle.dump(NoIndToCMInd,fd);
                    fd.close()

                if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
                if Time+1.e-6>=TimeNo:
                    WriteNodalData( f3, Time, NodeList, VecU, VecB)                 # write nodal data
                if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
                if f5!=None:
                    WriteNodes( f5, WrNodes, Time, VecU, VecB, VecP, i, counter)
                    f5.flush()
        Echo(f"total comp time {process_time()-stime:.0f} seconds", f6)
        f2.close()
        f3.close()
        if f5!=None: f5.close()
        f6.close()
        RC = FinishAllStuff(PloF, Name, ElList, NodeList,NoIndToCMInd, MatList, f2, VecU, WrNodes, ResType, ElPlotTimes)
#        StiffAlloc(MatK)
        return RC
 
if __name__ == "__main__":
#    numpy.seterr(all='raise')
    Name="../DataExamples/E08/E6-02"                                                # linear elastic plate
#    Name="../DataExamples/E09/E7-01"                                               # input data name - linear elastic slab
    Name="../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load"                  #  "Shell"  "bridge_el05m"
#    Name=str(raw_input('Filename without extension: '))                                 # input data name
    ConSimFem_ = ConSimFem()
    if Linalg_possible: LinAlgFlag = True                                           # True: store matrices as vector and use own c routines for LU-decomposition and backward substitution 
    else:               LinAlgFlag = False
    Plot = True                                   
    ConSimFem_.Run(Name, Plot, LinAlgFlag, "elemout") # flag for plot

