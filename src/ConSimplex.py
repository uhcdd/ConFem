# ConSimplex -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import *
from numpy import *
from scipy.linalg import *
from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import splu
from scipy.sparse.linalg import aslinearoperator
import matplotlib.pyplot as plt
from os import path as pth

from ConFemBasics import *
import ConFemMat
import ConFemElem
import ConFemSteps
from ConFemInOut import *
from ConLinAlg import *

def SetSimplexStage( N, NodeList, ElList, MatList, VecI, BoundL, CLoadL, ff):
    def WriteD( ff, BB):
        ff.write('   ')
        for j in range(BB.shape[1]): ff.write('%8i'%(j)) # 8.4f
        ff.write('\n')
        for i in range(BB.shape[0]):
            ff.write('%2i '%(i))
            for j in range(BB.shape[1]): ff.write('%8.4f'%(BB[i,j])) # 8.4f
            ff.write('\n')
        ff.write('\n')
        ff.flush()
        return 0
    def R2Update(i,j):                                  # eliminates column j using pivot in row i
        piv=1/BB[i,j]
        ra = 2*m+2
        xxx=piv*BB[i,0:ra]
        for k in range(nr+1):                           #m+n+2):      ###
            x = BB[k,j] 
            BB[k,0:ra] = BB[k,0:ra]-x*xxx
        BB[i,0:ra] = xxx
        return 0
    m = len(ElList)                                     # number of elements
    nc=2*m+1                                            # total number of columns
    MemF = zeros((m), dtype=double)                     # member forces
    MemU = zeros((m), dtype=double)                     # member limit force
    MemA = zeros((m), dtype=double)                     # member cross sectioo areas
    FInd = ones((m), dtype=int)                         # sign index for member forces
    Actv = zeros((nc), dtype=int)                       # active members in simplex alg
    ActR = array([-1 for i in range(nc)])              # active column to row / member to dof
    B = zeros((N,m), dtype=double)                      # equilibrium constraint matrix
    for i in range(m):                                 # loop over all elements
        Elem = ElList[i]                                # current element
        if not isinstance( Elem, ConFemElem.T2D2): raise NameError ("SimSim: only T2D2 element allowed")
        if not isinstance( MatList[Elem.MatN], ConFemMat.Mises): raise NameError ("SimSim: only Mises material allowed")
        A = Elem.Geom[1,0]                              # cross section area 
        MemA[i] = A                                     # cross section area vector
        MemU[i] = A*MatList[Elem.MatN].sigY             # fill vector for member limit load
        i0 = Elem.Inzi[0]
        i1 = Elem.Inzi[1]
        L = sqrt( (NodeList[i1].XCo-NodeList[i0].XCo)**2 + (NodeList[i1].YCo-NodeList[i0].YCo)**2 )# element length
        cosA = (NodeList[i1].XCo-NodeList[i0].XCo)/L        # direction cosine
        sinA = (NodeList[i1].YCo-NodeList[i0].YCo)/L        # direction sine
        B[Elem.DofI[0,0],i] = -cosA                     # fill equilibrium matrix
        B[Elem.DofI[0,1],i] = -sinA
        B[Elem.DofI[1,0],i] =  cosA
        B[Elem.DofI[1,1],i] =  sinA
        MemF[i] = A*Elem.Data[0][1]                     # member forces
        if MemF[i]<0: FInd[i]=-1                        # sign index for member forces
    ### consider boundary and load conditions
    BInd = []                                           # indices of dofs with boundary conditions
    UInd = []                                           # BB index to dof index 
    PInd = {}                                           # dicitionary for loaded nodes
    for i in range(len(BoundL)):                       # loop over all boundary conditions of step
# uhc        nI = BoundL[i].Index                            # node index of bc
        nI = FindIndexByLabel( NodeList, BoundL[i].NodeLabel) # node index of bc
        iterator = iter(NodeList[nI].DofT)              # iterator for set
        for j in range(len(NodeList[nI].DofT)):        # loop over all dofs of node
            if next(iterator)==BoundL[i].Dof: break    # dof type equals prescribed dof type
        k = NodeList[nI].GlobDofStart + j               # constrained dof global index
        BInd += [k]
    valm = 0                                            # largest load value
    LoadV = zeros((N))                                  # load vector all dofs for equilibrium control
    for i in range(len(CLoadL)):                       # loop over all nodal loads
        val = CLoadL[i].Val                             # prescribed load value
        if abs(val)>valm: valm=abs(val)
# uhc        nI = CLoadL[i].Index                            # node index of concentrated load
        nI = FindIndexByLabel( NodeList, CLoadL[i].NodeLabel) # node index of concentrated load
        iterator = iter(NodeList[nI].DofT)              # iterator for set
        for j in range(len(NodeList[nI].DofT)):        # loop over all dofs of node
            j0 = next(iterator)                        # dof type
            if j0==CLoadL[i].Dof: break                 # dof type equals prescribed dof type
        k = NodeList[nI].GlobDofStart + j               # loaded dof global index
        PInd[k] = val
        LoadV[k] = val                                  # nominal load vector
    for i in list(PInd.keys()): PInd[i] = PInd[i]/valm        # scale load values with maximum amount of 1 
    for i in range(N):
        if i not in BInd: UInd += [i]                   # append indices of active dof
    n = len(UInd)                                       # number of equilibrium conditions
    print('Total dofs ',N,' Elements ',m,'Active dofs ',n)

    ### set stage
#    ff=open( "xxx.txt", 'w')
    nr = m+n                                            # number of restraints (m number of members, n number of active dofs)
    ActI = array([-1 for i in range(nr)])              # index of active column in row / adof to member
    BB = zeros((nr+1,nc+1),dtype=double)                # set stage for simplex tableau nc=2*m+1 
    for i in range(m):                                 # first m rows and nc columns of simplex tableau
        BB[i,i] = 1                                    #       member reserve
        BB[i,i+m] = 1                                  #       member force
        BB[i,nc] = MemU[i]                             #       last column of simplex tableau, member limit load
    for i in range(n):                                 # final n rows of simplex tableau
        for j in range(m):  BB[m+i,m+j] =  FInd[j]*B[UInd[i],j]
        if UInd[i] in PInd: BB[m+i,2*m] = -PInd[UInd[i]]               # unit load vector
    BB[nr,2*m] = -1                                     # target function coefficient (last row for target function)
#    WriteD(ff, BB)
            
    ### create canonical form Luenberger p. 
    for i in range(nr):                                 # loop over rows / Adofs of BB but last to create canonical form
        for j in range(nc):                             # loop over columns / members of BB but last                       
            if abs(BB[i,j])<ZeroD or Actv[j]==1: continue
            else:
                Actv[j] = 1                             # make this member active
                ActR[j] = i                             # Member To Adof (active degree of freedom)
                ActI[i] = j                             # Adof To Member
                R2Update(i,j)                           # rank two update for new basic variable
                break
#    WriteD(ff, BB)

    ### simplex algorithm
    Tol = -1.e-12
    Bland = False
    sol = zeros((nc), dtype=double)
    for i in range(10):                                 # simplex algorithm loop
        ic = -1                                          # initialize index for variable / column to go in and become basic
        ir = -1
        rC = 0                                           # initialize relative cost
        rD = 9999
        if Bland:                                           # apply Bland's rule
            pass
#            for j in xrange(nc):
#                if BB[nr,j]<Tol:
#                   rC = BB[nr,j]
#                   ic=j
#                   break
#           print(ic,rC)
#           if ic<0: break                              # no decrease of target function
#           rD = 9999
#           rL = []
#           rL_= []
#           rV = []
#           for j in xrange(nr):                            # loop over rows / Adofs to find variable to go out of basic variables
#               if BB[j,ic]>ZeroD and BB[j,nc]/BB[j,ic]<rD-Tol: rD = BB[j,nc]/BB[j,ic] 
#           for j in xrange(nr):                            # loop over rows / Adofs to find variable to go out of basic variables
#               if BB[j,ic]>ZeroD and BB[j,nc]/BB[j,ic]<rD-Tol: 
#                   rL += [j]
#                   rL_+= [ActI[j]]
#                   rV += [BB[j,nc]/BB[j,ic]]                    
#           ir = min(rL)
        else:
            for j in range(nc):                        # loop over all variables / columns
                if BB[nr,j]<rC:                         # find cost factor
                    rC = BB[nr,j]                       #
                    ic = j                              # current column/member to become basic
            if rC>=Tol: break                           # no decrease of target function, presumbaly got solution
            for j in range(nr):                        # loop over rows / Adofs to find variable to go out of basic variables
                if BB[j,ic]>ZeroD:
                    xxx = BB[j,nc]/BB[j,ic]
                    if xxx < (rD-Tol):                  # Luenberger P. 37 
                        rD = xxx 
                        ir = j                          # current variable to go out of basic variables
        iC = ActI[ir]                                   # Adof To Member / column belonging to row
#        if Bland: print(rL,rL_,rV,iC)
#        else: print(rC,ic,ir,iC)
        Actv[ic] = 1                                    # make this column active
        Actv[iC] = 0                                    # make this column inactive
        ActR[ic] = ir                                   # Member/column to come in To Adof/row
        ActI[ir] = ic                                   # Adof To Member
        R2Update(ir,ic)                                 # rank two update for new basic variable
#        WriteD(ff, BB)
        for j in range(nc):                            # condense current solution
            if Actv[j]==1: sol[j]=BB[ActR[j],nc]
            else:          sol[j]=0
#        res = sol[0:m]
#        sss = sol[m:2*m]
#        sum = sss+res
#        print(res,sss,sum,BB[nr,nc])
        print('Rank two update; out, in:',iC,ic,'Load factor:',BB[nr,nc])

    # post processing
    sig = array([FInd[i]*sol[m+i] for i in range(m)])  # member forces of solution
    resi = dot(B,sig) - BB[nr,nc]/valm * LoadV
    for i in BInd: resi[i]=0
    print('Equilibrium control ',norm(resi))                                # equilibrium control
    plt.figure()
    plt.axis('equal')
    plt.grid()
    plt.title('ideal plastic solution')
    FonSizTi='x-large'                                                      # fontsize title x-large
    FonSizAx='x-large'                                                      # fontsize axis
    ff.write(" element       force\n")
    for i in range(m):                                                      # loop over all elements for plot
        Elem = ElList[i]                                                    # current element
        xN = [NodeList[Elem.Inzi[0]].XCo,NodeList[Elem.Inzi[1]].XCo]
        yN = [NodeList[Elem.Inzi[0]].YCo,NodeList[Elem.Inzi[1]].YCo]
        plt.plot(xN,yN, 'b--')
        plt.plot(xN,yN, 'ro')
        XX = 0.5*(NodeList[Elem.Inzi[0]].XCo+NodeList[Elem.Inzi[1]].XCo)
        YY = 0.5*(NodeList[Elem.Inzi[0]].YCo+NodeList[Elem.Inzi[1]].YCo)
        Angle = arcsin(Elem.Geom[1,2]/sqrt(Elem.Geom[1,3]))*180/3.141593 
        if sig[i]>0: plt.text(XX,YY,format(sig[i]/MemA[i],".2f"),ha='center',va='bottom',rotation=Angle,color='red',fontsize=FonSizAx)
        else: plt.text(XX,YY,format(sig[i]/MemA[i],".2f"),ha='center',va='bottom',rotation=Angle,color='green',fontsize=FonSizAx)
        ff.write('%8i%  12.4f\n'%(Elem.Label,sig[i]/MemA[i]))
#    ff.write('   %14.6f\n'%(BB[nr,nc]))
    plt.yticks(fontsize=FonSizAx) 
    plt.xticks(fontsize=FonSizAx) 
    return 0

class ConSimplex:
    def __init__(self):
        pass
    def Run(self, Name, LinAlgFlag, PloF, ElPlotTimes=['1.0000']):
        f1 =open( Name+".in.txt", 'r')
        f2 =open( Name+".elemout.txt", 'w')
        f3 =open( Name+".nodeout.txt", 'w')
        f6 =open( Name+".protocol.txt", 'w')
        f7=None
        print("ConSim: ", Name)
        NodeList, ElList, MatList, SecDic, NoLabToNoInd, StepList, ResultTypes, Header = DataInput( f1, f6, False)
        f1.close()
        NoIndToCMInd = [i for i in range(len(NodeList))]                    # maps identical before Cuthill McKee
        N, Skyline, SDiag, SLen, _ = AssignGlobalDof( NodeList, ElList, MatList, SecDic, NoIndToCMInd)              # assign degrees of freedom (dof) to nodes and elements -> see above
        # Initializations
        VecU = zeros((N),dtype=double)                                      # current displacement vector
        VecI = zeros((N),dtype=double)                                      # internal nodal forces vector
        VecA = zeros((N),dtype=double)                                      # nodal forces vector from prescribed constrained dofs
        VecP = zeros((N),dtype=double)                                      # load vector
        VecP0= zeros((N),dtype=double)                                      # dummy load vector
        VecZ = zeros((N),dtype=double)                                      # zero vector
        BCIn = ones((N),dtype=int)                                          # for compatibility with ConFem
        BCIi = zeros((N),dtype=int)                                         # for compatibility with SimFem
        # FE Calculation
        stime = process_time()
        Time = 1                                                            # set time
        if LinAlgFlag:
            KVecU = zeros(SLen, dtype=float)                                # Initialize upper right part of stiffness vector
            KVecL = zeros(SLen, dtype=float)                                # Initialize lower left part of stiffness vector
            IntForces( MatList,ElList, Time, VecZ,VecU,VecZ,VecZ,VecZ, VecI, None,None,KVecU,KVecL,SDiag, 2, None,False,False,False, 0)# internal nodal forces / stiffness matrix
        else:
            MatK = sparse.lil_matrix((N, N))                    # sparse stiffness matrix initialization
            IntForces( MatList,ElList, Time, VecZ,VecU,VecZ,VecZ,VecZ, VecI, MatK,None,None,None,None,    2, None,False,False,False, 0)# internal nodal forces / stiffness matrix
        StepList[0].NodalLoads( N, Time, Time, NodeList,NoLabToNoInd,NoIndToCMInd, VecP, VecP0)# introduce concentrated loads into system -> in ConFemSteps.py
        if LinAlgFlag:
            StepList[0].BoundCond( N, Time, 0, Time, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, None,KVecU,KVecL,Skyline,SDiag, 2, False, f6)# introduce boundary conditions
        else:
            StepList[0].BoundCond( N, Time, 0, Time, NodeList, VecU, VecI, VecP, VecP0, BCIn, BCIi, MatK,[],None,None,None, 2, False, f6)# introduce boundary conditions
        VecR = VecP + VecA - VecI                                           # residual vector
        print(f"time target {StepList[0].TimeTarg:f}, actual time {Time:f}, sum load vector {norm(VecR):e}")
        if LinAlgFlag:
            LinAlg2.sim0_lu(      KVecU, KVecL, SDiag, Skyline, N)
            LinAlg2.sim0_so(VecR, KVecU, KVecL, SDiag, Skyline, N)
            VecU[:] = VecR[:] #copy(VecR)
            IntForces( MatList,ElList, Time, VecZ,VecU,VecZ,VecZ,VecZ, VecI, None,None,KVecU,KVecL,SDiag, 1, None,False,False,False, 0)# internal nodal forces
        else:
#            K_LU = linsolve.splu(MatK.tocsc(),permc_spec=3)     #triangulization of stiffness matrix
            K_LU = splu(MatK.tocsc(),permc_spec=3)     #triangulization of stiffness matrix
            VecU = K_LU.solve( VecR )                           # solution of K*u=R -> displacement increment
            IntForces( MatList,ElList, Time, VecZ,VecU,VecZ,VecZ,VecZ, VecI, MatK,None,None,None,None,    1, None, False,False,False, 0)# internal nodal forces
        DataOutStress(Name+".elemout_.txt", ElList, NodeList, NoIndToCMInd, 1.0, "w", None)
        WriteElemData( f2, f7, Time, ElList, NodeList,NoIndToCMInd, MatList, [], ResultTypes)          # write element data
        if LinAlgFlag: NodeList.sort(key=lambda t: t.Label)
        WriteNodalData( f3, Time, NodeList, VecU, VecU)     # write nodal data
        DataOut(Name+".nodeout_.txt", NodeList, VecU, VecI, VecR, 1.0, "w")
        if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)
        f2.close()
        f3.close()
        f6.close()
        print(f"finished prior FEM calculation -- one step with nominal load")
        # post processing
        if PloF:
            from ConFemPost import ConFemPost 
            MaPlLib, VTK, VTK3D = True, False, False
            StepCounter = 0
            ConFemPost_ = ConFemPost()
            ConFemPost_.Run( "", "", Name,StepCounter, MaPlLib, VTK,VTK3D, 1, [])
        print(f"\ngo on with solving dual optimization problem with simplex method\nloading vector is scaled to unit loading (max entry 1)")
        # simplex
        f4=open( Name+".rigidplastic.txt", 'w')
        SetSimplexStage( N, NodeList, ElList, MatList, VecI, StepList[0].BoundList, StepList[0].CLoadList, f4)
        f4.close()
        print("computation time",(process_time()-stime))
        if PloF:
            plt.show()
            return 0
        else:
            import hashlib
            mmm = hashlib.md5()
            fp= open( Name+".rigidplastic.txt", "r")
            while True:
                data= fp.read(65536)
                if not data: break
                mmm.update(data.encode())
            fp.close()
            RC = mmm.hexdigest()
            return RC

if __name__ == "__main__":
    Name="../DataExamples/E05/E5-03"                         # input data name
    Name="../DataExamples/E05/E5-02"                         # input data name
    ConSimplex_ = ConSimplex()
    LinAlgFlag = True                                   # True: store matrices as vector and use own c routines for LU-decomposition and backward substitution 
    PlotFlag = True
    RC = ConSimplex_.Run(Name, LinAlgFlag, PlotFlag)
    print("finished")
