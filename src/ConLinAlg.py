# ConLinAlg -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
#from matplotlib.pyplot import matshow, show, grid, minorticks
import matplotlib.pyplot as plt
from matplotlib.pylab import *
from numpy import zeros, loadtxt, dot, transpose #, array, sqrt, double, fabs, ma
from ConFemBasics import FindIndexByLabel
import pickle
import imp

try:
    import LinAlg2
#    Linalg_possible=True
except ImportError:
    pass
#    Linalg_possible=False

#Plot Allocation of stiffness matrix
def StiffAlloc(MatK):
    MatK_plot=MatK.todense()
    for r in range(len(MatK_plot)):
        for c in range(len(MatK_plot)):
            if MatK_plot[r,c] != 0: MatK_plot[r,c]=0.5
            if r==c: MatK_plot[r,c]=1.0
    plt.matshow(abs(MatK_plot),cmap=cm.binary)
#    plt.minorticks_on()
    if MatK_plot.shape[0]<50:
        plt.xticks(arange(MatK_plot.shape[0]))
        plt.yticks(arange(MatK_plot.shape[1]))
    plt.grid(b=True, which='major', linestyle='-')
##    plt.grid(b=True, which='minor', linestyle='-')
    plt.show()
    
# Write Matrix stored as vectors
def SkyWriteFull(NN, SLen, VecU, VecL, Skyline, SDiag, VecR, VecD, MatK):
    # assemble stiffness matrix
    N = 20
    KK = zeros((N,N), dtype = float)
    if VecU!=None:
        for i in range(N):
            DiagAdr = SDiag[i]              # 
            SkylAdr = Skyline[i]
            for j in range(SkylAdr):
                jj = i-j
                KK[jj,i] = VecU[DiagAdr+j]
            for j in range(SkylAdr-1):
                jj = i-j-1
                KK[i,jj] = VecL[DiagAdr+1+j]

    if MatK!=None: MatK_=MatK.todense()
    ff=open( "matrix.txt", 'w')
    for i in range(N):
        for j in range(N): ff.write('%12.4e'%(KK[i,j]-MatK_[i,j]))
#        for j in xrange(N): ff.write('%12.4e'%(MatK_[i,j]))
        ff.write('\n')
    ff.close()

    # write vectors
#    ff=open( "NN.txt", 'w')
#    ff.write('%10i\n%10i\n'%(N,SLen))
#    ff.close()
#    ff=open( "Skyl.txt", 'w')
#    for i in xrange(N): ff.write('%10i\n'%(Skyline[i]))
#    ff.close()
#    ff=open( "SDiag.txt", 'w')
#    for i in xrange(N): ff.write('%10i\n'%(SDiag[i]))
#    ff.close()
#    ff=open( "VecU.txt", 'w')
#    for i in xrange(SLen): ff.write('%20.12e\n'%(VecU[i]))
#    ff.close()
#    ff=open( "VecL.txt", 'w')
#    for i in xrange(SLen): ff.write('%20.12e\n'%(VecL[i]))
#    ff.close()
    if VecR!=None:
        ff=open( "VecR.txt", 'w')
        for i in range(NN): ff.write('%20.12e\n'%(VecR[i]))
        ff.close()
    if VecD!=None:
        ff=open( "VecD.txt", 'w')
        for i in range(NN): ff.write('%20.12e\n'%(VecD[i]))
        ff.close()
#    Re = VecD[0:N]
#    MatK__ = MatK_[0:N,0:N]
#    print Re, MatK__
#    
#    Le = dot(MatK__,transpose(Re))
#    print Le
    raise NameError ("Stop in SkyWriteFull")

def CompareLU(N, SLen, Skyline, SDiag):
    LUVecU = zeros(SLen, dtype=float)    # Initialize upper right part of stiffness vector
    ff=open( "C:/Users/uhc/Documents/Work/SoftW/LinAlg2/res1.txt", 'r')
    i = 0
    z1 = ff.readline()                                  # 1st input line
    while z1!="":
        LUVecU[i]=float(z1.strip())
        i = i+1
        z1 = ff.readline()
    ff.close()
    LUVecL = zeros(SLen, dtype=float)    # Initialize lower left part of stiffness vector
    ff=open( "C:/Users/uhc/Documents/Work/SoftW/LinAlg2/res2.txt", 'r')
    j = 0
    z1 = ff.readline()                                  # 1st input line
    while z1!="":
        LUVecU[j]=float(z1.strip())
        j = j+1
        z1 = ff.readline()
    ff.close()
    KK = zeros((N,N), dtype = float)
    for i in range(N):
        DiagAdr = SDiag[i]              # 
        SkylAdr = Skyline[i]
        for j in range(SkylAdr):
            jj = i-j
            KK[jj,i] = LUVecU[DiagAdr+j]
        for j in range(SkylAdr-1):
            jj = i-j-1
            KK[i,jj] = LUVecL[DiagAdr+1+j]
    print(i,j)
    ff=open( "matrixLU.txt", 'w')
    for i in range(N):
        for j in range(N): ff.write('%12.4e'%(KK[i,j]))
        ff.write('\n')
    ff.close()
    raise NameError ("Stop in CompareLU")

def fill_vec(fname):
  return np.loadtxt(fname)

def PickleSysMat(NodeList, ElList, Step, N, MatK, MatM):
    def BC(Kmat, Flag):
        for i in range(len(Step.BoundList)):           # loop over all boundary conditions of step
            Found = False
            nI = FindIndexByLabel( NodeList, Step.BoundList[i].NodeLabel) # node index of bc
            iterator = iter(NodeList[nI].DofT)          # iterator for set
            for j in range(len(NodeList[nI].DofT)):    # loop over all dofs of node
                jj=next(iterator)
                if jj==Step.BoundList[i].Dof:           # dof type equals prescribed dof type 
                    Found = True
                    break 
            if not Found: raise NameError ("ConLinAlg::PickleSysMat:missing correspondence for dof type")
            k = NodeList[nI].GlobDofStart + j           # constrained dof global index
            for j in range(N):                         # loop over rows and columns simul
                Kmat[j,k] = 0.                          # modification stiffness matrix
                Kmat[k,j] = 0.
            if Flag: Kmat[k,k] = 1.                     # modification stiffness matrix
    Name = 'Sysmat'
    fd = open(Name+'.pkl', 'w') # Serialize data and store for restart
    ppp = pickle.Pickler(fd)
    ppp.dump(N)
    ppp.dump(NodeList)
    ppp.dump(ElList)
    ppp.dump(Step)
    if MatK!=None:
        BC(MatK, True)
        ppp.dump(MatK)
    if MatM!=None:
        BC(MatM, True)
        ppp.dump(MatM)
    fd.close()
    raise NameError("Exit from Pickle SysMat")

    
