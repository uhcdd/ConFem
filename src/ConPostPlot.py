# ConPostPlot -- 2014-01-13
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from time import *
from numpy import *
#from scipy.linalg import *
#from scipy import sparse
#from scipy.sparse.linalg.dsolve import linsolve
#from scipy.sparse.linalg import aslinearoperator
import matplotlib.pyplot as plt
from os import path as pth
import pickle

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

class ConPostPlot:
    def __init__(self):
        pass
    def Run(self, Name, ElPlotTimes):
        f6=open( Name+".protocol.txt", 'w')
        f1=open( Name+".in.txt", 'r')
        NodeList, ElList, MatList, StepList, NoLabToInd, SecDic = ReadInputFile(f1, f6, False) # read input file # print NodeList.__dict__
        f1.close()
#        N, Mask, Skyline, SDiag, SLen = AssignGlobalDof( NodeList, ElList, MatList, NoLabToInd)          # assign degrees of freedom (dof) to nodes and elements -> see above
        if pth.isfile(Name+".pkl"):
            fd = open(Name+'.pkl', 'rb')
            NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);FlElasticLT=pickle.load(fd);\
                VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
                VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
                TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);Step=pickle.load(fd);Mask=pickle.load(fd);Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);\
                NoLabToNoInd=pickle.load(fd);NoIndToCMInd=pickle.load(fd);ContinuumNodes=pickle.load(fd);CoNoToNoLi=pickle.load(fd);SecDic=pickle.load(fd);LinAlgFlag=pickle.load(fd);
            fd.close()
#            if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes ) # for search purposes, e.g. for EFG or aggregates or embedded truss elements
#            else:                     CoorTree = None
        else: 
            raise NameError ("PickeLoad: cannot read data")
        #
        if pth.isfile(Name+".opt.txt"):                     # read options file if there is any 
            f4=open( Name+".opt.txt", 'r')
            WrNodes, Lines, ReDes, MaxType = ReadOptionsFile(f4, NodeList, NoLabToInd, NoIndToCMInd)
            f4.close()
    
#        if pth.isfile(Name+".pkl"):                 # read restart file if there is one
#            fd = open(Name+'.pkl', 'rb')
#            NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);Flag=pickle.load(fd);\
#                VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
#                VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
#                TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);i=pickle.load(fd);Mask=pickle.load(fd);Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);
#            fd.close()
#            fd.close()
#        else: raise NameError ("PostFem: cannot read data")
    
        f2=open( Name+".elemout.txt", 'r')                  #
#        RC = FinishAllStuff(True, Name, ElList, NodeList, MatList, f2, VecU, WrNodes, None, ElPlotTimes)
        RC = FinishAllStuff(True, Name, ElList, NodeList,NoIndToCMInd, MatList, f2, VecU, WrNodes, None, ElPlotTimes)
        f2.close()
        f6.close

if __name__ == "__main__":
    Case = 0
    if Case == 0:
        ElPlotTimes =  ['1.0000','2.0','3.0','4.0']
        ElPlotTimes =  ['20.0']
        ElPlotTimes =  ["1.0"]
        Name="../_DataBeams/E3-02_Tens"                         # input data name
#        Name = "../DataSpecimen/One3D/WillamsTest"
        Name="../DataExamples/E10-02/ElasticLT/E8-01du02"                
#        Name="../DataExamples/E09/E7-04"                    
#        Name="../DataExamples/E09/E9-05c"                    
#        Name="../DataExamples/E08-03/MicroPl/E6-03"
        Name="C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1"
        Name="../_DataTmp/Deep_beam_AacNL"
        PostFem_ = ConPostPlot()
        PostFem_.Run(Name, ElPlotTimes)
    elif Case == 1:
        LogName="../LogFiles"                             # to log temporary data
        PlotResiduals( LogName, "log", 240, 99, 1.0e+0)
    elif Case == 2:
        LogName="../LogFiles"                             # to log temporary data
        WriteResiduals( LogName, 1, 8, 100)
    print('finish')
