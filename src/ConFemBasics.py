# ConFemBasics -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
import sys
import os
import numpy as np
from numpy.linalg import norm, det, cond, eigvals # ,inv
from math import pi, sqrt
try:
    from ConFemElemC import *
    ConFemElemCFlag = True    
except ImportError:
    ConFemElemCFlag = False    

_BeamsBern_         = ['B23',  'B23E', 'BAX23']
_BeamsBernEmbedded_ = ['B23I', 'B23EI','BAX23I','BAX23EI']
_BeamsBernAll_      = _BeamsBern_ + _BeamsBernEmbedded_
_BeamsTimo_         = ['B21', 'B21E', 'BAX21', 'BAX21E']
_BeamsTimoEmbedded_ = ['BAX21EI']
_BeamsTimoAll_      = _BeamsTimo_ + _BeamsTimoEmbedded_
_BeamsAll_          = _BeamsBernAll_ + _BeamsTimoAll_

_Truss2D_           = [ "T2D2",  "T2D3",  "TAX2", "TAX3"]
_Truss2DEmbedded_   = [ "T2D2I", "T2D3I", "TAX2I", "TAX3I" ]
_Truss2DAll_        = _Truss2D_ + _Truss2DEmbedded_

_Length2D_          = _Truss2DAll_ + _BeamsAll_
_Length1D_          = ["T1D2"]
_Bond2D_            = ["Bond2D2", "Bond2D3", "BondAX2", "BondAX3"]
_FourNodes2D_     = ['CPS4','CPE4','CAX4','CPS4S','CPE4S']
_Shells3D_        = ['SH4','SH3']
_Slabs2D_         = ['SB3']
_TriNodes_        = ["CPS3", 'CPS3S', "CPE3", 'CPE3S', 'CPS6', 'CPE6', 'CPS6S', 'CPE6S', "SB3", "SH3"]
_SpringsAll_      = ["S1D2","S2D6"]

__BeamsBern__ = [[i] for i in _BeamsBern_]
__BeamsTimo__ = [[i] for i in _BeamsTimo_]
__Length2D__   = [[i] for i in _Length2D_]
__Length1D__   = [[i] for i in _Length1D_]
__SpringsAll__ = [[i] for i in _SpringsAll_]
__TriNodes__   = [[i] for i in _TriNodes_]

SamplePoints={}
SampleWeight={}
# Sample points and weighting factors for Gauss Quadrature 1D
SamplePoints[0,0,0]=[0.,0.,0.]
SampleWeight[0,0,0]= 2.
SamplePoints[0,1,0]=[-0.577350269189626, 0., 0.]
SamplePoints[0,1,1]=[ 0.577350269189626, 0., 0.]
SampleWeight[0,1,0]= 1.
SampleWeight[0,1,1]= 1.
SamplePoints[0,2,0]=[-0.774596669241483, 0., 0.]
SamplePoints[0,2,1]=[ 0., 0., 0.]
SamplePoints[0,2,2]=[ 0.774596669241483, 0., 0.]
SampleWeight[0,2,0]= 0.555555555555556
SampleWeight[0,2,1]= 0.888888888888888
SampleWeight[0,2,2]= 0.555555555555556
# Sample points and weighting factors for Gauss Quadrature 2D
SamplePoints[1,0,0]=[0.,0.,0.]
SampleWeight[1,0,0]= 4. # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SamplePoints[1,1,0]=[-0.577350269189626,-0.577350269189626, 0.]
SamplePoints[1,1,1]=[-0.577350269189626, 0.577350269189626, 0.]
SamplePoints[1,1,2]=[ 0.577350269189626,-0.577350269189626, 0.]
SamplePoints[1,1,3]=[ 0.577350269189626, 0.577350269189626, 0.]
SampleWeight[1,1,0]= 1.
SampleWeight[1,1,1]= 1.
SampleWeight[1,1,2]= 1.
SampleWeight[1,1,3]= 1.
# Sample points and weighting factors for Gauss Quadrature 3D
SamplePoints[2,0,0]=[0.,0.,0.]
SampleWeight[2,0,0]= 8. # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SamplePoints[2,1,0]=[-0.577350269189626,-0.577350269189626, -0.577350269189626]
SamplePoints[2,1,1]=[-0.577350269189626,-0.577350269189626,  0.577350269189626]
SamplePoints[2,1,2]=[-0.577350269189626, 0.577350269189626, -0.577350269189626]
SamplePoints[2,1,3]=[-0.577350269189626, 0.577350269189626,  0.577350269189626]
SamplePoints[2,1,4]=[ 0.577350269189626,-0.577350269189626, -0.577350269189626]
SamplePoints[2,1,5]=[ 0.577350269189626,-0.577350269189626,  0.577350269189626]
SamplePoints[2,1,6]=[ 0.577350269189626, 0.577350269189626, -0.577350269189626]
SamplePoints[2,1,7]=[ 0.577350269189626, 0.577350269189626,  0.577350269189626]
SampleWeight[2,1,0]= 1.
SampleWeight[2,1,1]= 1.
SampleWeight[2,1,2]= 1.
SampleWeight[2,1,3]= 1.
SampleWeight[2,1,4]= 1.
SampleWeight[2,1,5]= 1.
SampleWeight[2,1,6]= 1.
SampleWeight[2,1,7]= 1.
# Sample points and weighting factors for triangular -elements, vgl. Zienk I, Table 8.2
SamplePoints[3,0,0]=[ 0.333333333333333, 0.333333333333333, 0.333333333333333]
SampleWeight[3,0,0]= 1.
SamplePoints[3,1,0]=[ 0.5,               0.5,               0.0]
SamplePoints[3,1,1]=[ 0.0,               0.5,               0.5]
SamplePoints[3,1,2]=[ 0.5,               0.0,               0.5]
SampleWeight[3,1,0]= 0.333333333333333
SampleWeight[3,1,1]= 0.333333333333333
SampleWeight[3,1,2]= 0.333333333333333
SamplePoints[3,2,0]=[ 0.333333333333333, 0.333333333333333, 0.333333333333333]
SamplePoints[3,2,1]=[ 0.2,               0.6,               0.2]
SamplePoints[3,2,2]=[ 0.2,               0.2,               0.6]
SamplePoints[3,2,3]=[ 0.6,               0.2,               0.2]
SampleWeight[3,2,0]=-0.562500000000000
SampleWeight[3,2,1]= 0.520833333333333
SampleWeight[3,2,2]= 0.520833333333333
SampleWeight[3,2,3]= 0.520833333333333
# specht, int. j. num. meth. eng., 1988, 705
SamplePoints[3,3,0]=[ 0.666666666666667, 0.166666666666667, 0.166666666666667]
SamplePoints[3,3,1]=[ 0.166666666666667, 0.666666666666667, 0.166666666666667]
SamplePoints[3,3,2]=[ 0.166666666666667, 0.166666666666667, 0.666666666666667]
SampleWeight[3,3,0]= 0.333333333333333
SampleWeight[3,3,1]= 0.333333333333333
SampleWeight[3,3,2]= 0.333333333333333
# Sample points and weighting factors for shells SH4
SamplePoints[4,0,0] =[-0.577350269189626,-0.577350269189626, 0.]
SamplePoints[4,0,1] =[-0.577350269189626, 0.577350269189626, 0.]
SamplePoints[4,0,2] =[ 0.577350269189626,-0.577350269189626, 0.] 
SamplePoints[4,0,3] =[ 0.577350269189626, 0.577350269189626, 0.] 
SampleWeight[4,0,0] = 2.
SampleWeight[4,0,1] = 2. 
SampleWeight[4,0,2] = 2.
SampleWeight[4,0,3] = 2.
SamplePoints[4,1,0] =[-0.577350269189626,-0.577350269189626, -0.861136311594053] #-0.577350269189626]
SamplePoints[4,1,1] =[-0.577350269189626,-0.577350269189626, -0.339981043584856] #0.577350269189626]
SamplePoints[4,1,2] =[-0.577350269189626,-0.577350269189626,  0.339981043584856] #-0.577350269189626]
SamplePoints[4,1,3] =[-0.577350269189626,-0.577350269189626,  0.861136311594053] #0.577350269189626]
SamplePoints[4,1,4] =[-0.577350269189626, 0.577350269189626, -0.861136311594053] #-0.577350269189626]
SamplePoints[4,1,5] =[-0.577350269189626, 0.577350269189626, -0.339981043584856] #0.577350269189626]
SamplePoints[4,1,6] =[-0.577350269189626, 0.577350269189626,  0.339981043584856] #-0.577350269189626]
SamplePoints[4,1,7] =[-0.577350269189626, 0.577350269189626,  0.861136311594053] #0.577350269189626]
SamplePoints[4,1,8] =[ 0.577350269189626,-0.577350269189626, -0.861136311594053] #-0.577350269189626]
SamplePoints[4,1,9] =[ 0.577350269189626,-0.577350269189626, -0.339981043584856] #0.577350269189626]
SamplePoints[4,1,10]=[ 0.577350269189626,-0.577350269189626,  0.339981043584856] #-0.577350269189626]
SamplePoints[4,1,11]=[ 0.577350269189626,-0.577350269189626,  0.861136311594053] #0.577350269189626]
SamplePoints[4,1,12]=[ 0.577350269189626, 0.577350269189626, -0.861136311594053] #-0.577350269189626]
SamplePoints[4,1,13]=[ 0.577350269189626, 0.577350269189626, -0.339981043584856] #0.577350269189626]
SamplePoints[4,1,14]=[ 0.577350269189626, 0.577350269189626,  0.339981043584856] #-0.577350269189626]
SamplePoints[4,1,15]=[ 0.577350269189626, 0.577350269189626,  0.861136311594053] #0.577350269189626]
SampleWeight[4,1,0] = 0.347854845137454
SampleWeight[4,1,1] = 0.652145154862546 
SampleWeight[4,1,2] = 0.652145154862546
SampleWeight[4,1,3] = 0.347854845137454
SampleWeight[4,1,4] = 0.347854845137454
SampleWeight[4,1,5] = 0.652145154862546
SampleWeight[4,1,6] = 0.652145154862546
SampleWeight[4,1,7] = 0.347854845137454
SampleWeight[4,1,8] = 0.347854845137454
SampleWeight[4,1,9] = 0.652145154862546
SampleWeight[4,1,10]= 0.652145154862546
SampleWeight[4,1,11]= 0.347854845137454
SampleWeight[4,1,12]= 0.347854845137454
SampleWeight[4,1,13]= 0.652145154862546
SampleWeight[4,1,14]= 0.652145154862546
SampleWeight[4,1,15]= 0.347854845137454
SamplePoints[4,4,0] =[-0.577350269189626,-0.577350269189626, -0.9061798459386640]
SamplePoints[4,4,1] =[-0.577350269189626,-0.577350269189626, -0.5384693101056831]
SamplePoints[4,4,2] =[-0.577350269189626,-0.577350269189626,  0.0]
SamplePoints[4,4,3] =[-0.577350269189626,-0.577350269189626,  0.5384693101056831]
SamplePoints[4,4,4] =[-0.577350269189626,-0.577350269189626,  0.9061798459386640]
SamplePoints[4,4,5] =[-0.577350269189626, 0.577350269189626, -0.9061798459386640]
SamplePoints[4,4,6] =[-0.577350269189626, 0.577350269189626, -0.5384693101056831]
SamplePoints[4,4,7] =[-0.577350269189626, 0.577350269189626,  0.0]
SamplePoints[4,4,8] =[-0.577350269189626, 0.577350269189626,  0.5384693101056831]
SamplePoints[4,4,9] =[-0.577350269189626, 0.577350269189626,  0.9061798459386640]
SamplePoints[4,4,10]=[ 0.577350269189626,-0.577350269189626, -0.9061798459386640]
SamplePoints[4,4,11]=[ 0.577350269189626,-0.577350269189626, -0.5384693101056831]
SamplePoints[4,4,12]=[ 0.577350269189626,-0.577350269189626,  0.0]
SamplePoints[4,4,13]=[ 0.577350269189626,-0.577350269189626,  0.5384693101056831]
SamplePoints[4,4,14]=[ 0.577350269189626,-0.577350269189626,  0.9061798459386640]
SamplePoints[4,4,15]=[ 0.577350269189626, 0.577350269189626, -0.9061798459386640]
SamplePoints[4,4,16]=[ 0.577350269189626, 0.577350269189626, -0.5384693101056831]
SamplePoints[4,4,17]=[ 0.577350269189626, 0.577350269189626,  0.0]
SamplePoints[4,4,18]=[ 0.577350269189626, 0.577350269189626,  0.5384693101056831]
SamplePoints[4,4,19]=[ 0.577350269189626, 0.577350269189626,  0.9061798459386640]
SampleWeight[4,4,0] = 0.2369268850561891
SampleWeight[4,4,1] = 0.4786286704993665 
SampleWeight[4,4,2] = 0.5688888888888889
SampleWeight[4,4,3] = 0.4786286704993665
SampleWeight[4,4,4] = 0.2369268850561891
SampleWeight[4,4,5] = 0.2369268850561891
SampleWeight[4,4,6] = 0.4786286704993665
SampleWeight[4,4,7] = 0.5688888888888889
SampleWeight[4,4,8] = 0.4786286704993665
SampleWeight[4,4,9] = 0.2369268850561891
SampleWeight[4,4,10]= 0.2369268850561891
SampleWeight[4,4,11]= 0.4786286704993665
SampleWeight[4,4,12]= 0.5688888888888889
SampleWeight[4,4,13]= 0.4786286704993665
SampleWeight[4,4,14]= 0.2369268850561891
SampleWeight[4,4,15]= 0.2369268850561891
SampleWeight[4,4,16]= 0.4786286704993665
SampleWeight[4,4,17]= 0.5688888888888889
SampleWeight[4,4,18]= 0.4786286704993665
SampleWeight[4,4,19]= 0.2369268850561891
# Sample points and weighting factors for shells SH3
SamplePoints[5,1,0] =[.16666666666666666667, .16666666666666666667, -.906179845938664]
SamplePoints[5,1,1] =[.16666666666666666667, .16666666666666666667, -.538469310105683]
SamplePoints[5,1,2] =[.16666666666666666667, .16666666666666666667, 0.]
SamplePoints[5,1,3] =[.16666666666666666667, .16666666666666666667,  .538469310105683]
SamplePoints[5,1,4] =[.16666666666666666667, .16666666666666666667,  .906179845938664]
SamplePoints[5,1,5] =[.66666666666666666667, .16666666666666666667, -.906179845938664]
SamplePoints[5,1,6] =[.66666666666666666667, .16666666666666666667, -.538469310105683]
SamplePoints[5,1,7] =[.66666666666666666667, .16666666666666666667, 0.]
SamplePoints[5,1,8] =[.66666666666666666667, .16666666666666666667,  .538469310105683]
SamplePoints[5,1,9] =[.66666666666666666667, .16666666666666666667,  .906179845938664]
SamplePoints[5,1,10]=[.16666666666666666667, .66666666666666666667, -.906179845938664]
SamplePoints[5,1,11]=[.16666666666666666667, .66666666666666666667, -.538469310105683]
SamplePoints[5,1,12]=[.16666666666666666667, .66666666666666666667, 0.]
SamplePoints[5,1,13]=[.16666666666666666667, .66666666666666666667,  .538469310105683]
SamplePoints[5,1,14]=[.16666666666666666667, .66666666666666666667,  .906179845938664]
# divided by 8 to make unit volume to 1/2, see SH4::         # uhc220827 sums up to 1 which should be correct for unit element with base area 1/2 and t -> [-1..+1]
                                                             # furthermore, weighting factors of SH4 are multiplied by 4/3 to compensate for 3 base points instead of 4
SampleWeight[5,1,0] = 0.31590251340825200000/8.
SampleWeight[5,1,1] = 0.63817156066582133333/8.
SampleWeight[5,1,2] = 0.75851851851851866666/8.
SampleWeight[5,1,3] = 0.63817156066582133333/8.
SampleWeight[5,1,4] = 0.31590251340825200000/8.
SampleWeight[5,1,5] = 0.31590251340825200000/8.
SampleWeight[5,1,6] = 0.63817156066582133333/8.
SampleWeight[5,1,7] = 0.75851851851851866666/8.
SampleWeight[5,1,8] = 0.63817156066582133333/8.
SampleWeight[5,1,9] = 0.31590251340825200000/8.
SampleWeight[5,1,10]= 0.31590251340825200000/8.
SampleWeight[5,1,11]= 0.63817156066582133333/8.
SampleWeight[5,1,12]= 0.75851851851851866666/8.
SampleWeight[5,1,13]= 0.63817156066582133333/8.
SampleWeight[5,1,14]= 0.31590251340825200000/8.
# for reinforcement of shells, see SH4.__init__
SamplePointsRCShell={}
SampleWeightRCShell={}
# Microplane
I21Points=np.array([[1.,             0.,             0.            ],
                 [0.,             1.,             0.            ],
                 [0.,             0.,             1.            ],
                 [0.707106781187, 0.707106781187, 0.            ],
                 [0.707106781187,-0.707106781187, 0.            ],
                 [0.707106781187, 0.            , 0.707106781187],
                 [0.707106781187, 0.            ,-0.707106781187],
                 [0.            , 0.707106781187, 0.707106781187],
                 [0.            , 0.707106781187,-0.707106781187],
                 [0.387907304067, 0.387907304067, 0.836095596749],
                 [0.387907304067, 0.387907304067,-0.836095596749],
                 [0.387907304067,-0.387907304067, 0.836095596749],
                 [0.387907304067,-0.387907304067,-0.836095596749],
                 [0.387907304067, 0.836095596749, 0.387907304067],
                 [0.387907304067, 0.836095596749,-0.387907304067],
                 [0.387907304067,-0.836095596749, 0.387907304067],
                 [0.387907304067,-0.836095596749,-0.387907304067],
                 [0.836095596749, 0.387907304067, 0.387907304067],
                 [0.836095596749, 0.387907304067,-0.387907304067],
                 [0.836095596749,-0.387907304067, 0.387907304067],
                 [0.836095596749,-0.387907304067,-0.387907304067]])
I21Weights=np.array([0.0265214244093,0.0265214244093,0.0265214244093,
                  0.0199301476312,0.0199301476312,0.0199301476312,0.0199301476312,0.0199301476312,0.0199301476312,
                  0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487,0.0250712367487])

ZeroD = 1.e-9                                                           # Smallest float for division by Zero
StrNum0 = 3.*sqrt(3.)/2.                                                # strange number used by ConFemMat::IsoDamage

def Echo(*args):
    print(args[0])
    if len(args)>1: print(args[0], file=args[1])
#def EchoM( *args ):                                                     # formatted output of matrix -- 0: title, 1: matrix, 2: format length
def EchoM( M, Label=" ", L=12, ff=sys.stdout ):
    n = M.shape[0]                                                      # number of rows
    m = M.shape[1]                                                      # number of columns
    fI = '%'+str(L)+'i'
    fE = '%'+str(L)+'.'+str(L-8)+'e'
    ff.write('\n%s\n'%Label)
    ff.write(L*' ')
    for j in range(m): ff.write(fI%(j))
    ff.write('\n')
    for i in range(n):
        ff.write(fI%i)
        for j in range(m): ff.write(fE%(M[i,j]))
        ff.write('\n')
    ff.write('\n')
def EchoSys( f6, SymSys, LinAlgFlag, ConFemMatCFlag, ConFemElemCFlag, N, ElList, ActiveNodes, SLen, MatList, Element, ElemDataAll):
        Echo(f"symmetric system {SymSys}, used LinAlg2 {LinAlgFlag}, used ConFemMatC {ConFemMatCFlag}, used CaeFemElemC {ConFemElemCFlag}", f6)
        Echo(f"ndofs {N:d}, elems {len(ElList):d}, active nodes {ActiveNodes:d}, abs.size {1.*SLen/1024:.2f} MItems, rel.size {SLen/(N**2):.4f} (abs.size/(ndofs**2)", f6)
#        if len(Element.CharLengthCPE3)>0: 
#            Element.CharLengthCPE3.sort(key=lambda t: t[0]) 
#            Echo(f'used element CP*3 [char. length, scaling type, scaling factor] from el {Element.CharLengthCPE3[0][1]:d},{Element.CharLengthCPE3[0][0]:.3f}, to el {Element.CharLengthCPE3[-1][1]:d},{Element.CharLengthCPE3[-1][0]:.3f}', f6)
        if "CPS3" in ElemDataAll:
            ElDat =  ElemDataAll["CPS3"]
            ElDat.sort(key=lambda t: t[0])
            Echo(f'used element CP*3 [char. length] from el {ElDat[0][1]:d},{ElDat[0][0]:.3f}, to el {ElDat[-1][1]:d},{ElDat[-1][0]:.3f}', f6)
#        if len(Element.CharLengthCPE4)>0:                                   # see CPE4:Ini2
#            Element.CharLengthCPE4.sort(key=lambda t: t[0]) 
#            Echo(f'used element CP*4 [char. length, scaling type, scaling factor] span from el {Element.CharLengthCPE4[0][0]:d},{Element.CharLengthCPE4[0][1]:.3f},{Element.CharLengthCPE4[ 0][2]:d},{Element.CharLengthCPE4[0][3]:.3f} to el {Element.CharLengthCPE4[-1][0]:d},{Element.CharLengthCPE4[-1][1]:.3f},{Element.CharLengthCPE4[-1][2]:d},{Element.CharLengthCPE4[-1][3]:.3f}', f6)
        if "CPS4" in ElemDataAll:
            ElDat =  ElemDataAll["CPS4"]
            ElDat.sort(key=lambda t: t[0])
            Echo(f'used element CP*4 [char. length, scaling type, scaling factor] span from el {ElDat[0][0]:d},{ElDat[0][1]:.3f},{ElDat[0][2]:d},{ElDat[0][3]:.3f} to el {ElDat[-1][0]:d},{ElDat[-1][1]:.3f},{ElDat[-1][2]:d},{ElDat[-1][3]:.3f}', f6)
        for m in MatList:
            m_ = MatList[m]
            if m_.Used:
                if m_.Type in ['ISODAMAGE']: #,'MICRODAMAGE']:
                    Echo(f'used material {m:s}, spec. inertial mass {m_.Density:7.2e}, type {m_.Type:s}, fct {m_.fct:.2f}, DStrength {m_.DStrength:.4f}, reg.type {m_.RType:d}, spec.cr.en {m_.SpecCrEn:7.2e}, cr.bw {m_.bw:.2f}', f6)
                else: 
                    Echo(f'used material {m:s}, spec. inertial mass {m_.Density:7.2e}, type {m_.Type:s}', f6)
                if m_.Type in ['ISODAMAGE'] and m_.RegPar == 2:
                    Echo(f'IsoDam interpolated crack band reg data small el {m_.SmallElReg[0]:12.6e}, {m_.SmallElReg[1]:12.6e}, {m_.SmallElReg[2]:12.6e}, {m_.SmallElReg[3]:12.6e}, {m_.SmallElReg[4]:12.6e}', f6) # SmallElReg = [min(CrL),max(CrL),sep,a[0],a[1]]
                    Echo(f'IsoDam interpolated crack band reg data large el {m_.LargeElReg[0]:12.6e}, {m_.LargeElReg[1]:12.6e}, {m_.LargeElReg[2]:12.6e}, {m_.LargeElReg[3]:12.6e}, {m_.LargeElReg[4]:12.6e}, {m_.LargeElReg[5]:12.6e}', f6) # LargeElReg = [min(CrL),max(CrL),sep,a[0],a[1],a[2]]

def DirFilFromName( Name ):
    FilName = os.path.basename( Name )
    FilName = FilName.split(".")
    DirName = os.path.dirname( Name )
    return DirName, FilName[0]

def FindIndexByLabel(NodeList, Key):                                    # find index from label -- must not be NodeList
    if len(NodeList)==0: raise NameError ("no item defined")
    for i, Node in enumerate(NodeList):                                 # loop over all nodes
        if Node.Label == Key: return i
    if i==(len(NodeList)-1):  return -1                                 # no node found

def PrinC( xx, yy, xy):                                                 # calculation of principal values
    if ZeroD < abs(xy):
        h1 = xx+yy
        h2 = np.sqrt((xx-yy)**2+4*xy**2);
        s1 = .5*h1+.5*h2
        s2 = .5*h1-.5*h2
        h = (s1-xx)/xy
        L = np.sqrt(h**2+1)
        n11 = 1./L
        n12 = h/L
        h = (s2-xx)/xy
        L = np.sqrt(h**2+1)
        n21 = 1./L
        n22 = h/L 
    else:
        s1 = xx
        n11 = 1.
        n12 = 0
        s2 = yy
        n21 = 0
        n22 = 1. 
    return( [s1, n11, n12, s2, n21, n22] ) 

#def PrinCLT(v1, v2, v3):                                                # princial values, largest, direction
#    pp = PrinC( v1, v2, v3)                                             # principal values
#    if pp[0]>pp[3]:
#        if abs(pp[1])>ZeroD: phi = np.arctan(pp[2]/pp[1])  # direction of larger tensile principal stress
#        else: phi = pi/2
#        pig = pp[0]                                     # larger principal value
#        pig_= pp[3]                                     # lower principal value
#    else:
#        if abs(pp[4])>ZeroD: phi = np.arctan(pp[5]/pp[4])
#        else: phi = pi/2
#        pig = pp[3]                                     # larger principal value
#        pig_= pp[0]                                     # lower principal value
#    return pig, phi, pig_
def PrinCLT_(vx, vy, vxy):                                
    a = 0.5*(vx-vy)
    b = np.sqrt(a**2+vxy**2)
    if b>ZeroD:
        c = 0.5*(vx+vy)
        v1 = c + b
        v2 = c - b
        if vxy<0.: si = -1.
        else:      si = 1.
        ang = 0.5*si*np.arccos(a/b)
    else:
        v1 = vx
        v2 = vy
        ang = 0.
    if v1>v2: return v1, ang, v2
    else:     return v2, ang+0.5*pi, v1

def AssignGlobalDof( NodeList, ElList, MatList, SecDict, NoIndToCMInd):     # assign dof indices -- called by ConFem
    ElemDataAll = {}                                                        # collection informative data for all element types
    for El in ElList:
        if El.Type not in ElemDataAll: ElemDataAll[El.Type] = []
    for El in ElList:
        ElemData = El.Ini2( NodeList,NoIndToCMInd, MatList, SecDict)        # initialization of data depending on Sequence of nodes in NodeList which has been determined by Cuthill McKee
        El.ElemDimData( NodeList,NoIndToCMInd )
        El.Ini3( NodeList,NoIndToCMInd )
        ElemDataAll[El.Type] += [ElemData]
    Index = 0
    for No in NodeList:                                                     # loop over nodes
        No.GlobDofStart = Index                                             #
        Index = Index + len(No.DofT)                                        # Node.DofT filled during data input / element initialization
    for El in ElList:                                                       # loop over all elements to assign global dof index to element table if element shares this dof (which is not mandatory)
        for j in range(El.nNod):                                            # loop over nodes of element
            No = NodeList[ NoIndToCMInd[El.Inzi[j]] ]                       # global node of element
            set1 = El.DofT[j]                                               # set of dof types of element node
            set2 = No.DofT                                                  # set of dof types of global node - not the same as element node due to contributions of other elements
            for k, k0 in enumerate(set1):                                   # loop over number of dofs of element node
                for l, l0 in enumerate(set2):                               # loop over number of dofs of global node
                    if k0==l0:                                              # element node dof found as global node dof
                        El.DofI[j,k] = No.GlobDofStart + l                  # assign global dof index to element table
                        break
#    Skyline = np.zeros((Index), dtype=np.int)                               # switch back for win10, otherwise this will not work for large fields / large number of dofs (est. > 10**5)
#    SDiag = np.zeros((Index), dtype=np.int64)
    Skyline = np.zeros((Index), dtype=int)                               # switch back for win10, otherwise this will not work for large fields / large number of dofs (est. > 10**5)
    SDiag = np.zeros((Index), dtype=int)
    for El in ElList:                                                       # upper right skyline of system matrix
        minL = []
        for j in El.DofI:                                                   # global dof index for each element node j -> list for each node
            try:    jc = j.compressed()                                     # for masked array, see e.g.B23E
            except: jc = j
            minL += [min( set(jc).difference(set([-1])) )]                  # subtracts -1 from jc, if it there (to get rid of -1 which might be there for some element types) -> minimum dof index
        min_ = min(minL)                                                    # minimum global dof index for element
        for j in El.DofI:
            for k in j:                                                     # find skyline for every global dof
                if k>=0:
                    sl1 = k - min_
                    sl2 = min(k,sl1)+1
                    sl3 = Skyline[k]
                    if sl2>sl3: Skyline[k] = sl2
    for i in range(Index-1): 
        SDiag[i+1]=SDiag[i]+Skyline[i]                                      # addresses of diagonal entries of system matrix in system vector
    SLen = SDiag[Index-1]+Skyline[Index-1]                                  # total length of system matrix stored as skyline vector
    return Index, Skyline, SDiag, SLen, ElemDataAll                         # return last global dof index

def IPMassMatrix( Elem, mm, j, r,s,t,f, mmat, B23Flag):
    det = Elem.JacoD(r, s, t)
    if B23Flag and Elem.Type in ["B23E","BAX23EI"]:                          # ordinary approach yields negative mass matrix entries for this element type, see contimech.tex -- Elementtechnologie, B23X-Element
        if j == 0:                                                          # following analytically integrated
            L = 2. * Elem.Geom[0, 0]
#            df_ = det * mm[0, 0] * Elem.Geom[0, 0]  # Elem.Geom[0, 0] --> L/2
            df_ = f * mm[0,0]                                                # f includes Jacobian for integration, integration weighting factor, radius in case of axisymmetry, nRebarEl
            mmat += np.array([[0.3333333333 * df_, 0., 0., 0., 0., 0., 0.],
                             [0., df_, 0., 0., 0., 0., 0.],
                             [0., 0., df_ * 0.4761904760e-02 * L * L, 0., 0., 0., 0.],
                             [0., 0., 0., 1.3333333333 * df_, 0., 0., 0.],
                             [0., 0., 0., 0., 0.3333333333 * df_, 0., 0.],
                             [0., 0., 0., 0., 0., df_, 0.],
                             [0., 0., 0., 0., 0., 0., df_ * 0.4761904760e-02 * L * L]])
        else:
            pass
    else:
        N = Elem.FormN(r, s, t)  # shape function
        mmat += det*f*Elem.Geom[1,0] * np.dot(np.transpose(N), np.dot(mm, N))  # element mass matrix
    return 0

def RSTf( CalcType,Elem, InT,nint, j):
    if Elem.ShellRCFlag and j >= Elem.nIntLi:  # for reinforcement layers
        eSet = Elem.Set
        r = SamplePointsRCShell[eSet, InT, nint - 1, j][0]
        s = SamplePointsRCShell[eSet, InT, nint - 1, j][1]
        t = SamplePointsRCShell[eSet, InT, nint - 1, j][2]
        f = Elem.Geom[0, 0] * SampleWeightRCShell[eSet, InT, nint - 1, j]  # weighting factor
    else:
        r = SamplePoints[InT, nint - 1, j][0]
        s = SamplePoints[InT, nint - 1, j][1]
        t = SamplePoints[InT, nint - 1, j][2]
        f = Elem.Geom[0, 0] * SampleWeight[InT, nint - 1, j]  # weighting factor; Geom[0,0] might hold Jacobian
    if Elem.dim == 4:  # axisymmetric CAX4
        x = Elem.X[j]
        f = 2. * x * pi * SampleWeight[InT, nint - 1, j]  # Jacobian for integration comes into play with det
    if Elem.dim in [5, 12, 13, 96]:  # axisymmetric TAX, BAX23, BAX21, Bond
        x = Elem.X[j]
        if CalcType in [10, 11] and Elem.Type in ['BAX21E', 'BAX21EI', 'BAX23E', 'BAX23EI']:  # for BAX21E
            x = np.mean(Elem.X)  # mass vector will have zero elements otherwise ???
        w = SampleWeight[InT, nint - 1, j]
        f = Elem.Geom[0, 0] * 2. * x * pi * w  # Geom[0,0] -> Jacobian for integration
    #
    f = f * Elem.nRebarEl  # nRebar (default 1) scaling for embedded beams, trusses for number of bars embedded
    return r, s, t, f

#@profile
def IntForces( MatList, ElList, Dt, VecC,VecU,VecD,VecS,VecT, VecI, MatK,MatG,KVecU,KVecL,SDiag, CalcType, ff, NLGeom, SymS, Buckl, eqiter, nTmp=0):
                                # CalcType 0: check system 1: internal forces only 2: internal forces and tangential stiffness matrix 10: mass matrix only
    VecI[:] = 0.
    sysMass = 0                      
    for i, Elem in enumerate(ElList):
        if Elem.Active:
            if NLGeom and not Elem.NLGeomI: raise NameError("ConFem::ConFemBasics:Intforces: NLGEOM not yet implemented for this element typ",Elem.Type)
            # Initialize variables
            nIntL= Elem.nIntL
            nint = Elem.nInt                                            # integration order
            nnod = Elem.nNod                                                # number of nodes per element
            matN = Elem.MatN                                                # name of material
            InT  = Elem.IntT                                                # integration type of element
            uvec = np.zeros((Elem.DofE), dtype=np.double)                   # element displacements  # DofE: number of dofs for whole element
            dvec = np.zeros((Elem.DofE), dtype=np.double)                   # element displacement increments from last time step
            xvec = np.zeros((Elem.DofE), dtype=np.double)                   # element displacement iteration increments
            tvec = np.zeros((Elem.DofEini), dtype=np.double)                # element temperatures
            svec = np.zeros((Elem.DofEini), dtype=np.double)                # element temperature increments from last time step
            CFlag, CFlag1 = False, False
            if ConFemElemCFlag and Elem.Type in ['SH4']:         CFlag = True
            if ConFemElemCFlag and Elem.Type in ['CPS4','CPE4']: CFlag1 = True
            if nTmp==0:
                if Elem.Type in _BeamsAll_: tempI, dtmpI = np.zeros((2), dtype=float), np.zeros((2), dtype=float)
                else:                                                                       tempI, dtmpI = 0., 0.

            # Determination of internal forces and element stiffness ans mass matrix
            ndof0, ndof0_ = 0, 0
            for j in range(nnod):                                           # node loop
                for k in range(Elem.DofN[j]):                               # dof loop
                    kk = Elem.DofI[j,k]                                     # global index of local dof k of local node j
                    uvec[ndof0+k] = VecU[kk]                                # element displacements from global displacements
                    dvec[ndof0+k] = VecU[kk]-VecC[kk]                       # element displacement time step increments
                    xvec[ndof0+k] = VecD[kk]                                # element displacement iteration increments
                ndof0 = ndof0 + Elem.DofN[j]                                # update entry index for element displacements
                if nTmp>0:
                    for k in range(Elem.DofNini[j]):                        # dof loop
                        kk = Elem.DofI[j,k]                                 # global index of local dof k of local node j
                        tvec[ndof0_+k]= VecT[kk]                            # element node temperatures
                        svec[ndof0_+k]= VecT[kk]-VecS[kk]                   # element node temperature time step increments
                    ndof0_= ndof0_+ Elem.DofNini[j]                         # update entry index for element displacements, initial values
            if NLGeom and Elem.NLGeomI and not Elem.RotG: Elem.UpdateCoord(uvec, dvec) # update current element coordinates for all elements but sh4
            if Elem.Rot: 
                uvec = np.dot(Elem.Trans,uvec)                              # transform global to local displacements
                dvec = np.dot(Elem.Trans,dvec)                              # transform global to local displacement increments
            if NLGeom and Elem.NLGeomI and Elem.RotG: 
                Elem.UpdateCoord(uvec, dvec)                                # update current element coordinates for shell element SH4 (must be done here due to variable number of dofs per node before)
                Elem.UpdateElemData()
            Elem.Nodisp[:] = uvec[:]

            # SDA internal equilibrium integration point loop
            if Elem.Type in ["CPS4S", "CPE4S", "CPS3S", "CPE3S", "C3D8S", "CPE6S", "CPS6S"]: SDAFlag, NWiter = True, 6 # NWiter = 3, 7 iteration limiter  is currently trial and error
            else:                                                                            SDAFlag, NWiter = False, 1
            for witer in range(NWiter):                                     # iteration for internal equilibrium of SDA;
                rvec = np.zeros((Elem.DofEini), dtype=np.double)            # element internal forces # dofEini: as DofE may be subject to change, e.g. for sh4
                kmat = np.zeros((Elem.DofEini, Elem.DofEini), dtype=np.double)  # element stiffness
                if NLGeom: gmat = np.zeros((Elem.DofEini, Elem.DofEini), dtype=np.double)# geometric stiffness
                #
                if SDAFlag: #Elem.Type in ["CPS4S","CPE4S","CPS3S","CPE3S","C3D8S","CPE6S","CPS6S"]:
                    fw, Kww, Kuw, Kwu = Elem.CrIni( eqiter, witer, xvec, ff)
                # integration point loop
                # integration jacobian is either from Geom[0,0] which goes into f or from det from FormB
                for j in range(nIntL):                                      # build element stiffness with integration loop
                    if False:
                       if Elem.ShellRCFlag and j>=Elem.nIntLi:                     # for reinforcement layers
                            eSet = Elem.Set
                            r = SamplePointsRCShell[eSet,InT,nint-1,j][0]
                            s = SamplePointsRCShell[eSet,InT,nint-1,j][1]
                            t = SamplePointsRCShell[eSet,InT,nint-1,j][2]
                            f = Elem.Geom[0,0]*SampleWeightRCShell[eSet,InT,nint-1,j] # weighting factor
                       else:
                            r = SamplePoints[InT,nint-1,j][0]
                            s = SamplePoints[InT,nint-1,j][1]
                            t = SamplePoints[InT,nint-1,j][2]
                            f = Elem.Geom[0,0]*SampleWeight[InT,nint-1,j]           # weighting factor; Geom[0,0] might hold Jacobian
                       if Elem.dim==4:                                             # axisymmetric CAX4
                            x = Elem.X[j]
                            f = 2. * x * pi *  SampleWeight[InT,nint-1,j]           # Jacobian for integration comes into play with det
                       if Elem.dim in [5,12,13,96]:                                # axisymmetric TAX, BAX23, BAX21, Bond
                            x = Elem.X[j]
                            if CalcType in [10,11] and Elem.Type in ['BAX21E','BAX21EI','BAX23E', 'BAX23EI']:  # for BAX21E
                                x = np.mean(Elem.X)                                 # mass vector will have zero elements otherwise ???
                            w = SampleWeight[InT,nint-1,j]
                            f = Elem.Geom[0,0] * 2.*x*pi * w                        # Geom[0,0] -> Jacobian for integration

                       f = f*Elem.nRebarEl                                         # nRebar (default 1) scaling for embedded beams, trusses for number of bars embedded
                    else:
                        r,s,t,f = RSTf( CalcType,Elem, InT,nint, j)
                    #
                    if CalcType in [10,11]:                                      # mass matrix -- extra contributions to continuum in case of ShellRCFlag - index 11 for explicit
                        mm  = MatList[matN].Mass(Elem)                           # element mass in integration point
                        IPMassMatrix(Elem, mm, j,r,s,t,f, kmat, (CalcType==11) )
                        continue

                    # preliminaries
                    if Elem.RegType in [1,4]:                                   # regularization types 0 no, 1 gradient, 2 crack band width, 3 SDA, 4 phase field
                        B, BN, det, TM = Elem.FormB(r,s,t, NLGeom)              # shape function derivative for regularization 
                    else:
                        if CFlag: 
                            B = np.zeros((6,20), dtype=float)
                            Data_ = np.zeros((1), dtype=float)
                            TM = np.zeros((6,6), dtype=float)
                            _ = SH4FormBC( r, s, t, B, Data_, Elem.XX, Elem.a, Elem.Vn, Elem.EdgeDir, Elem.gg[0], Elem.gg[1], TM )
                            det = Data_[0]
                        elif CFlag1:
                            B = np.zeros((3,8), dtype=float)
                            det_ = np.zeros((1), dtype=float)
                            Data_ = np.zeros((1), dtype=float)
                            _ = CPS4FormBC( Elem.X0,Elem.Y0,Elem.X1,Elem.Y1,Elem.X2,Elem.Y2,Elem.X3,Elem.Y3, r, s, t, B, det_, Data_ )
                            det = det_[0]
                        else:
                            B, det, TM = Elem.FormB(r,s,t, NLGeom)              # shape function derivative, Jacobian for integration
                    if det<= 0.:
                        pass
                        raise NameError("ConFemBasics::IntForces: element with <= 0 determinant", Elem.Label, det)
                    epsI = np.dot( B, uvec)                                     # integration point strains
                    dpsI = np.dot( B, dvec)                                     # integration point strain increments

                    if SDAFlag:
                        if not Elem.SDANew: Bw = Elem.FormBw( r, s, t )         # Mw comes here into play
                        else:               Bw = Elem.Bw[j]
                        epsI = epsI - np.dot(Bw,Elem.ww)
                        dpsI = dpsI - np.dot(Bw,Elem.dw)

                    if nTmp>0:
                        T = Elem.FormT(r,s,t)                                   # shape function for temperature interpolation
                        tempI= np.dot( T, tvec)                                 # integration point temperatures
                        dtmpI= np.dot( T, svec)                                 # integration point temperature increments
                    # stresses, material stiffness
                    if Elem.RegType in [1,4]: epsR = np.dot( BN, uvec)          # integration point nonlocal equivalent strains
                    if Elem.RotM:
                        epsI = np.dot(TM,epsI)
                        dpsI = np.dot(TM,dpsI)
                    if Elem.RegType in [1,4]:                                   # 1: gradient damage, 4: phase field
                        sigI, C, sigR, CR, Data = MatList[matN].Sig( ff, CalcType, Dt,i,j, Elem, dpsI, epsI, dtmpI, tempI, epsR)# stress, stress incr, tangential stiffness 
                    else:         
                        sigI, C, Data = MatList[matN].Sig( ff, CalcType, Dt,i,j, Elem, dpsI, epsI, dtmpI, tempI, [])# stress, stress incr, tangential stiffness
                    if CalcType==0:
                        continue                                                # next integration point in case of system check

                    # nodal forces
                    if Elem.RotM:
                        C = np.dot(np.transpose(TM), np.dot(C,TM))
                        sigI = np.dot(np.transpose(TM),sigI)
                    for k in range(len(Data)): Elem.Data[j, k] = Data[k]        # element data (strains, stresses etc.)
                    if Elem.RegType in [1,4]:
                        rvec = rvec + det*f*Elem.Geom[1,0]*(np.dot(np.transpose(B), sigI)+np.dot(np.transpose(BN), sigR)) # element internal forces
                    else:
                        if ConFemElemCFlag:
                            fact = det*f * Elem.Geom[1,0]
                            rc = BTransSig( rvec, fact, B, sigI)
                            if rc!=0: raise NameError("ConFemBasics::Intforces:BTransSig",rc)
                        else:
                            rvec = rvec + det*f * Elem.Geom[1,0] * np.matmul(np.transpose(B),sigI)  # element internal forces

                    if not SDAFlag and CalcType==1:
                        continue                                                # nodal forces only -- will not work with SDA initial convergence control

                    # stiffness matrix
                    f_ = det*f*Elem.Geom[1,0]
                    if Elem.RegType in [1,4]: kmat = kmat + f_ *(np.dot(np.transpose(B), np.dot(C,B))+np.dot(np.transpose(BN), np.dot(CR,BN)))# element stiffness
                    else:
                        if ConFemElemCFlag:
                            X = np.zeros((1,1), dtype=float)                          # currently dummy
                            Data_ = np.zeros((1), dtype=float)
                            rc = BTxCxB( kmat, f_, B, C, Data_, X)
                            if rc != 0: raise NameError("ConFemBasics::Intforces:BTxCxB", rc)
                        else:
                            kmat = kmat + f_ * np.dot(np.transpose(B), np.dot(C,B)) # element stiffness
                    if SDAFlag: #Elem.Type in ["CPS4S","CPE4S","CPS3S","CPE3S","C3D8S","CPE6S","CPS6S"]:
                        Kww = Kww + f_*np.dot(np.transpose(Bw),np.dot(C,Bw))
                        Kuw = Kuw + f_*np.dot(np.transpose(B), np.dot(C,Bw))
                        Kwu = Kwu + f_*np.dot(np.transpose(Bw),np.dot(C,B))
                        fw  = fw  + f_*np.dot(np.transpose(Bw),sigI)        # fw: integrated discontinuity forces from bulk stresses
                    # geometric stiffness matrix
                    if NLGeom:
                        if CFlag:                                           # ConFemElemCFlag and Elem.Type=='SH4':
                            GeomSti = np.zeros((20,20), dtype=float)
                            DataX = np.zeros((1), dtype =float)
                            _ = SH4FormGeomC( r, s, t, GeomSti, DataX, sigI, Elem.gg[0], Elem.gg[1] )
                        else: 
                            GeomSti = Elem.GeomStiff(r,s,t,sigI)
                        if Elem.NLGeomCase==0: gmat = gmat + det*f*Elem.Geom[1,0]                                             *GeomSti# Elem.GeomStiff(r,s,t,sigI)
                        else:                  gmat = gmat + Elem.Geom[1,0]*(SampleWeight[InT,nint-1,j]/SampleWeight[InT,0,0])*GeomSti# Elem.GeomStiff(r,s,t,sigI)
                    # end of integration point loop
                    
                if SDAFlag: #Elem.Type in ["CPS4S","CPE4S","CPS3S","CPE3S","C3D8S","CPE6S","CPS6S"]:
                    Elem.CrFin( fw, Kww, Dt, ff )                           # Elem.fwd, Elem.Dd are built here    
                    kmat = kmat - np.dot( Kuw, np.dot(Elem.Dd,Kwu) )        # modification of bulk stiffness matrix
#                    rvec = rvec + np.dot( Kuw, np.dot(Elem.Dd,Elem.fwd))            # uhc belongs here to strict theory (?) but practically disturbes convergence behavior
                    Elem.Kwu = Kwu
                    if eqiter==0 or norm(Elem.fwd)<1.0e-6: break            # internal SDA equilibrium reached; !!! threshold is currently trial and error / ddw = 0 for eqiter==0
                                                                            # NWiter = 1 if not SDA
            # end if intetration point loop
            if CalcType==0:
                continue                                                    # next element in case of system check
            elif CalcType  in [10,11]:
                sysMass += sum(sum(kmat))

            # Transform local to global values
            if Elem.Rot: 
                rvec = np.dot(Elem.Trans.transpose(),rvec)# transform local to global forces
            if Elem.Rot and CalcType==2: 
                kmat = np.dot(np.dot(Elem.Trans.transpose(),kmat),Elem.Trans)# transform local to global element stiffness
                if NLGeom and Elem.RotG: gmat = np.dot(np.dot(Elem.Trans.transpose(),gmat),Elem.Trans)
                
            # some OPTIONAL hard coded output for error search ####################################################
#            if CalcType==10: # and Elem.Type in ['BondAX3']: #
#            if CalcType==11 and Elem.Label == 633:
#                print('xxx',Elem.Label,Elem.Type)
#                EchoM(kmat,L=12)
#                exit()
#            if Elem.Type in ['B23EI','B23I']:
#                if Elem.Type == 'B23EI': kmat[3,3]=1.0
#                print('XXX\n',kmat)
#                exit()
#                kmat_ = kmat[0:2,0:2]
#                con = "%12.4e"%cond(kmat_)
#                eiV = eigvals(kmat_)
##                EchoM( kmat_, Label='kmat_ '+str(Elem.Label)+' '+con+' '+str(eiV))
#                EchoM( kmat_, Label='kmat_ '+str(Elem.Label)+' '+con+' '+str(eiV), ff=ff)
#            if CalcType==2 and Elem.Type=="Bond2D3":
#                EchoM(kmat,L=16)
#                exit()

#            if CalcType==10:
#                for k in range(kmat.shape[0]):
#                    xxx = 0.
#                    for kk in range(kmat.shape[0]):
#                        xxx += kmat[k,kk]
#                    if xxx<ZeroD:
#                        print('aaa',f)
#                        print("xxx",Elem.Label,Elem.Type,k,xxx)
#                        print('\n',kmat)
#                        raise NameError("xxx")

            #######################################################################################################

            # Fill system vectors (internal forces) and matrices (stiffness matrix)
            ndof0 = 0
            for j0 in range(nnod):                                          # assemble rows
                for k in range(Elem.DofN[j0]):                              # loop over element row dofs
                    kk = Elem.DofI[j0,k]                                    # global dof row index
                    VecI[kk] = VecI[kk] + rvec[ndof0+k]
                ndof0 = ndof0 + Elem.DofN[j0]
            if CalcType>=2:                                                 # asssemble system matrices
                ndof0 = 0
                if MatK!=None:                                              # case sparse matrix
                    for j0 in range(nnod):                                  # assemble rows
                        ndof1 = 0
                        for j1 in range(nnod):                              # assemble columns
                            for k in range(Elem.DofN[j0]):                  # loop over element row dofs
                                kk = Elem.DofI[j0,k]                        # global dof row index
                                for l in range(Elem.DofN[j1]):              # loop over element column dofs
                                    ll = Elem.DofI[j1,l]                    # global dof column index
                                    MatK[kk,ll] = MatK[kk,ll] + kmat[ndof0+k,ndof1+l]
                                    if NLGeom: 
                                        if not Buckl: MatK[kk,ll] = MatK[kk,ll] + gmat[ndof0+k,ndof1+l]
                                        else:         MatG[kk,ll] = MatG[kk,ll] + gmat[ndof0+k,ndof1+l]
                            ndof1 = ndof1 + Elem.DofN[j1]                   # update entry for element dof column index
                        ndof0 = ndof0 + Elem.DofN[j0]
                else:                                                       # case skyline matrix
                    for j0 in range(nnod):                                  # assemble rows
                        ndof1 = 0
                        for j1 in range(nnod):                              # assemble columns
                            for k in range(Elem.DofN[j0]):                  # loop over element row dofs
                                kk = Elem.DofI[j0,k]                        # global dof row index
                                for l in range(Elem.DofN[j1]):              # loop over element column dofs
                                    ll = Elem.DofI[j1,l]                    # global dof column index
                                    if ll>=kk:                              # upper right part column >= row
                                        jj = SDiag[ll] + ll - kk            # index in stiffness vector
                                        KVecU[jj] = KVecU[jj] + kmat[ndof0+k,ndof1+l] # update stiffness vector
                                    elif not SymS:                          # lower left part column < row
                                        jj = SDiag[kk] + kk - ll            # index in stiffness vector
                                        KVecL[jj] = KVecL[jj] + kmat[ndof0+k,ndof1+l] # update stiffness vector
                                    if NLGeom:
                                        if ll>=kk:                          # upper right part column >= row
                                            jj = SDiag[ll] + ll - kk        # index in stiffness vector
                                            KVecU[jj] = KVecU[jj] + gmat[ndof0+k,ndof1+l] # update stiffness vector
                                        elif not SymS:                      # lower left part column < row
                                            jj = SDiag[kk] + kk - ll        # index in stiffness vector
                                            KVecL[jj] = KVecL[jj] + gmat[ndof0+k,ndof1+l] # update stiffness vector
                            ndof1 = ndof1 + Elem.DofN[j1]                   # update entry for element dof column index
                        ndof0 = ndof0 + Elem.DofN[j0]
    
    #    if Buckl:
    #        sys.stdout.write('  ')
    #        for kk in xrange(MatG.shape[1]): sys.stdout.write('%10i'%(kk))
    #        sys.stdout.write('\n')
    #        for k in xrange(MatG.shape[0]):
    #            sys.stdout.write('%2i'%(k))
    #            for kk in xrange(MatG.shape[1]): sys.stdout.write('%10.2f'%(MatG[k,kk]))
    #            sys.stdout.write('\n')
    return sysMass

def ArcLength(gamma, UI, dUII, DU, DY, Mask):                           # arc length control, arclength measure, load displ, residuum displ, displ step incr,
    aa = MaskedP( UI,UI,Mask)                                           # parameter a
#    if aa<ZeroD: raise NameError("Arc length control: exit 1")
    if aa<0.01*ZeroD: raise NameError("Arc length control: exit 1",aa)
    bb =   MaskedP( DU,UI,Mask)   + MaskedP( dUII,UI,Mask)              # parameter b
    cc = 2*MaskedP( DU,dUII,Mask) + MaskedP( dUII,dUII,Mask)            # parameter c
    dd =   MaskedP( DU,DU,Mask)                                         # parameter d
    cp = cc+dd-gamma**2
    qq = bb**2-aa*cp
#    
#    print('XXX',qq,sqrt(abs(qq)),'_',-bb/aa+sqrt(abs(qq))/aa,-bb/aa-sqrt(abs(qq))/aa)
#    
    if qq<0: 
        return -bb/(2*aa)
    else:
        L1 = (-bb-np.sqrt(qq))/aa                                       # arc length solution 1
        L2 = (-bb+np.sqrt(qq))/aa                                       # arc length solution 2
        x1 = MaskedP( DY,(DY+dUII),Mask)
        x2 = MaskedP( DY,UI,Mask)
        g1 = x1 + L1*x2
        g2 = x1 + L2*x2
        if g1>g2: return L1
        else:     return L2
def MaskedP( V1, V2, Mask):
    xx = 0.
    for v1, v2, ma in zip( V1, V2, Mask): xx = xx +v1*v2*ma
    return xx

def EvNodeEl( NodeList, ElList, NoLabToInd):                                # evaluate element indices which each node is part of required for cuthill-mckee
    for i, el in enumerate(ElList):
        elType = el.Type
        for j in el.Inzi:
            no = NodeList[j]
            no.NodeEl += [i]
            if elType in ["CPS3","CPE3","CPS4","CPE4","CAX4"]:
                if   no.Type=="":     no.Type = "C2D"
                elif no.Type=="C2D":  pass
                else: raise NameError("EvNodeEl: exit 1",elType, no.Type)
            elif elType in ["C3D8"]:
                if   no.Type=="":     no.Type = "C3D"
                elif no.Type=="C3D":  pass
                else: raise NameError("EvNodeEl: exit 2",elType, no.Type)
            elif elType in ["T3D2I","T3D3I"]:
                if   no.Type=="":     no.Type = "T3D"
                elif no.Type=="T3D": pass
                else: raise NameError("EvNodeEl: exit 3",elType, no.Type)
            elif elType in ["T2D2I","T2D3I"]:
                if   no.Type=="":     no.Type = "T2DE"
                elif no.Type=="T2DE": pass
                else: raise NameError("EvNodeEl: exit 4",elType, no.Type)
            elif elType in ["TAX2I", "TAX3E"]:                              # ?????????????????
                if no.Type == "":     no.Type = "TAXE"
                elif no.Type=="TAXE": pass
                else: raise NameError("EvNodeEl: exit 5", elType, no.Type)
    for No in NodeList:
        No.NodeEl = list(set(No.NodeEl))                                    # to remove double entries

def EvNodeC2D( NodeList):                                                   # build coordinate list for nodes belonging to continuum elements - mainly used for search tree
    NodeCoorList2D = []
    NodeCoorList2DtoNodeList = [] 
    for i in range(len(NodeList)):
        j = NodeList[i]
        if j.Type in ["C2D","C3D"]: 
            NodeCoorList2D += [[ j.XCo,  j.YCo, j.ZCo]]
            NodeCoorList2DtoNodeList += [i]
    return NodeCoorList2D, np.array(NodeCoorList2DtoNodeList)

def FinishEquilibIteration( MatList, ElList, NodeList,NoIndToCMInd, ff, NLGeom, TimeTargetActiveFlag):
    Flag = False
    for Elem in ElList:                                                     # loop over all elements
        if len(np.shape(Elem.Data)) != 0:                                   # Data is array
            Elem.DataP[:,:] = Elem.Data[:,:]
        else:                                                               # data is dictionary
            for k in Elem.Data.keys():
                Elem.DataP[k][:,:] = Elem.Data[k][:,:]
        if len(Elem.StateVar) > 0:
            Mat = MatList[Elem.MatN] 
            if Mat.Update:
                Flag2 = Mat.UpdateStateVar(Elem, ff)
                if not Flag and Flag2: Flag = True                          # Flag triggers UpdateStat2Var
                # if Flag2: Flag = True # ? 
            else:
                Elem.StateVar[:,:] = Elem.StateVarN[:,:]
    for Elem in ElList:                                                     # loop over all elements ???
        if len(Elem.StateVar) > 0:                                          # ???
            Mat = MatList[Elem.MatN] 
            if Mat.Updat2:                                                  # currently seems to be for ElastLT only
                Mat.UpdateStat2Var(Elem, ff, Flag, TimeTargetActiveFlag)
    for Elem in ElList:
        if Elem.Type in["CPE4S","CPE3S","CPS4S","CPS3S", "C3D8S","CPE6S","CPS6S"]:
            Elem.Update( MatList[Elem.MatN], NodeList, NoIndToCMInd, ff, 0) 
    return 0

def CuthillMckee_(NodeList, ElList, NoLabToInd ):
    mapNoInd = np.zeros((len(NodeList)), dtype=int)                     # to map from original to new index in NodeList
    # for each node create information regarding connected nodes
#    for el in ElList:
#        el.Inzi_ = [ NoLabToInd[i] for i in el.InzList ]
    for i, node in enumerate(NodeList):
        ConNodes = set()
        for j in node.NodeEl:
            elem = ElList[j]
            ConNodes = ConNodes.union(set(elem.Inzi_))
        node.ListOfConnectedNodes_ = list(ConNodes.difference(set([i])))
        node.ConnectionDegree_ = len(node.ListOfConnectedNodes_)        # grade of connection
        node.reOrdered = False                                          # flag whether this got a reordered index
        node.Index = i
    # search nodes for the smallest grade of connection
    NodeWithLowestConDeg_ = 100                                         # initial large value
    for node in NodeList:
        if node.ConnectionDegree_ <= NodeWithLowestConDeg_ and node.ConnectionDegree_>0:
            NodeWithLowestConDeg_ = node.ConnectionDegree_
    # collect all node indices which have minimum grade of connection
    NodeListWithLowestConDeg_ = []
    for i, node in enumerate(NodeList):
        if node.ConnectionDegree_ <= NodeWithLowestConDeg_ and node.ConnectionDegree_>0: # consider less equal
            NodeListWithLowestConDeg_.append(i)                         # collects all with minimum grade of connection 
    # CM level 0 node
    StartNodeIndex = NodeListWithLowestConDeg_[0]                       # this choice is NOT mandatory
    startNode = NodeList[StartNodeIndex]
    startNode.CMIndex =  0
    startNode.reOrdered = True
    usedNodes = [StartNodeIndex]                                        # initialize list of nodes which have already been processed for reordering
    # prepare for reordering
    currentIndex = 1
    levelNodes = [StartNodeIndex]                                       # initialize list of nodes of current level
    # recursive scheme to reorder, hopefully comes to an end
    def Relabel(currentIndex, levelNodes, usedNodes):
        curlevNodes = set()
        for i in levelNodes:                                    
            curlevNodes = curlevNodes.union(set(NodeList[i].ListOfConnectedNodes_)) # indices of nodes of current CM level
        curlevNodes = curlevNodes.difference(set(usedNodes))            # disregard nodes which already have been reordered
        curlevNodes = list(curlevNodes)                                 # transform set into list
        curlevNodes =  sorted(curlevNodes, key=lambda t: NodeList[t].ConnectionDegree_) # sort
        # reorder according to sequence in curlevNodes
        for u in curlevNodes:
            currentNode = NodeList[u]
            if not currentNode.reOrdered:
                currentNode.CMIndex = currentIndex
                currentNode.reOrdered = True
                usedNodes += [ u ]
                currentIndex += 1
        return curlevNodes, currentIndex
    while len(levelNodes)>0:
        levelNodes, currentIndex = Relabel(currentIndex, levelNodes, usedNodes)
    # assign reordered index
    for i, node in enumerate(NodeList):
        if not node.reOrdered: node.CMIndex = len(NodeList)             # e.g. node is read in but not used
        try:    mapNoInd[i] = node.CMIndex
        except: raise NameError("CaeFemMaint::CuthillMcKee: node not mapped: "+str(node.Label))
        del node.ListOfConnectedNodes_
        del node.ConnectionDegree_
#    for el in ElList: del el.Inzi_
    return mapNoInd

def CuthillMckee(NodeList, ElList):
    NoIndToCMInd = np.zeros((len(NodeList)), dtype=int)                        # to map from original to new index in NodeList
#    CMIndToNoInd = full((len(NodeList)), -1, dtype=int)                    # inverse Cuthill-McKee
    # for each node create information regarding connected nodes
    for i in range(len(NodeList)):
        node = NodeList[i]
        ConNodes = set()
        #
        for j in node.NodeEl:
            elem = ElList[j]
            ConNodes = ConNodes.union(set(elem.Inzi))
#        for elem in ElList:
#            if i in elem.Inzi:
#                InzSet = set(elem.Inzi)
#                ConNodes=ConNodes.union(InzSet)
#        print 'XXX', node.Label, ConNodes_, ConNodes
        #
        node.ConnectedNodes = list(ConNodes.difference(set([i])))
        node.ConnectionDegree = len(node.ConnectedNodes)# grade of connection
        node.reOrdered = False                                  # flag whether this got a reordered index
        node.Index = i
    # search nodes for the smallest grade of connection
    LowestConDeg = 100                                 # initial large value
    for node in NodeList:
        if node.ConnectionDegree > 0 and node.ConnectionDegree <= LowestConDeg:
            LowestConDeg = node.ConnectionDegree
    # collect all node indices which have minimum grade of connection
    NodeListWithLowestConDeg_ = []
    for i in range(len(NodeList)):
        node = NodeList[i]   
        if node.ConnectionDegree > 0 and node.ConnectionDegree <= LowestConDeg: # consider less equal
            NodeListWithLowestConDeg_.append(i)                 # collects all with minimum grade of connection 
    # CM level 0 node
#    StartNodeIndex = NodeListWithLowestConDeg_[0]               # this choice is NOT mandatory
    StartNodeIndex = NodeListWithLowestConDeg_[0]               # this choice is NOT mandatory
#    print('CM',NodeListWithLowestConDeg_,StartNodeIndex)
    startNode = NodeList[StartNodeIndex]
    startNode.CMIndex =  0
    startNode.reOrdered = True
    usedNodes = [StartNodeIndex]                                # initialize list of nodes which have already been processed for reordering
    # prepare for reordering
    currentIndex = 1
    levelNodes = [StartNodeIndex]                               # initialize list of nodes of current level
    # recursive scheme to reorder, hopefully comes to an end
    def Relabel(currentIndex, levelNodes, usedNodes):
        curlevNodes = set()
        for i in levelNodes:                                    
            node = NodeList[ i]
            curlevNodes = curlevNodes.union(set(node.ConnectedNodes)) # indices of nodes connected of previous level / nodes in levelnodes
        curlevNodes = curlevNodes.difference(set(usedNodes))    # disregard nodes which already have been reordered
        curlevNodes = list(curlevNodes)                         # transform set into list
        curlevNodes =  sorted(curlevNodes, key=lambda t: NodeList[t].ConnectionDegree) # sort current level with respect to increasing connection degree
        # reorder according to sequence in curlevNodes
        for u in curlevNodes:
            currentNode = NodeList[u]
            if not currentNode.reOrdered:
                currentNode.CMIndex = currentIndex
                currentNode.reOrdered = True
                usedNodes += [ u ]
                currentIndex += 1
        return curlevNodes, currentIndex
    # recursively takes "actual" level as input and generates new level
#    levelCounter, maxLevelLength, levelAcc = 0, 0, 0                                                            # auxiliary
    while len(levelNodes)>0:
#        levelCounter +=1                                                                    # auxiliary
        levelNodes, currentIndex = Relabel(currentIndex, levelNodes, usedNodes)
#        if len(levelNodes)>maxLevelLength: maxLevelLength = len(levelNodes)                             # auxiliary
#        levelAcc +=   len(levelNodes)
#        print('CM',levelCounter,currentIndex,'_',len(levelNodes),levelAcc,'_',maxLevelLength,) 
    # assign reordered index
    numberOfUnusedNodes = 0
    for i in range(len(NodeList)):
        node = NodeList[i]
        if not node.reOrdered: numberOfUnusedNodes += 1
    numberOfUsedNodes = len(NodeList)-numberOfUnusedNodes
    for i in range(len(NodeList)):
        node = NodeList[i]
        if not node.reOrdered: node.CMIndex = len(NodeList)     # e.g. node is read in but not used
        # reverse CM
#        if node.reOrdered: node.CMIndex = numberOfUsedNodes-1-node.CMIndex
        #
        try:    NoIndToCMInd[i] = node.CMIndex
        except: raise NameError("CaeFemMain::CuthillMcKee: node not mapped: "+str(node.Label))
#        print(i,NoIndToCMInd[i],node.Label,node.NodeEl,'__',node.ConnectionDegree,node.ConnectedNodes,node.CMIndex,'_',numberOfUsedNodes)
        del node.ConnectedNodes
        del node.ConnectionDegree
        del node.reOrdered
    return NoIndToCMInd #, CMIndToNoInd

def Sloan(NodeList, ElList):
    # recursive scheme to generate level structure, hopefully comes to an end
    def CreateLevelStructure( CurrentLevel, UsedNodes):
        CurLevNodes = set()
        for i in CurrentLevel:                                    
            node = NodeList[ i]
            CurLevNodes = CurLevNodes.union(set(node.ConnectedNodes)) # indices of nodes connected of previous level / nodes in levelnodes
        CurLevNodes = CurLevNodes.difference(set(UsedNodes))    # disregard nodes which already have been reordered
        CurLevNodes = list(CurLevNodes)                         # transform set into list
        CurLevNodes =  sorted(CurLevNodes, key=lambda t: NodeList[t].ConnectionDegree) # sort current level with respect to increasing connection degree
        for i in CurLevNodes: 
            UsedNodes.append(i)
        return CurLevNodes
    def CreateLevelStructureWrapper( RootNodeInd ):
        LevelStructure, CurrentLevel, UsedNodes, LevelDepth, LevelWidth = [], [RootNodeInd], [RootNodeInd], 1, 0
        while len(CurrentLevel)>0:
            CurrentLevel = CreateLevelStructure( CurrentLevel, UsedNodes)
            CurrentLevelL = len(CurrentLevel)
            if CurrentLevelL>0: 
                LevelDepth += 1
                LevelStructure += [CurrentLevel]
                if CurrentLevelL>LevelWidth: LevelWidth= CurrentLevelL
        return LevelDepth, LevelWidth, LevelStructure
    def ShrinkLevel( LastLevel ):
        LastLevelS = [LastLevel[0]]
        for i in range(len(LastLevel)-1):
            i0, i1 = LastLevel[i], LastLevel[i+1]
            if NodeList[i1].ConnectionDegree>NodeList[i0].ConnectionDegree: LastLevelS += [i1]
        return LastLevelS
    def DiaMtr( s, rlist, LevelDepth ):     # r: index root node; rlist: shrinked list of nodes of last level of levellist of r, LevelDepth: depth of levelist 
        LevelDepthMax, LevelWidthMin, e = LevelDepth, 1000000, -1     
        for i in rlist:
            LevelDepth_, LevelWidth_, LevelStructure_ = CreateLevelStructureWrapper( i )
            if LevelDepth_ > LevelDepth:
#                raise NameError("check this case")
                return True,  i, LevelDepth_, None, LevelStructure_
            if LevelWidth_ < LevelWidthMin: 
                e, eLevelStructure, LevelWidthMin = i, LevelStructure_, LevelWidth_
        return         False, s, None,       e,     eLevelStructure
    
    NoIndToCMInd = np.zeros((len(NodeList)), dtype=int)                        # to map from original to new index in NodeList
    # connection lists for each node
    for i in range(len(NodeList)):
        node = NodeList[i]
        ConNodes = set()
        for j in node.NodeEl:
            elem = ElList[j]
            ConNodes = ConNodes.union(set(elem.Inzi))
        node.ConnectedNodes = list(ConNodes.difference(set([i])))           # list of connected nodes
        node.ConnectionDegree = len(node.ConnectedNodes)                    # grade of connection
        node.reOrdered = False                                  # flag whether this got a reordered index
        node.Index = i
#        if node.ConnectionDegree < 2: print("Sloan:: node with connection degree < 2, be cautious",node.Label,node.ConnectionDegree)
    # nodes with smallest grade of connection
    LowestConDeg = 100                                                      # initial large value
    for node in NodeList:
        if node.ConnectionDegree>0 and node.ConnectionDegree<=LowestConDeg: 
            LowestConDeg = node.ConnectionDegree
    NodesWithLowestConDeg = []
    for i in range(len(NodeList)):
        node = NodeList[i]   
        if node.ConnectionDegree>0 and node.ConnectionDegree<=LowestConDeg: # consider less equal
            NodesWithLowestConDeg.append(i)                                 # collects all with minimum grade of connection
    # rooted level structure
    LevelDepth, LevelWidth, LevelStructure = CreateLevelStructureWrapper( NodesWithLowestConDeg[0] )
    # look for pseudo peripheral nodes
    Flag, s, LevelDepth_ = True, NodesWithLowestConDeg[0], LevelDepth      # initial values for DiaMtr
    while Flag:
        Flag, s, LevelDepth_, e, LevelStructure = DiaMtr( s, ShrinkLevel(LevelStructure[-1]), LevelDepth_ )  # s: root node (index in NodeList); e: end node of peripheral pair;
                                                                                                             # LevelStructure finally: level structure with e as root
                                                                                                             # take care: LevelStructure is formerly used for first root node
    # compute distance / priority ..... of nodes to e
    W1, W2 = 1, 2
    W1, W2 = 2, 1
    for i in range(len(LevelStructure)):
        for j in LevelStructure[i]:
            node = NodeList[j]
            node.LevelDist = i+1                                            # Level Distance  from each node to end node
            node.Priority  = W1*node.LevelDist - W2*node.ConnectionDegree
            node.Status    = -3                                             # >=0: post active; -1: active; -2: pre active; -3: inactive 
    # what about e
    node = NodeList[e] 
    node.LevelDist = 0
    node.Priority  = W1*node.LevelDist - W2*node.ConnectionDegree
    node.Status    = -3
    # initialize with start node s
    LSTNUM, QQ = 0, [s]
    NodeList[s].Status = -2                                                 # pre-active
    while len(QQ)>0:
        # scan queue for node with max priority
        ADDRES = 0
        MAXPRT = NodeList[QQ[0]].Priority
        for i in range(1,len(QQ)):
            PRTY = NodeList[QQ[i]].Priority
            if PRTY>MAXPRT:
                ADDRES = i
                MAXPRT = PRTY
        # NEXT is the node to be numbered next
        NEXT = QQ[ADDRES]
        # delete node NEXT from queue
        del QQ[ADDRES]
        if NodeList[NEXT].Status == -2:
            neighbours = NodeList[NEXT].ConnectedNodes
            for i in neighbours:
                # decrease current degree of neighbor by -1 (???)
                NodeList[i].Priority += W2                                  # ???
                # add neighbor to queue if it is inactive, assign it a pre-active status
                if NodeList[i].Status == -3:
                    QQ.append(i)
                    NodeList[i].Status = -2                                 # status pre-active 
        # store new node number, status of node NEXT is now post-active
        NoIndToCMInd[NEXT]  = LSTNUM
        NodeList[NEXT].CMIndex = LSTNUM
        NodeList[NEXT].Status  = LSTNUM 
        NodeList[NEXT].reOrdered = True
        LSTNUM += 1
        # search for pre-active neighbors of node NEXT
        for i in neighbours:                                                # 80
            if NodeList[i].Status == -2:
                # decrease current degree of preactive neighbor by -1, assign neighbor an active status
                NodeList[i].Priority += W2                                  # ???
                NodeList[i].Status = -1                                     # status active
                # loop over the nodes adjacent to preactive neighbor
                neighboursPre = NodeList[i].ConnectedNodes
                for j in neighboursPre:                                     # 60
                    # decrease current degree of adjacent node by -1
                    NodeList[j].Priority += W2
                    if NodeList[j].Status == -3:
                        NodeList[j].Status = -2
                        QQ.append(j)
        # 80
        QQ = list(set(QQ))
        
#        print('BBB',NEXT,N,'__',QQ)
#        postactive, active, preactive, inactive = [], [], [], []
#        for i in range(len(NodeList)):
#            node = NodeList[i]
#            if   node.Status >= 0:  postactive += [i]
#            elif node.Status == -1: active     += [i]
#            elif node.Status == -2: preactive  += [i]
#            elif node.Status == -3: inactive   += [i]
#            else: raise NameError("error")
#        print('CCC\n','post',postactive,'\n','acti',active,'\n','prea',preactive,'\n','inac',inactive)
#        
    for i in range(len(NodeList)):
        node = NodeList[i]
        if not node.reOrdered: node.CMIndex = len(NodeList)                 # e.g. node is read in but not used
    # housekeeping
    for i in range(len(NodeList)):
        node = NodeList[i]
        del node.ConnectedNodes
        del node.ConnectionDegree
        del node.reOrdered
        try:
            del node.Status
            del node.LevelDist
            del node.Priority
        except:
            pass
    return NoIndToCMInd

def WithoutBandwidthOpti( NodeList ):
    NoIndToCMInd = np.zeros((len(NodeList)), dtype=int)                        # to map from original to new index in NodeList
    j = 0
    jj= len(NodeList)
    for i in range(len(NodeList)):
        node = NodeList[i]
        node.Index   = i
        if len(node.NodeEl)>0:
            NoIndToCMInd[i] = j
            node.CMIndex    = j
            j += 1
        else:
            NoIndToCMInd[i] = jj
            node.CMIndex    = jj
    return NoIndToCMInd

def IsSymSys( MatList ):
    for i in list(MatList.values()):                                    # if there is at least one unsymmetric material the whole system is unsymmetric
        if (not i.Symmetric) and (i.Used):
            return False
    return True 

def PickleSysMat(NodeList, ElList, Step, N, MatK, MatM):                # where is this used, also defined in CpnLinAlg
    import pickle
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

#@profile
def Eigen( ne, lim, N, NodeList, ElList, Step, Kmat, Mmat, NoLabToNoInd,NoIndToCMInd, PloF, ff ): # NoLabToInd, mapNoInd, PloF, ff ):
    from scipy.sparse.linalg.dsolve import linsolve
    from scipy.sparse.linalg import aslinearoperator
    from numpy import zeros, dot, transpose
    from numpy.linalg import norm
    from scipy.linalg import eigh

    Echo(f"Eigenvalue analysis: Matrices have been computed", ff)
    BCIndex, QList = zeros((N),dtype=int), []                           # QList to sort ratios K/M
    for bl in Step.BoundList:                                           # loop over all boundary conditions of step
        Found = False
#        nI = NoLabToNoInd[bl.NodeLabel]
        nI = FindIndexByLabel(NodeList, bl.NodeLabel)  # node index of bc

        for j, jj in enumerate(NodeList[nI].DofT):                      # loop over all dofs of node
            if jj==bl.Dof:                                              # dof type equals prescribed dof type 
                Found = True
                break 
        if not Found:
            raise NameError ("ConFemBasics::Eigen:missing correspondence for dof type", bl.NodeLabel,bl.Dof,NodeList[nI].DofT)
        k = NodeList[nI].GlobDofStart + j                               # constrained dof global index
        BCIndex[k] = 1                                                  #  mark dofs with boundary conditions to build start vectors
        for j in range(N):                                              # loop over rows and columns simul
            Kmat[j,k] = 0.                                              # modification stiffness matrix
            Kmat[k,j] = 0.
            Mmat[j,k] = 0.                                              # modification geometric stiffness matrix
            Mmat[k,j] = 0.
#        Kmat[:,k] = 0.
#        Kmat[k,:] = 0.
        Kmat[k,k] = 1.
#        Mmat[:,k] = 0.
#        Mmat[k,:] = 0.                                                  # modification stiffness matrix
        Mmat[k,k] = 1.
    Echo(f"Eigenvalue analysis: kinematic boundary conditions considered", ff)

    K_LU = linsolve.splu(Kmat.tocsc(),permc_spec=3)                     # triangularization of stiffness matrix
    np_ = min(2*ne, ne+8, N)                                            # dimension of subspace with ne number of required eigenpairs, Bathe p. 963
    KK_, MM_, YY, XX = zeros((np_,np_),dtype=float), zeros((np_,np_),dtype=float), zeros((N,np_),dtype=float),zeros((N,np_),dtype = float)
    # starting vectors YY for iteration
    for i in range(N): 
        if BCIndex[i]==0: YY[i,0] = Mmat[i,i]                           # first column of auxiliary vector YY
    QList = [ [Mmat[i,i]/Kmat[i,i],i] for i in range(N)]                # build list of "lazyness"-values 
    QList.sort( reverse=True)                                           # sort for most "lazy" dofs, largest values M/K first 
    i, j = 0, 1
    while j<np_:                                                        # 0 column of YY already occupied  
        index = QList[i][1]                                             # extract index for most lazy dofs 
        if BCIndex[index] == 0:                                         # should not be dof with boundary condition
            YY[index,j] = 1.                                            # further columns of auxiliary vector YY
            j += 1
        i +=1
    Echo(f"Eigenvalue analysis: prereqs for iterations ready", ff)
    
    # iteration loop for eigenform evaluation
    rr, tol, iiter, YY_ = [1.0], 1.e-6, 0, zeros((N,np_), dtype=float)
    while max(rr)>tol and iiter<lim:                                    
        for i in range(np_):     XX[:,i] = K_LU.solve( YY[:,i] )
        for i in range(np_): 
            for j in range(np_): KK_[i,j] = dot(XX[:,i],YY[:,j])        # Stiffness transformed to subspace
        for i in range(np_):     YY[:,i] = aslinearoperator(Mmat).matvec(XX[:,i]) # dot(MM,XX[:,i])
        for i in range(np_): 
            for j in range(np_): MM_[i,j] = dot(XX[:,i],YY[:,j])        # Mass transformed to subspace
        evals, evecs = eigh(KK_, MM_)                                   # eigenvalues, eigenvectors in subspace
        rr = [sqrt(abs(1.-evals[i]**2/dot(transpose(evecs[:,i]),evecs[:,i]))) for i in range(ne)] # to control convergence, see Bathe 1996, 11.6.4
        for i in range(np_):
            YY_[:,i] = 0.
            for j in range(np_): 
                YY_[:,i] = YY_[:,i] + YY[:,j]*evecs[j,i]                # update auxiliary vector
        YY[:] = YY_[:]
        Echo(f"iter {iiter:d}, residual {norm(rr):f}", ff)
        iiter += 1
    # write data
    # postprocessing and control
    EV = zeros((N,np_), dtype=float)
    for i in range(ne):                                   
        for j in range(np_):                               
            EV[:,i] = EV[:,i] + XX[:,j]*evecs[j,i]                      # derive eigenvectors
        xx = EV[:,i]
        re = aslinearoperator(Kmat).matvec(xx)
        yy = aslinearoperator(Kmat).matvec(xx)-evals[i]*aslinearoperator(Mmat).matvec(xx) # eigenform condition
        Echo(f"eigenform {i:d}, value {evals[i]:f}, rel. acc {norm(yy)/norm(re):e}", ff)

        if PloF: 
            from ConFemPost import PostNode3D 
            PostNode3D( ElList,NodeList,NoIndToCMInd, None,xx, evals[i], 1.0) # xx is eigenvector
    Echo(f"Eigenvalue analysis: finished", ff)



