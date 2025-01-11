# ConFemInOut -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
#import matplotlib           as mpl
#import matplotlib.pyplot    as plt
from matplotlib.pyplot import figure, plot, grid, title, text, contour, clabel, show, axis, xticks, yticks, ylim, annotate
#import matplotlib.colors    as colors
#import mpl_toolkits.mplot3d as a3
#import pylab                as pl
#from pylab import arange, griddata
import numpy                as np
#from numpy        import meshgrid, pi, arcsin, fabs, sqrt
#from scipy.linalg import norm
from os           import path, makedirs
import copy

from ConFemBasics import _BeamsAll_, _BeamsBernAll_, _BeamsTimoAll_, _Truss2DEmbedded_
from ConFemMat import *
from ConFemElem import *
from ConFemSteps import *
#
from ConFemRandom import RandomField_Routines
from numpy.linalg import lstsq

#Colors = ['tab:green','tab:brown','magenta','darkviolet','blue','green','yellow']
#Colors = ['royalblue','magenta','darkviolet','darkcyan','blue','green','yellow']
#FonSizTi='x-large'   # fontsize title x-large
#FonSizAx='x-large'     # fontsize axis
#FonSizAx='medium'     # fontsize axis
#LiWiCrv=3            # line width curve in pt
#LiWiCrvMid=2            # line width curve in pt
#LiWiCrvThin=1.0        # line width curve in pt

def ReadPlotOptionsFile(f1, SecDic):                            # read plot options
    SN, SN2, PE2DFlag, PE3DFlag  = 1.0, 0.0, True, False
    PlotTimes,  ScaleStress, Post1DFlag, PostNode, ShellL, ScaleShellL, Contour2D = [], [], True, True, False, 0., {}
    z1 = f1.readline()                                  # 1st input line
    z2 = z1.split(",")
    zLen = len(z2)
    while z1!="":
        if z2[0]=="" or z2[0].find("**")>-1:
            pass
        else:
            for i in range(zLen): z2[i] = z2[i].strip()
            if z2[0].upper()=="*PLOTTIMES":  
                for i in range(1,len(z2)):          PlotTimes += [float(z2[i])]
            elif z2[0].upper()=="*SCALESTRESS":
                n = int((len(z2)-1)/2)                                      # two data - elset, scaling factor - expected per item  
                for j in range(n):
                    ScaleStress += [[z2[2*j+1].upper(),float(z2[2*j+2])]]
            elif z2[0].upper()=="*SCALEDISPL":      SN  =float(z2[1])
            elif z2[0].upper()=="*SCALEDISPL2D":    SN2 =float(z2[1])       # for deformed 2D elements with stresses
            elif z2[0].upper()=="*NOPOSTNODE":      PostNode = False
            elif z2[0].upper()=="*NOPOSTELEM1D":    Post1DFlag = False
            elif z2[0].upper()=="*NOPOSTELEM2D":    PE2DFlag = False
            elif z2[0].upper()=="*POSTELEM3D":      PE3DFlag = True
            elif z2[0].upper()=="*POSTSHELLLAYERS": ShellL, ScaleShellL = True, float(z2[1])
            elif z2[0].upper()=="*POSTELEM2DCONTOUR": # four data - elset, label, index for float data (0 is starter for 1st float) for data position in *.elemout_.txt, length for triangle masking, indicator 'H' for histogram - expected per item
                r = 5
                n = int((len(z2)-1)/r)
                for j in range(n):
                    key = z2[r*j+1].upper()+" "+z2[r*j+2]
                    Contour2D[key] = [ int(z2[r*j+3]), z2[r*j+4] , z2[r*j+5]]
            else: raise NameError("ConFemInOut::ReadPlotOptionsFile: unknown keyword "+z2[0])
        z1 = f1.readline()                                          # read next line
        z2 = z1.split(",")
        zLen = len(z2)
    for s in ScaleStress:
        elset = s[0]
        if elset in SecDic: pass
        else:               print("ConFemInOut::ReadPlotOptionsFile: unknown elset for stress scaling "+elset)
    return SN, SN2, PE2DFlag, PE3DFlag, PlotTimes, Post1DFlag, ScaleStress, PostNode, ShellL, ScaleShellL, Contour2D

def ReadOptionsFile(f1, NodeList, NoLabToNoInd, NoIndToCMInd):              # read options
    LS = [1,1.]                                                             # default values for line search
    R, L, MaxType, MaxEquiIter, StabSys,StabTolF, SoftSys,SoftRed = [],[],[],None, False,None, False,None
    with ( f1 ) as ff:
        for z1 in ff:
            z2 = z1.split(",")
            zLen = len(z2)
            for i in range(zLen): z2[i] = z2[i].strip()
            if z2[0].upper()=="WRITENODES":                                 # write single nodal displacement, force during whole course of computation
                n = int((len(z2)-1)/3)                                      # three data - node, dof, outtype - expected per item
                for j in range(n):
                    iP = int(z2[3*j+1])                                     # node
                    iD = int(z2[3*j+2])                                     # dof
                    ty = z2[3*j+3]                                          # outtypes p-u, p-t, u-t, b-u, b-t
                    try: node = NodeList[ NoIndToCMInd[NoLabToNoInd[iP]] ] # !!! might not correctly work for unused nodes where  NoIndToCMInd[NoLabToNoInd[iP]] presumably points to zero - don't refer to them in opt-file
                    except: raise NameError ("ConFemInOut::ReadOptionsFile: no such node available in current dataset "+str(iP))
                    DofTL = list(node.DofT)                                 # DofT is set of markers for dof types for current node
                    iD_ = -1
                    for i in range(len(DofTL)):
                        if DofTL[i]==iD:
                            iD_ = i                                         # iD_ and iD differ!
                    if iD_==-1: raise NameError ("ConFemInOut::ReadOptionsFile: invalid dof type indicator ")
                    ii = node.GlobDofStart+iD_                              # global index
                    L += [[iP,iD_,ii,ty]]
            elif z2[0].upper()=="LINESEARCH":                               # parameter for line search iteration method
                if int(z2[1])>1:
                    LS[0] = int(z2[1])
                    LS[1] = float(z2[2])
            elif z2[0].upper()=="REINFORCEMENTDESIGN":                      # parameter reinforcement design - ConPlad
                if int(z2[1])>1:
                    R += [float(z2[1])]                                     # design value reinforcement strength
                    R += [float(z2[2])]                                     # design value concrete compressive strength
                    R += [float(z2[3])]                                     # minimum reinforcement ratio x (%), minimum reinforcement ratio labs   ( formerly geometric height )
                    R += [float(z2[4])]                                     # minimum reinforcement ratio y (%)   ( formerly structural height)
                    R += [float(z2[5])]                                     # geometric height
                    R += [float(z2[6])]                                     # geometric structural height
                    R += [float(z2[7])]                                     # specific weight reinforcing steel
            elif z2[0].upper()=="MAXVALUES":
                MaxType = z2[1:]
            elif z2[0].upper()=="MAX_EQUI_ITERATIONS":                      # overrides EquiFailedMax in ConFem
                MaxEquiIter = int(z2[1])
            elif z2[0].upper()=="STABILITY_TOLERANCE_FACTOR":               #
                StabSys  = True
                try:    StabTolF = float(z2[1])
                except: raise NameError("ConFemInOut::ReadOptionsFile: a value must be defined for "+z2[0])
            elif z2[0].upper()=="TERMINATE_SOFTENING_SYSTEM":               # "
                SoftSys = True
                try:    SoftRed = float(z2[1])
                except: raise NameError("ConFemInOut::ReadOptionsFile: a value must be defined for "+z2[0])
            elif z2[0].find("**") > -1:
                pass
            else:
                raise NameError("ConFemInOut::ReadOptionsFile: unknown keyword "+z2[0])

    return L,R, MaxType,MaxEquiIter, StabSys,StabTolF, SoftSys,SoftRed

def WriteNodes( f5, WrNodes, Time, VecU,VecB,VecP, step,counter, queues,ndeq,maxWriteNodes):   # write single nodal displacement, force during whole course of computation into file
    for j in range(len(WrNodes)):
        ii = WrNodes[j][2]
        if WrNodes[j][3]  =='u-t': f5.write('%8.4f%12.6f' %(Time,VecU[ii]))
        elif WrNodes[j][3]=='t-u': f5.write('%12.6f%8.4f' %(fabs(VecU[ii]),Time,))
        elif WrNodes[j][3]=='b-u':
            f5.write('%14.8f%14.6e'%(VecU[ii],VecB[ii]))
            if j==0:
                queues[0].append(VecU[ii])
                queues[1].append(VecB[ii])
                queues[2].append(counter)
                qVecB = np.mean(list(queues[1]))
                if len(queues[0])==ndeq and abs(qVecB) > maxWriteNodes[1]:
                    maxWriteNodes[0] = abs(np.mean(list(queues[0])))
                    maxWriteNodes[1] = abs(qVecB)
                    maxWriteNodes[2] = counter
        elif WrNodes[j][3]=='b-t': f5.write('%8.4f%12.6f'%(Time,VecB[ii]))
        elif WrNodes[j][3]=='p-u': f5.write('%14.8f%14.8f'%(VecU[ii],VecP[ii]))
        elif WrNodes[j][3]=='p-t': f5.write('%8.4f%12.6f'%(Time,VecP[ii]))
        else: raise NameError ("no type defined for writeNodes in *.opt.txt")
    f5.write('%3d%8d\n'%(step,counter))
    return 0

def MatInput( z3, zLen, MatList, matName, IType2, mass, LineCounter, L1, RotCrack, PrinStrain, ShearRetFac, S2ndCrack, S2ndCrackA, Visco, PhaseField, SDA, TraPre, LatBonTeFac, RedInt, sval, ff ):
    #
    if z3[0].upper()=="*DENSITY":    
        IType2="density"
    elif   z3[0].upper()=="*ELASTIC":    
        for i in z3:
            if i.upper().find("PHASE_FIELD")>-1:            PhaseField = True
        for i in z3:
            if i.upper().find("VISCO")>-1:                  Visco = True
        IType2="elastic"
    elif z3[0].upper()=="*ELASTICR1D":
        IType2 = "elastic1dr"
    elif z3[0].upper()=="*ELASTICC1D":
        IType2 = "elasticc1d"
    elif z3[0].upper()=="*ELASTICLT":  
        IType2="elasticlt"
        for i in z3:
            if i.upper().find("REDUCED_INTEGRATION")>-1:    RedInt = True
            if i.upper().find("ROTATING_CRACK")>-1:         RotCrack = True
            if i.upper().find("PRINCIPAL_STRESS")>-1:       PrinStrain = False
            if i.upper().find("2NDCRACK")>-1:               
                S2ndCrack = True
                S2ndCrackA = i.split("=")[1]                                                # minimum deviation of 2nd crack from 1st crack in degree
            if i.upper().find("SHEAR_RETENTION_FACTOR")>-1: ShearRetFac=i.split("=")[1]     # shear retention factor
            if i.upper().find("SDA")>-1:                    SDA = True
        for i in z3:
            if SDA and i.upper().find("FULLINT")>-1:        RedInt = False                  # e.g. for one element systems
            else:                                           RedInt = True                   # default approach with reduced integration for SDA elements
    elif z3[0].upper()=="*ELASTIC_SDA":
#        IType2="elastic_sda"  be careful to activate original ELASTIC_SDA, input format changed for last items
        SDA = True
        for i in z3:
            if i.upper().find("REDUCED_INTEGRATION")>-1:    RedInt = True
            if i.upper().find("ROTATING_CRACK")>-1:         RotCrack = True
            if i.upper().find("PRINCIPAL_STRESS")>-1:       PrinStrain = False
            if i.upper().find("2NDCRACK")>-1:               
                S2ndCrack = True
                S2ndCrackA = i.split("=")[1]                                                # minimum deviation of 2nd crack from 1st crack in degree
            if i.upper().find("SHEAR_RETENTION_FACTOR")>-1: ShearRetFac=float(i.split("=")[1])     # shear retention factor
    elif z3[0].upper()=="*ELASTICORTHO":        
        IType2="elasticOR"
    elif z3[0].upper()=="*MISES":               
        IType2="mises"
    elif z3[0].upper()=="*ISODAMAGE":  
        IType2="isodam"
        for i in z3:
            if i.upper().find("REDUCED_INTEGRATION")>-1:    RedInt = True
            if i.upper().find("ROTATING_CRACK")>-1:         RotCrack = True
            if i.upper().find("PRINCIPAL_STRESS")>-1:       PrinStrain = False
            if i.upper().find("2NDCRACK")>-1:               
                S2ndCrack = True
                S2ndCrackA = i.split("=")[1]                    # minimum deviation of 2nd crack from 1st crack in degree
            if i.upper().find("SHEAR_RETENTION_FACTOR")>-1: ShearRetFac=float(i.split("=")[1])     # shear retention factor
    elif z3[0].upper()=="*BOND":  
        IType2="bond"
        for i in z3:                    # transverse pressure option - Ahmad
            if i.upper().find("TRANSVERSE_PRESSURE")>-1:    TraPre = True
            else :                                          TraPre = False     
            if i.upper().find("LAT_BOND_TEN_FACTOR")>-1:    LatBonTeFac=i.split("=")[1]     # factor for lateral bond tension stiffness
            else:                                           LatBonTeFac=1.0
#    elif z3[0].upper()=="*MICROPLANE_DAMAGE1":
#        IType2="mp_dam1"
#        for i in z3:
#            if i.upper().find("REDUCED_INTEGRATION")>-1:    RedInt = True
#            if i.upper().find("ROTATING_CRACK")>-1:         RotCrack = True
#            if i.upper().find("PRINCIPAL_STRESS")>-1:       PrinStrain = False
#            if i.upper().find("2NDCRACK")>-1:               
#                S2ndCrack = True
#                S2ndCrackA = i.split("=")[1]                    # minimum deviation of 2nd crack from 1st crack in degree
#            if i.upper().find("SHEAR_RETENTION_FACTOR")>-1: ShearRetFac=float(i.split("=")[1])     # shear retention factor
    elif z3[0].upper()=="*MICRODAMAGE":         IType2="mp_dam_"
    elif z3[0].upper()=="*NLSLAB":              IType2="nlslab"
    elif z3[0].upper()=="*RESHEET":             IType2="misesre"                                     # reinforcement membrane
#    elif z3[0].upper()=="*MICROPLANE_DAMAGE1":  IType2="micropl"                             # microplane
    elif z3[0].upper()=="*LUBLINER":            IType2="lubliner"
    elif z3[0].upper()=="*SPRING":              IType2="spring"                                        # nonlinear spring - currently 1D
    elif z3[0].upper()=="*RCSHELL":             IType2="rcshell"
    elif z3[0].upper()=="*RCBEAM":                                      # sub-items MATERIAL
        IType2="rcBeam"
        LineCounter = 0                                         # more than one material data line
    elif z3[0].upper()=="*POLYBEAM":
        IType2="polyBeam"
        LineCounter = 0                                         # more than one material line
    elif z3[0].upper() in ["*CONCRETE DAMAGED PLASTICITY","*CONCRETE COMPRESSION HARDENING","*CONCRETE TENSION STIFFENING"]:
        IType2 = "AbaConc"
        sval = z3[0]
    elif z3[0].upper() == "*MISESBEAM":         IType2="misesbeam"
    #
    elif IType2=="density":
        mass = float(z3[0])
    elif IType2=="elastic":
        if PhaseField and Visco: raise NameError ("ConFemInOut::DataInput: *ELASTIC options PHASE_FIELD and VISCO not allowed simultaneously",matName)           
        if PhaseField:              
            if zLen<6: raise NameError ("ConFemInOut::DataInput: not enough data for *ELASTIC with PHASEFIELD option")
            MatList[matName] = Elastic( [float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), 0., 0., mass])
        elif Visco:
            if zLen<5: raise NameError ("ConFemInOut::DataInput: not enough data for *ELASTIC with VISCO option")
#           MatList[matName] = ElasticLT([float(z3[0]),float(z3[1]),1.,1.,1.,float(z3[3]),1./float(z3[4]),0.,0.,0., mass] )
            MatList[matName] = Elastic( [float(z3[0]), float(z3[1]), float(z3[2]), 0., 0., 0., float(z3[3]), float(z3[4]), mass])
        else: 
            if zLen==3: MatList[matName] = Elastic( [float(z3[0]), float(z3[1]), float(z3[2]), 0., 0., 0., 0., 0., mass])
            else:       MatList[matName] = Elastic( [float(z3[0]), float(z3[1]), 0.,           0., 0., 0., 0., 0., mass])
    elif IType2=="elastic1dr":
        MatList[matName] =              ElasticR1D([float(z3[0]), float(z3[1]), float(z3[2]), mass])
    elif IType2 =="elasticc1d":
        MatList[matName] =              ElasticC1D( [float(z3[0]), float(z3[1]), float(z3[2]), mass])
    elif IType2=="elasticlt":
        if SDA:
            if zLen<8: raise NameError ("ConFemInOut::DataInput: not enough data for *ELASTIC option SDA",matName)
#           MatList[matName] = ElasticSDA( [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]), mass],RedInt,RotCrack,PrinStrain, ShearRetFac, S2ndCrack,float(S2ndCrackA) )
            if float(z3[5])>0.: S2ndCrack = True
            else:               S2ndCrack = False
            PrinStrain = False
            MatList[matName] = ElasticSDA( [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[6]),float(z3[7]), mass],RedInt,RotCrack,PrinStrain, float(z3[4]),S2ndCrack,float(z3[5]) )
        else:
            if zLen<7: raise NameError ("ConFemInOut::DataInput: not enough data for *ELASTICLT",matName)
#           MatList[matName] = ElasticLT([float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]),int(z3[7]),float(z3[8]),float(z3[9]), mass] )
            MatList[matName] = ElasticLT([float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]), mass] )
    elif IType2=="elasticOR":                # data MATERIAL elastic
        MatList[matName] = ElasticOrtho( [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]),float(z3[7]),float(z3[8]), mass])
    elif IType2=="elastic_sda" and SDA:                     # not addressed anymore above - replaced by ELASTICLT, SDA
        MatList[matName] = ElasticSDA( [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]), mass],RedInt,RotCrack,PrinStrain, ShearRetFac, S2ndCrack,float(S2ndCrackA) )
    elif IType2=="mises":
        if zLen<7: raise NameError ("ConFemInOut::DataInput: not enough data for *MISES",matName,z3)
        MatList[matName] =     Mises(         [float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6]), mass], [0, 0] )
        MatList[matName+'1'] = MisesUniaxial( [float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6]), mass], [0, 0] ) # provision for uniaxial behavior
        MatList[matName+'2'] = MisesBeam2D_([float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6]), mass], [0, 0]) # provision for uniaxial behavior with beams
    elif IType2=="misesre":
        if   len(z3)==6: val = 0.
        elif len(z3)==7: val = float(z3[6])
        else: raise NameError ("ConFemInOut::DataInput: number of mises reinforcement membrane data")
        MatList[matName] = WraMisesReMem( [float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), 0., val, mass], [0, 0] )
    elif IType2=="isodam":
        if zLen<12: raise NameError ("ConFemInOut::DataInput: not enough data for *ISODAMAGE",zLen)
        iList = []
        iList=[float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),int(z3[4]),float(z3[5]),float(z3[6]),float(z3[7]),int(z3[8]),float(z3[9]),float(z3[10]),float(z3[11]),0.]
        iList.append(mass)
        MatList[matName] = IsoDamage( iList, RotCrack,PrinStrain,ShearRetFac, S2ndCrack,float(S2ndCrackA))
    elif IType2=="mp_dam_":
        if zLen<11: raise NameError ("ConFemInOut::DataInput: not enough data for *MICRODAMAGE")
        iList = []
        iList=[float(z3[0]),float(z3[1]),int(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]),int(z3[7]),float(z3[8]),float(z3[9]),float(z3[10]), 0.]
        iList.append(mass)
        MatList[matName] = MicroPlaneDam( iList, None,None,None, None,None, ff)
    elif IType2=="lubliner":
        if   zLen==13: pass
        else: raise NameError ("ConFemInOut::ReadInputFile: number of Lubliner data")
        MatList[matName] = Lubliner( [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]),float(z3[7]),float(z3[8]),float(z3[9]),float(z3[10]),float(z3[11]),float(z3[12]), mass])
    elif IType2=="rcBeam":                   # data MATERIAL rcBeam -- also notice tcbeam 
        if LineCounter==0:                  # material data concrete
            if    zLen<10: raise NameError ("ConFemInOut::DataInput: too less concrete material data",z3)
            L1 = [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),int(z3[5]),float(z3[6]),float(z3[7]),float(z3[8]),mass,float(z3[9])] # material data concrete
            LineCounter = LineCounter+1
        else:                               # material data reinforcement of mises type
            if    zLen<7: raise NameError ("ConFemInOut::DataInput: too less reinforcement material data")
            MatList[matName] = RCBeam( [L1, [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]), float(z3[6]), mass]])  # material data reinforcement
    elif IType2=="nlslab":
        IList =                [float(z3[0]), float(z3[1]),  float(z3[2]),  float(z3[3]),  float(z3[4]),  float(z3[5]),  float(z3[6]),  float(z3[7]), float(z3[8])]
        if len(z3)>9: IList += [float(z3[9]), float(z3[10]), float(z3[11]), float(z3[12]), float(z3[13]), float(z3[14]), float(z3[15]), float(z3[16])]
        MatList[matName] = NLSlab( IList )
    elif IType2=="rcshell":                  # data MATERIAL wrapper for reinforced shell
        if   zLen==7: val, matn, matT = 0., z3[6].strip(), None
        elif zLen==8: val, matn, matT = float(z3[6]), z3[7].strip(), None
        elif zLen==9: val, matn, matT = float(z3[6]), z3[7].strip(), z3[8]   # to allow for a selection of reinforcement type: TR -> TextileR, otherwise -> MisesUniaxial
        if matn in MatList: 
            L1 = MatList[matn].PropMat   
        else: raise NameError ("ConFemInOut::DataInput: unknown bulk for RC shell wrapper - maybe wrong sequence of materials in input",val,matn)
        L2 = [float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), val, mass] # see mises:__init__ for meaning of parameters
        MatList[matName] = WraRCShell([L1, L2], MatList[matn], matT, ff)
    elif IType2=="spring":
        MatList[matName] = Spring( [float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6])] )
    elif IType2=="bond":
        if LineCounter == 0:
            L1 = [float(z3[0]), float(z3[1]), float(LatBonTeFac)]
            LineCounter += 1
        elif LineCounter == 1:
            L2 = [float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6])]
            MatList[matName] = Bond( L1, L2, TraPre )                       # adding transverse pressure option to bond material - Ahmad
    elif IType2=="polyBeam":
        if LineCounter == 0:
            L1 = [[float(z3[0]), float(z3[1])]]
            LineCounter += 1
        elif LineCounter == 1:
            L1 += [[float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5])]]
            LineCounter += 1
        elif LineCounter == 2:
            L1 += [[float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5])]]
            MatList[matName] = PolyBeam( matName, L1[0][0],L1[0][1], L1[1],L1[2])
    elif IType2=="AbaConc":
        Echo(f"ignored {sval:s}", ff)
    elif IType2 == "misesbeam":
        MatList[matName] = MisesBeam2D([float(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6]), mass], [0, 0])
    else: 
        raise NameError("ConFemInOut::DataInput: material not defined for name", matName, IType2)
    return IType2, mass, LineCounter, L1, RotCrack, PrinStrain, ShearRetFac, S2ndCrack, S2ndCrackA, Visco, PhaseField, SDA, TraPre, LatBonTeFac, RedInt, sval 

def DataInput( f1, ff, Restart):
    # first check input for correct keywords
    z1 = f1.readline()
    z2 = z1.strip()
    z3 = z2.split(',')
    while z1!="":
        if   z3[0].find("**")>-1 :pass
        elif z3[0].find("*")>-1 and z3[0].upper() not in ["*HEADING","*NODE","*ELEMENT","*MATERIAL","*BEAM REINF","*BEAM SECTION","*SOLID SECTION","*SHELL SECTION","*STEP","*END STEP",
                                   "*DENSITY","*RCBEAM","*TCTBAR","*ELASTIC","*ELASTICORTHO","*ELASTICLT","*RCSHELL","*ISODAMAGE","*SPRING","*MISES","*LUBLINER","*RESHEET","*NLSLAB","*MICRODAMAGE","*BOND","*ELASTIC_SDA","*ELASTICR1D","*ELASTICC1D",
                                   "*CONTROLS","*STATIC","*DYNAMIC", "*DAMPING","*BUCKLING","*BOUNDARY","*DLOAD","*CLOAD","*TEMPERATURE","*EL FILE","*NODE FILE","*SOLUTION TECHNIQUE","*AMPLITUDE","*PRESTRESS","*ELSET","*NSET",
                                   "*CONCRETE DAMAGED PLASTICITY","*CONCRETE COMPRESSION HARDENING","*CONCRETE TENSION STIFFENING","*CONCRETE COMPRESSION DAMAGE","*CONCRETE TENSION DAMAGE","*DEPVAR","*USER MATERIAL","*MISESBEAM",
                                   "*PREPRINT","*PART","*END PART","*ASSEMBLY","*INSTANCE","*END INSTANCE","*END ASSEMBLY","*RESTART","*OUTPUT","*PREPRINT","*SYSTEM","*NODE OUTPUT","*ELEMENT OUTPUT","*MONITOR","*POLYBEAM","*RESTART FILE",
                                                          "*EIGENMODES"]:
            raise NameError ("ConFemInOut::DataInput: unknown keyword", z3[0])
        z1 = f1.readline()
        z2 = z1.strip()
        z3 = z2.split(',')
    f1.seek(0)
    #
    Header  = []
    MatList = {}                                                            # dictionary for materials
    SecDic  = {}                                                            # local, dictionary for sections
    BeamSec = {}                                                            # dictionary for sections, for beam section only
    BreiDic = {}                                                            # dictionary for beam reinforcement, local use only, key is a label, data is list of lists
    BreiDic["DEFAULT"] = []                                                 # default value empty dic / no reinforcement
    ElemList  = []
    ElemList_ = []
    NodeList  = []
    StepList  = []
    ElSetDic  = {}
    NSetDic   = {}
    AmpDict = {}                                                            # local, also defined for steps where it is persistent -- both are combined in following post section  
    NoLabToNoInd = {}
    ResultTypes = {}                                                        # dictionary for each element type
    RandomData = {}                                                         # local -- random data for material property definition per element
    NoCounter = 0
    StepCounter = 0
    z1 = f1.readline()
    z2 = z1.strip()
    z3 = z2.split(',')
    zLen = len(z3)
    IType  = ""
    IType2 = ""                                                     # initialization of key
    while z1!="":
        if   z3[0].upper()=="*HEADING":
            IType = "heading"
        elif z3[0].upper()=="*NODE":        
            IType = "node"
        elif z3[0].upper()=="*NSET":        
            IType = "nset"
            generate = False
            for i in z3:
                if i.upper().find("NSET=")>-1:    nset =i.split("=")[1]
                if i.upper().find("GENERATE")>-1: generate = True
                L1 = []
        #    
        elif z3[0].upper()=="*ELEMENT":
            IType = "element"
            elset = None
            bondl = None
#            bStiff= None
            for i in z3:                                                    # eltype must be determined first  
                if i.upper().find("TYPE")>-1:      eltype= i.split("=")[1].strip().upper()
                elif i.upper().find("ELSET")>-1:   elset = i.split("=")[1].strip()
                elif i.upper().find("BONDLAW")>-1: bondl = i.split("=")[1].strip().upper()
#                elif i.upper().find("BSTIFF")>-1:  bStiff = float(i.split("=")[1])
        elif z3[0].upper()=="*ELSET":
            IType = "elset"
            generate, internal = False, False
            for i in z3:
                if i.upper().find("ELSET=")>-1:   elset =i.split("=")[1]  
                if i.upper().find("GENERATE")>-1: generate = True
                if i.upper().find("INTERNAL")>-1: internal = True           # currently not used
            L1 = []
        #
        elif z3[0].upper()=="*MATERIAL":
            IType = "mat"
            for i in z3:
                if i.upper().find("NAME")>-1:     
                    matName=i.split("=")[1].upper()
                    if MatList.get(matName)!=None: raise NameError ("ConFemInOut::DataInput: material with same name already exists")
            IType2= " "                                                     # IType2 will become active ???
            mass  = 0.
            LineCounter = 0                                                 # option may have more than one data line with different format
            RedInt, RotCrack, PrinStrain, S2ndCrack, S2ndCrackA, Visco, PhaseField, SDA, sval = False, False, True, False, 0., False, False, False, "" 
            ShearRetFac=-1                                                  # see default value in init
            L1 = []
            LatBonTeFac = 1.0
            TraPre = False
        #
        elif z3[0].upper()=="*BEAM REINF":                                  # key BEAM REINFORCEMENT - another option for input of reinforcement geometry presumably if beam section is not RECT
            IType = "beamreinforcement"
            for i in range(1,zLen):                                         # loop starts with 1
                if z3[i].find("NAME")>-1:
                    breiName=z3[i].split("=")[1]                            # set current reinforcement name
                    if BreiDic.get(breiName)!=None: raise NameError ("ConFemInOut::DataInput: beam reinf with same name")
            linecounter = 0
        #
        elif z3[0].upper()=="*BEAM SECTION":                                # key BEAM SECTION
            IType = "beamsection"
            reinf = "DEFAULT"
            nRebar = 1.0   #None
            bStiff_ = None
            for i in range(1,len(z3)):                                      # loop starts with 1
                if   z3[i].upper().find("SECTION")>-1:  beamSecType=z3[i].split("=")[1]
                elif z3[i].upper().find("ELSET")>-1:    elset=z3[i].split("=")[1]
                elif z3[i].upper().find("MATERIAL")>-1: mat=z3[i].split("=")[1].upper()
                elif z3[i].upper().find("REINF")>-1:    reinf=z3[i].split("=")[1]
                elif z3[i].upper().find("NREBAR") >-1:  nRebar = float(z3[i].split("=")[1])
                elif z3[i].upper().find("BSTIFF") >-1:  bStiff_ = float(z3[i].split("=")[1])
            linecounter = 0
        elif z3[0].upper()=="*SOLID SECTION":                               # key SOLID SECTION
            IType = "solidsection"
            nRebar = None
            for i in z3:
                if   i.upper().find("ELSET")>-1:        elset   =i.split("=")[1]
                elif i.upper().find("MATERIAL")>-1:     mat     =i.split("=")[1].upper()
                elif i.upper().find("NREBAR") >-1:      nRebar = float(i.split("=")[1])
        elif z3[0].upper()=="*SHELL SECTION":                               # key SHELL SECTION
            IType = "shellsection"
            for i in range(1,zLen):                                         # loop starts with 1
                if   z3[i].find("ELSET")>-1:    elset=z3[i].split("=")[1]
                elif z3[i].find("MATERIAL")>-1: mat=z3[i].split("=")[1].upper()
            linecounter = 0
        #
        elif IType!= 'step' and z3[0].upper()=="*AMPLITUDE":                # amplitude definition out of step
            IType = "amplitude" 
            for i in z3:  
                if i.upper().find("NAME")>-1:      AmpName= i.split("=")[1] # name of amplitude
        #
        elif z3[0].upper()=="*STEP":
            ### build final list for elements #############################################################
            if len(ElemList_)>0:
                for i in range(len(ElemList_)):
                    eltype = ElemList_[i][0]
                    #                    B23__init__(self, Label,          SetLabel,       InzList,        MatName, NodeList, BeamSec, ReinfD, StateV, NData, ff, NoLabToNoInd, Material)
                    if   eltype=='B23':  ElemList +=[B23(  ElemList_[i][1],ElemList_[i][3],ElemList_[i][2],None,    NodeList, BeamSec, None,   None,   None,  ff, NoLabToNoInd, None)]
                    elif eltype=='B23I':                                    # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat 
                        #B23I__init__(self,Label,         InzList,        BeamSec, NoList,  elSet, MatName, BondLaw,                          NoLabToNoInd,Material,    NoIndToCMInd)
                        ElemList +=[B23I( ElemList_[i][1],ElemList_[i][2],BeamSec, NodeList,elset, mat,     ElemList_[i][4],NoLabToNoInd,MatList[mat],[])]
                        #               B23E__init__(self, Label,          SetLabel,       InzList,        MatName, NodeList, BeamSec, ReinfD, StateV, NData, ff, NoLabToNoInd, Material)
                    elif eltype=='B23E': ElemList +=[B23E( ElemList_[i][1],ElemList_[i][3],ElemList_[i][2],None,    NodeList, BeamSec, None,   None,   None,  ff, NoLabToNoInd, None)]
                    elif eltype=='B23EI':                                   # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat 
                        #B23EI__init__(self,Label,        InzList,        BeamSec, NoList,  elSet, MatName, BondLaw,                          NoLabToNoInd, Material,    NoIndToCMInd)
                        ElemList +=[B23EI(ElemList_[i][1],ElemList_[i][2],BeamSec, NodeList,elset, mat,     ElemList_[i][4],NoLabToNoInd, MatList[mat],[])]
                    elif eltype=='BAX23': ElemList += [BAX23(ElemList_[i][1],ElemList_[i][3],ElemList_[i][2], NodeList,NoLabToNoInd)]
                    elif eltype=='BAX23I':        #          label,           elset            InzList          BondLaw
                        ElemList += [                 BAX23I(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], ElemList_[i][4], NodeList, NoLabToNoInd)]
                    elif eltype=='BAX23E':                 # label,           elset            InzList          BondLaw
                        ElemList += [                 BAX23EI(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], None,            NodeList, NoLabToNoInd)]
                    elif eltype=='BAX23EI':                 # label,           elset            InzList          BondLaw
                        ElemList += [                 BAX23EI(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], ElemList_[i][4], NodeList, NoLabToNoInd)]
                    elif eltype=='B21':   ElemList +=[B21(   ElemList_[i][1],ElemList_[i][2],ElemList_[i][3], NodeList,NoLabToNoInd)]
                    elif eltype=='BAX21': ElemList +=[BAX21( ElemList_[i][1],ElemList_[i][2],ElemList_[i][3], NodeList,NoLabToNoInd)]
#                   elif eltype=='BAX21I':ElemList +=[BAX21I(ElemList_[i][1],ElemList_[i][2],ElemList_[i][3], NodeList, NoLabToNoInd)]
                    elif eltype=='B21E':  ElemList +=[B21E(  ElemList_[i][1],ElemList_[i][2],ElemList_[i][3], NodeList,NoLabToNoInd)]
                    elif eltype=='BAX21E':ElemList +=[BAX21E(ElemList_[i][1],ElemList_[i][2],ElemList_[i][3], NodeList,NoLabToNoInd)]
                    elif eltype=='BAX21EI':
                        ElemList += [BAX21EI( ElemList_[i][1], ElemList_[i][2],ElemList_[i][3],ElemList_[i][4], NodeList,NoLabToNoInd)]
                    elif eltype=='S2D6': ElemList +=[S2D6( ElemList_[i][1],ElemList_[i][2],ElemList_[i][3],           ElemList_[i][4], ElemList_[i][5], SecDic[ElemList_[i][2]], MatList[ElemList_[i][4]].StateVar, MatList[ElemList_[i][4]].NData, NoLabToNoInd)]
                    #                   T2D2__init__(self, Label,         SetLabel,        InzList,         MatName, NodeList, SolSecDic, StateV, NData, Val,                  NoLabToNoInd)
                    elif eltype=='T2D2': ElemList +=[T2D2( ElemList_[i][1],ElemList_[i][3],ElemList_[i][2], None,    NodeList, SecDic,    None,   None,  ElemList_[i][4], NoLabToNoInd)]
                    elif eltype=='T2D2I':                                   # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat
                        # T2D2I__init__(self,Label,         InzList,        SolSecDic,Val,            NoList,  elSet,          MatName,  BondLaw,       NoLabToNoInd, Material, NoIndToCMInd)
                        ElemList += [T2D2I(ElemList_[i][1], ElemList_[i][2], SecDic, ElemList_[i][4], NodeList, ElemList_[i][3], mat, ElemList_[i][5], NoLabToNoInd, MatList[mat], [])]
                    elif eltype=='T2D3':
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat
                        #T2D3__init__(self, Label,           InzList,         AA,              NoList,   elSet,           MatName, NoLabToNoInd, Material, NoIndToCMInd):
                        ElemList += [ T2D3( ElemList_[i][1], ElemList_[i][2],        ElemList_[i][4],NodeList,SecDic, elset,          mat,                                      NoLabToNoInd, MatList[mat], [])]
                    elif eltype=='T2D3I':                                   # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat
                        # T2D3I__init__(self,Label,         InzList,        SolSecDic,Val,             NodeList, elSet, MatName, BondLaw,                         NoLabToNoInd, Material, NoIndToCMInd):
                        ElemList += [T2D3I(ElemList_[i][1], ElemList_[i][2], SecDic, ElemList_[i][4], NodeList, elset, mat, ElemList_[i][5], NoLabToNoInd, MatList[mat], [])]
                    elif eltype=='T3D2': ElemList += [T3D2(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], None, NodeList, None, None, None, ElemList_[i][4], NoLabToNoInd)]
                    elif eltype=='T3D2I':                                   # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat
                        ElemList += [T3D2I(ElemList_[i][1], ElemList_[i][2], SecDic, ElemList_[i][4], NodeList, elset, mat, ElemList_[i][5], NoLabToNoInd, MatList[mat], [])]
                    elif eltype=='T3D3':
                        # T3D3__init__(self,Label,         SetLabel,        InzList,         MatName, NodeList, SolSecDic, StateV, NData, Val,             NoLabToNoInd)
                        ElemList += [T3D3(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], None,    NodeList, None,      None,   None,  ElemList_[i][4], NoLabToNoInd)]
                    elif eltype=='T3D3I':                                   # full data needed for creation of bond elements
                        elset = ElemList_[i][3]
                        mat   = SecDic[elset].Mat
                        # T3D3I__init__(self,Label,         InzList,        SolSecDic,Val,             NodeList, elSet, MatName, BondLaw,         NoLabToNoInd, Material, NoIndToCMInd)
                        ElemList += [T3D3I(ElemList_[i][1], ElemList_[i][2], SecDic, ElemList_[i][4], NodeList, elset, mat, ElemList_[i][5], NoLabToNoInd, MatList[mat], [])]
                    elif eltype=='TAX2':
                        # T2D2__init__(self,Label,           SetLabel,        InzList,         NodeList, Val,             NoLabToNoInd)
                        ElemList += [TAX2(  ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], NodeList, SecDic, ElemList_[i][4], NoLabToNoInd)]
                    elif eltype=='TAX2I':
                        ElemList += [TAX2I(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], NodeList, SecDic, ElemList_[i][4], NoLabToNoInd, ElemList_[i][5])]
                    elif eltype=='TAX3':
                        ElemList += [TAX3( ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], NodeList, SecDic, ElemList_[i][4], NoLabToNoInd)]
                    elif eltype=='TAX3I':
                        ElemList += [TAX3I(ElemList_[i][1], ElemList_[i][3], ElemList_[i][2], NodeList, SecDic, ElemList_[i][4], NoLabToNoInd, ElemList_[i][5])]

                ElemList_ = []
            # finally, complete assignment of elset to elements for considering *ELSET keyword                    
            for elset in ElSetDic:                                          # loop over elsets defined with *ELSET, ElSetDic is local
                    for i in ElSetDic[elset]:                               # loop over elements
                        j = FindIndexByLabel(ElemList, i)
                        el = ElemList[j] 
                        if el.Set != None: 
                            Echo(f"element {el.Label:d} has already set directly assigned, skipped in related assignment", ff)
                        else:
                            el.Set = elset
            # finally assign element types to SecDic    
            for i, e in enumerate(ElemList):
                elset  = e.Set
                eltype = e.Type
                if elset in SecDic:                                         # SecDic is global
                    s = SecDic[elset]
                    if eltype not in s.ElemTypes:
                        s.ElemTypes += [eltype]
                else: raise NameError("ConFemInOut::DataInput: section not found for element set",elset,e.Label)
            for s in SecDic:
                if len(SecDic[s].ElemTypes)==1: 
                    pass
                else:
                    for t in SecDic[s].ElemTypes:
                        if t in ['CPS3','CPE3','CPS4','CPE4']: pass # CP*3, CP*4 may be in same set
                        else: raise NameError("ConFemInOut::DataInput: element set should have one element type assigned",SecDic[s].ElemTypes)
            ### end of build final list for elements #######################################################
            ###
#            from ConFem import CheckInz
#            CheckInz(ElemList, 681)
            ###
            # initialize step
            IType = "step"
            IType2=" "                                                  # IType2 will become active ???
            StepCLoad  = []
            StepBoundC = []
            StepDLoad  = []
            StepTemp   = []
            StepCounter += 1
            if len(StepList) == 0:
                StepList += [Step()]                                        # initializing step
                for i in range(1,zLen):                                     # loop starts with 1
                    if z3[i].upper().find("NLGEOM")>-1:
                        if z3[i].split("=")[1]=="YES": StepList[len(StepList)-1].NLGeom = True # acts on itself
                    if z3[i].upper().find("SCALEMASS") > -1:
                        StepList[len(StepList) - 1].ScaleMass = float(z3[i].split("=")[1])   # acts on itself
            else:
                StepList += [copy.deepcopy(StepList[-1])]                   # take over everything for additional step from previous step
                StepList[-1].PrestList = []                                 # prestressing data may have different formats in subsequent step, are fully initialized for each step  
                for i in range(1,zLen):                                     # loop starts with 1
                    if z3[i].find("NLGEOM")>-1:                             
                        if z3[i].split("=")[1]=="YES": StepList[len(StepList)-1].NLGeom = True # acts on last

            Step_ = StepList[len(StepList)-1]                               # Step_ is a pointer to last entry of StepList -- in case of 1st step an initialized instance of object Step
            # end initialize step
        elif z3[0].upper() in ['*PREPRINT','*PART','*END PART','*ASSEMBLY','*INSTANCE','*END INSTANCE','*END ASSEMBLY','*RESTART','*OUTPUT','*SYSTEM']:
            Echo(f"ignored {z3[0]:s}", ff)                                  # Abaqus stuff
        elif z3[0]=="" or z3[0].find("**")>-1: # or (z3[0][0]=="*" and IType2==""): 
            pass                                                            # blank line / comment line / no secondary keyword
        elif z3[0].find("*")>-1 and IType!="mat" and IType!="step":     # unknown key, comment, connected line
            raise NameError("ConFemInOut::DataInput: unknown keyword",z3[0])
        #
        elif IType!="":
            if   IType=="heading":
                Header += [z3[0]]
            elif IType=="node":
                if zLen<3: raise NameError ("not enough nodal data")
                if zLen==3:   NodeList += [Node( int(z3[0]), float(z3[1]), float(z3[2]), 0., 0.,0.,0. )]
                elif zLen==4: NodeList += [Node( int(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), 0.,0.,0. )]
                elif zLen==7: NodeList += [Node( int(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]), float(z3[4]), float(z3[5]), float(z3[6]) )]
                NoLabToNoInd[ int(z3[0]) ] = NoCounter
                NoCounter += 1
            #
            elif IType=="nset":
                if generate:
                    L1 = [ i for i in range( int(z3[0]), int(z3[1])+1, int(z3[2]) ) ]
                else:
                    if len(z3[-1])==0: nL = len(z3)-1                       # for data like .... 632, 
                    else:              nL = len(z3)
                    L1 = L1 + [int(z3[i]) for i in range(nL)]                                            # 
                NSetDic[nset] = L1                                          # dic only used locally, brute repeated updating, not a smart way
            #
            elif IType=="element":
                # what follows refers to ElList_
                if eltype=='B23':
                    if bondl!=None: ElemList_+=[['B23I', int(z3[0]),[int(z3[1]),int(z3[2])], elset,bondl]]
                    else:           ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2])], elset]]
                elif eltype=='B23E':
                    if len(z3)==6:          
                        '''accept 2 randomly distributed variables'''
                        RandomField_Routines().add_elset(elset)
                        RandomField_Routines().add_propertyEntry(float(z3[4]),elset,int(z3[0]))
                        RandomField_Routines().add_propertyEntry(float(z3[5]),elset,int(z3[0]))
#                    if bondl!=None:  ElemList_+=[['B23EI',int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])],elset,bondl]]
#                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])],elset      ]]
                    if bondl!=None:  ElemList_+=[['B23EI',int(z3[0]),[int(z3[1]),int(z3[3]),int(z3[2])],elset,bondl]] # enhanced node moves last in input, right node now in center of incidence
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[3]),int(z3[2])],elset      ]]
                elif eltype == 'BAX23':
                    if bondl!=None:  ElemList_ += [["BAX23I", int(z3[0]), [int(z3[1]), int(z3[2])], elset, bondl]]
                    else:            ElemList_ += [[eltype, int(z3[0]), [int(z3[1]), int(z3[2])], elset]]
                elif eltype == 'BAX23E':
                    if bondl!=None:  ElemList_ += [["BAX23EI", int(z3[0]),[int(z3[1]),int(z3[3]),int(z3[2])],elset,bondl]]
                    else:            ElemList_ += [[eltype,    int(z3[0]),[int(z3[1]),int(z3[3]),int(z3[2])], elset]]
                elif eltype in ['B21','BAX21']:
                    if bondl!=None:  ElemList_+=[[eltype, int(z3[0]),elset,[int(z3[1]),int(z3[2])],bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),elset,[int(z3[1]),int(z3[2])]]]
                elif eltype=='B21E':
                    ElemList_+=[[eltype, int(z3[0]),elset,[int(z3[1]),int(z3[2]),int(z3[3])]]]
                elif eltype=='BAX21E':
                    if bondl!=None:  ElemList_+=[['BAX21EI', int(z3[0]),elset,[int(z3[1]),int(z3[2]),int(z3[3])],bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),elset,[int(z3[1]),int(z3[2]),int(z3[3])]]]
                elif eltype=='S2D6': ElemList_+=[[eltype, int(z3[0]),elset,[int(z3[1]),int(z3[2])],           None,NodeList,None]]
                elif eltype=='T2D2':
                    if zLen==3: Val = 1.0
                    else:       Val = float(z3[3])                          # modifying factor for cross section
                    if bondl!=None: ElemList_+=[['T2D2I',int(z3[0]),[int(z3[1]),int(z3[2])], elset, Val, bondl]]
                    else:           ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2])], elset, Val]]
                elif eltype=='T2D3':
                    if zLen==4: Val = 1.0
                    else:       Val = float(z3[4])                          # modifying factor for cross section
                    if bondl!=None:  ElemList_+=[['T2D3I',int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, Val, bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, Val]]
                elif eltype=='T3D2':
                    if zLen==3: Val = 1.0
                    else:       Val = float(z3[3])                          # modifying factor for cross section
                    if bondl!=None:  ElemList_+=[['T3D2I',int(z3[0]),[int(z3[1]),int(z3[2])], elset, Val, bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2])], elset, Val]]
                elif eltype=='T3D3':
                    if zLen==4: Val = 1.0
                    else:       Val = float(z3[4])                          # modifying factor for cross section
                    if bondl!=None:  ElemList_+=[['T3D3I',int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, Val, bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, Val]]
                elif eltype=='TAX2':
                    if bondl!= None: ElemList_+=[["TAX2I",int(z3[0]),[int(z3[1]),int(z3[2])], elset, 1.0, bondl]]
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2])], elset, 1.0 ]]
                elif eltype=='TAX3':
                    if bondl!= None: ElemList_+=[["TAX3I",int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, 1.0, bondl, None]]
                    else:            ElemList_+=[[eltype, int(z3[0]),[int(z3[1]),int(z3[2]),int(z3[3])], elset, 1.0 ]]

                # what follows refers to ElemList
                elif eltype in ['CPS4','CPE4']:
                    if eltype=='CPS4': PlStrFl = True                       # flag for plane stress (True->plane stress, False->plane strain)
                    else:              PlStrFl = False
                    l3 = len(z3)
                    if z3[-1].strip()=="": l3=l3-1                          # white spaces after comma
                    if   l3==5:                                             # quad element
                        ElemList += [CPE4( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], NodeList, PlStrFl, 2, NoLabToNoInd)]
                    elif l3==9:                                             # quad element -- with further random ip parameter
                        ElemList += [CPE4( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], NodeList, PlStrFl, 2, NoLabToNoInd)]
                        RandomData[int(z3[0])] = [float(z3[5]),float(z3[6]),float(z3[7]),float(z3[8])]
                    elif l3==4:                                             # tri element
                        ElemList += [CPE3( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3])], NodeList, PlStrFl, NoLabToNoInd)]
                    else: raise NameError("ConFemInOut::ReadInputFile: wrong number of nodes for CPXn elements")
                elif eltype in ['CPS4R']:
                    PlStrFl = True
                    ElemList += [CPE4( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], None,None,        NodeList, None,          None,                  None,               PlStrFl,  1, None,               ff, NoLabToNoInd)]
                elif eltype in ['CPS3','CPE3']:
                    if eltype=='CPS3': PlStrFl = True                       # flag for plane stress (True->plane stress, False->plane strain)
                    else:              PlStrFl = False
                    ElemList += [CPE3( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3])], NodeList, PlStrFl, NoLabToNoInd)]
                elif eltype=='CAX4':
                    ElemList += [CAX4(int(z3[0]), elset, [int(z3[1]), int(z3[2]), int(z3[3]), int(z3[4])], NodeList, 2,NoLabToNoInd)]
                elif eltype=='CPS6':
                    ElemList += [CPE6( int(z2[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4]),int(z3[5]),int(z3[6])], None,None, NodeList, None, None, None, True, NoLabToNoInd)]
                elif eltype=='CPE6': 
                    ElemList += [CPE6( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4]),int(z3[5]),int(z3[6])], None,None, NodeList, None, None, None, False,NoLabToNoInd)]
#                elif eltype=='C3D4':
#                    if Flag:
##                        SectionList += [ SolidSection( elset, mat, elsetVal, "FEM") ]
#                        ElemList += [C3D4( int(z3[0]), [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], elsetVal, NodeList, elset, mat, NoLabToNoInd, None, [])]
#                        Flag = False
#                    else:
#                        ElemList += [C3D4( int(z3[0]), [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], elsetVal, NodeList, elset, mat, NoLabToNoInd, None, [])]
                elif eltype=='C3D8':
#                   ElemList += [C3D8( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4]),int(z3[5]),int(z3[6]),int(z3[7]),int(z3[8])], mat,MatList[mat], NodeList, MatList[mat].StateVar, MatList[mat].NData,    MatList[mat].RType, NoLabToNoInd)]
                    ElemList += [C3D8( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4]),int(z3[5]),int(z3[6]),int(z3[7]),int(z3[8])], None,None,        NodeList, None,                  None,               2, None,               NoLabToNoInd)]# for 2 -->nI for default integration
                elif eltype=='C3D8R':
                    ElemList += [C3D8( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4]),int(z3[5]),int(z3[6]),int(z3[7]),int(z3[8])], None,None,        NodeList, None,                  None,               1, None,               NoLabToNoInd)]# for 2 -->nI for default integration
                elif eltype=='S1D2':
                    if len(z3)==3: LenFactor = 1.0
                    else:          LenFactor = float(z3[3])
                    ElemList += [S1D2( int(z3[0]), elset, [int(z3[1]),int(z3[2])], None, NodeList, None, None,None, NoLabToNoInd, LenFactor)]
                elif eltype=='T1D2': 
                    if len(z3)==4:
                        '''accept 1 randomly distributed variables'''
                        RandomField_Routines().add_elset(elset)
                        RandomField_Routines().add_propertyEntry(float(z3[3]),elset,int(z3[0]))
                    ElemList += [T1D2( int(z3[0]), elset, [int(z3[1]),int(z3[2])], None,None, NodeList, None, None, None, None, NoLabToNoInd)]
                elif eltype=='SB3':
                    ElemList+= [SB3( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3])], None,NodeList, None, None, None, NoLabToNoInd)]
                elif eltype=='SH4':
                    ElemList+= [SH4( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3]),int(z3[4])], None,None, NodeList, None, None, None, None, NoLabToNoInd)]
                elif eltype=='SH3':
#                    NData = MatList[mat].NData
#                    NStateVar = MatList[mat].StateVar
#                    ElemList+= [SH3( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3])],            mat,MatList[mat], NodeList, SecDic[elset], NStateVar, NData, False, NoLabToNoInd)]
                    ElemList+= [SH3( int(z3[0]), elset, [int(z3[1]),int(z3[2]),int(z3[3])],            None,None, NodeList, None, None, None, None, NoLabToNoInd)]
                else:
                    raise NameError("unknown element type "+eltype)
            elif IType=="elset":
                if generate:
                    L1 = [ i for i in range( int(z3[0]), int(z3[1])+1, int(z3[2]) ) ]
                else:
                    if len(z3[-1])==0: nL = len(z3)-1                       # for data like .... 632, 
                    else:              nL = len(z3)
                    L1 = L1 +  [int(z3[i]) for i in range(nL)]              # 
                ElSetDic[elset] = L1                                        # ElSetDic only used locally, brute repeated updating, not a smart way 
            # materials
            elif IType=="mat": 
                IType2, mass, LineCounter, L1, RotCrack, PrinStrain, ShearRetFac, S2ndCrack, S2ndCrackA, Visco, PhaseField, SDA, TraPre, LatBonTeFac, RedInt, sval \
                        = MatInput( z3, zLen, MatList, matName, IType2, mass, LineCounter, L1, RotCrack, PrinStrain, ShearRetFac, S2ndCrack, S2ndCrackA, Visco, PhaseField, SDA, TraPre, LatBonTeFac, RedInt, sval, ff )
            # reinforcement
            elif IType=="beamreinforcement":
                if zLen<4: raise NameError("ConFemInOut::DataInput: too less reinforcement data")
                if zLen>4 and z3[4][0:2]=='TC': RType='TC'                  # textile or other reinforcement than ordinary rebars
                else:                           RType='RC'
                if linecounter==0:
                    BreiDic[breiName] = [[float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),RType]] # reinforcement cross section and tension stiffening data
                    linecounter = linecounter+1
                else:
                    BreiDic[breiName] += [[float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),RType]] # reinforcement cross section and tension stiffening data
            # sections
            elif IType=="solidsection":                                     # data SOLID SECTION
                SecDic[elset] = SolidSection( elset, mat,nRebar, float(z3[0]))
            elif IType=="beamsection":                                      # data BEAM SECTION 
                if beamSecType=="RECT":
                    if linecounter==0:
                        SecDic[elset] = BeamSectionRect( beamSecType,elset,mat,nRebar,bStiff_, float(z3[0]), float(z3[1]), [], ReinfName=reinf) # concrete cross section data
                        linecounter = linecounter+1
                    else:
                        if zLen<4: raise NameError ("ConFemInOut::DataInput: too less reinforcement data")
                        if zLen>4 and z3[4][0:2]=='TC': RType='TC' # textile or other reinforcement than ordinary rebars
                        else:                           RType='RC'
                        SecDic[elset].Reinf += [[float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),RType]] # reinforcement cross section and tension stiffening data
                elif beamSecType.upper()=="CIRCLE":
                    SecDic[elset] = BeamSectionCircle( beamSecType,elset,mat,nRebar,bStiff_, float(z3[0]), ReinfName=reinf) # concrete cross section data
                elif beamSecType.upper()=="POLYLINE":
                    if linecounter==0:
                        SecDic[elset] = BeamSectionPoly( beamSecType,elset,mat,nRebar,bStiff_, [], ReinfName=reinf) # concrete cross section data
                        SecDic[elset].AddPoint(float(z3[0]), float(z3[1]))
                        linecounter = linecounter+1
                    else:
                        SecDic[elset].AddPoint(float(z3[0]), float(z3[1]))
                else: raise NameError ("ConFemInOut::DataInput: unknown beam section type", beamSecType)
            elif IType=="shellsection":                                     # data SHELL SECTION
                if linecounter==0:
                    SecDic[elset] = ShellSection( elset, mat, float(z3[0]), [] )
                    linecounter = linecounter+1
                else:
                    SecDic[elset].Reinf += [[float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4])]] # reinforcement cross section and tension stiffening data, see also ConFemEelem::SH4.__init__
            #  amplitude definition out of step
            elif IType == 'amplitude':
                AList_ = []
                for i in range(0,len(z3),2): 
                    try:     AList_ += [[float(z3[i]),float(z3[i+1])]]
                    except:  break 
                AmpDict[AmpName] = AList_                                   # AmpDict is a local dictionary
            # step ##########################################################
            elif IType=="step":
                if z3[0].find("*")>-1 and IType2=="stepPrestress":
                    if PreName not in Step.PreSDict:                        # underlined defined below 
                        Step.PreSDict[PreName] = PreStress(PreType, L1, PList, AmpLbl) # PreSDict is common for all steps -- underlined defined below
                                                                            # construct should work for multiple *prestress sections, current *PTRESTRESS is closed by new  
#                Step_ = StepList[len(StepList)-1]                           # Step_ is a pointer to last entry of StepList -- in case of 1st step an initialized instance of object Step
                if   z3[0].upper()=="*STATIC":
                    IType2 = "static"
                    Step_.Dyn        = False
                    Step_.Damp       = False
                    Step_.Eigenmodes = False
                    Step_.Buckl      = False
                    for i in z3:
                        if i.upper().find("ARCLENGTH")>-1: 
                            Step_.ArcLen = True
                            if i.upper().find("=")>-1:
                                Step_.ArcLenV = float(i.split("=")[1])
                            else:
                                Step_.ArcLenV = -1                          # should be computed from time step with first time steo increment
                            Step_.ArcLenNodes = []
                            if len(z3)>2:                                   # nodes for arc length selected
                                for n in z3[2:]:
                                    Step_.ArcLenNodes += [ int(n) ]
                elif z3[0].upper()=="*EIGENMODES":
                    Step_.Eigenmodes = True
                    Step_.Dyn = True
                    IType2 = "eigenmodes"
                elif z3[0].upper()=="*DYNAMIC":    
                    Step_.Dyn = True
                    Step_.varTimeSteps = False
                    Step_.Eigenmodes = False
                    IType2 = "dynamic"
                    for i in z3:
                        if i.upper().find("EIGENMODE2TIMESTEP") > -1:
                            Step_.Eigenmode2TS = True
                            try:
                                i.find("=") > -1
                                Step_.Eigenmode2TSrel = float(i.split("=")[1])
                            except:
                                raise NameError("ConFemInOut::DataInput: value required for ",i)
                elif z3[0].upper()=="*SOLUTION TECHNIQUE":                  # key STEP SOLUTION TECHNIQUE
                    for i in range(1,zLen):
                        if z3[i].upper().find("TYPE")>-1: SolType=z3[i].split("=")[1].upper()
                    if SolType=="QUASI-NEWTON": Step_.SolType = 'BFGS'
                    elif SolType=="MODIFIED-NR": Step_.SolType = 'MNR'
                elif z3[0]=="*BUCKLING":                                    # key STEP BUCKLE
                    Step_.Buckl = True
                elif z3[0]=="*CONTROLS" and z3[1].strip().upper()=="PARAMETERS=FIELD":          # key STEP CONTROLS FIELD
                    IType2 = "stepControls0"
                elif z3[0]=="*CONTROLS" and z3[1].strip().upper()=="PARAMETERS=TIME INCREMENTATION":# key STEP CONTROLS TIME INCREMENTATION
                    IType2 = "stepControls1"                
                elif z3[0].upper()=="*CONTROLS":  
                    for i in z3:  
                        if i.upper().find("ITOL")>-1:  Step_.IterTol = float(i.split("=")[1])             # tolerance for equilibrium residuum
                        if i.upper().find("NITER")>-1: Step_.IterNum =  int(i.split("=")[1])   # maximum number of equilibrium iterations allowed
                elif z3[0].upper()=="*AMPLITUDE":
                    IType2 = "amplitude" 
                    for i in z3:  
                        if i.upper().find("NAME")>-1:      AmpName= i.split("=")[1]                     # name of amplitude
                    if AmpName in Step_.AmpDict: raise NameError ("ConFemInOut::DataInput: amplitude with same name already exists",AmpName)
                elif z3[0].upper()=="*BOUNDARY":   
                    IType2 = "boundary"
                    BoundAmp = "Default"                                                        # default amplitude initialized in Step.__init__
                    AddVal = False                                                              # default                                          
                    for i in z3:  
                        if i.upper().find("AMPLITUDE")>-1: BoundAmp = i.split("=")[1]                   # amplitude of current *BOUNDARY set
                    for i in z3:
                        if i.upper().find("OP")>-1:
                            if i.split("=")[1].upper()=="ADD": AddVal = True                    # add given value to final value of last step to get actual prescribed value
                elif z3[0].upper()=="*CLOAD":      
                    IType2 = "cload"
                    CLoadAmp = "Default"                                                        # default amplitude initialized in Step.__init__
                    for i in z3:  
                        if i.upper().find("AMPLITUDE")>-1: CLoadAmp = i.split("=")[1]           # amplitude of current *CLOAD set
                elif z3[0].upper()=="*DLOAD":                   # key STEP DLOAD
                    IType2 = "stepDLoad"
                    DLoadAmp = "Default"
                    for i in z3:
                        if i.upper().find("AMPLITUDE")>-1: DLoadAmp=i.split("=")[1]
                elif z3[0].upper()=="*TEMPERATURE":             # key STEP TEMPERATURE
                    IType2 = "stepTemp"
                    TempAmp="Default"
                    for i in range(1,zLen):
                        if z3[i].find("AMPLITUDE")>-1: TempAmp=z3[i].split("=")[1]
                elif z3[0].upper()=="*PRESTRESS":              # key STEP PRESTRESS
                    IType2 = "stepPrestress"
                    AmpLbl="Default"
                    for i in range(1,zLen):
                        if z3[i].upper().find("AMPLITUDE")>-1: AmpLbl=z3[i].split("=")[1]
                        if z3[i].upper().find("TYPE")>-1: PreType=(z3[i].split("=")[1])         # used above
                        if z3[i].upper().find("NAME")>-1: PreName=z3[i].split("=")[1]
                    Step_.PrestList += [PreStreSt(PreName, AmpLbl)]                             # collects all prestressing labels and amplitudes for current step
                    PList = []
                    linecounter = 0
                elif z3[0].upper()=="*EL FILE":   
                    IType2 = "stepElFile"
                    for i in z3:  
                        if i.find("FREQUENCY")>-1:
                            Step_.ElFilList += [float(i.split("=")[1])]                         # element output interval
                elif z3[0].upper()=="*NODE FILE":
                    IType2 = "stepNoFile"
                    for i in z3:  
                        if i.find("FREQUENCY")>-1:
                            Step_.NoFilList += [float(i.split("=")[1])]                         # node output interval
                elif z3[0].upper() == "*RESTART FILE":
                    for i in z3:
                        if i.find("FREQUENCY") > -1:
                            Step_.ReFilList += [float(i.split("=")[1])]                         # restart output interval
                elif z3[0].upper()=="*NODE OUTPUT":
                    IType2 = "stepNoOut"
                elif z3[0].upper()=="*ELEMENT OUTPUT":
                    IType2 = "stepElOut"
                elif z3[0].upper()=="*DAMPING":                                                 # key STEP DAMPING
                    Step_.Damp = True
                    for i in z3:
                        if i.upper().find("ALPHA")>-1: Step_.RaAlph=float(i.split("=")[1])
                    for i in z3:
                        if i.upper().find("BETA")>-1:  Step_.RaBeta=float(i.split("=")[1])
                    for i in z3:
                        if i.upper().find("EIGENVAL2BETA") > -1:
                            Step_.EigenVal2Beta = True
                            Step_.Zeta   = float(i.split("=")[1])                               # damping related to critical damping
                #
                # IType -- finalizing current step
                elif z3[0].upper().strip().replace(" ","")=="*ENDSTEP":
                    if len(StepBoundC)> 0:                                  # already taken from previous step in case StepCounter>1 and len==0
                        bcKeys = []
                        if StepCounter>1:                                   
                            for bc in StepBoundC: bcKeys += [str(bc.NodeLabel)+str(bc.Dof)]
                            for bc in StepList[-2].BoundList:               # check whether identical bc (node, dof) already defined in previous step -- if not append 
                                if str(bc.NodeLabel)+str(bc.Dof) not in bcKeys: StepBoundC += [bc]
                        else:
                            StepBoundC_ = []
                            for bc in StepBoundC:                           # check whether identical bc (node, dof) already defined in current step -- if not append
                                if str(bc.NodeLabel)+str(bc.Dof) not in bcKeys: 
                                    StepBoundC_ += [bc]
                                    bcKeys += [str(bc.NodeLabel)+str(bc.Dof)]
                                else: 
                                    Echo(f"Boundary condition Node {bc.NodeLabel:d}, dof {bc.Dof:d} has already been defined for this step -- used the former one, skipped this current", ff)
                            StepBoundC = StepBoundC_
                        Step_.BoundList = StepBoundC                        # Step_ is current step
                    if len(StepCLoad) > 0: Step_.CLoadList = StepCLoad
                    if len(StepDLoad) > 0: Step_.DLoadList = StepDLoad
                    if len(StepTemp)  > 0: Step_.TempList  = StepTemp
                # end finalizing current step
                #
                # all IType2 for IType "step"
                elif IType2 in ["static"]:
                    Step_.TimeStepVar, Step_.TimeTargVar =  [float(z3[0])], [float(z3[1])]
#                    for i in range( 2, int(len(z3)/2) ):                    # following target decreases
                    for i in range(1, int(len(z3) / 2)):                    # this loop should be active for more than one pair
#                        if float(z3[2*i+1]) < Step_.TimeTargVar[-1]:
                        if float(z3[2 * i + 1]) > Step_.TimeTargVar[-1]:
                            Step_.TimeTargVar += [float(z3[2*i+1])]
                            Step_.TimeStepVar += [float(z3[2*i  ])]
                        else: Echo(f"skipped time increment / target {z3[2*i]:s}/{z3[2*i+1]:s} as time target decreases", ff)
                    if len(Step_.TimeStepVar) != len(Step_.TimeTargVar): raise NameError ("ConFemInOut::DataInput: fault in time step, time target definition")
                    if len(Step_.TimeStepVar) == 1:
                        Step_.TimeStep = Step_.TimeStepVar[0]
                        Step_.TimeTarg = Step_.TimeTargVar[0]
                    else:
                        Step_.varTimeSteps = True
                elif IType2 == "eigenmodes":
                    Step_.EigenmodesN = int(z3[0])
                elif IType2=="dynamic":                                     # variable time step not allowed within step as mass matrix is pre-determined
                    Step_.TimeStep = float(z3[0])
                    Step_.TimeTarg = float(z3[1])
                elif IType2=="stepControls0":                              # data STEP CONTROLS
                    Step_.IterTol = float(z3[0])
                elif IType2=="stepControls1":                               # data STEP CONTROLS
                    Step_.IterNum = int(z3[0])
                elif IType2=="amplitude":
                    AList_ = []
                    for i in range(0,len(z3),2): 
                        try:     AList_ += [[float(z3[i]),float(z3[i+1])]]
                        except:  break 
                    Step_.AmpDict[AmpName] = AList_                         # AmpDict is a dictionary with key! More are allowed for this step.
                elif IType2=="boundary":
                    def AddToStepBoundC( z_, StepBoundC):
                        i0 = FindIndexByLabel(NodeList, int(z_))
                        if i0>-1:
                            nDofRange =  int(z3[2])-int(z3[1])+1
                            if   zLen==3:                                   # thats for omitted value indicating 0 but there should be no comments following
                                for i in range(nDofRange): StepBoundC += [Boundary( int(z_), int(z3[1])+i, 0.,           NodeList, BoundAmp, AddVal)]
                            else: 
                                for i in range(nDofRange): StepBoundC += [Boundary( int(z_), int(z3[1])+i, float(z3[3]), NodeList, BoundAmp, AddVal)]
                        else: raise NameError ("ConFemInOut::DataInput: boundary condition node not found")
                    try:                                                    # for node number given
                        AddToStepBoundC( int(z3[0]), StepBoundC)
                    except:                                                 # for node set label given
                        z4 = z3[0].split('.')                               # Abaqus may extend set labels xy with Part_a_b.xy
                        if len(z4)>1: L = NSetDic[z4[1]]
                        else:         L = NSetDic[z3[0]]                    # z3[0] is node set label
                        for z in L: AddToStepBoundC( z, StepBoundC)
                elif IType2=="cload":
                    def AddToStepCLoad( z_, StepCLoad): 
                        i0 = FindIndexByLabel(NodeList, int(z_))   # int(z3[0])
                        if i0>-1: 
                            StepCLoad += [CLoad( z_, int(z3[1]), float(z3[2]), NodeList, CLoadAmp)]
                        else:     raise NameError ("ConFemInOut::DataInput: concentrated load node not found")
                    try:
                        AddToStepCLoad( int(z3[0]), StepCLoad)
                    except: 
                        z4 = z3[0].split('.')                               # Abaqus may extend set labels xy with Part_a_b.xy
                        if len(z4)>1: L = NSetDic[z4[1]]
                        else:         L = NSetDic[z3[0]]                    # z3[0] is node set label
                        for z in L: AddToStepCLoad( z, StepCLoad)
                elif IType2=="stepDLoad":                                   # data STEP DLOAD
                    elset = z3[0]
                    DLoadFlag = False
                    for i in SecDic:
                        if i == elset:
                            DLoadFlag = True
                            break
                    if not DLoadFlag: raise NameError ("ConFemInOut::DataInput: element set not found for *DLOAD definition",elset)
                    StepDLoad += [DLoad( elset, int(z3[1]), float(z3[2]), DLoadAmp)]
                elif IType2=="stepTemp":                 # data STEP TEMPERATURE
                    i0 = FindIndexByLabel(NodeList, int(z3[0]))
                    if i0>-1: 
                        StepTemp += [Temperature( int(z3[0]), [float(z3[i]) for i in range(1,len(z3))], NodeList, TempAmp)]                
                    else:     raise NameError ("ConFemInOut::DataInput: temperature load node not found")
                elif IType2=="stepPrestress":            # data STEP PRESTRESS
                    if linecounter == 0:
                        if   zLen==7: val = 0.
                        elif zLen==8: val = float(z3[7])
                        else: raise NameError ("ConFemInOut::ReadInputFile: number of reinforcement data")
                        L1 = [float(z3[0]),float(z3[1]),float(z3[2]),float(z3[3]),float(z3[4]),float(z3[5]),float(z3[6]),val,mass] # prestressing parameter
                        linecounter = linecounter +1
                    else: PList += [[ FindIndexByLabel(ElemList, int(z3[0])), [float(z3[i]) for i in range(1,len(z3))] ]]
                elif IType2=="stepElFile":               # data STEP EL File  -- for compatibility reasons
                    pass
                elif IType2=="stepNoFile":               # data STEP NODE File -- for compatibility reasons
                    pass
                elif IType2=="stepNoOut":               # data STEP NODE File -- for compatibility reasons
                    Echo(f"ignored *Node Output", ff)
                elif IType2=="stepElOut":               # data STEP NODE File -- for compatibility reasons
                    Echo(f"ignored *Element Output", ff)
            # end step ######################################################
        #
        z1 = f1.readline()
        z2 = z1.strip()
        if z2==",": 
            z2 = "1.0"                                                      # might happen with Abaqus input  
        z3 = z2.split(',')
        zLen = len(z3)
        if z3[-1].strip() == "": zLen=zLen-1                                # for blanks after last comma
    f1.close()
    #
    # post -- input data checks and completion
    for elset in SecDic:
        # checks for materials defined with sections
        mat = SecDic[elset].Mat
        try:    Material = MatList[mat]
        except: raise NameError("ConFemInOut::DataInput: Material not defined",mat)
        Material.Used = True
        if SecDic[elset].bStiff != None:
            Material.matbStiff = SecDic[elset].bStiff
        if len(SecDic[elset].ElemTypes)>0: eltype = SecDic[elset].ElemTypes[0] 
        else:
            raise NameError("ConFemInOut::DataInput: seems to be no element for element set",elset)
        # assigns result types to element sets
        if eltype in _BeamsBernAll_:
            if bondl!=None:                          ResultTypes[elset] = ('longitudinal strain','curvature','normal force','moment','sig_bottom','sig_top') # 4 more bond items in DataOutStress
            elif isinstance( Material, RCBeam):      ResultTypes[elset] = ('longitudinal strain','curvature','normal force','moment','max reinf strain','min conc strain')
            elif isinstance( Material, Elastic):     ResultTypes[elset] = ('longitudinal strain','curvature','normal force','moment')
            elif isinstance( Material, Mises):       ResultTypes[elset] = ('longitudinal strain','curvature','normal force','moment')
            elif isinstance(Material, MisesBeam2D):  ResultTypes[elset] = ('longitudinal strain', 'curvature', 'normal force', 'moment')
            elif isinstance( Material, PolyBeam):    ResultTypes[elset] = ('longitudinal strain','curvature','normal force','moment')
        elif eltype in _BeamsTimoAll_:
            if isinstance( Material, RCBeam):        ResultTypes[elset]=('longitudinal strain','curvature','shear def','normal force','moment','shear force','max reinf strain','min conc strain')
            elif isinstance( Material, Elastic):     ResultTypes[elset]=('longitudinal strain','curvature','shear def','normal force','moment','shear force')
            elif isinstance( Material, MisesBeam2D): ResultTypes[elset]=('longitudinal strain', 'curvature','shear def', 'normal force', 'moment','shear force')
        elif eltype in ['T1D2','T2D2','T3D2','T2D3','T3D3','T2D2I','T2D3I','T3D2I','T3D3I']:
            ResultTypes[elset]=('strain','stress')
            if isinstance( Material, Elastic):
                if Material.PhaseField:              ResultTypes[elset]=('strain','stress','damage')
        elif eltype in ['TAX2','TAX3','TAX2I','TAX3I']:
            ResultTypes[elset] = ('eps_xx','eps_zz','sig_xx', 'sig_zz')
        elif eltype=='S1D2':                         ResultTypes[elset]=('displ','stress')
        elif eltype in ['CPE4','CPE4R','CPE3','CPE6','CPS4','CPS4R','CPS3','CPS6']: ResultTypes[elset]=('strain','stress')
        elif eltype=='SB3':                          ResultTypes[elset]=('curv x','curv y','curv xy','mom x','mom y','mom xy')
        elif eltype=='SH4':                          ResultTypes[elset]=('n_x','n_y','n_xy','m_x','m_y','m_xy','q_x','q_y')
        elif eltype=='SH3':                          ResultTypes[elset]=('n_x','n_y','n_xy','m_x','m_y','m_xy','q_x','q_y')
        elif eltype=='C3D8':                         ResultTypes[elset]=('strain','stress')
        elif eltype=='CAX4':                         ResultTypes[elset] = ('strain', 'stress')
        elif eltype=='S2D6':                         ResultTypes[elset]=('displ','force')
        elif eltype in ["Bond2D2","Bond2D3","Bond3D2","Bond3D3"]: ResultTypes[elset]=('slip','bond stress') # ??? are never in input data
        else: raise NameError ("no result types for this element defined", eltype)
    #
    for el in ElemList:
        # assign material to elements
        elset = el.Set
        mat = SecDic[elset].Mat
        eltype = el.Type
        if eltype in ['T1D2','T2D2','T3D2','T2D3','T3D3','T2D2I','T2D3I','T3D2I','T3D3I']:
            if isinstance( MatList[mat], Mises): mat = mat+'1'              # becomes uniaxial Mises, mat --> element definition --> Elem.MatN: seems to be a label only
        if eltype in _BeamsBernAll_:                                             # mises not appropriate for B21,'B21E due to shear component
            if isinstance( MatList[mat], Mises): mat = mat+'2'              # becomes Bernoulli beam Mises, mat --> element definition --> Elem.MatN: seems to be a label only, shows up in written output
        el.MatN = mat 
        el.Material = MatList[mat]
        elLabel     = el.Label
        Material    = el.Material                                           # material set here  !!!
        NData       = Material.NData
        NStateVar   = Material.StateVar
        RandDat     = []
        # init storage for data and state variables
        if eltype in _BeamsAll_:
            BeamSec = SecDic[elset]
            ReinfName = BeamSec.ReinfName
            if el.dim == 10:                                                # Bernoulli beam - NData size of filed Data, DatP
                if   Material.Type in ['RCBEAM','MISES']:                           NData = 6
                elif Material.Type in ['ELASTIC']:                                  NData = 6
            else:                                                           # Timoshenko beam
                if   Material.Type in ['RCBEAM']:                                   NData = 8
                elif Material.Type in ['ELASTIC']:                                  NData = 6
#                else: raise NameError("ConFemInOut::DataInput: material type not implemented for Timoshenko beam", elset,Material.Type)
            el.IniBeam( BeamSec, BreiDic[ReinfName], NData, NStateVar, ff, Material) # defines also Data, StateVar
        else:
            if el.Type in ["T1D2","T2D2","T2D2I","T2D3","T2D3I","T3D2","T3D2I","T3D3","T3D3I"]:
                if   Material.Type in ['ELASTICLT','ELASTIC','ELASTIC_VISCO','ELASTICC1D']:      NData = 2
                elif Material.Type in ['ISODAMAGE','MISES','ELASTIC_PHASEFIELD']:   NData = 3
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for truss element", elset,Material.Type)
            elif eltype in ['TAX2','TAX3',"TAX2I","TAX3I","TAX3E"]:
                if Material.Type in ['ELASTIC','ELASTIC1DR','ELASTICC1D']:          NData = 4
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for axisym truss element", elset,Material.Type)
            elif eltype in ['CPS4','CPE4','CPS3','CPE3','CPS6','CPE6','CAX4']:
                if   Material.Type in ['ELASTIC']:                                  NData = 8
                elif Material.Type in ['MISES','ISODAMAGE','MICRODAMAGE']:          NData = 10 # 9 # gradient damage ???
                if Material.Type in ['ISODAMAGE'] and elLabel in RandomData:# prepare data for element initialization of random data
                    RandData = RandomData[elLabel]                          # RandomData is local dictionary and holds random values for each integration point
                    if len(RandData) != el.nIntL: raise NameError("ConFemInOut::DataInput: number of ip mismatch with random data",elset)
                    for r in RandData:
                        _,_,_,_,cc0,cc1,cc2,cc3 = Material.EvalPar_(Material.Emod, Material.nu,Material.fc,r,Material.beta, Material.ga,Material.a3)
                        RandDat += [[cc0,cc1,cc2,cc3,r]]                    # RandDat is used for initialization of element data, see el.InitData
            elif eltype in ['CPS4S','CPE4S','CPS3S','CPE3S','CPS6S','CPE6S']:   # strong discontinuity elements
                pass                                                        # cannot yet be addressed 
            elif eltype in ['C3D8']: 
                if   Material.Type in ['ELASTIC','ELASTICLT_SDA']:                  NData = 12
                elif Material.Type in ['MISES','ISODAMAGE','MICRODAMAGE']:          NData = 14  # gradient damage ???
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for C3D* element", elset,Material.Type)
#                NData = NData+2                        # ???
            elif eltype in ['SB3']: 
                if   Material.Type in ['ELASTIC','NLSLAB']:                         NData = 8   # storage for element data, extra 2 for shear forces computed from moment derivatives
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for SB3 element", elset,Material.Type)
            elif eltype in ['SH4','SH3']:
                el.ShellRCFlag = False
                if   Material.Type in ['ELASTIC','MISES']:                          NData = 6
                elif isinstance( Material, WraRCShell):
                    if eltype=='SH3': raise NameError('not yet implemented')
                    nre = len(SecDic[elset].Reinf) 
                    if nre>0: el.nIntL = el.nIntL+4*nre                     # for reinforcement integration points in base area
                    NData = NData+6
                    el.ShellRCFlag = True
                    if Material.StateVar==None: NStateVar = 5               # for mises reinforcement, TR included (-> 1)
                    else:                       NStateVar = max(5,NStateVar)# 1st for mises reinforcement, 2nd for bulk material -- both use different ip points, see SH4.__init__
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for SH* element", elset,Material.Type)
            elif eltype in ['S1D2']: 
                if   Material.Type in ['SPRING','ELASTIC']:                         NData = 2
                else: raise NameError("ConFemInOut::DataInput: material type not implemented for S1D2 element", elset,Material.Type)
            elif eltype in ['S2D6']: 
                raise NameError("ConFemInOut::DataInput: S2D6 still under construction", elset,Material.Type)
            # assign storage for data and state variables
            el.InitData(el.nIntL, NData, NStateVar, RandDat)
    # assigns list of elements to each Section
    for i, e in enumerate(ElemList):
        SecDic[e.Set].Elems += [i]
    # assigns globally defined amplitudes to steps
    for a in AmpDict:
        AmpFlag = True
        for s in StepList:
            for sa in s.AmpDict:
                if a == sa: 
                    AmpFlag = False
                    break
            if AmpFlag: 
                s.AmpDict[a] = AmpDict[a]
            else:       
                Echo(f"global amplitude {a:s} overridden by step amplitude", ff)
    # post ready
    if Restart: return MatList, StepList
    else:       return NodeList, ElemList, MatList, SecDic, NoLabToNoInd, StepList, ResultTypes, Header

def WriteElemData( f2, f7, Time, ElList, NodeList,NoIndToCMInd, MatList, MaxType, ResultTypes):
    f2.write('%s%8.4f\n'%('Time ',Time))
    MaxOut = False
    if f7!=None and len(MaxType)>0:                                 # to extract maximum and minimum values of element data
        MaxOut = True 
        f7.write('%s%8.4f\n'%('Time ',Time))
        MaxVal, MinVal = {}, {}
        for i in MaxType: MaxVal[i], MinVal[i] = [0, 0.,0.,0., 0.], [0, 0.,0.,0., 0.] # element number, coordinates, value
    for Elem in ElList:
        if Elem.Type=='SB3':
            for j in Elem.Inzi:
                no = NodeList[ NoIndToCMInd[j] ]
                no.mx, no.my, no.mxy, no.C = 0,0,0,0                         # seems to be not needed
    for Elem in ElList:
        if not Elem.Active: continue
        xN, yN, zN = [], [], []
        f2.write('%s%7i%8s%40s%40s%20s'%('El',Elem.Label,Elem.Type,Elem.Set,Elem.MatN,MatList[Elem.MatN].Type))
        for j in Elem.Inzi:
            no = NodeList[NoIndToCMInd[j]]
            f2.write('%7i'%(no.Label))
            xN += [no.XCo]                                                  # nodal coordinates
            yN += [no.YCo]
            zN += [no.ZCo]
        if Elem.Type not in ['SH4','SH3']: f2.write('\n')
        if Elem.Type=='SB3': Elem.ShearForces( NodeList, NoIndToCMInd, True)# post processing kirchhoff slab for shear forces
        #
        if Elem.Type=='SH4' or Elem.Type=='SH3':                            # post processing continuum based shell
            IntF = {}                                                       # internal forces dictionary initialization
            if Elem.ShellRCFlag: nRe, marker = Elem.Geom.shape[0]-2, 'C'    # number of reinforcement layers
            else:           nRe, marker = 0, 'U'
            if   Elem.nInt==2: nC = 4                                       # number of concrete layers
            elif Elem.nInt==5: nC = 5
            else:              nC = 0
            f2.write('     %7i%7i\n'%(nC,nRe))
            Lis = Elem.Lists1()                                             # indices for start integration points in base area
            if   Elem.Type=='SH4': Corr, CorM =0.5, 0.5             # to compensate for local iso-parametric coordinate t in range [-1..1] --> 2
            elif Elem.Type=='SH3': Corr, CorM =3.0, 0.5             # " + Corr -->  0.5 * 8 *3/4 = 3, see ConFemBasic::SampleWeights SH3
            for j in Lis:                                                   # loop over base area
                Lis2, Lis3 = Elem.Lists2( nRe, j)                           # integration point indices along cross section height specific for current start point in base
                                                                            # Lis2 all, Lis3 reinforcement layers only in case of RC
                r = SamplePoints[Elem.IntT,Elem.nInt-1,j][0]                # intT: integration type, nInt: integration order
                s = SamplePoints[Elem.IntT,Elem.nInt-1,j][1]
                xI = dot( Elem.FormX(r,s,0), xN)
                yI = dot( Elem.FormX(r,s,0), yN)
                zI = dot( Elem.FormX(r,s,0), zN)
                aa = dot( Elem.FormX(r,s,0), Elem.a)                        # interpolated shell thickness from node thickness
                for k in ResultTypes[Elem.Set]: IntF[k] = 0. 
                for k in Lis2:                                              # loop over cross section height indices
                    jj = j+k                                                # j is base index
                    if (k in Lis3):                                         # reinforcement only
                        t  =         SamplePointsRCShell[Elem.Set,Elem.IntT,Elem.nInt-1,jj][2]
                        ff = Corr*aa*SampleWeightRCShell[Elem.Set,Elem.IntT,Elem.nInt-1,jj]  
                    else:           
                        t  =                SamplePoints[         Elem.IntT,Elem.nInt-1,jj][2]
                        ff = Corr*aa*       SampleWeight[         Elem.IntT,Elem.nInt-1,jj]  
                    ppp = Elem.Data[jj]                                     # content ruled through the material methods returning Data for Elem.dim==21
                    # isodamage: sig[0],sig[1],sig[2],sig[3],sig[4],sig[5], Eps__[0],Eps__[1],Eps__[5], sig[0],sig[1],sig[5], D, svs
                    IntF['n_x']  = IntF['n_x']  + ff*ppp[0]  
                    IntF['n_y']  = IntF['n_y']  + ff*ppp[1]
                    IntF['q_y']  = IntF['q_y']  + ff*ppp[3]
                    IntF['q_x']  = IntF['q_x']  + ff*ppp[4]
                    IntF['n_xy'] = IntF['n_xy'] + ff*ppp[5]
                    IntF['m_x']  = IntF['m_x']  - CorM*aa*t*ff*ppp[0]
                    IntF['m_y']  = IntF['m_y']  - CorM*aa*t*ff*ppp[1]
                    IntF['m_xy'] = IntF['m_xy'] - CorM*aa*t*ff*ppp[5]
                    la_,v_ = eigh( [[ppp[0],ppp[5],ppp[4]] , [ppp[5],ppp[1],ppp[3]] , [ppp[4],ppp[3],ppp[2]]])  # principal values/eigenvalues -- in ascending order -- and eigenvectors of stress tensor
#                    D =   ppp[-2]                                           # for damage of isodam
#                    kap = ppp[-1]                                           # for equivalent damage strain of isodam
                    # case R: reinforcement data          
                    if (k in Lis3):                  #    R                 # Lis3 reinforcement layers only
                        f2.write('%10.4f%10.4f%10.4f%7.3f R%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e                          '\
                                       %(xI,yI,zI,t,ppp[0],ppp[1],ppp[2],ppp[3],ppp[4],ppp[5])) 
                        f2.write('%9.2f  %7.3f%7.3f%7.3f%9.2f  %7.3f%7.3f%7.3f%9.2f  %7.3f%7.3f%7.3f'%(la_[0],v_[0,0],v_[1,0],v_[2,0],la_[1],v_[0,1],v_[1,1],v_[2,1],la_[2],v_[0,2],v_[1,2],v_[2,2]))
                        f2.write('\n')
                    # case marker == C but not in Lis3: concrete Data
                    elif marker=='C':               #     C
                        f2.write('%10.4f%10.4f%10.4f%7.3f C%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e'\
                                       %(xI,yI,zI,t,ppp[0+6],ppp[1+6],ppp[2+6],ppp[3+6],ppp[4+6],ppp[5+6],ppp[6+6],ppp[7+6]))
                        # eigenvalues + eigenvectors in ascending order
                        f2.write('%9.2f  %7.3f%7.3f%7.3f%9.2f  %7.3f%7.3f%7.3f%9.2f  %7.3f%7.3f%7.3f'%(la_[0],v_[0,0],v_[1,0],v_[2,0],la_[1],v_[0,1],v_[1,1],v_[2,1],la_[2],v_[0,2],v_[1,2],v_[2,2]))
                        f2.write('\n')
                _, _, _, _, _, vv = Elem.Basics( r, s, 0.)                  # vv: orientation of local coordinate system 1st aligned to first edge, 2nd column perp in-plane, 3rd column director, see also SH4:Basics
                # case Z                            Z
                f2.write('%10.4f%10.4f%10.4f        Z%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e   %8.4f%8.4f%8.4f %8.4f%8.4f%8.4f\n'\
                                       %(xI,yI,zI,IntF['n_x'],IntF['n_y'],IntF['n_xy'],IntF['m_x'],IntF['m_y'],IntF['m_xy'],IntF['q_x'],IntF['q_y'],vv[0,0],vv[1,0],vv[2,0],vv[0,1],vv[1,1],vv[2,1])) # integrated internal forces
                if MaxOut:
                    for k in MaxType:
                        if IntF[k]>MaxVal[k][4]: MaxVal[k][0], MaxVal[k][1],MaxVal[k][2],MaxVal[k][3], MaxVal[k][4] = Elem.Label, xI,yI,zI, IntF[k]                
                        if IntF[k]<MinVal[k][4]: MinVal[k][0], MinVal[k][1],MinVal[k][2],MinVal[k][3], MinVal[k][4] = Elem.Label, xI,yI,zI, IntF[k]                
        # 3D
        elif Elem.Type in ['C3D8','C3D8S','T3D2','T3D3','T3D2I','T3D3I',"Bond3D2","Bond3D3"]:
            for j in range(Elem.nIntL):                                     # write sample point coordinates and element data into file
                r = SamplePoints[Elem.IntT,Elem.nInt-1,j][0]
                s = SamplePoints[Elem.IntT,Elem.nInt-1,j][1]
                t = SamplePoints[Elem.IntT,Elem.nInt-1,j][2]
                xI = dot( Elem.FormX(r,s,t), xN)
                yI = dot( Elem.FormX(r,s,t), yN)
                zI = dot( Elem.FormX(r,s,t), zN)
                f2.write('%9.4f%9.4f%9.4f'%(xI,yI,zI))
                Data_ = Elem.Data[j]
                for D_ in Data_: f2.write('%13.4e'%(D_))
                f2.write('\n')
            if Elem.Type in ['C3D8S']:                                      # discontinuity data for SDA element - not yet finished
                f2.write('%9.4f%9.4f%9.4f           C1'%(Elem.CrXCo,Elem.CrYCo,Elem.CrZCo)) 
                for i in range(Elem.n_vertices): f2.write('%2i%11.4e%7.4f%11.4e%7.4f%11.4e%7.4f'%(i,Elem.wwL[i][0],Elem.wwT[i][0],Elem.wwL[i][1],Elem.wwT[i][1],Elem.wwL[i][2],Elem.wwT[i][2]))          
                f2.write('%12.6f%12.6f%12.6f\n'%(Elem.CrNx,Elem.CrNy,Elem.CrNz))           
        # else than shells and 3D
        else:
            for j in range(Elem.nIntL):                                     # write sample point coordinates and element data into file
                r = SamplePoints[Elem.IntT,Elem.nInt-1,j][0]
                s = SamplePoints[Elem.IntT,Elem.nInt-1,j][1]
                t = SamplePoints[Elem.IntT,Elem.nInt-1,j][2]
                xI = dot( Elem.FormX(r,s,t), xN)
                yI = dot( Elem.FormX(r,s,t), yN)
                f2.write('%9.4f%9.4f'%(xI,yI))
                for k in range(Elem.Data.shape[1]): f2.write('%13.4e'%(Elem.Data[j,k]))
                f2.write('\n')
                if Elem.Type in ['CPS4S','CPE4S','CPS3S','CPE3S','CPS6S','CP63S']:  # discontinuity data for SDA elements
                    f2.write('%9.4f%9.4f  C1%9.4f%9.4f%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e\n'%(Elem.CrXCo,Elem.CrYCo,Elem.CrN[0,0],Elem.CrN[0,1],Elem.wwL[0,0,0],Elem.wwL[0,0,1],Elem.wwL[0,1,0],Elem.wwL[0,1,1],Elem.wwT[0,0,0],Elem.wwT[0,0,1],Elem.wwT[0,1,0],Elem.wwT[0,1,1]))          
                    f2.write('%9.4f%9.4f  C2%9.4f%9.4f%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e%13.4e\n'%(Elem.CrXCo,Elem.CrYCo,Elem.CrN[1,0],Elem.CrN[1,1],Elem.wwL[1,0,0],Elem.wwL[1,0,1],Elem.wwL[1,1,0],Elem.wwL[1,1,1],Elem.wwT[1,0,0],Elem.wwT[1,0,1],Elem.wwT[1,1,0],Elem.wwT[1,1,1]))          
    #
    if MaxOut: 
        for k in MaxType:
            f7.write('max%6s%8i%10.4f%10.4f%10.4f%13.4e\n'%(k,MaxVal[k][0],MaxVal[k][1],MaxVal[k][2],MaxVal[k][3],MaxVal[k][4]))
            f7.write('min%6s%8i%10.4f%10.4f%10.4f%13.4e\n'%(k,MinVal[k][0],MinVal[k][1],MinVal[k][2],MinVal[k][3],MinVal[k][4]))
    return 0

def WriteNodalData( f3, Time, NodeList, VecU, VecB):
    f3.write('%8.4f\n'%(Time))
    for i in range(len(NodeList)):
        Node = NodeList[i]
        iS = Node.GlobDofStart
        if len(Node.DofT)>0:
            f3.write('%5i%9.4f%9.4f%9.4f'%(Node.Label,Node.XCo,Node.YCo,Node.ZCo))
            for j in range(len(Node.DofT)): f3.write('%12.4e'%(VecU[iS+j]))                            # DofT might be 5 or 6 for SH$, the latter for presumably folded edges / unfavorable director
            if Node.c>0: f3.write('  %12.4e%12.4e%12.4e  '%(Node.mx/Node.c,Node.my/Node.c,Node.mxy/Node.c)) # only for element type SB3
            for j in range(len(Node.DofT)): f3.write('%12.4e'%(VecB[iS+j]))
            f3.write('\n')
    return 0

def FinishAllStuff(PloF, DirName, FilName, Name, ResType, VTK):
    if PloF:
        from ConFemPost import ConFemPost
        MaPlLib = True
        StepCounter = 0
        ConFemPost_ = ConFemPost()
        ConFemPost_.Run( DirName,FilName, Name,StepCounter, MaPlLib, VTK, 1, [])
        show()
    else:
        import hashlib
        mmm = hashlib.md5()
        fp= open( Name+"."+ResType+".txt", "r")
        while True:
            data= fp.read(65536)
            if not data: break
            mmm.update(data.encode())
        fp.close()
        RC = mmm.hexdigest()
        return RC

def LogResiduals(LogName, counter, i, NodeList, VecR, VecU):
    if not path.exists(LogName):
        makedirs(LogName)
    ff = open( LogName+"/log-"+str(counter)+"-"+str(i)+".txt", 'w')              #
    f1 = open( LogName+"/logD-"+str(counter)+"-"+str(i)+".txt", 'w')              #
    for Node in NodeList:
        iS = Node.GlobDofStart
        if len(Node.DofT)>0:
            ff.write('%5i%12.4f%12.4f%12.4f%4i'%(Node.Label,Node.XCo,Node.YCo,Node.ZCo,len(Node.DofT)))
            f1.write('%5i%12.4f%12.4f%12.4f%4i'%(Node.Label,Node.XCo,Node.YCo,Node.ZCo,len(Node.DofT)))
            xxx = 0.
            for j in range(len(Node.DofT)):
                xxx = xxx +VecR[iS+j]**2
                ff.write('%12.4e'%(VecR[iS+j]))
                f1.write('%12.4e'%(VecU[iS+j]))
            ff.write('   %12.4e\n'%(sqrt(xxx)))
            f1.write('   %12.4e\n'%(sqrt(xxx)))
    ff.close()
    f1.close()
    return

#def PlotResiduals( LogName, Name, counter, i, scale):
#    def Stuff( Pl, Name, scale, scale2):
#        ff, xx, yy, rr = open( Name, 'r'), [], [], 0.
#        z1 = ff.readline()                                  # 1st input line
#        z2 = z1.split()
#        while z1!="":
#            x, y, rx, ry = float(z2[1]),float(z2[2]), float(z2[5]), float(z2[6])
##            x, y, rx, ry, ax, ay = float(z2[1]),float(z2[2]), float(z2[4]), float(z2[5]), float(z2[7]), float(z2[8])
#            rr = rr + rx**2 + ry**2
#            xx += [x]
#            yy += [y]
#            Pl.plot([x,x+scale*rx],[y,y+scale*ry],color='red')
##            Pl.plot([x,x+scale*ry],[y,y+scale*rx],color='red')
#    #        Pl.plot([x,x+scale2*ax],[y,y+scale2*ay],color='blue')
#            z1 = ff.readline()                                  # 1st input line
#            z2 = z1.split()
#        ff.close()
#        print(sqrt(rr))
#        Pl.plot(xx,yy,'x')
#        Pl.set_aspect('equal')
#        Pl.grid()
#        return
#    Name = LogName+"/log-"+str(counter)+"-"+str(i)+".txt"
#    Pl1 = plt.figure().add_subplot(111,title='Resid '+str(counter)+'/'+str(i))
#    Stuff( Pl1, Name, scale, 1.0e1) # 1e1
##    Stuff( Pl1, Name, 2.0e+3, 1.0e3) # 1e1
##    Name = LogName+"/logD-"+str(counter)+"-"+str(i)+".txt"
##    Pl2 = plt.figure().add_subplot(111,title='Displ '+str(counter)+'/'+str(i))
##    Stuff( Pl2, Name, 1.0e2, 1.0e1)
#    show()
#    return

#def WriteResiduals( LogName, counter, i, N):
#    class Residual():                                  # 1D Truss 2 nodes
#        def __init__(self, Label, X, Y, Z, Res):
#            self.Label = Label
#            self.X   = X
#            self.Y   = Y
#            self.Z   = Z
#            self.Res = Res
#    Name = LogName+"/log-"+str(counter)+"-"+str(i)+".txt"
#    Nam1 = LogName+"/logWr-"+str(counter)+"-"+str(i)+".txt"
#    ResList = []
#    ff  = open( Name, 'r')
#    f1  = open( Nam1, 'w')
#    z1 = ff.readline()
#    z2 = z1.split()
#    while z1!="":
#        nL, x, y, z, _, rr  = int(z2[0]), float(z2[1]), float(z2[2]), float(z2[3]), int(z2[4]), float(z2[-1])
#        ResList += [Residual( nL, x, y, z, rr)]
#        f1.write("%5i ; %8.2f ; %8.2f ; %8.2f ; %16.6e\n"%(nL,x,y,z,rr))
#        z1 = ff.readline()
#        z2 = z1.split()
#    ff.close()
#    f1.close()
#    n = len(ResList)
#    ResList.sort(key=lambda t: t.Res)
#    print(Name)
#    for i in range(N):
#        r = ResList[n-1-i]
#        print(i+1, r.Label, r.Res)
        
def DataOut(fileName, NodeList, uu, ff, rL, Time, WriteType ):
    fp = open(fileName,WriteType)
    fp.write('Time, %8.4f\n'%Time)
    for no in NodeList:
        fp.write('%7i, %8.4f, %8.4f, %8.4f,   %2i,  '%(no.Label,no.XCo,no.YCo,no.ZCo,len(no.DofT)))
        s = no.GlobDofStart
        for j in range(len(no.DofT)): fp.write('%12.4e,'%(uu[s+j]))         # displacements
        fp.write('    ')
        for j in range(len(no.DofT)): fp.write('%12.4e,'%(ff[s+j]))         # residual / reaction force
        fp.write('    ')
        for j in range(len(no.DofT)): fp.write('%12.4e,'%(rL[s+j]))         # residual / reaction force
        fp.write('\n')
    fp.flush()
    fp.close()
  
def DataOutStress(fileName, ElemList, NodeList, NoIndToCMInd, Time, WriteType, ff):
    fp = open(fileName, WriteType)
    fp.write('Time, %8.4f\n'%Time)
    counter = 0
    for i in ElemList :
        if i.Active:
            Type  = i.Type
            Label = i.Label
            ElSet = i.Set
            Mat   = i.MatN
            xyN   = i.NodalCoordinates( NodeList, NoIndToCMInd)
            if i.Type in ['SH4','SH3']:
                offset, nRe = 0, 0                                      # currently not clear for what offset is, does this reinforcement work ???  uhc 220621
                if i.ShellRCFlag: nRe = i.Geom.shape[0]-2               # number of reinforcement layers
                else:             nRe = 0
                Lis = i.Lists1()                                        # indices for first integration points in base area
                for j in Lis:
                    r = SamplePoints[i.IntT,i.nInt-1,j][0]              # intT: integration type, nInt: integration order
                    s = SamplePoints[i.IntT,i.nInt-1,j][1]
                    aa = dot(i.FormX(r, s, 0), i.a)                     # interpolated shell thickness from node thicknesses
                    XX  = i.FormX_( r, s, 0)
                    xyP = dot(XX,array(xyN))                            # global integration point coordinates
                    fp.write('%6i, %4s, %10s, %6s, %3i, '%(Label,ElSet,Mat,Type,j))
                    for x in xyP:               fp.write('%10.4f,'%(x))         # write nodal coordinates
                    nx, ny, nxy, qx, qy, mx, my, mxy, _ = i.StressIntegration( j, offset, nRe, None, ff) # None --> uses stresses from data to derive internal forces --> StressIntegration is method shellss
                    fp.write('%14.6e,%14.6e,%14.6e,%14.6e,%14.6e,%14.6e,%14.6e,%14.6e,%14.6e'%(nx, ny, nxy, mx, my, mxy, qx, qy, aa))
                    fp.write('\n')
            else:
                if i.Type == 'SB3':
                    i.ShearForces( NodeList,NoIndToCMInd, False)        # post processing for shear forces
                for k in range(i.nIntL):
                    Sig   = i.Data[k]                                   # consider: material type assigned to integration might be different for one element with side effect for actual length of Data
                    r = SamplePoints[i.IntT,i.nInt-1,k][0]
                    s = SamplePoints[i.IntT,i.nInt-1,k][1]
                    t = SamplePoints[i.IntT,i.nInt-1,k][2]
                    XX  = i.FormX_( r, s, t)
                    xyP = dot(XX,array(xyN))                            # global integration point coordinates
                    fp.write('%6i, %4s, %10s, %6s, %3i, '%(Label,ElSet,Mat,Type,k))
                    for x in xyP: fp.write('%10.4f,'%(x))               # write integration point coordinates
                    # write element stresses, strains
                    if i.Type in ["T3D2I","T3D3I"]:                     # for truss elements, required for bond data in VTK
                        fp.write('%14.6e,%14.6e,'%(Sig.tolist()[0],Sig.tolist()[1])) # looks like longituinal strain, stress
                    else:
                        for x in Sig.tolist():  fp.write('%14.6e,'%(x))
                    # append bond data for embedded elements - 4 items for 2D, 6 items for 3D: 2 bond, 4 lateral
                    if i.Type in ['T2D2I','T2D3I',"TAX2I","TAX3I",'T3D2I','T3D3I',"B23I","BAX23I","BAX23EI","B23EI","BAX21EI"] and i.BondLaw!=None:
                        bel  = ElemList[i.BondElInd[k]]                 # corresponding bond element
                        if counter!=bel.iT or k!=bel.ElInt: raise NameError("CaeFemInOut::DataOutStress: wrong assignment of embedded truss to bond element ")
                        for kk in range(len(bel.Data[0])): fp.write('%14.6e,'%(bel.Data[0,kk])) # bond element should have only one integration point
                    fp.write('\n')
                if i.Type in ['CPS4S','CPE4S','CPS3S','CPE3S','CPS6S','CP63S']:
                    fp.write('%6i, %4s, %10s, %6s, Cra, %10.4f,%10.4f,%12.6f,%12.6f,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e\n'%(Label,ElSet,Mat,Type,i.CrXCo,i.CrYCo,i.CrN[0,0],i.CrN[0,1],i.wwL[0,0,0],i.wwL[0,0,1],i.wwL[0,1,0],i.wwL[0,1,1],i.wwT[0,0,0],i.wwT[0,0,1],i.wwT[0,1,0],i.wwT[0,1,1]))
                    fp.write('%6i, %4s, %10s, %6s,Cra2, %10.4f,%10.4f,%12.6f,%12.6f,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e,%12.4e\n'%(Label,ElSet,Mat,Type,i.CrXCo,i.CrYCo,i.CrN[1,0],i.CrN[1,1],i.wwL[1,0,0],i.wwL[1,0,1],i.wwL[1,1,0],i.wwL[1,1,1],i.wwT[1,0,0],i.wwT[1,0,1],i.wwT[1,1,0],i.wwT[1,1,1]))
                if i.Type in ['C3D8S']:
                    ab = lstsq(i.wwIntPCoor, i.wwTn, rcond=None)[0]         # linear regression - should be four conditions (-> wwTn) and three parameters (-> coordinate coefficients)
                    vT = []
                    for c in i.xyzC: vT += [ dot( ab, c) ]
                    ab = lstsq(i.wwIntPCoor, i.wwLn, rcond=None)[0]         # linear regression
                    vL = []
                    for c in i.xyzC: vL += [ dot( ab, c) ]
                    for j, c in enumerate(i.xyzC):                          # loop over vertices of crack face
                        fp.write('%6i, %4s, %10s, %6s, Cra, %10.4f,%10.4f,%10.4f,%12.4e,%12.4e\n'%(Label,ElSet,Mat,Type,c[0],c[1],c[2],vL[j],vT[j]))
        counter+=1
    fp.flush()
    fp.close()

def SelOutIni( FilName, Type, OutList):                 # initialize files for single element / node output
    FileHandles = []
    for i in OutList:
        if len(i)>0:
            StrName = FilName+Type+str(i[0])+"_"+str(i[1])+".txt".strip()
            fx = open( StrName, 'w')
            FileHandles += [fx]
    return FileHandles
def SelOutWrite( Type, Outlist, FileHandles, Time, List, NoIndToCMInd, NoLabToNoInd, uu, step, cc): # output for single elements / nodes over all time steps 
    counter = 0
    for i in Outlist:
        if Type=="elem":
            ind = FindIndexByLabel(List, i[0])
            if ind==-1: raise NameError("CAEFem::Main: SelOut unknown element / node",i) 
            el = List[ind]
            if el.Active:
                fx = FileHandles[counter]
                fx.write("%14.6e%8i%3i"%(Time,el.Label,i[1]))
                for j in range(len(el.Data[i[1]])): fx.write("%14.5e"%(el.Data[i[1],j]))
                if el.Type in ['CPS4S','CPE4S','CPS3S','CPE3S']:
                    fx.write('  %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e'%(el.wwL[0,0,0],el.wwL[0,0,1],el.wwL[0,1,0],el.wwL[0,1,1],el.wwT[0,0,0],el.wwT[0,0,1],el.wwT[0,1,0],el.wwT[0,1,1]))          
                    fx.write('  %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e'%(el.wwL[1,0,0],el.wwL[1,0,1],el.wwL[1,1,0],el.wwL[1,1,1],el.wwT[1,0,0],el.wwT[1,0,1],el.wwT[1,1,0],el.wwT[1,1,1]))          
                fx.write('%3d%5d\n'%(step,cc))
        elif Type=="node":
            i1 = NoLabToNoInd[i[1]]
            i2 = NoIndToCMInd[i1]
            no = List[i2]
            s = no.GlobDofStart
            fx = FileHandles[counter]
            fx.write("%14.6e%8i"%(Time,no.Label))
            for j in range(len(no.DofT)): fx.write("%14.5e"%(uu[s+j]))
            fx.write("\n")
#            fx.flush()
        counter += 1
