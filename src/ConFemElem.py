# ConFemElem -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
"""Element library for ConFem"""

from numpy import zeros, array, sqrt, double, fabs, ma, dot, transpose, full, cross, shape
from numpy.linalg import norm, inv, det
from bisect import bisect_left #, bisect_right
from math import acos, cos, sin, pi
#import sys 

from ConFemBasics import FindIndexByLabel, SamplePoints, SampleWeight, SamplePointsRCShell, SampleWeightRCShell, ZeroD, Echo, EchoM
from ConFemMat import Mises

def pClosest0(list, p): 
    list.sort(key = lambda a: sqrt((a[0]-p[0])**2 + (a[1]-p[1])**2 + (a[2]-p[2])**2))
    LL = [list[1], list[0], list[2]]
    return LL
def pClosest1(list, p): 
    list.sort(key = lambda a: sqrt((a[0]-p[0])**2 + (a[1]-p[1])**2 + (a[2]-p[2])**2)) 
    LL = list[1]
    return LL

class Node(object):
    
    ActiveNodes = 0
  
    def __init__(self, Label, XCo, YCo, ZCo, XDir, YDir, ZDir):
        self.DofT = set([])
        self.GlobDofStart = 0
        self.c = 0
        self.Label = Label
        self.XCo = XCo
        self.YCo = YCo
        self.ZCo = ZCo
        self.XDi = XDir
        self.YDi = YDir
        self.ZDi = ZDir
        self.CMLabel_ = 0                                           # modified label from application or cuthill mckee
        self.NodeEl = []                                            # element indices which node is part of
        self.Type = ""                                              # C2D for nodes belonging to 2D continuum elements, C3D, T3De
        self.CMIndex = 0
        self.used = False
        # used for VTK -> DataOutVtk
#        self.BufSig = zeros((6),dtype=float)                        # buffer for nodal stress values used for VTK
#        self.BufDis = zeros((3),dtype=float)                        # buffer for nodal displacement
#        self.BufBond = 0.                                           # buffer for nodal bond values
#        self.BufCount = 0
#        self.BufBondCount = 0

class Section(object):
    def __init__(self, SetLabel, Mat):
        self.Set = SetLabel
        self.Mat = Mat
        self.Elems = []                                                     # list for elements of this section
#        self.Elems_= []
        self.ElemTypes = []
#        self.ElemTypes_= []

class SolidSection(Section):
    def __init__(self, SetLabel, Mat, nRebar, Val):
        Section.__init__(self,  SetLabel, Mat)
        self.nRebar = nRebar
        self.Val = Val
        self.bStiff = None

class BeamSectionRect(Section):
    Echo = []
    def __init__(self, Type,SetLabel,Mat,nRebar,bStiff, Width, Height, Reinf, ReinfName=None):
        Section.__init__(self,  SetLabel, Mat)
        self.Type = Type                                                    # Shape, Type
        self.Width = Width                                                  # Width
        self.Height = Height                                                # Height
        self.Reinf = Reinf                                                  # Reinforcement
        self.ReinfName = ReinfName
        self.nRebar = nRebar
        self.bStiff = bStiff
class BeamSectionPoly(Section):
    Echo = []
    def __init__(self, Type,SetLabel,Mat,nRebar,bStiff, Reinf, ReinfName=None):
        Section.__init__(self,  SetLabel, Mat)
        self.Type = Type                                                    # Shape, Type
        self.Reinf = Reinf                                                  # Reinforcement
        self.Poly = []                                                      # Description of cross section by polyline
        self.ReinfName = ReinfName
        self.nRebar = nRebar
        self.bStiff = bStiff
    def AddPoint(self, zz, bb):
        nn = len(self.Poly)
        if nn>0: 
            lz=self.Poly[nn-1][0]
            if lz>zz: raise NameError ("ConFemElements::AddPoint: wrong data order")
        self.Poly += [[zz,bb]]
class BeamSectionCircle(Section):
    Echo = []
    def __init__(self, Type,SetLabel,Mat,nRebar,bStiff, Diam, ReinfName=None):
        Section.__init__(self,  SetLabel, Mat)
        self.Type = Type                                                    # Shape, Type
        self.Diam= Diam                                                     # diameter
        self.Reinf = []
        self.ReinfName = ReinfName
        self.nRebar = nRebar
        self.bStiff = bStiff

class ShellSection(Section):
    def __init__(self, SetLabel, Mat, Height, Reinf):
        Section.__init__(self,  SetLabel, Mat)
        self.Height = Height                                                # Value
        self.Reinf = Reinf                                                  # Reinforcement
        self.bStiff = None

class Element(object):
    """ holds element parameters.
    """
    BarLength = 0.                                                  # some auxiliary value for control purposes
    BarMass = 0.                                                    # some auxiliary value for control purposes
    LastLabel = 0
    def __init__(self, TypeVal,InzList, nFieVal, IntTVal,nIntVal,nIntLVal, DofTVal, dimVal,NLGeomIVal,Label,elSet,SpaceDimVal, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd):
        self.Rot = False                                            # Flag for elements with local coordinate system
        self.RotM= False                                            # Flag for elements with local coordinate system for materials
        self.RotG= False                                            # Flag for elements with local coordinate system for geometric stiffness
        self.ZeroD = 1.e-9                                          # Smallest float for division by Zero
        self.Type = TypeVal                                         # element type
        self.Label= Label
        try: self.nNod = len(DofTVal)                                    # number of nodes
        except: pass
        self.nFie = nFieVal                                         # degrees of freedom for external loading
                                                                    # corresponds to largest load-DofT used by element and its dimension/structure of shape matrix N 
        self.IntT = IntTVal                                         # Integration type 0: 1D Gaussian, 1: 2D Gaussian, 2: 3D Gaussian, 3: 2D triangular, 4: shells
        self.nInt = nIntVal                                         # integration order (--> index for points and weights = nInt-1)
        self.nIntL= nIntLVal                                        # total number of integration points 
                                                                    # may be increased dynamically for some element/material type combinations, nInt is probably used for distinction of cases 
        self.nIntLi = nIntLVal                                      # initial value for number of integration points 
                                                                    # is not increased dynamically, used for distinction of bulk and reinforcement
        self.DofT = DofTVal                                         # tuple, type of dof for every node of this element
                                                                    # 1 -> u_x, 2->u_y, 3->u_z, 4->phi around x, 5->phi around y, 6->phi around z, 7->gradient field 
                                                                    # 4->phi around x, 5->phi around y, 6->phi around z (currently not used)
                                                                    # may be subject to change for FEM: 7, 8, 9 -> discontinuity
        self.dim  = dimVal                                          # material type index
        self.NLGeomI = NLGeomIVal                                   # Flag whether large deformations are already implemented for this element
        self.RegType = 0                                            # default value no regularization
        self.ElemUpdate = False                                     # default for update of element data
        self.NLGeomCase = 0                                         # determines how geometric stiffness is built from element contributions
        self.ShellRCFlag = False                                    # indicator for combination of reinforcement with bulk. Currently used for postprocessing of SH4
        self.CrBwS = 1.                                             # scaling factor for crack regularization
        self.Inzi = []
        for i in InzList:
            try: 
                self.Inzi += [ NoLabToNoInd[i] ]                    # indices of nodes belonging to this element, maybe later subject to change, e.g. for EFG
            except: 
                raise NameError("Element.__ini__:", self.Label, "node not found",i)
        self.SDim = SpaceDimVal                                     # dimension in space
        self.Label= Label                                           # element number in input
        self.Set  = elSet                                           # label of corresponding element set
        try:    self.MatN  = MatName.upper()                        # name of default material of element
        except: pass
        self.PlSt = None
        self.Active = True
        self.TensStiff = False                                       # flag for tension stiffening
        self.nRebarEl  = 1.                                          # is adressed for SOLID SECTION (embedded truss elements) and BEAM SECTION, CIRCLE (embedded beam elements)
        self.FlexBonded= False                                       # embedded truss and beam elements may have rigig bond with shared nodes
        if Label>Element.LastLabel: Element.LastLabel = Label
    
    def InitData(self, nIntLVal, NData, StateV, RanDat, xIp=[], yIp=[]):    # is called in post-phase of InOut::DataInput -- might be problem for elements which are created on the fly
        self.Data = zeros((nIntLVal, NData),dtype=float)                    # ??? ! actually used value may be subject to change in case of individual ip materials !
        self.DataP= zeros((nIntLVal, NData),dtype=float)                    # ??? ! currently only required for "mises" !
        self.StateVar = []                                                  # default no state variables
        self.RandomData = []
        self.StateUpdate = False
        if StateV!=None:
            self.StateVar = zeros((nIntLVal,StateV), dtype=float)
            self.StateVarN= zeros((nIntLVal,StateV), dtype=float)
            self.StateUpdate = True
        if len(RanDat)>0:
            nr = max([ len(r) for r in RanDat ])
            self.RandomData = zeros((nIntLVal,nr), dtype=float)
            for i, r in enumerate(RanDat):
                for j in range(len(r)):
                    self.RandomData[i,j] = r[j]
        if len(xIp)>0 and len(xIp)==len(yIp):
            self.xIp = xIp
            self.yIp = yIp

    def CrossSecGeom(self, zc1, zc2):                               # cross section area, statical moment, inertial moment for 2D cross section defined with polyline within z limits
        pList = self.CrSecPolyL
        n = len(pList)
        if n<2: raise NameError ("ConFemElements::CrossSecGeom: not enough data")
        AA, SS, JJ = 0, 0, 0
        On6, On12 = 1./6, 1./12
        for i in range(n-1):
            z1, b1 = max(zc1,pList[i][0]),   pList[i][1]
            z2, b2 = min(zc2,pList[i+1][0]), pList[i+1][1]
            if z1>z2: z1=z2                                         # Section with tension only
            AA = AA+0.5*(b1+b2)*(z2-z1)                             # cross section area
            S1, S2, S3 =(2*b2+b1), (b1-b2), (2*b1+b2)
            SS = SS + On6*(S1*z2**2+S2*z1*z2-S3*z1**2)              # cross section statical moment
            J1, J2, J3 = (3*b2+b1), (b1-b2), -(3*b1+b2) 
            JJ = JJ + On12*(J1*z2**3+J2*(z1*z2**2+z1**2*z2)+J3*z1**3) # cross section inertial moment
        return AA, SS, JJ, pList[n-1][0]-pList[0][0]
    def WidthOfPolyline(self, zz):                                  # width of 2D cross section defined by polyline depending on z, used by RCBeam.Sig
        pList = self.CrSecPolyL
        n = len(pList)
        for i in range(n-1):
            z1, b1 = pList[i][0], pList[i][1]
            z2, b2 = pList[i+1][0], pList[i+1][1]
#            if zz>=z1: # and (zz-z2)<=ZeroD:
            if zz <= (z2+1.e-6): # and (zz-z2)<=ZeroD:
                if (z2-z1)>ZeroD: bb = ((b2-b1)*zz+z2*b1-b2*z1)/(z2-z1)
                else:             bb = b1
                return bb
    def CrossSecGeomA(self, zc1, zc2, s1, s2):                      # matrix AZ, geometric tangential stiffness for 2D polyline cross section 
        pList = self.CrSecPolyL
        n = len(pList)
        AZ = zeros((2,2))
        for i in range(n-1):
            z1, b1, z2, b2 = pList[i][0], pList[i][1], pList[i+1][0], pList[i+1][1]
            if zc1>=z1:
                if (z2-z1)>ZeroD: 
                    AZ[0,0] = -1./6*s1*b2-1./3*s1*b1-1./3*s2*b2-1./6*s2*b1
                    AZ[1,0] = 1./6*(b2*zc1+zc2*b2+zc1*b1)*s2+1./6*(b2*zc1-zc2*b1+3.*zc1*b1)*s1
                else:             AZ[0,0] = 0; AZ[1,0] = 0
                break
        for i in range(n-2,0,-1):
            z1, b1, z2, b2 = pList[i][0], pList[i][1], pList[i+1][0], pList[i+1][1]
            if zc2<=z2:
                if (z2-z1)>ZeroD: 
                    AZ[0,1] = 1./6*s1*b2+1./3*s1*b1+1./3*s2*b2+1./6*s2*b1;
                    AZ[1,1] = 1./6*(b2*zc1-3.*zc2*b2-zc2*b1)*s2-(1./6*zc1*b1+1./6*zc2*b2+1./6*zc2*b1)*s1
                else:             AZ[0,1] = 0; AZ[1,1] = 0
                break
        return AZ
    def ElemDimData(self, NoList,NoIndToCMInd ):                            # called by ConFemBasics:AssignGlobalDof
        if self.nNod != len(self.Inzi):
            raise NameError ("element node number inconsistent",self.Label,self.nNod,self.Inzi)
        self.DofN = []                                               
        for i in self.DofT: self.DofN += [len(i)]                           # list of number of dofs per node ( -> len(DofT[i]) ) -- maybe subject to later change
        self.DofE = sum(self.DofN)                                          # number of dofs for whole element -- maybe subject to later change, presumably only for CPS4_XFEM 
        self.DofNini = self.DofN.copy() #copy(self.DofN)                                      # remember initial value as DofN may be subject to change, e.g. for sh4
        self.DofEini = sum(self.DofN) #copy(self.DofE)
        self.Nodisp = zeros((self.DofE), dtype=float) # uhc
        i_ = 0
        for i in self.Inzi:
            if len(NoIndToCMInd)==0: j_ = i 
            else:                    j_ = NoIndToCMInd[i]
            No = NoList[j_]
#            NoList[j_].DofT = NoList[j_].DofT.union(self.DofT[i_])          # carry dof types from element to nodes/NodeList belonging to this element
            No.DofT = No.DofT.union(self.DofT[i_])                          # carry dof types from element to nodes/NodeList belonging to this element
            i_ += 1
        self.DofI = full( (self.nNod,max(self.DofN)), -1, dtype=int)        # initialize table for global indices of dofs of every element, self.DofN is an array
    def NodalCoordinates(self, NodeList, NoIndToCMInd):
        xy = []
        for j in self.Inzi:                                                 # loop over nodes defining element geometry
            j_ = NoIndToCMInd[j]
            node = NodeList[j_]
            if   self.SDim==1: xy += [node.XCo]
            elif self.SDim==2: xy += [node.XCo,node.YCo]                    # nodal coordinates
            elif self.SDim==3: xy += [node.XCo,node.YCo,node.ZCo]
        return array(xy)
    def CrBScaleType(self):                                                 # determining scaling type and factor for crack band regularization, called by Ini2
        self.ScaleType, InPol  = -1, False                                  # InPol flag for lin/quad interpolation of scaling factors, False is default
        if self.Material.RType == 2:
            if self.Material.Type in ['ISODAMAGE',"MicroPl",'MICRODAMAGE']:
                pass
            elif self.Material.Type in ['RCSHELL'] and self.Material.ConcS in ['ISODAMAGE',"MicroPl",'MICRODAMAGE']:
                pass
            else:
                raise NameError("ConFemElem::CrBScaleType: regularizaton type 2 not implenmented for material type ",self.Material.Type,self.Material.ConcS )
            L = self.Lch_
            if L > self.Material.bw*self.Material.ElCharLenBound:           # Material.ElCharLenBound material-spcific typically 0.278
                self.ScaleType = 1                                          # type 1 of crack band scaling for this element
            else:
                self.ScaleType = 2                                          # type 2 of crack band scaling for this element
            self.CrBwS = self.Material.CrackBandScale(self.Label, L, self.ScaleType) # element specific scaling factor
            if InPol:                                                       # currently better not use
                if self.ScaleType == 1:
                    self.CrBwS = self.Material.LargeElReg[3]*L*L + self.Material.LargeElReg[4]*L +  self.Material.LargeElReg[5] # quadratic interpolation
                if self.ScaleType == 2:
                    self.CrBwS = self.Material.SmallElReg[3] * L + self.Material.SmallElReg[4]  # linear interpolation
        else:
            self.CrBwS=1.
    #
    def Ini3(self, NoList, NoIndToCMInd):
        pass
    # initialization for all 2D beam elements applied after Element.__init__ in Post section of ConFemInOut::DataInput and may perform redefinitions
    def IniBeam(self, BeamSec, ReinfD, NData, StateV, ff, Material):
        self.TensStiff = True                                       # flag for tension stiffening
        self.NLGeomCase = 1                                         # determines how geometric stiffness is built from element contributions, see ConFemBasics:Intforces
        nRe = len(ReinfD)                                           # ReinfD comes from *BEAM REINF input, required in case of POLYLINE due to distinguished input format
        if nRe>0: ReiP = ReinfD
        else:     nRe, ReiP = len(BeamSec.Reinf), BeamSec.Reinf     # comes from *BEAM SSECTION / RECT input
        #
        self.Geom = zeros( (2+nRe,6), dtype=double)                 # all geometry data for concrete and reinforcement, tension stiffening
        self.Geom[0,0] = 1.                                         # Jacobian for integration set in Ini2
        self.Geom[1,0] = 1                                          # for compatibility reasons
        if self.Type in ["B23I","B23EI","BAX21EI","BAX23","BAX23EI"]:
            if BeamSec.Type != "CIRCLE":
                raise NameError("ConFemElem::IniBeam: option SECTION=CIRCULAR allowed for element types B23I, B23EI, BAX21EI, BAX23EI only",self.Type)
        if BeamSec.bStiff == None:
            bStiff = 1.                                             # used in the following
        else:
            bStiff = BeamSec.bStiff                                 # to modify bending stiffness to consider variations of cross section shape
        if not self.FlexBonded:                                     #
            self.nRebarEl = BeamSec.nRebar                          # nrebar introduced also with TxDxI for flexible bonded elements
        #                                                             nrebar should include thickness of underlying continuum with this usage -- 1 in case of axisymmetry
        if BeamSec.Type=="POLYLINE":
            self.CrossSecType = "POLYLINE"
            self.CrSecPolyL = BeamSec.Poly                          # Poly holds geometry data
            AA, SS, JJ, _ = self.CrossSecGeom(self.CrSecPolyL[0][0],self.CrSecPolyL[len(self.CrSecPolyL)-1][0])# area, statical moment, inertial moment
            if not BeamSec.Set in BeamSec.Echo:
                Echo(f"{BeamSec.Set:s}, Area {AA:f}, 1st moment area {SS:e}, 2nd moment area {JJ:e} ", ff)
                ff.flush()
                BeamSec.Echo += [BeamSec.Set]
            self.Geom[1,1] = self.CrSecPolyL[0][0]                  # lower coordinate
            self.Geom[1,2] = self.CrSecPolyL[-1][0]                 # upper coordinate
            self.zLow      = self.CrSecPolyL[0][0]                  # found no other place for lower coordinate to make this compatible with other cross section types
            self.zUpp      = self.CrSecPolyL[-1][0]                 # found no other place for upper coordinate to make this compatible with other cross section types
        elif BeamSec.Type=="CIRCLE":                                # nRebar is applied in CreateBond
            self.CrossSecType = "CIRCLE"
            Diam = BeamSec.Diam
            AA = 0.25*Diam**2*pi
            SS = 0.
            JJ = bStiff * pi*(Diam**4)/64.                          # bStiff not part of public input -> docu, default 1
            self.Geom[1,2] = Diam
            self.zLow      = -0.5*Diam                              # found no other place for lower coordinate
            self.zUpp      =  0.5*Diam                              # found no other place for upper coordinate
        elif BeamSec.Type == "RECT":
            self.CrossSecType = "RECT"
            if self.Type in ["BAX23"]: Width = 1.                   # corresponds to unit length of circumference
            else:                      Width = BeamSec.Width
            AA = Width*BeamSec.Height
            SS = 0
            JJ = bStiff * Width*BeamSec.Height**3/12.
            self.Geom[1,1] = Width
            self.Geom[1,2] = BeamSec.Height
            self.zLow      = -0.5*BeamSec.Height                    # found no other place for lower coordinate
            self.zUpp      =  0.5*BeamSec.Height                    # found no other place for upper coordinate
        else:
            raise NameError("ConFemElem::IniBeam: *BEAM SECTION option SECTION must be one out of RECT, CIRCLE, POLYLINE")
        self.Geom[1,3] = AA                                         # might be modified by nRebar in TxDx::CreateBond
        self.Geom[1,4] = SS
        self.Geom[1,5] = JJ
        self.RType = []                                             # Type of reinforcement, RC ordinary rebar, TC textile or else
        for j in range(nRe):
            self.Geom[2+j,0] = ReiP[j][0]                           # reinforcement cross section
            self.Geom[2+j,1] = ReiP[j][1]                           # " lever arm
            self.Geom[2+j,2] = ReiP[j][2]                           # effective reinforcement ratio
            self.Geom[2+j,3] = ReiP[j][3]                           # tension stiffening parameter betat
            self.RType += [ReiP[j][4]]                              # RC / TRC
        #
        self.InitData( self.nIntL, NData, Material.StateVar, [])
        #
        # redefinition of StateVar only
        if Material.Type.upper() != "MISES":                        # adjust StateVar field for reinforcement
            if StateV!=None:
                self.StateVar = zeros((self.nIntL*(nRe+1),StateV), dtype=float)
                self.StateVarN= zeros((self.nIntL*(nRe+1),StateV), dtype=float)
                if Material.Type.upper() == "MISESBEAM":
                    for j in range(self.nIntL):
                        self.StateVar[j][0]  = Material.fY         # attention: position defined in __init__, may be subject to change from time to time
                        self.StateVar[j][6]  = Material.fY         # attention: position defined in __init__, may be subject to change from time to time
                        self.StateVar[j][4] = -999.         # "
                        self.StateVar[j][10] = 999         # "
        # redefinition of StateVar and introducing Datai instead of Data
        else:                                                       # mises beam
            if nRe!=0: raise NameError("ConFemElements::IniBeam: reinf not possible with mises beam")
            nI = Material.nI
            self.Datai, self.DataPi = {}, {}
            for d in range(self.nIntL):
                self.Datai[d]  = zeros((nI, NData), dtype=float)
                self.DataPi[d] = zeros((nI, NData), dtype=float)
            if StateV!=None:
                self.StateVar, self.StateVarN = {}, {}
                for d in range(self.nIntL):
                    self.StateVar[d]  = zeros((nI, StateV), dtype=float)
                    self.StateVarN[d] = zeros((nI, StateV), dtype=float)

class S1D2(Element):                                                        # 1D Spring 2 nodes
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, SolSecDic, StateV, NData, NoLabToNoInd, LenFactor):
#        Element.__init__(self,"S1D2",InzList, 1, 0,1,1, (set([1]),set([1])), 99,False,Label,SetLabel,1, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        Element.__init__(self,"S1D2",InzList, 1, 0,1,1, (set([1]),set([1])), 99,False,Label,SetLabel,1, None,None,None, NoList,NoLabToNoInd,[])
        self.LenFactor = LenFactor                                          # to account for length of T1D2-elements (bond!) connected by S1D2
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        self.DofI = zeros( (self.nNod,1), dtype=int)                        # indices of global dofs per node
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 0.5                                                # Jacobi determinant value
        self.Geom[1,0] = SecDict[self.Set].Val*self.LenFactor
        mat = SecDict[self.Set].Mat
        self.MatN = mat
        self.Material = MaList[mat]
        self.Lch_ = 0                                                       # dummy
        return []
    def FormX(self, r, s, t):
        X = array([ 0.5, 0.5])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5, 0.5],
                   [ 0.0, 0.0]])
        return X
    def FormN(self, r, s, t):
        N = array([[ -1., 1.], [ 0, 0]])
        return N
    def FormB(self, r, s, t, NLg):
        B = array([[ -1., 1.], [ 0, 0]])
        return B, 1, 0
    def FormT(self, r, s, t):                                               # interpolation on temperature
        T = array([ 0., 0.])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1
class S2D6(Element):                                                # 2D rotational Spring 2 nodes
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, SolSecDic, StateV, NData, NoLabToNoInd):
        Element.__init__(self,"S2D6",InzList, 3, 0,1,1, (set([1,2,6]),set([1,2,6])), 98,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        self.DofI = zeros( (self.nNod,3), dtype=int)                # indices of global dofs per node
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 0.5                                        # Jacobi determinant value
        self.Geom[1,0] = SolSecDic.Val
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        self.Lch_ = 0                                               # dummy
        return []
    def FormX(self, r, s, t):
        X = array([ 0.5, 0.5])
        return X
    def FormN(self, r, s, t):
        N = array([[ -1.,  0.,  0., 1., 0., 0.],
                   [  0., -1.,  0., 0., 1., 0.],
                   [  0.,  0., -1., 0., 0., 1.]])
        return N
    def FormB(self, r, s, t, NLg):
        B = array([[ -1.,  0.,  0., 1., 0., 0.],
                   [  0., -1.,  0., 0., 1., 0.],
                   [  0.,  0., -1., 0., 0., 1.]])
        return B, 1, 0
    def FormT(self, r, s, t):                                       # interpolation on temperature -- here dummy, actually not interpolated
        T = array([[ 0., 0., 0., 0., 0., 0.], [ 0., 0., 0., 0., 0., 0.], [ 0., 0., 0., 0., 0., 0.]])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

class T1D2(Element):                                                        # 1D Truss 2 nodes
    def __init__(self, Label, SetLabel, InzList, MatName,Material, NoList, SolSecDic, StateV, NData, RegType, NoLabToNoInd):
#        if RegType in [1,4]: Element.__init__(self,"T1D2",InzList, 1, 0,1,1, (set([1,7]),set([1,7])), 1,False,Label,SetLabel,1, MatName,StateV,NData, NoList,NoLabToNoInd,[])
#        else:                        Element.__init__(self,"T1D2",InzList, 1, 0,1,1, (set([1]),set([1])),     1,False,Label,SetLabel,1, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        Element.__init__(self,"T1D2",InzList, 1, 0,1,1, None, 1,False,Label,SetLabel,1, None,None,None, NoList,NoLabToNoInd,[])
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        L = sqrt( (NoList[i1].XCo-NoList[i0].XCo)**2 + (NoList[i1].YCo-NoList[i0].YCo)**2 )# element length
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 0.5*L                                              # Jacobi determinant value / length measure
        self.Geom[1,0] = SecDict[self.Set].Val                              # cross sectional area 
        self.Lch_ = L                                                       # length in undeformed configuration !
        mat = SecDict[self.Set].Mat
        if isinstance( MaList[mat], Mises): self.MatN = mat+'1'
        else:                               self.MatN = mat
        self.Material = MaList[mat]
        self.RegType  = self.Material.RType                                 # regularization, 0 no, 1 gradient, 2 crack band width, 3 SDA, 4 phase field
        if self.RegType in [1,4]: self.DofT = (set([1,7]),set([1,7]))       # tuple, type of dof for every node of this element: 1 -> u_x, 7->gradient field
        else:                     self.DofT = (set([1]),  set([1]))
        self.nNod = len(self.DofT)                                          # number of nodes
        if self.RegType in [1,4]: self.DofI = zeros( (self.nNod,2), dtype=int)
        else:                     self.DofI = zeros( (self.nNod,1), dtype=int)         # indices of global dofs per node
        self.CrBwS = 1.
        self.CrBScaleType()
        return []
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0.5*(1+r)],
                   [ 0.,        0.]])
        return X
    def FormN(self, r, s, t):
        N = array([[ 0.5*(1-r), 0.5*(1+r)],[ 0, 0]])
        return N
    def FormB(self, r, s, t, NLg):                                  # includes Jacobian
        L = 2.*self.Geom[0,0]
        if self.RegType in [1,4]: 
            B = array([[ -1./L, 0, 1./L, 0], [ 0, -1./L, 0, 1./L]])
            BN= array([[ -1./L, 0, 1./L, 0], [ 0, 0.5*(1-r), 0, 0.5*(1+r)]])
            return B, BN, 1, 0
        else:               
            B = array([[ -1./L, 1./L], [ 0, 0]])
            return B, 1, 0
    def FormT(self, r, s, t):
        if self.RegType in [1,4]: T = array([ 0.5*(1-r), 0, 0.5*(1+r), 0])
        else:                     T = array([ 0.5*(1-r), 0.5*(1+r)])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

class T2D2( Element ):                                                      # 2D Truss 2 nodes
    def __init__(self, Label, elset, InzList, MatName, NoList, SolSecDic, StateV, NData, Val, NoLabToNoInd):
        Element.__init__(self,"T2D2",InzList, 1, 0,1,1, (set([1,2]),set([1,2])), 1,True,Label,elset,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        self.Geom = zeros((2, 4), dtype=double)                             #
        self.Geom[1, 0] = SolSecDic[elset].Val * Val                        #
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        self.CoordRef = array( [NoList[i0].XCo, NoList[i0].YCo, NoList[i1].XCo, NoList[i1].YCo]) # Nodal coordinates 
        self.CoordDis = array( [NoList[i0].XCo, NoList[i0].YCo, NoList[i1].XCo, NoList[i1].YCo]) # Nodal coordinates displaced
        dx = self.CoordRef[2]-self.CoordRef[0]
        dy = self.CoordRef[3]-self.CoordRef[1]
        L = sqrt(dx**2 + dy**2)
        self.Geom[0,0] = 0.5*L                                              # Jacobian
        self.Geom[1,1] = dx
        self.Geom[1,2] = dy
        self.Geom[1,3] = L*L
        self.Lch_ = L                                                       # length in undeformed configuration !
        self.LL = L
        self.Lch= L
        self.sinA = (NoList[i1].YCo-NoList[i0].YCo)/L                       # not used for T2D2 but for derived classes
        self.cosA = (NoList[i1].XCo-NoList[i0].XCo)/L
        Element.BarLength += L
        return []
    def UpdateCoord(self, dis, ddis ):
        self.CoordDis[0] = self.CoordRef[0]+dis[0]
        self.CoordDis[1] = self.CoordRef[1]+dis[1]
        self.CoordDis[2] = self.CoordRef[2]+dis[2]
        self.CoordDis[3] = self.CoordRef[3]+dis[3]
        self.Geom[1,1] = self.CoordDis[2]-self.CoordDis[0]                  # updated dx used for FormB
        self.Geom[1,2] = self.CoordDis[3]-self.CoordDis[1]                  # updated dy
        return
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((self.DofE, self.DofE), dtype=float)
        val = sig[0]/self.Geom[1,3]                                         # Geom[1,3] undeformed length, see Basics.tex 17.4.2.3
        GeomK[0,0] = val                                                    #     val    0      -val    0
        GeomK[0,2] = -val                                                   #     0      val      0     -val
        GeomK[1,1] = val                                                    #    -val    0       val    0
        GeomK[1,3] = -val                                                   #     0      -val     0     val
        GeomK[2,0] = -val
        GeomK[2,2] = val
        GeomK[3,1] = -val
        GeomK[3,3] = val
        return GeomK
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X
    def FormN(self, r, s, t):
        N = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return N
    def FormB(self, r, s, t, NLg):
        cc = self.Geom[1,1]/self.Geom[1,3]                                  # continuously updated
        ss = self.Geom[1,2]/self.Geom[1,3]                                  #
        B = array([[ -cc, -ss, cc, ss], [ 0, 0, 0, 0]])                     # includes strain-Jacobian in cc, ss by Geom[1,3]
        return B, 1, 0                                                      # integration-Jacobian in Geom[0,0] in IntForces
    def FormT(self, r, s, t):
        T = array([ 0.5*(1-r), 0.5*(1+r), 0, 0])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1
# 2D truss element  with additional longitudinal degree of freedom -- seems to be only used to derive T2D3I
class T2D3(T2D2):
    # Ini2 from T2D2
    def __init__(self, Label, InzList, AA, NoList, SolSecDic, elset, MatName, NoLabToNoInd, Material, NoIndToCMInd):
        Element.__init__(self,"T2D3",InzList, 2, 0,2,2, (set([1,2]),set([1,2]),set([1])), 1,False,Label,elset,   2, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,[]) # integration order 2 required to have at least 2 bond elements per truss element, otherwise the latter might become unstable in case of embedding
        self.Geom = zeros((2, 4), dtype=double)  #
        self.Geom[1, 0] = SolSecDic[elset].Val * AA  #
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r), 0.])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0.,        0.5*(1+r), 0.,        0., 0.],
                   [ 0.,        0.5*(1-r), 0.,        0.5*(1+r), 0., 0.]])          # last two entries for internal node with coordinates [0,0]
        return X
    def FormN(self, r, s, t):
        cc = self.cosA
        ss = self.sinA
        c2 = self.cosA**2
        s2 = self.sinA**2
        cs = self.cosA*self.sinA
        r2 = r**2
        N  = 0.5*array([[(-r+r2)*c2+(1.-r)*s2, (-r+r2)*cs-(1.-r)*cs, (r+r2)*c2+(1.+r)*s2, (r+r2)*cs-(1.+r)*cs, 2.*cc*(1.-r2)*(cc+ss)],
                        [(-r+r2)*cs-(1.-r)*cs, (-r+r2)*s2+(1.-r)*c2, (r+r2)*cs-(1.+r)*cs, (r+r2)*s2+(1.+r)*c2, 2.*ss*(1.-r2)*(cc+ss)]])
                        # last factor for compensation that this dof is referred to local system, relevant, e.g., for mass matrix
        return N
    def FormB(self, r, s, t, NLGeom):
        L = self.LL
        Li = 1./L                                                           # strain-jacobian, 2 is in following B
        cc = self.cosA
        ss = self.sinA
        r2 = 2.*r
        B = array([[(-1.+r2)*cc*Li, (-1.+r2)*ss*Li, (1.+r2)*cc*Li, (1.+r2)*ss*Li, -2.*r2*Li],[0,0,0,0,0]])
#        return B, L/2., 0                                                   # returns integration jacobian -- Geom[0,0 shold be set set to 1 in __init__ or Ini2
        return B, 1., 0                                                     #
    def FormT(self, r, s, t):
        T = array([ 0.5*(1-r), 0.5*(1+r), 0, 0, 0])                         # remains to be controlled / implemented
        return T

class TAX2( T2D2 ):                                                         # axisymmetric Truss 2 nodes
    # Ini2 from T2D2
    def __init__(self, Label, elset, InzList, NoList, SolSecDic, Val, NoLabToNoInd):
#       Element.__init__(self,"TAX2",InzList, 1, 0,1,1, (set([1,2]),set([1,2])), 5,False,Label,SetLabel,2, None,None,None, NoList,NoLabToNoInd,[])
        Element.__init__(self,"TAX2",InzList, 1, 0,2,2, (set([1,2]),set([1,2])), 5,False,Label,elset,2, None,None,None, NoList,NoLabToNoInd,[])
        self.Geom = zeros((2, 4), dtype=double)  #
        self.Geom[1, 0] = SolSecDic[elset].Val * Val  #
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT,self.nInt-1,j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
    def FormB(self, r, s, t, NLg):
        cc = self.Geom[1, 1] / self.Geom[1, 3]                              # presumably includes strain jacobian witj Geom[1,3]
        ss = self.Geom[1, 2] / self.Geom[1, 3]
        x  = 0.5*( (1-r)*self.X0 + (1+r)*self.X1 )                          # 0.5 for interpolation of x-coordinate
        x_ = 0.5/x                                                          # 0.5 for interpolation of displacement
        B = array([[-cc, -ss, cc, ss], [x_*(1-r), 0, x_*(1+r), 0]])         # includes strain-Jacobian in cc, ss by Geom[1, 3], Jacobian for 2nd row already included in x_
        return B, 1, 0                                                      # integration-Jacobian in Geom[0,0] in IntForces
class TAX3( T2D3 ):                                                         # axisymmetric Truss 2 + 1 nodes
    # Ini2 from T2D3 --> T2D2
    def __init__(self, Label, elset, InzList, NoList, SolSecDic, Val, NoLabToNoInd):
        Element.__init__(self,"TAX3",InzList, 1, 0,2,2, (set([1,2]),set([1,2]),set([1])), 5,False,Label,elset,2, None,None,None, NoList,NoLabToNoInd,[])
        self.Geom = zeros((2, 4), dtype=double)  #
        self.Geom[1, 0] = SolSecDic[elset].Val * Val  #
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT, self.nInt - 1, j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
    def FormB(self, r, s, t, NLGeom):
        L = self.LL
        Li = 1./L                                                           # strain-Jacobian, , 2 is in following B 1st row
        cc = self.cosA
        ss = self.sinA
        c2 = self.cosA**2
        s2 = self.sinA**2
        r2 = 2.*r
        x  = 0.5*( (1-r)*self.X0 + (1+r)*self.X1 )                          # 0.5 for interploation of coordinate
        x_ = 0.5/x                                                          # 0.5 for interpolation of displacement in following B 2nd row
        B = array([[         (-1.+r2)*cc*Li,(-1.+r2)*ss*Li,         (1.+r2)*cc*Li, (1.+r2)*ss*Li,                -2.*r2*Li],
                   [x_*(-r+r2)*c2+(1.-r)*s2,             0,x_*(r+r2)*c2+(1.+r)*s2,             0, x_*2.*cc*(1.-r2)*(cc+ss)]])
        return B, 1, 0                                                      # returns dummy for integration jacobian; integration-jacobian in Geom[0,0] from T2D2:Ini2

class T3D2( Element ):                                                      # 3D Truss 2 nodes
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, SolSecDic, StateV, NData, Val, NoLabToNoInd):
        Element.__init__(self,"T3D2",InzList, 3, 0,1,1, (set([1,2,3]),set([1,2,3])), 1,True,Label,SetLabel,3, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        self.CrossSecVal = Val                                              # correction factor for cross section
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        n0, n1 = NoList[i0], NoList[i1]
        self.CoordRef = array( [n0.XCo, n0.YCo, n0.ZCo, n1.XCo, n1.YCo, n1.ZCo]) # Nodal coordinates 
        self.CoordDis = array( [n0.XCo, n0.YCo, n0.ZCo, n1.XCo, n1.YCo, n1.ZCo]) # Nodal coordinates displaced 
        dx = self.CoordRef[3]-self.CoordRef[0]
        dy = self.CoordRef[4]-self.CoordRef[1]
        dz = self.CoordRef[5]-self.CoordRef[2]
        L = sqrt(dx**2 + dy**2 + dz**2)
        self.Lch_ = L                                                       # length in undeformed configuration !
        self.LL = L
        elset = self.Set
        self.Geom = zeros( (2,5), dtype=double)
        self.Geom[1,0] = SecDict[elset].Val*self.CrossSecVal                     # cross section area
        self.Geom[0,0] = 0.5*L                                              # Jacobian
        self.Geom[1,1] = dx
        self.Geom[1,2] = dy
        self.Geom[1,3] = dz
        self.Geom[1,4] = L*L
        vL = array([ dx/L, dy/L, dz/L ])                                # longitudinal orientation
        # create local coordinate system
        v3 = cross(vL,array([1., 0., 0.]))
        if norm(v3)<1.0e-3: v3 = cross(vL,array([0., 1., 0.]))
        L3 = sqrt( v3[0]**2 + v3[1]**2 + v3[2]**2 )
        v3 = v3/L3
        v2 = cross(v3, vL)
        self.RotMat = array([[vL[0],vL[1],vL[2]],
                            [v2[0],v2[1],v2[2]],
                            [v3[0],v3[1],v3[2]]])
        Element.BarLength += L
        return []
    def UpdateCoord(self, dis, ddis ):
        self.CoordDis[0] = self.CoordRef[0]+dis[0]
        self.CoordDis[1] = self.CoordRef[1]+dis[1]
        self.CoordDis[2] = self.CoordRef[2]+dis[2]
        self.CoordDis[3] = self.CoordRef[3]+dis[3]
        self.CoordDis[4] = self.CoordRef[4]+dis[4]
        self.CoordDis[5] = self.CoordRef[5]+dis[5]
        self.Geom[1,1] = self.CoordDis[3]-self.CoordDis[0]
        self.Geom[1,2] = self.CoordDis[4]-self.CoordDis[1]
        self.Geom[1,3] = self.CoordDis[5]-self.CoordDis[2]
        return
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((self.DofE, self.DofE), dtype=float)
        val = sig[0]/self.Geom[1,4]
        GeomK[0,0] = val
        GeomK[0,3] = -val
        GeomK[1,1] = val
        GeomK[1,4] = -val
        GeomK[2,2] = val
        GeomK[2,5] = -val
        GeomK[3,0] = -val
        GeomK[3,3] = val
        GeomK[4,1] = -val
        GeomK[4,4] = val
        GeomK[5,2] = -val
        GeomK[5,5] = val 
        return GeomK
    def FormN(self, r, s, t):
        N = array([[ 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0.],
                   [ 0., 0., 0.5*(1-r), 0., 0., 0.5*(1+r)]])
        return N
    def FormB(self, r, s, t, NLg):
        ll = self.Geom[1,1]/self.Geom[1,4]                          # Jacobian incorporated
        mm = self.Geom[1,2]/self.Geom[1,4]
        nn = self.Geom[1,3]/self.Geom[1,4]
        B = array([[ -ll, -mm, -nn, ll, mm, nn], [ 0, 0, 0, 0, 0, 0]])
        return B, 1, 0
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0.],
                   [ 0., 0., 0.5*(1-r), 0., 0., 0.5*(1+r)]])                  # last three entries for internal node with coordinates [0,0] 
        return X
    def FormT(self, r, s, t):
        T = array([ 0.5*(1-r),0.5*(1-r),0.5*(1-r), 0.5*(1+r),0.5*(1+r),0.5*(1+r)])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

# 3D truss element   with additional longitudinal degree of freedom
class T3D3(T3D2):                                  # 3D truss 2 nodes
    # Ini2 from T3D2
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, SolSecDic, StateV, NData, Val, NoLabToNoInd):
        Element.__init__(self,"T3D3", InzList, 3, 0,2,2, (set([1,2,3]),set([1,2,3]),set([1])), 1,False,Label,SetLabel,3, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        self.CrossSecVal = Val                                              # correction factor for cross section
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r), 0.])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0., 0., 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0., 0., 0.],
                   [ 0., 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0., 0.]])                  # last three entries for internal node with coordinates [0,0] 
        return X
    def FormN(self, r, s, t):
        g00 = self.RotMat[0,0]
        g01 = self.RotMat[0,1]
        g02 = self.RotMat[0,2]
        g10 = self.RotMat[1,0]
        g11 = self.RotMat[1,1]
        g12 = self.RotMat[1,2]
        g20 = self.RotMat[2,0]
        g21 = self.RotMat[2,1]
        g22 = self.RotMat[2,2]
        r2 = r**2
        N00 = (-r+r2)*g00*g00 + (1.-r)*g10*g10 + (1.-r)*g20*g20
        N01 = (-r+r2)*g00*g01 + (1.-r)*g10*g11 + (1.-r)*g20*g21
        N02 = (-r+r2)*g00*g02 + (1.-r)*g10*g12 + (1.-r)*g20*g22
        N03 = ( r+r2)*g00*g00 + (1.+r)*g10*g10 + (1.+r)*g20*g20
        N04 = ( r+r2)*g00*g01 + (1.+r)*g10*g11 + (1.+r)*g20*g21
        N05 = ( r+r2)*g00*g02 + (1.+r)*g10*g12 + (1.+r)*g20*g22
        N06 = 2.*g00*(1-r2)
        N10 = (-r+r2)*g00*g01 + (1.-r)*g11*g10 + (1.-r)*g20*g21
        N11 = (-r+r2)*g01*g01 + (1.-r)*g11*g11 + (1.-r)*g21*g21
        N12 = (-r+r2)*g02*g01 + (1.-r)*g11*g12 + (1.-r)*g22*g21
        N13 = ( r+r2)*g00*g01 + (1.+r)*g11*g10 + (1.+r)*g20*g21
        N14 = ( r+r2)*g01*g01 + (1.+r)*g11*g11 + (1.+r)*g21*g21
        N15 = ( r+r2)*g02*g01 + (1.+r)*g11*g12 + (1.+r)*g22*g21
        N16 = 2.*g01*(1-r2)
        N20 = (-r+r2)*g00*g02 + (1.-r)*g12*g10 + (1.-r)*g20*g22
        N21 = (-r+r2)*g01*g02 + (1.-r)*g12*g11 + (1.-r)*g21*g22
        N22 = (-r+r2)*g02*g02 + (1.-r)*g12*g12 + (1.-r)*g22*g22
        N23 = ( r+r2)*g00*g02 + (1.+r)*g12*g10 + (1.+r)*g20*g22
        N24 = ( r+r2)*g01*g02 + (1.+r)*g12*g11 + (1.+r)*g21*g22
        N25 = ( r+r2)*g02*g02 + (1.+r)*g12*g12 + (1.+r)*g22*g22
        N26 = 2.*g02*(1-r2)
        N  = 0.5*array([[N00, N01, N02, N03, N04, N05, N06],[N10, N11, N12, N13, N14, N15, N16],[N20, N21, N22, N23, N24, N25, N26]]) # requires presumably modification similar to T2D3
        return N
    def FormB(self, r, s, t, NLg):
        L = self.LL
        Li = 1./L
        g00 = self.RotMat[0,0]
        g01 = self.RotMat[0,1]
        g02 = self.RotMat[0,2]
        r2 = 2.*r
        B = array([[(-1.+r2)*Li*g00, (-1.+r2)*Li*g01, (-1.+r2)*Li*g02, (1.+r2)*Li*g00, (1.+r2)*Li*g01, (1.+r2)*Li*g02, -2.*r2*Li],[0,0,0,0,0,0,0]])
        return B, 1., 0
    def FormT(self, r, s, t):
        T = array([ 0.5*(1-r), 0.5*(1+r), 0, 0, 0, 0, 0])           # remains to be controlled / implemented
        return T

# embedded truss elements -- created with BONDLAW-option  for elements TxDx
class TxDxI(object):
    def __init__(self, NoList):
        self.BondElInd  = zeros((self.nIntL), dtype=int)
        self.FlexBonded = True
    # determine CxD elements where TxDxI integration points lie in and create bond elements
    def CreateBond(self, ElList,NodeList,SecDict, CoorTree_,NoCoLitoNoLi, MatList, NoLabToNoInd,NoIndToCMInd): # called by ["T2D2I","T2D3I","T3D2I","T3D3I","B23I","B23EI"]
        iself = FindIndexByLabel( ElList, self.Label)                       # index in element list
        xy = self.NodalCoordinates( NodeList, NoIndToCMInd )                # nodal coordinates
        for i in range(self.nIntL): # continuum element to be an input in the bond element class - Ahmad
            k, rs, R, IndR, rst, intW = self.IntPointInBulk( i, xy, CoorTree_,NoCoLitoNoLi, NodeList, ElList, NoIndToCMInd)
            # k: index of continuum element, rs: local coor continuum element, R: form function values continuum nodes, indR: index continuum nodes,
            # rst: local coor truss, intW: integration weight of T2D2I element
            el = ElList[k]
            #
            if el.Type in ["CPE3","CPE3S","CPE6","CPE6S","CPE4","CPE4S"  ,  "CPS3","CPS3S","CPS6","CPS6S","CPS4","CPS4S"]:
                ContThick = SecDict[el.Set].Val                             # thickness of underlying 2D continuum -- might depend on integration point of embedded element
            elif el.Type in ["CAX4"]:
                ContThick = 1.0                                             # unit thickness -- circumferential direction -- for axissymmetric continuum
            elif el.Type in ["C3D8","C3D8S"]:
                pass
            else:
                raise NameError("ComFemElem::CreateBond: unknown continuum element type", el.Type)
            #
            if i==0:
                if SecDict[self.Set].nRebar != None and self.Type != "TAX2I": # presumably because TAX2I is defined as sheet
                    self.nRebarEl = SecDict[self.Set].nRebar * ContThick
                else:
                    self.nRebarEl = 1.
            #
            NData = 4                                                       # for 2D only
            if self.Type in   ["T2D2I"]:
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label] # subject to change in the following to include nodes of continuum
                ElList += [Bond2D2(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]],               self.Geom[1,0], self.BondLaw,self.nRebarEl, MatList[self.BondLaw].StateVar, NData, NoLabToNoInd, NoIndToCMInd, el, self)]
            elif self.Type in  ["T2D3I"]:
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[self.Inzi[2]].Label] # subject to change in the following to include nodes of continuum
                ElList += [Bond2D3(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]],               self.Geom[1,0], self.BondLaw,self.nRebarEl, MatList[self.BondLaw].StateVar, NData, NoLabToNoInd, NoIndToCMInd, el, self)]
            elif self.Type in ["TAX2I"]:
                nrebar =  SecDict[self.Set].nRebar
                AA = self.Geom[1,0]/nrebar                                  # cross section of single rebar
                xx = 0.5*((1.-rst[0])*xy[0] + (1.+rst[0])*xy[2])            # radial distance from center
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [BondAX2(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]], xx,           AA,self.BondLaw, nrebar, NData, NoLabToNoInd, el, self)]
            elif self.Type in ["TAX3I"]:
                nrebar =  SecDict[self.Set].nRebar
                AA = self.Geom[1,0]/nrebar                                  # cross section of single rebar
                xx = 0.5*((1.-rst[0])*xy[0] + (1.+rst[0])*xy[2])            # radial distance from center
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[self.Inzi[2]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [BondAX3(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]], xx,           AA,self.BondLaw, nrebar, NData, NoLabToNoInd, el, self)]
            #
            elif self.Type in ["T3D2I"]:
                InzList = [NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [Bond3D2(0, InzList, NodeList, self.Set, self.BondLaw, IndR, R, intW, k, rs, iself, i, [rst[0]], self.Geom[1,0], NoLabToNoInd, MatList[self.BondLaw].StateVar, MatList[self.BondLaw].NData, NoIndToCMInd, el, self)]
            elif self.Type in ["T3D3I"]:
                InzList = [NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[ self.Inzi[2]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [Bond3D3(0, InzList, NodeList, self.Set, self.BondLaw, IndR, R, intW, k, rs, iself, i, [rst[0]], self.Geom[1,0], NoLabToNoInd, MatList[self.BondLaw].StateVar, MatList[self.BondLaw].NData, NoIndToCMInd, el, self)]
            #
            elif self.Type in ["B23I"]:
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label] # subject to change in the following to include nodes of continuum
                ElList += [Bond2D2(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]],              self.Geom[1,3], self.BondLaw,self.nRebarEl, MatList[self.BondLaw].StateVar, NData, NoLabToNoInd, NoIndToCMInd, el, self)]
            elif self.Type in ["B23EI"]:
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[self.Inzi[2]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [Bond2D3(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]],              self.Geom[1,3], self.BondLaw,self.nRebarEl, MatList[self.BondLaw].StateVar, NData, NoLabToNoInd, NoIndToCMInd, el, self)]
            elif self.Type in ["BAX23I"]:
                xx = 0.5*((1.-rst[0])*xy[0] + (1.+rst[0])*xy[2])            # radial distance from center -- x, y in one array
                InzList = [NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label]                                 # subject to change in the following to include nodes of continuum
                ElList += [BondAX2(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]], xx,          self.Geom[1, 3], self.BondLaw, self.nRebarEl, NData, NoLabToNoInd, el, self)]
            elif self.Type in ["BAX23EI"]:
                xx = 0.5*((1.-rst[0])*xy[0] + (1.+rst[0])*xy[4])            # radial distance from center -- x, y in one array -- different node ordering compared to BAX21EI
                InzList = [NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[self.Inzi[2]].Label]   # subject to change in the following to include nodes of continuum
                ElList += [BondAX3(0, InzList, NodeList, self.Set, IndR, R, intW, k, rs, iself, i, [rst[0]], xx,          self.Geom[1, 3], self.BondLaw, self.nRebarEl, NData, NoLabToNoInd, el, self)]
            elif self.Type in ["BAX21EI"]:
                xx = 0.5*((1.-rst[0])*xy[0] + (1.+rst[0])*xy[2])            # radial distance from center
                InzList = [ NodeList[self.Inzi[0]].Label, NodeList[self.Inzi[1]].Label, NodeList[self.Inzi[2]].Label]  # subject to change in the following to include nodes of continuum
                ElList += [BondAX3(0, InzList, NodeList, self.Set, IndR, R,intW,k,rs,iself,i,[rst[0]],xx,                 self.Geom[1,3], self.BondLaw,self.nRebarEl, NData, NoLabToNoInd, el, self)]
            else:
                raise NameError("ComFemElem::CreateBond: unknown bond element type", self.Type)
            self.BondElInd[i] = len(ElList)-1                               # index of current corresponding bond element
    # find background bulk element, local coordinates of respective point and weighting factors
    def IntPointInBulk(self, i, xy, CoorTree,NoCoLitoNoLi, NodeList, ElList, NoIndToCMInd):
        rst = SamplePoints[self.IntT,self.nInt-1,i]                         # local integration point coordinates of bonded element
        intW= SampleWeight[self.IntT,self.nInt-1,i]                         # integration weight of bonded element
        XX  = self.FormX_( rst[0], 0., 0.)
        xyP = dot(XX,array(xy))                                             # global integration point coordinates
        if len(xyP)==2: res = CoorTree.query([array([xyP[0],xyP[1],0.])],k=2) # 2D, k nearest points, 2 for coinciding points
        else:           res = CoorTree.query([xyP],k=2)                     # 3D,       "
        dd = res[0][0]                                                      # distances to xyP, second index obviously dummy
        jj = NoCoLitoNoLi[res[1][0]]                                        # list of length k, NoCoLiToNoLi maps indices of list of continuum nodes (special list for CoorTree) to those in NodeList
        Finis = False
        # list of elements which have this node
        if abs(dd[0]-dd[1])<ZeroD: NodeEl = list(set(NodeList[jj[0]].NodeEl+NodeList[jj[1]].NodeEl))
        else:                      NodeEl = NodeList[jj[0]].NodeEl
        for k in NodeEl:                                                    # loop over respective elements, k is index for ElList
            el  = ElList[k]
            if el.Type in ["CPS3","CPS3S","CPE3","CPE3S","CPS4","CPS4S","CPE4","CPE4S","C3D8","C3D4","CAX4"]:  # continuum elements relevant only
                # rs: local coordinates of bond point in continuum element, R: form function values from continuum nodes; IndR: indices of corresponding nodes
                Finis, rs, R, IndR = el.LocatePointInElement( NodeList, xyP, NoIndToCMInd) # --> ElementC2D
                if Finis: return k, rs, R, IndR, rst, intW                  # terminates with 1st admissable continuum element found
        if not Finis: raise NameError ("ComFemElem::CreateBond: Continuum Element not found for truss integration point ",self.Label,i,rst,xyP)

# embedded truss element 2D
class T2D2I(T2D2, TxDxI):
    # Ini2 from T2D2
    # uses LocatPointInEelement, IniBond, FormX, FormN, FormB from T2D2
    # uses CreateBond, IntPointInBulk from TxDxI
    def __init__(self, Label, InzList, SolSecDic,Val, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
        Element.__init__(self,"T2D2I", InzList,2, 0,2,2, (set([1,2]),set([1,2])), 1,False,Label,elSet,2, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
        self.Geom = zeros( (2,4), dtype=double)                             # ? should be there again, see below T2D3I, skipped below with ini2 and deferred definition in DataIn
        self.Geom[1,0] = SolSecDic[elSet].Val*Val                           # cross section area - SolSecDic available in ConFemInOut only
        self.BondLaw = BondLaw
    def FormB(self, r, s, t, val):                                          # not taken from T2D2 which has large displacments
        L = self.LL
        B = array([[ -self.cosA/L, -self.sinA/L, self.cosA/L, self.sinA/L],[0,0,0,0]])
        return B, 1., 0
# embedded truss element with additional longitudinal degree of freedom
#class T2D3I(T2D3, T2D2I):  # sequence is significant as both parents implement FormX, FormN, FormB; first should be used (inheritance from left to right, see Lutz/Ascher p. 405)
class T2D3I(T2D3, T2D2I):  # sequence is significant as both parents implement FormX, FormN, FormB; first should be used (inheritance from left to right, see Lutz/Ascher p. 405)
    # Ini2 from T2D3 -> T2D2
    # uses FormX, FormN, FormB from T2D3
    # uses CreateBond from TxDxI
    def __init__(self, Label, InzList, SolSecDic,Val, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
        Element.__init__(self,"T2D3I", InzList,2, 0,2,2, (set([1,2]),set([1,2]),set([1])), 1,False,Label,elSet,2, MatName,None,             None,           NoList,NoLabToNoInd,NoIndToCMInd)
#       Element.__init__(self,"T2D3I", InzList, nFieVal,IntTVal,nIntVal,nIntLVal, DofTVal, dimVal,NLGeomIVal,Label,elSet,SpaceDimVal MatName,StateV,NData,  NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
        self.Geom = zeros( (2,4), dtype=double)
        self.Geom[1,0] = SolSecDic[elSet].Val*Val                           # cross section area - SolSecDic available in ConFemInOut only
        self.BondLaw = BondLaw

class TAX2I(TAX2, TxDxI):                                                   # embedded axisymmetric Truss 2 nodes
    # generates BondAX2
    # Ini2 from T2D2
    def __init__(self, Label, elset, InzList, NoList, SolSecDic, Val, NoLabToNoInd, BondLaw):
        Element.__init__(self,"TAX2I",InzList, 1, 0,2,2, (set([1,2]),set([1,2])), 5,False,Label,elset,2, None,None,None, NoList,NoLabToNoInd,[])
        TxDxI.__init__(self, NoList)
        self.Geom = zeros( (2,4), dtype=double)                             #
        self.Geom[1,0] = SolSecDic[elset].Val*Val                           #
        self.BondLaw = BondLaw                                              # transferred to bond elements
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT, self.nInt - 1, j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
class TAX3I( TAX3, TxDxI):                                                 # embeddes axisymmetric Truss 2 + 1 nodes
    # generates BondAX2
    # Ini2 from T2D2
    def __init__(self, Label, elset, InzList, NoList, SolSecDic, Val, NoLabToNoInd, BondLaw):
        Element.__init__(self,"TAX3I",InzList, 1, 0,2,2, (set([1,2]),set([1,2]),set([1])), 5,False,Label,elset,2, None,None,None, NoList,NoLabToNoInd,[])
        TxDxI.__init__(self, NoList)
        self.Geom = zeros((2, 4), dtype=double)  #
        self.Geom[1, 0] = SolSecDic[elset].Val * Val  #
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        self.BondLaw = BondLaw                                              # transferred to bond elements
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT, self.nInt - 1, j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)

# embedded truss element 3D
class T3D2I(T3D2, TxDxI):
    # uses IniBond, LocatePointInElement, FormX, FormN, FormB from  T3D2
    # uses CreateBond, IntPointInBulk from TxDxI
    def __init__(self, Label, InzList,    SolSecDic,          Val, NoList,    elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
#       Element.__init__(self,"T3D2I",InzList, 3, 0,2,2, (set([1,2,3]),set([1,2,3])), 1,False,Label,elSet,3, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,NoIndToCMInd)
        Element.__init__(self,"T3D2I",InzList, 3, 0,2,2, (set([1,2,3]),set([1,2,3])), 1,False,Label,elSet,3, MatName,None,             None,           NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
#        self.IniBond( AA, NoList)
        self.CrossSecVal = Val
        self.Geom = zeros( (2,5), dtype=double)
        self.Geom[0,0] = 1.
        self.Geom[1,0] = SolSecDic[elSet].Val*Val                           # cross section area - SolSecDic available in ConFemInOut only
        self.BondLaw = BondLaw
    def FormNI(self):                                                       # N-inverse for gaussian 3D order 2x2x2 integration, see N_Inverse.py, to get nodal values from integration point values
        return array([ [ 1.366025404, -0.3660254039],                       # extrapolation to left node, see ExtrapolateIntpointsToNodesT3D3E_mws.pdf
                       [-0.366025404,  1.366025404]])                       # extrapolation to right node
# embedded truss element 3D with additional longitudinal degree of freedom
class T3D3I(T3D3, TxDxI):
    # uses IniBond, LocatePointInElement, FormX from  T3D2
    # uses CreateBond, IntPointInBulk from TxDxI
#   def __init__(self, Label, InzList, AA, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
    def __init__(self, Label, InzList, SolSecDic,Val, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
#       Element.__init__(self,"T3D3I",InzList, 3, 0,2,2, (set([1,2,3]),set([1,2,3]),set([1])), 1,False,Label,elSet,3, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,NoIndToCMInd)
        Element.__init__(self,"T3D3I",InzList, 3, 0,2,2, (set([1,2,3]),set([1,2,3]),set([1])), 1,False,Label,elSet,3, MatName,None,             None,           NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
        self.CrossSecVal = Val                                              # correction factor for cross section - needed for Ini2
        self.Geom = zeros( (2,5), dtype=double)
        self.Geom[0,0] = 1.
        self.Geom[1,0] = SolSecDic[elSet].Val*Val                           # will be repeated in Ini2 which is for all T3D*
        self.BondLaw = BondLaw
    def FormNI(self):                                                       # N-inverse for gaussian 3D order 2x2x2 integration, see N_Inverse.py, to get nodal values from integration point values
        return array([[  1.366025404,  -0.3660254039],                      # extrapolation to left node, see ExtrapolateIntpointsToNodesT3D3E_mws.pdf
                      [ -0.366025404 ,  1.366025404],                       # extrapolation to right node
                      [  0.5,           0.5 ]])                             # interpolation to center node

class B21(Element):                                                 # Timoshenko Beam 2D 2 nodes, linear shape
    # IniBeam in InOut-Post
    def __init__(self, Label, SetLabel, InzList, NoList, NoLabToNoInd):
        Element.__init__(self,"B21",InzList, 3, 0,1,1, (set([1, 2, 6]),         set([1, 2, 6])), 11,False,Label,SetLabel,2, None,None,None, NoList,NoLabToNoInd,[])
        self.TensStiff = True                                       # flag for tension stiffening
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        L = sqrt( (NoList[i1].XCo-NoList[i0].XCo)**2 + (NoList[i1].YCo-NoList[i0].YCo)**2 )# element length
        self.Geom[0,0] = 0.5*L                                      # Jacobi determinant value
        self.Lch_ = L                                               # characteristic length
        cosA = (NoList[i1].XCo-NoList[i0].XCo)/L                    # direction cosine
        sinA = (NoList[i1].YCo-NoList[i0].YCo)/L                    # direction sine
        if fabs(sinA) >= self.ZeroD:                                # element direction rotated to global axis
            self.Rot = True
            self.Trans = zeros((self.DofE, self.DofE), dtype=float) # coordinate transformation matrix
            for i in range(self.DofE): self.Trans[i,i] = cosA       # fill transformation axis
            if self.Type in ["B21"]:
                self.Trans[2,2] = 1.
                self.Trans[5,5] = 1.
                self.Trans[0,1] = sinA
                self.Trans[1,0] = -sinA
                self.Trans[3,4] = sinA
                self.Trans[4,3] = -sinA
            elif self.Type in ["B21E"]:
                self.Trans[0, 1] = sinA
                self.Trans[1, 0] = -sinA
                self.Trans[2, 2] = 1.
                self.Trans[3, 4] = sinA
                self.Trans[4, 3] = -sinA
                self.Trans[5, 5] = 1.
                self.Trans[6, 7] = sinA
                self.Trans[7, 6] = -sinA
                self.Trans[8, 8] = 1.
            else: raise NameError("XXX")
        return []
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X
    def FormN(self, r, s, t):
        N = array([[0.5*(1-r),0.,0.,0.5*(1+r),0.,0.],
                   [0.,0.5*(1-r),0.,0.,0.5*(1+r),0.],
                   [0.,0.,0.5*(1-r),0.,0.,0.5*(1+r)]])
        return N
    def FormB(self, r, s, t, NLg):
        L = 2.*self.Geom[0,0]
        B = array([[-1./L,0,0,1./L,0,0],
                   [0,0,-1./L,0,0,1./L],
                   [0,-1./L,-0.5*(1-r),0,1./L,-0.5*(1+r)]])
        return B, 1, 0
    def FormT(self, r, s, t):                                       # for temperatures
        T = array([[ 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0.]])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1
class B21E( B21 ):                                                   # Timoshenko Beam 2D 3 nodes, quadratic shape for displacements
    # dim  = 11                                                       # material type index, timoshenko beam
    # NLGeomI = False                                                 # Flag whether large deformations are already implemented for this element
    # IniBeam in InOut-Post
    # Ini2 from B21
    def __init__(self, Label, SetLabel, InzList, NoList, NoLabToNoInd):
        Element.__init__(self,"B21E",InzList, 3,0,2,2, (set([1, 2, 6]),set([1, 2, 6]),set([1, 2, 6])), 11,False,Label,SetLabel,2, None,None,None, NoList, NoLabToNoInd, [])
        self.TensStiff = True                                       # flag for tension stiffening
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r), 0.])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0., 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r), 0., 0.]])
        return X
    def FormN(self, r, s, t):
        N = array([[0.5*(r**2-r),0.,0.,  0.5*(r**2+r),0.,0.,   1.-r**2,0.,0.],
                   [0.,0.5*(r**2-r),0.,  0.,0.5*(r**2+r),0.,   0.,1.-r**2,0.],
                   [0.,0.,0.5*(r**2-r),  0.,0.,0.5*(r**2+r),   0.,0.,1.-r**2]])
        return N
    def FormB(self, r, s, t, NLg):
        L = 2.*self.Geom[0,0]                                       # Geom[0,0] is jacobi determinant
        B = array([[(2*r-1)/L,0.,0.,              (2*r+1)/L,0.,0.,              -4*r/L,0.,    0.],
                   [ 0.,0.,(2*r-1)/L,             0.,0.,(2*r+1)/L,               0.,   0.,   -4*r/L],
                   [ 0.,(2*r-1)/L,-0.5*(r**2-r),  0.,(2*r+1)/L,-0.5*(r**2+r),    0.,  -4*r/L,-(1-r**2)]])
        return B, 1, 0
    def FormT(self, r, s, t):                                               # for temperatures
        T = array([[ 0.5*(1-r),0.,0.,  0.5*(1+r),0.,0.,  0.,0.,0.],
                   [ 0.,0.5*(1-r),0.,  0.,0.5*(1+r),0.,  0.,0.,0.]])
        return T
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

class B23(Element):                                                 # Bernoulli Beam 2D 2 nodes, cubic shape
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, BeamSec, ReinfD, StateV, NData, ff, NoLabToNoInd, Material):
        self.Material = Material
#       Element.__init__(self,"B23",InzList, 2, 0,2,2, (set([1, 2, 6]),         set([1, 2, 6])),  10,True,Label,SetLabel,2, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,[])
        Element.__init__(self,"B23",InzList, 2, 0,2,2, (set([1, 2, 6]),         set([1, 2, 6])),  10,True,Label,SetLabel,2, None,   None,             None,           NoList,NoLabToNoInd,[])
        self.TensStiff = True                                               # flag for tension stiffening
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):                   # called by ConFemBasics:AssignGlobalDof
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        self.CoordRef = array( [NoList[i0].XCo, NoList[i0].YCo, NoList[i1].XCo, NoList[i1].YCo]) # Nodal coordinates
        self.CoordDis = array( [NoList[i0].XCo, NoList[i0].YCo, NoList[i1].XCo, NoList[i1].YCo]) # Nodal coordinates displaced
        L =       sqrt( (NoList[i1].XCo-NoList[i0].XCo)**2 + (NoList[i1].YCo-NoList[i0].YCo)**2 )# element length
        self.Geom[0,0] = 0.5*L                                              # Jacobi determinant value
        self.LL = L
        self.Lch_ = L
        self.cosA = (NoList[i1].XCo-NoList[i0].XCo)/L                       # direction cosine
        self.sinA = (NoList[i1].YCo-NoList[i0].YCo)/L                       # direction sine
        if fabs(self.sinA) >= self.ZeroD:                                   # element direction rotated to global axis
            self.Rot = True
            C, S = self.cosA, self.sinA
            self.Trans = array([[ C, S, 0, 0, 0, 0],                        # coordinate transformation matrix
                              [-S, C, 0, 0, 0, 0],
                              [ 0, 0, 1, 0, 0, 0],
                              [ 0, 0, 0, C, S, 0],
                              [ 0, 0, 0,-S, C, 0],
                              [ 0, 0, 0, 0, 0, 1]])
        return []
    def UpdateCoord(self, dis, ddis ):
        if not self.Rot: self.Trans = zeros((self.DofE, self.DofE), dtype=float)# coordinate transformation matrix
        self.Rot = True
        self.CoordDis[0] = self.CoordRef[0]+dis[0]
        self.CoordDis[1] = self.CoordRef[1]+dis[1]
        self.CoordDis[2] = self.CoordRef[2]+dis[3]
        self.CoordDis[3] = self.CoordRef[3]+dis[4]
        L = sqrt( (self.CoordDis[2]-self.CoordDis[0])**2 + (self.CoordDis[3]-self.CoordDis[1])**2 )# element length
        cosA = (self.CoordDis[2]-self.CoordDis[0])/L                        # direction cosine
        sinA = (self.CoordDis[3]-self.CoordDis[1])/L                        # direction sine
        for i in range(self.DofE): self.Trans[i,i] = cosA                   # fill transformation axis
        self.Trans[2,2] = 1.
        self.Trans[5,5] = 1.
        self.Trans[0,1] = sinA
        self.Trans[1,0] = -sinA
        self.Trans[3,4] = sinA
        self.Trans[4,3] = -sinA
        self.Geom[0,0] = 0.5*L
        self.sinA = sinA
        self.cosA = cosA
        return
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((self.DofE, self.DofE), dtype=float)
        sinA = self.sinA
        cosA = self.cosA
        L = 2.*self.Geom[0,0]
        NI = -sig[0]/L
        VI =  sig[1]*6*r/L**2
        AI = sinA*NI + cosA*VI
        BI = cosA*NI - sinA*VI
        GeomK[0,0] = -AI*sinA
        GeomK[0,1] =  AI*cosA
        GeomK[0,3] =  AI*sinA
        GeomK[0,4] = -AI*cosA
        GeomK[1,0] =  BI*sinA
        GeomK[1,1] = -BI*cosA
        GeomK[1,3] = -BI*sinA
        GeomK[1,4] =  BI*cosA

        GeomK[3,0] = -GeomK[0,0]
        GeomK[3,1] = -GeomK[0,1]
        GeomK[3,3] = -GeomK[0,3]
        GeomK[3,4] = -GeomK[0,4]
        GeomK[4,0] = -GeomK[1,0]
        GeomK[4,1] = -GeomK[1,1]
        GeomK[4,3] = -GeomK[1,3]
        GeomK[4,4] = -GeomK[1,4]
        return GeomK
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X
    def FormN(self, r, s, t):
        L = self.LL
        N = array([[0.5*(1-r),0.,0.,0.5*(1+r),0.,0.],
                   [0, 0.25*(r**3-3*r+2), L*0.125*(r**3-r**2-r+1),0, 0.25*(-r**3+3*r+2), L*0.125*(r**3+r**2-r-1)]])
        return N
    def FormB(self, r, s, t, NLg):
        L = self.LL
        B = array([[-1./L,0,0,1./L,0,0],
                   [0, 6*r/L**2, (3*r-1)/L, 0, -6*r/L**2, (3*r+1)/L]])
        return B, 1, 0
    def FormT(self, r, s, t):                                               # for temperatures
        T = array([[ 0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0.5*(1+r), 0.]])
        return T
    def FormP(self, r, s):                                                  # prestressing interpolation, 1st row: force, 2nd row: local vertical coordinate, 3rd row: local inclination
        L = self.LL
        P = array([[0.5*(1-r), 0., 0., 0.5*(1+r), 0., 0.],
                   [0, 0.25*(r**3-3*r+2), L*0.125*(r**3-r**2-r+1),0, 0.25*(-r**3+3*r+2), L*0.125*(r**3+r**2-r-1)],
                   [0, 0.5* (3*r**2-3)/L, 0.25*(3*r**2-2*r-1),    0, 0.5* (-3*r**2+3)/L, 0.25*(3*r**2+2*r-1)]])
        return P
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

class B23E(Element):
    """ Bernoulli Beam 2D 2 nodes, cubic shape for transverse displacements.
    """
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, BeamSec, ReinfD, StateV, NData, ff, NoLabToNoInd, Material):
        Element.__init__(self,"B23E",InzList,2, 0,3,3, (set([1, 2, 6]),set([1]),set([1, 2, 6])), 10,True,Label,SetLabel,2, None,   None,  None, NoList,NoLabToNoInd,[])
#        Element.__init__(self,"B23E",InzList,2, 0,3,3, (set([1, 2, 6]),set([1, 2, 6]),set([1])), 10,True,Label,SetLabel,2, None,   None,  None, NoList,NoLabToNoInd,[])
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
#        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        n0, n2 = NoList[i0], NoList[i2]
        self.CoordRef = array( [n0.XCo, n0.YCo, n2.XCo, n2.YCo])            # Nodal coordinates
        self.CoordDis = array( [n0.XCo, n0.YCo, n2.XCo, n2.YCo])            # Nodal coordinates displaced
        L = sqrt( (n2.XCo-n0.XCo)**2 + (n2.YCo-n0.YCo)**2 )                 # element length
        self.LL = L
        self.Geom[0,0] = 0.5*L                                              # Jacobi determinant value
        self.Lch_ = L                                                       # characteristic length
        cosA = (n2.XCo-n0.XCo)/L                                            # direction cosine
        sinA = (n2.YCo-n0.YCo)/L                                            # direction sine
        self.cosA = cosA
        self.sinA = sinA
        if fabs(sinA) >= self.ZeroD:                                        # element direction rotated to global axis
            self.Rot = True
            self.Trans = array([[ cosA, sinA, 0, 0, 0,    0,    0],
                                [-sinA, cosA, 0, 0, 0,    0,    0],
                                [ 0,    0,    1, 0, 0,    0,    0],
                                [ 0,    0,    0, 1, 0,    0,    0],
                                [ 0,    0,    0, 0, cosA, sinA, 0],
                                [ 0,    0,    0, 0,-sinA, cosA, 0],
                                [ 0,    0,    0, 0, 0,    0,    1]])
        return []
    def UpdateCoord(self, dis, ddis ):
        if not self.Rot: self.Trans = zeros((self.DofE, self.DofE), dtype=float)# coordinate transformation matrix
        self.Rot = True
        self.CoordDis[0] = self.CoordRef[0]+dis[0]
        self.CoordDis[1] = self.CoordRef[1]+dis[1]
        self.CoordDis[2] = self.CoordRef[2]+dis[4]
        self.CoordDis[3] = self.CoordRef[3]+dis[5]
        L = sqrt( (self.CoordDis[2]-self.CoordDis[0])**2 + (self.CoordDis[3]-self.CoordDis[1])**2 )# element length
        cosA = (self.CoordDis[2]-self.CoordDis[0])/L                        # direction cosine
        sinA = (self.CoordDis[3]-self.CoordDis[1])/L                        # direction sine
        for i in range(self.DofE): self.Trans[i,i] = cosA                   # fill transformation axis
        self.Trans[2,2] = 1.
        self.Trans[3,3] = 1.                                                # longitudinal displ of center node remains local!
        self.Trans[6,6] = 1.
        self.Trans[0,1] = sinA
        self.Trans[1,0] = -sinA
        self.Trans[4,5] = sinA
        self.Trans[5,4] = -sinA
        self.Geom[0,0] = 0.5*L
        self.sinA = sinA
        self.cosA = cosA
        return
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((self.DofE, self.DofE), dtype=float)
        L = 2.*self.Geom[0,0]
        sinA = self.sinA
        cosA = self.cosA
        NI = -sig[0]/L
        VI =  sig[1]*6*r/L**2
        AI = sinA*NI + cosA*VI
        BI = cosA*NI - sinA*VI
        GeomK[0,0] = -AI*sinA
        GeomK[0,1] =  AI*cosA
        GeomK[0,4] =  AI*sinA
        GeomK[0,5] = -AI*cosA
        GeomK[1,0] =  BI*sinA
        GeomK[1,1] = -BI*cosA
        GeomK[1,4] = -BI*sinA
        GeomK[1,5] =  BI*cosA
        GeomK[4,0] = -GeomK[0,0]
        GeomK[4,1] = -GeomK[0,1]
        GeomK[4,4] = -GeomK[0,4]
        GeomK[4,5] = -GeomK[0,5]
        GeomK[5,0] = -GeomK[1,0]
        GeomK[5,1] = -GeomK[1,1]
        GeomK[5,4] = -GeomK[1,4]
        GeomK[5,5] = -GeomK[1,5]
        return GeomK
    def FormX(self, r, s, t):
        X = array([ 0.5*(1-r), 0., 0.5*(1+r)])
        return X
    def FormX_(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0., 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0., 0.5*(1+r)]])
        return X
    def FormN(self, r, s, t):
        L = 2.*self.Geom[0,0]
        N = array([[0.5*r*(r-1),0.,               0.,                     1-r**2,0.5*r*(1+r),0.,                0.],
                   [0.,         0.25*(r**3-3*r+2),L*0.125*(r**3-r**2-r+1),0.,    0.,         0.25*(-r**3+3*r+2),L*0.125*(r**3+r**2-r-1)]])
        return N
    def FormB(self, r, s, t, NLg):
#        L = self.LL                                                        # will not work for large displacements as LL is not updated during deformation
        L = 2.*self.Geom[0,0]
        B = array([[(2*r-1)/L, 0., 0., -4.*r/L, (2*r+1)/L, 0., 0.],
                   [0., 6*r/L**2, (3*r-1)/L, 0., 0., -6*r/L**2, (3*r+1)/L]])
        return B, 1, 0
    def FormT(self, r, s, t):                                               # for temperatures
        T = array([[ 0.5*(1-r), 0., 0., 0., 0.5*(1+r), 0., 0.],
                   [ 0., 0.5*(1-r), 0., 0., 0., 0.5*(1+r), 0.]])
        return T
    def FormP(self, r, s):
        """ Form functions for prestressing interpolation.
        1st row: force, 3rd row: local vertical coordinate, 4th row: local inclination
        """
        L = 2.*self.Geom[0,0]
        P = array([[0.5*(1-r),0.,                0.,                     0., 0.5* (1+r), 0.,                 0.],
                   [0.,       0.25*(r**3-3*r+2), L*0.125*(r**3-r**2-r+1),0., 0.,         0.25*(-r**3+3*r+2), L*0.125*(r**3+r**2-r-1)],
                   [0.,       0.5* (3*r**2-3)/L, 0.25*(3*r**2-2*r-1),    0., 0.,         0.5* (-3*r**2+3)/L, 0.25*(3*r**2+2*r-1)]])
        return P
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1

# embedded beam element 2D
class B23I(B23, TxDxI):
    # Ini2 from B23
    # uses LocatPointInEelement, IniBond, FormX, FormN, FormB from B23
    # uses CreateBond, IntPointInBulk from TxDxI
    def __init__(self, Label, InzList, BeamSec, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
        self.Material = Material
        Element.__init__(self,"B23I", InzList,2,      0,2,2, (set([1,2,6]),set([1,2,6])), 10,    False,     Label,elSet,2,           MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
        self.BondLaw = BondLaw
# embedded beam element 2D
class B23EI(B23E, TxDxI):
    # Ini2 from B23E
    # uses LocatPointInEelement, IniBond, FormX, FormN, FormB from B23
    # uses CreateBond, IntPointInBulk from TxDxI
    def __init__(self, Label, InzList, BeamSec, NoList, elSet, MatName, BondLaw, NoLabToNoInd, Material, NoIndToCMInd):
        self.Material = Material
#       Element.__init__(self,"B23EI",InzList,2, 0,3,3, (set([1,2,6]),set([1]),set([1,2,6])), 10,False,Label,elSet,2, MatName,Material.StateVar,Material.NData, NoList,NoLabToNoInd,NoIndToCMInd)
        Element.__init__(self,"B23EI",InzList,2, 0,2,2, (set([1,2,6]),set([1]),set([1,2,6])), 10,False,Label,elSet,2, None,   None,             None,           NoList,NoLabToNoInd,NoIndToCMInd)
#       Element.__init__(self,"B23EI",InzList,2, 0,3,3, (set([1,2,6]),set([1]),set([1,2,6])), 10,False,Label,elSet,2, None,   None,             None,           NoList,NoLabToNoInd,NoIndToCMInd)
        TxDxI.__init__(self, NoList)
        self.BondLaw = BondLaw
    def NodalCoordinates(self, NodeList, NoIndToCMInd):                     # supersedes NodalCoordinates from Element base class
        n0 = NoIndToCMInd[self.Inzi[0]]
        n2 = NoIndToCMInd[self.Inzi[2]]                                     # loop over nodes defining element geometry
        no0 = NodeList[n0]
        no2 = NodeList[n2]
        return array([no0.XCo,no0.YCo,0.,0.,no2.XCo,no2.YCo])

class BAX21(B21):                                                           # axisymmetric Timoshenko Beam 2D 2 nodes, linear shape
    # Ini2 from B21
    def __init__(self, Label, SetLabel, InzList, NoList, NoLabToNoInd):
        Element.__init__(self,"BAX21",InzList, 3, 0,1,1, (set([1, 2, 6]),         set([1, 2, 6])), 13,False,Label,SetLabel,2, None,None,None, NoList,NoLabToNoInd,[])
        self.TensStiff = True                                               # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X  = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT,self.nInt-1,j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
class BAX21E( B21E ):                                                       # axisymmetric Timoshenko Beam 2D 3 nodes, quadratic shape for displacements
    # Ini2 from B21
    # dim  = 13                                                             # material type index, timoshenko beam axisymmetric
    # NLGeomI = False                                                       # Flag whether large deformations are already implemented for this element
    def __init__(self, Label, SetLabel, InzList, NoList, NoLabToNoInd):
#        Element.__init__(self,"BAX21E",InzList, 3,0,2,2, (set([1, 2, 6]),set([1, 2, 6]),set([1, 2, 6])), 13,False,Label,SetLabel,2, None,None,None, NoList, NoLabToNoInd, [])
        Element.__init__(self,"BAX21E",InzList, 3,0,3,3, (set([1, 2, 6]),set([1, 2, 6]),set([1, 2, 6])), 13,False,Label,SetLabel,2, None,None,None, NoList, NoLabToNoInd, [])
        self.TensStiff = True                                       # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X  = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT,self.nInt-1,j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
class BAX21EI(BAX21E, TxDxI):                                               # embedded axisymmetric Timoshenko Beam 2D 3 nodes, quadratic shape for displacements
    # Ini2 from B21
    # dim  = 13                                                             # material type index, timoshenko beam axisymmetric
    # NLGeomI = False                                                       # Flag whether large deformations are already implemented for this element
    def __init__(self, Label,SetLabel,InzList, BondLaw, NoList,NoLabToNoInd):
#        Element.__init__(self,"BAX21EI",InzList, 3,0,2,2, (set([1, 2, 6]),set([1, 2, 6]),set([1, 2, 6])), 13,False,Label,SetLabel,2, None,None,None, NoList, NoLabToNoInd, [])
        Element.__init__(self,"BAX21EI",InzList, 3,0,3,3, (set([1, 2, 6]),set([1, 2, 6]),set([1, 2, 6])), 13,False,Label,SetLabel,2, None,None,None, NoList, NoLabToNoInd, [])
        self.TensStiff = True                                               # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X  = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT,self.nInt-1,j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
        TxDxI.__init__(self, NoList)
        self.BondLaw = BondLaw
#        self.Y0 = NoList[NoLabToNoInd[InzList[0]]].YCo                      # needed for FormB
#        self.Y1 = NoList[NoLabToNoInd[InzList[1]]].YCo
#        dx = self.X1 - self.X0
#        dy = self.Y1 - self.Y0
#        L = sqrt(dx**2 + dy**2)
#        self.sinA = dy / L
#        self.cosA = dx / L

class BAX23(B23):                                                           # axisymmetric Bernoulli Beam 2D 2 nodes, cubic shape
    # Ini2 form B23
    def __init__(self, Label, SetLabel, InzList, NoList, NoLabToNoInd):
        Element.__init__(self,"BAX23",InzList, 2, 0,2,2, (set([1, 2, 6]),   set([1, 2, 6])),  12,False,Label,SetLabel,2, None,   None, None, NoList,NoLabToNoInd,[])
        self.TensStiff = True                                               # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo                      # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X  = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):                                         # global x-coordinates of integration points
            r = SamplePoints[self.IntT,self.nInt-1,j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
class BAX23I(B23, TxDxI):                                                   # embdedded axisymmetric Bernoulli Beam 2D 2 nodes, cubic shape
    # Ini2 from B23
    def __init__(self, Label,SetLabel,InzList, BondLaw, NoList,NoLabToNoInd):
        Element.__init__(self, "BAX23I", InzList, 2,0,2,2, (set([1, 2, 6]), set([1, 2, 6])), 12, False, Label, SetLabel, 2, None, None, None, NoList, NoLabToNoInd, [])
        TxDxI.__init__(self, NoList)
        self.BondLaw = BondLaw
        self.TensStiff = True  # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo  # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[1]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):  # global x-coordinates of integration points
            r = SamplePoints[self.IntT, self.nInt - 1, j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
class BAX23EI(B23E, TxDxI):                                                   # embdedded axisymmetric Bernoulli Beam 2D 2 nodes, cubic shape
    # Ini2 from B23E
    def __init__(self, Label,SetLabel,InzList, BondLaw, NoList,NoLabToNoInd):
        #def __init__(self, TypeVal, InzList, nFieVal, IntTVal,nIntVal,nIntLVal, DofTVal,             dimVal, NLGeomIVal, Label,elSet,SpaceDimVal, MatName, StateV, NData, NoList, NoLabToNoInd, NoIndToCMInd):
#        Element.__init__(self, "BAX23EI", InzList, 2,0,3,3, (set([1, 2, 6]),set([1]),set([1, 2, 6])), 12, False, Label,SetLabel, 2, None, None, None, NoList, NoLabToNoInd, [])
        Element.__init__(self, "BAX23EI", InzList, 2,0,2,2, (set([1, 2, 6]),set([1]),set([1, 2, 6])), 12, False, Label,SetLabel, 2, None, None, None, NoList, NoLabToNoInd, [])
        #Element.__init__(self,"B23E",    InzList, 2,0,3,3, (set([1, 2, 6]),set([1]),set([1, 2, 6])), 10, True,  Label,SetLabel, 2, None, None, None, NoList,NoLabToNoInd,[])
        TxDxI.__init__(self, NoList)
        self.TensStiff = True                                               # flag for tension stiffening
        self.X0 = NoList[NoLabToNoInd[InzList[0]]].XCo  # needed for FormB
        self.X1 = NoList[NoLabToNoInd[InzList[2]]].XCo
        self.X = zeros((self.nIntL), dtype=float)
        for j in range(self.nIntL):  # global x-coordinates of integration points
            r = SamplePoints[self.IntT, self.nInt - 1, j][0]
            self.X[j] = 0.5 * ((1 - r) * self.X0 + (1 + r) * self.X1)
        self.BondLaw = BondLaw

# bond elements
class BondXDX(object):
    # redefine some element dimensioning data
    def ElementDofs(self, NodeList, NoIndToCMInd):
        self.DofT = []
        for j in self.Inzi:
            node = NodeList[j]
            self.DofT += [node.DofT]                                        # redefine list of sets/tuples for nodal dofs
        self.ElemDimData( NodeList, NoIndToCMInd )
    def Ini2(self, NodeList,NoIndToCMInd, MatList, SecDict):
        self.Rot = False                                                    # already set in Element ???
        mat = self.MatN
        Material = MatList[mat]
        NData = self.NData
        self.InitData( self.nIntL, NData, Material.StateVar, [])            # must be here as bond elements are created after post phase in InOut
        return []
    def JacoD(self, r, s, t):
        """ Dummy for jacobian determinant.
        """
        return 1
# 2D bond element corresponding to T2D2I
class Bond2D2(Element, BondXDX):
    # uses ElementDofs, LocatePointInElement from BxDxE
    def __init__(self, Label, InzList, NoList,elSet, InziC,BondElWe,IntWeight,kC,rsC,kT,kI,rsT,AA, MatName,nRebar,StateV,NData, NoLabToNoInd,NoIndToCMInd, contelem,bar):
                                                                            # kC: index continuum element; rsC: local coor cont el; kT: index of truss element;
        # ElDofs collects all degrees of freedom of element - in case of bond element these are allo dofs of connected elements
        if   bar.Type == 'T2D2I': ElDofs = [set([1,2]),  set([1,2])]
        elif bar.Type == 'B23I':  ElDofs = [set([1,2,6]),set([1,2,6])]      # rotatianal dofs includede, see FormN below -- presumably not necessary
        for _ in range(len(InziC)): ElDofs += [set([1,2])]                  # InziC is list of nodes of underlying continuum element
        #                                                          dim
        Element.__init__(self,"Bond2D2", InzList,2, 0,1,1, ElDofs, 97,False,Label,elSet,2, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA,nRebar)
        self.contelem = contelem                                            # underlying continuum element - required for transverse stress
        self.NData = NData                                                  # to make this independent from general material NData -- InitData is in BxDxE.Ini2
    # also used from B2D3E
    def IniBond(self, NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA,nRebar):
        # inziC: indices of underlying continuum nodes, other see below
        for i in InziC: self.Inzi += [i] 
        if self.bar.Type in ['B23EI','BAX23EI']: L = sqrt( (NoList[self.Inzi[2]].XCo-NoList[self.Inzi[0]].XCo)**2 + (NoList[self.Inzi[2]].YCo-NoList[self.Inzi[0]].YCo)**2 ) # element length
        else:                                    L = sqrt( (NoList[self.Inzi[1]].XCo-NoList[self.Inzi[0]].XCo)**2 + (NoList[self.Inzi[1]].YCo-NoList[self.Inzi[0]].YCo)**2 ) # element length
        self.LL = L
        self.Lch_ = 0
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 1.                                                 # dummy for Area / Jacobi determinant used instead
#        F1 = 3.544907702 * sqrt( nRebar * AA )  # F1 circumference [length] from cross section of rebar - to scale from bond stress on a circular tube [MN/m^2] -- nRebar should have a value already set in CreateBond
        F1 = 3.544907702 * nRebar * sqrt(AA)  # F1 circumference [length] from cross section of SINGLE rebar - to scale from bond stress on a circular tube [MN/m^2] -- nRebar should have a value already set in CreateBond 3.5.. is 2*sqrt(pi)
        #        to bond stress flow for bar element [MN/m], 3.544907702 is 2 sqrt(pi): U = 2 sqrt(pi) sqrt(A)
                                                                            # in case of multiple rebars U = n * 2 sqrt(pi) sqrt(A/n)
        F2 = IntWeight/2.                                                   # adjust integration weighting factor by backdoor
                                                                            # IntWeight: integration weight of corresponding T2D2I element for respective integration point
                                                                            # /2. compensates for original 0-integration order integration weight 2.
        self.Geom[1,0] = F1 * F2
        if self.bar.Type in ['B23EI','BAX23EI']:
            self.sinA = (NoList[self.Inzi[2]].YCo-NoList[self.Inzi[0]].YCo)/L
            self.cosA = (NoList[self.Inzi[2]].XCo-NoList[self.Inzi[0]].XCo)/L
        else:
            self.sinA = (NoList[self.Inzi[1]].YCo-NoList[self.Inzi[0]].YCo)/L
            self.cosA = (NoList[self.Inzi[1]].XCo-NoList[self.Inzi[0]].XCo)/L
        self.angle = (acos(self.cosA))
        if   self.bar.Type in ["T2D2I","T2D3I","TAX2I","TAX3I"]:            # transformation for 2D slip in integration point of bonded element
            self.Trans = array([[self.cosA,self.sinA    ],[-self.sinA,self.cosA]])
        elif self.bar.Type in ["B23I","B23EI","BAX23I","BAX23EI",'BAX21EI']:
            self.Trans = array([[self.cosA,self.sinA    ],[-self.sinA,self.cosA]])
        else:
            pass
        self.BondElWe = BondElWe                                            # form function values from continuum nodes
        self.iC   = kC                                                      # index of continuum element
        self.rsC  = rsC                                                     # local coordinates of bond point in continuum element
        self.iT   = kT                                                      # index of corresponding truss element in ElemList
        self.ElInt= kI                                                      # integration point on corresponding truss element
        self.rsT  = rsT                                                     # local integration point coordinate (= bond point) of truss element
    def FormX(self, r, s, t):                                               # dummy
        nc = self.nNod                                                      # nc: number of connected nodes -> includes nodes of bars and neighbor nodes of continuum
        return zeros(( 1, 1*nc), dtype=float)
    def FormX_(self, r, s, t):                                              # used e.g. for output
        r = self.rsT[0]
        nc = self.nNod                                                      # nc: number of connected nodes -> includes nodes of bars and neighbor nodes of continuum
        X = zeros(( 2, 2*nc), dtype=float)
        X[0,0] = 0.5*(1-r)                                                  # the following four correspond to the bar element nodes
        X[0,2] = 0.5*(1+r)                                                  # first row is x-component
        X[1,1] = 0.5*(1-r)                                                  # second row is y-component
        X[1,3] = 0.5*(1+r) 
        return X
    def FormN(self, r, s, t):
        r  = self.rsT[0]
        nc = self.nNod                                                      # nc: number of connected nodes, DofE: number of dofs for whole element -> includes nodes of bars and neighbor nodes of continuum
        ne = self.DofE                                                      # number of dofs for whole element including connected nodes
        NN = zeros(( 2, ne), dtype=float)
        N_ = self.BondElWe                                                  # this continuum contribution might not be a partition of unity for the respective integration point, e.g. in case of XFEM with dof types 7,8 involved
        if self.bar.Type in ["T2D2I","TAX2I"]:
            NN = zeros(( 2, ne), dtype=float)
            NN[0,0] = 0.5*(1-r)                                             # the following four correspond to the bar elements which together form a partition of unity
            NN[0,2] = 0.5*(1+r)
            NN[1,1] = 0.5*(1-r) 
            NN[1,3] = 0.5*(1+r) 
            for j in range(nc-2):                                           # this corresponds to the neighbor continuum elements
                NN[0,4+2*j]   = -N_[j]
                NN[1,4+2*j+1] = -N_[j]
        elif self.bar.Type in ["B23I","BAX23I"]:
            NN = zeros(( 2, ne), dtype=float)
            NN[0,0] = 0.5*(1-r)                                             # the following four correspond to the bar elements which together form a partition of unity
            NN[0,3] = 0.5*(1+r)                                             # omits rotational dofs to derive translational beam displacement in integration points
            NN[1,1] = 0.5*(1-r)
            NN[1,4] = 0.5*(1+r)
            for j in range(nc-2):                                           # this corresponds to the neighbor continuum elements
                NN[0,6+2*j]   = -N_[j]                                      # 5 is highest index for beam nodal dofs
                NN[1,6+2*j+1] = -N_[j]
        return NN
    def FormB(self, r, s, t, NLGeom):
        L = self.LL
        BB = dot(self.Trans,self.FormN( r, s, t))                           # transformation / rotation into local system
        return BB, L/2., 0                                                  # Jacobian for integration
    def FormT(self, r, s, t):                                               # dummy
        ne = self.DofE                                                      # number of dofs for whole element including connected nodes
        return zeros(( 2, ne), dtype=float)
    def transverse_stress(self):                                            # Ahmad
        def NodalFromIntPoint( SI ):                                        # this relates to CPE4 and strongly depends on the sequence of integration points, 
                                                                            # see ConFemBasics in connection how N is organized --> ExtraploateStressesInQuad.pdf
            sigN = zeros((4), dtype=float)
            sigN[0] =  .1339745962*SI[0] - .5000000000*SI[2] + 1.866025404*SI[3] - .5000000000*SI[1]
            sigN[1] = -.5000000000*SI[0] + .1339745962*SI[2] - .5000000000*SI[3] + 1.866025404*SI[1]
            sigN[2] =  1.866025404*SI[0] - .5000000000*SI[2] + .1339745962*SI[3] - .5000000000*SI[1]
            sigN[3] = -.5000000000*SI[0] + 1.866025404*SI[2] - .5000000000*SI[3] + .1339745962*SI[1]
            return sigN
        sigcon = zeros((3),dtype = float)
        Data = self.contelem.Data
        sxxN = NodalFromIntPoint( [Data[0][4],Data[1][4],Data[2][4],Data[3][4]] )   # sigma_xx in nodes according to sequence in Inzi
        syyN = NodalFromIntPoint( [Data[0][5],Data[1][5],Data[2][5],Data[3][5]] )   # sigma_yy in nodes
        szzN = NodalFromIntPoint( [Data[0][6],Data[1][6],Data[2][6],Data[3][6]] )   # sigma_zz in nodes - !!! requires that lateral stress is provided mat material routine   
        sxyN = NodalFromIntPoint( [Data[0][7],Data[1][7],Data[2][7],Data[3][7]] )   # sigma_xy in nodes
        sigcon[0] = dot(self.contelem.FormX(self.rsC[0], self.rsC[1], 0.),sxxN)     # rsC should be local coordinates of bond element in underlying continuum element
        sigcon[1] = dot(self.contelem.FormX(self.rsC[0], self.rsC[1], 0.),syyN)
        sigcon[2] = dot(self.contelem.FormX(self.rsC[0], self.rsC[1], 0.),sxyN)
        sig_ZZ    = dot(self.contelem.FormX(self.rsC[0], self.rsC[1], 0.),szzN)
        ang = self.angle                                                    # angle of bond/truss element
        sig_mean = 0.5*( 0.5*(sigcon[0]+sigcon[1]) - 0.5*(sigcon[0]-sigcon[1])*cos(2.*ang) - sigcon[2]*sin(2.*ang) )   #stress of continuum element at bond intpt perpendicular on truss element  
        if not self.contelem.PlSt:                                          # plane strain
            sig_mean = sig_mean + 0.5*sig_ZZ  
        fc =  self.contelem.Material.fc                                     # fc of continuum element
        fct = self.contelem.Material.fct                                    # fct of continuum element
        return sig_mean, fc, fct
# 2D bond element corresponding to T2D3I
class Bond2D3(Bond2D2):
    def __init__(self, Label, InzList, NoList,elSet, InziC,BondElWe,IntWeight,kC,rsC,kT,kI,rsT,AA, MatName,nRebar,StateV,NData, NoLabToNoInd,NoIndToCMInd, contelem,bar ):
                                                                            # kC: index continuum element; rsC: local coor cont el; kT: index of truss element;
        if   bar.Type == 'T2D3I': ElDofs = [set([1,2]),  set([1,2]),  set([1])]
        elif bar.Type == 'B23EI': ElDofs = [set([1,2,6]),set([1]),set([1,2,6])]
        for _ in range(len(InziC)): ElDofs += [set([1,2])]                  # InziC is list of nodes of underlying continuum element
        #                                                          dim
        Element.__init__(self,"Bond2D3", InzList,2, 0,1,1, ElDofs, 97,False,Label,elSet,2, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA,nRebar)
        self.contelem = contelem                                            # underlying continuum element
        self.NData = NData                                                  # to make this independent from general material NData
    def FormX(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                      # nc: number of connected nodes -> includes nodes of bars and neighbor nodes of continuum
        X = zeros(( 1, 1*nc), dtype=float)
        if self.bar.Type in ['T2D3I','BAX21EI']:
            X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
            X[0,1] = 0.5*(1+r)
        elif self.bar.Type=="B23EI":
            X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
            X[0,2] = 0.5*(1+r)
        return X
    def FormX_(self, r, s, t):                                          # used e.g. for output
        r = self.rsT[0]
        nc = self.nNod                                                  # nc: number of connected nodes -> includes nodes of bars and neighbor nodes of continuum
        X = zeros(( 2, 2*nc), dtype=float)
        if self.bar.Type in ['T2D3I','BAX21EI']:
            X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
            X[0,2] = 0.5*(1+r)                                              # first row is x-component
            X[1,1] = 0.5*(1-r)                                              # second row is y-component
            X[1,3] = 0.5*(1+r) 
        elif self.bar.Type=="B23EI":
            X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
            X[0,4] = 0.5*(1+r)                                              # first row is x-component
            X[1,1] = 0.5*(1-r)                                              # second row is y-component
            X[1,5] = 0.5*(1+r) 
        return X
    def FormN(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                      # nc: number of connected nodes, DofE: number of dofs for whole element -> includes nodes of bars and neighbor nodes of continuum
        ne = self.DofE                                                      # number of dofs for whole element including connected nodes
        NN = zeros(( 2, ne), dtype=float)
        N_ = self.BondElWe                                                  # this continuum contribution might not be a partition of unity for the respective integration point, e.g. in case of XFEM with dof types 7,8 involved
        cc = self.cosA
        ss = self.sinA
        c2 = self.cosA**2
        s2 = self.sinA**2
        cs = self.cosA*self.sinA
        r2 = r**2
        r3 = r2*r
        if self.bar.Type in ["T2D3I","TAX3I"]:                              # see Bond_TxDLQ.mws -- presumably global --> local --> global due to enhancing 1-dof
            NN[0,0] = 0.5*((-r+r2)*c2+(1.-r)*s2)
            NN[0,1] = 0.5*((-r+r2)*cs-(1.-r)*cs)
            NN[0,2] = 0.5*(( r+r2)*c2+(1.+r)*s2)
            NN[0,3] = 0.5*(( r+r2)*cs-(1.+r)*cs) 
            NN[0,4] = cc*(1.-r2) 
            NN[1,0] = 0.5*((-r+r2)*cs-(1.-r)*cs)
            NN[1,1] = 0.5*((-r+r2)*s2+(1.-r)*c2)
            NN[1,2] = 0.5*(( r+r2)*cs-(1.+r)*cs)
            NN[1,3] = 0.5*(( r+r2)*s2+(1.+r)*c2) 
            NN[1,4] = ss*(1.-r2)
            for j in range(nc-3):                                           # corresponds to the neighbor continuum elements
                NN[0,5+2*j]   = -N_[j]
                NN[1,5+2*j+1] = -N_[j]
        elif self.bar.Type == 'BAX21EI':
            NN[0,0] = (.5*r2-.5*r)
            NN[0,3] = (.5*r2+.5*r)
            NN[0,6] = (1.-r2)
            NN[1,1] = (.5*r2-.5*r)
            NN[1,4] = (.5*r2+.5*r)
            NN[1,7] = (1.-r2)
            for j in range(nc-3):                                           # corresponds to the neighbor continuum elements
                k = 9+2*j
                NN[0,k]   = -N_[j]
                NN[1,k+1] = -N_[j]
        elif self.bar.Type in ["B23EI","BAX23EI"]:                                      #
            L = self.bar.LL
            # see Bond_BeamxDLQ.mws -- transforms beam nodal degrees of freedom into local, applies shape functions and transforms local translations back into global translations
            # the following correspond to the bar elements which together form a partition of unity
            NN[0,0]=  .5*c2*r*(r-1)+s2*(.25*r3-.75*r+.50)
            NN[0,1]=  .5*cc*r*(r-1)*ss-ss*(.25*r3-.75*r+.50)*cc
            NN[0,2]= -.125*ss*L*(r3-r2-r+1) 
            NN[0,3]=   cc*(1.-r2)
            NN[0,4]=  .5*c2*r*(1+r)+s2*(-.25*r3+.75*r+.50)
            NN[0,5]=  .5*cc*r*(1+r)*ss-ss*(-.25*r3+.75*r+.50)*cc
            NN[0,6]= -.125*ss*L*(r3+r2-r-1)
            NN[1,0]=  .5*cc*r*(r-1)*ss-ss*(.25*r3-.75*r+.50)*cc
            NN[1,1]=  .5*s2*r*(r-1)+(.25*r3-.75*r+.50)*c2
            NN[1,2]=  .125*cc*L*(r3-r2-r+1)
            NN[1,3]=   ss*(1.-r2)
            NN[1,4]=  .5*cc*r*(1+r)*ss-ss*(-.25*r3+.75*r+.50)*cc
            NN[1,5]=  .5*s2*r*(1+r)+(-.25*r3+.75*r+.50)*c2
            NN[1,6]=  .125*cc*L*(r3+r2-r-1)
            for j in range(nc-3):                                           # corresponds to the neighbor continuum elements
                NN[0,7+2*j]   = -N_[j]                                      # 6 is highest index for enhanced beam nodal dofs 
                NN[1,7+2*j+1] = -N_[j]
        return NN
# axisymmetric bond element corresponding to TAX2I
class BondAX2( Bond2D2 ):
    # uses ElementDofs, LocatePointInElement from BxDxE
    def __init__(self, Label, InzList, NoList,elSet, InziC,BondElWe,IntWeight,kC,rsC,kT,kI,rsT,xx, AA, MatName,nRebar, NData, NoLabToNoInd, contelem,bar):
                                                                            # kC: index continuum element; rsC: local coor cont el; kT: index of truss element;
        if   bar.Type in ["TAX2I"]:  ElDofs = [set([1,2]),  set([1,2])]
        elif bar.Type in ["BAX23I"]: ElDofs = [set([1,2,6]),set([1,2,6])]
        else:                       raise NameError("ComFemElem::BondAX2: not allowed for element", self.Type)
        for _ in range(len(InziC)):  ElDofs += [set([1,2])]                  # InziC is list of nodes of underlying continuum element
        #                                                          dim
        Element.__init__(self,"BondAX2", InzList,2, 0,1,1, ElDofs, 96,False,Label,elSet,2, MatName,None,None, NoList,NoLabToNoInd,[]) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
#        Element.__init__(self,"BondAX2", InzList,2, 0,1,1, ElDofs, 5,False,Label,elSet,2, MatName,None,None, NoList,NoLabToNoInd,[]) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA,nRebar) # from Bond2D2
        self.contelem = contelem                                            # underlying continuum element - required for transverse stress
        self.NData = NData                                                  # to make this independent from general material NData -- InitData is in BondXDXE.Ini2
        self.X = array([xx])
    # Ini2 from BondXDX
class BondAX3( Bond2D3 ):                                                   # looks basically like BonAX2 but inherits methods of Bon2D3
    # uses ElementDofs, LocatePointInElement from BxDxE
    def __init__(self, Label, InzList, NoList,elSet, InziC,BondElWe,IntWeight,kC,rsC,kT,kI,rsT,xx, AA, MatName,nRebar, NData, NoLabToNoInd, contelem,bar):
                                                                            # kC: index continuum element; rsC: local coor cont el; kT: index of truss element;
        if bar.Type   == 'TAX3I':   ElDofs = [set([1,2]),  set([1,2]),  set([1])]
        elif bar.Type == 'BAX21EI': ElDofs = [set([1,2,6]),set([1,2,6]),set([1,2,6])]
        elif bar.Type == 'BAX23EI': ElDofs = [set([1,2,6]),set([1]),    set([1,2,6])]
        else:                       raise NameError("ComFemElem::BondAX3: not allowed for element", bar.Type)
        for _ in range(len(InziC)): ElDofs += [set([1,2])]                  # InziC is list of nodes of underlying continuum element
        #                                                          dim
        Element.__init__(self,"BondAX3", InzList,2, 0,1,1, ElDofs, 96,False,Label,elSet,2, MatName,None,None, NoList,NoLabToNoInd,[]) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA,nRebar) # from Bond2D2
        self.contelem = contelem                                            # underlying continuum element - required for transverse stress
        self.NData = NData                                                  # to make this independent from general material NData -- InitData is in BondXDXE.Ini2
        self.X = array([xx])
    # Ini2 from BondXDX
# 3D bond element corresponding to T3D2I
class Bond3D2(Element, BondXDX):
    # uses ElementDofs, LocatePointInElement from BxDxE
    def __init__(self, Label, InzList, NoList, elSet, MatName, InziC,BondElWe,IntWeight, kC,rsC,kT,kI,rsT, AA, NoLabToNoInd,StateV,NData, NoIndToCMInd, contelem,bar ):
        ElDofs = [set([1,2,3]),set([1,2,3])]
        for _ in range(len(InziC)): ElDofs += [set([1,2,3])]
        #                                                          dim
        Element.__init__(self,"Bond3D2", InzList,2, 0,1,1, ElDofs, 95,False,Label,elSet,3, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
#        Element.__init__(self,"Bond3D2", InzList,2, 0,1,1, ElDofs, 1,False,Label,elSet,3, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA)
        self.contelem = contelem                                            # underlying continuum element
        self.NData = NData                                                  # to make this independent from general material NData
    def IniBond(self, NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA):
        for i in InziC: self.Inzi += [i]                                # add node indices from continuum element where current bond element is embedded in 
        xL = NoList[self.Inzi[1]].XCo-NoList[self.Inzi[0]].XCo
        yL = NoList[self.Inzi[1]].YCo-NoList[self.Inzi[0]].YCo
        zL = NoList[self.Inzi[1]].ZCo-NoList[self.Inzi[0]].ZCo
        L  = sqrt( xL**2 + yL**2 + zL**2 )                              # element length
        vL = array([ xL/L, yL/L, zL/L ])                                # longitudinal orientation
        # create local coordinate system
        v3 = cross(vL,array([0., 0., 1.]))
        if norm(v3)<1.e-3: v3 = cross(vL,array([0., 1., 0.]))
        L3 = sqrt( v3[0]**2 + v3[1]**2 + v3[2]**2 )
        v3 = v3/L3
        v2 = cross(v3, vL)
        self.RotMat = array([[vL[0],vL[1],vL[2]],
                             [v2[0],v2[1],v2[2]],
                             [v3[0],v3[1],v3[2]]])
        #
        self.LL = L
        self.Lch_ = 0
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 1.                                             # dummy for Area / Jacobi determinant used instead
        F1 = 3.544907702 * sqrt(AA)  # circumference [length] from cross section - to scale from bond stress on a circular tube [MN/m^2] to bond stress flow for bar element [MN/m] 
        F2 = IntWeight/2.                                               # adjust integration weighting factor by backdoor
                                                                        # IntWeight: integration weight of corresponding T2D2I element for respective integration point
                                                                        # /2. compensates for original 0-integration order integration weight 2. 
        self.Geom[1,0] = F1 * F2
        self.BondElWe = BondElWe                                        # form function values from continuum nodes
        self.iC   = kC                                                  # index of continuum element
        self.rsC  = rsC                                                 # local coordinates of bond point in continuum element
        self.iT   = kT                                                  # index of corresponding truss element
        self.rsT  = rsT                                                 # local integration point coordinate (= bond point) of truss element
        self.ElInt= kI                                                  # integration point on corresponding truss element  
    def FormX(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                  # nc: number of connected nodes -> includes nodes of bars and neighbor nodes of continuum
        X = zeros(( 1, 1*nc), dtype=float)
        X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
        X[0,1] = 0.5*(1+r)
        return X
    def FormX_(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                  # nc: number of connected nodes
        X = zeros(( 3, 3*nc), dtype=float)
        X[0,0] = 0.5*(1-r)                                              # the following four correspond to the bar element nodes
        X[1,1] = 0.5*(1-r) 
        X[2,2] = 0.5*(1-r) 
        X[0,3] = 0.5*(1+r)
        X[1,4] = 0.5*(1+r) 
        X[2,5] = 0.5*(1+r) 
        return X
    def FormN(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                  # nc: number of connected nodes
        ne = self.DofE                                                  # number of dofs for whole element including connected nodes
        NN = zeros(( 3, ne), dtype=float)
        N_ = self.BondElWe                                              # this continuum contribution might not be a partition of unity for the respective integration point, e.g. in case of XFEM with dof types 7,8 involved
        NN[0,0] = 0.5*(1-r)                                             # the following four correspond to the bar element nodes
        NN[1,1] = 0.5*(1-r) 
        NN[2,2] = 0.5*(1-r) 
        NN[0,3] = 0.5*(1+r)
        NN[1,4] = 0.5*(1+r) 
        NN[2,5] = 0.5*(1+r) 
#        for j in range(nc-3):                                          # this corresponds to the neighbor continuum elements
        for j in range(nc-2):                                          # this corresponds to the neighbor continuum elements
            NN[0,6+3*j]   = -N_[j]
            NN[1,6+3*j+1] = -N_[j]
            NN[2,6+3*j+2] = -N_[j]
        return NN
    def FormB(self, r, s, t, val):
        L = self.LL
        BB = dot(self.RotMat,self.FormN( r, s, t))
        return BB, L/2., 0
    def FormT(self, r, s, t):                                               # dummy !!!
        ne = self.DofE                                                      # number of dofs for whole element including connected nodes
        return zeros(( 3, ne), dtype=float)
    def transverse_stress(self): # Ahmad
        Bcon, _ = self.contelem.FormB(self.rsC[0], self.rsC[1], self.rsC[2], None)   #B matrix of continuum element at bond intpt
        epscon = dot(Bcon, self.contelem.Nodisp)                     #strain of continuum element at bond intpt
        Cb = self.contelem.Cbond                          
        sigcon = dot(Cb, epscon)           #stress of continuum element at bond intpt
        sigcon_ = array([[sigcon[0], sigcon[5], sigcon[4]],
                         [sigcon[5], sigcon[1], sigcon[3]],
                         [sigcon[4], sigcon[3], sigcon[2]]])
        sig_rot =  dot(self.Rot,dot(sigcon_,transpose(self.Rot)))
        sig_mean = (1/2)*sig_rot[1,1] + (1/2)*sig_rot[2,2] + (2/pi)*sig_rot[1,2]   #stress of continuum element at bond intpt perpendicular on truss element  
        fc =  self.contelem.Material.fc            # fc of continuum element
        fct = self.contelem.Material.fct           # fct of continuum element
        return sig_mean, fc, fct
# 3D bond element corresponding to T2D3I
class Bond3D3(Bond3D2):
    def __init__(self, Label, InzList, NoList, elSet, MatName, InziC,BondElWe,IntWeight, kC,rsC,kT,kI,rsT, AA, NoLabToNoInd,StateV,NData, NoIndToCMInd, contelem,bar ):
        ElDofs = [set([1,2,3]),set([1,2,3]),set([1])]
        for _ in range(len(InziC)): ElDofs += [set([1,2,3])]
        #                                                          dim
        Element.__init__(self,"Bond3D3", InzList,2, 0,1,1, ElDofs, 95,False,Label,elSet,3, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
#        Element.__init__(self,"Bond3D3", InzList,2, 0,1,1, ElDofs, 1,False,Label,elSet,3, MatName,StateV,NData, NoList,NoLabToNoInd,NoIndToCMInd) # integration order 1 is mandatory as every truss integration point leads to an B2D2E element
        self.bar = bar
        self.IniBond( NoList, InziC, BondElWe, IntWeight, kC, rsC, kT,kI, rsT, AA)
        self.contelem = contelem                                            # underlying continuum element
        self.NData = NData                                                  # to make this independent from general material NData
    def FormN(self, r, s, t):
        r = self.rsT[0]
        nc = self.nNod                                                  # nc: number of connected nodes, DofE: number of dofs for whole element -> includes nodes of bars and neighbor nodes of continuum
        ne = self.DofE                                                  # number of dofs for whole element including connected nodes
        NN = zeros(( 3, ne), dtype=float)
        N_ = self.BondElWe                                              # this continuum contribution might not be a partition of unity for the respective integration point, e.g. in case of XFEM with dof types 7,8 involved
        g00 = self.RotMat[0,0]                                             # corresponds to FormN from T3D3I
        g01 = self.RotMat[0,1]
        g02 = self.RotMat[0,2]
        g10 = self.RotMat[1,0]
        g11 = self.RotMat[1,1]
        g12 = self.RotMat[1,2]
        g20 = self.RotMat[2,0]
        g21 = self.RotMat[2,1]
        g22 = self.RotMat[2,2]
        r2 = r**2
        NN[0,0] = 0.5*((-r+r2)*g00*g00 + (1.-r)*g10*g10 + (1.-r)*g20*g20)
        NN[0,1] = 0.5*((-r+r2)*g00*g01 + (1.-r)*g10*g11 + (1.-r)*g20*g21)
        NN[0,2] = 0.5*((-r+r2)*g00*g02 + (1.-r)*g10*g12 + (1.-r)*g20*g22)
        NN[0,3] = 0.5*(( r+r2)*g00*g00 + (1.+r)*g10*g10 + (1.+r)*g20*g20)
        NN[0,4] = 0.5*(( r+r2)*g00*g01 + (1.+r)*g10*g11 + (1.+r)*g20*g21)
        NN[0,5] = 0.5*(( r+r2)*g00*g02 + (1.+r)*g10*g12 + (1.+r)*g20*g22)
        NN[0,6] = g00*(1-r2)
        NN[1,0] = 0.5*((-r+r2)*g00*g01 + (1.-r)*g11*g10 + (1.-r)*g20*g21)
        NN[1,1] = 0.5*((-r+r2)*g01*g01 + (1.-r)*g11*g11 + (1.-r)*g21*g21)
        NN[1,2] = 0.5*((-r+r2)*g02*g01 + (1.-r)*g11*g12 + (1.-r)*g22*g21)
        NN[1,3] = 0.5*(( r+r2)*g00*g01 + (1.+r)*g11*g10 + (1.+r)*g20*g21)
        NN[1,4] = 0.5*(( r+r2)*g01*g01 + (1.+r)*g11*g11 + (1.+r)*g21*g21)
        NN[1,5] = 0.5*(( r+r2)*g02*g01 + (1.+r)*g11*g12 + (1.+r)*g22*g21)
        NN[1,6] = g01*(1-r2)
        NN[2,0] = 0.5*((-r+r2)*g00*g02 + (1.-r)*g12*g10 + (1.-r)*g20*g22)
        NN[2,1] = 0.5*((-r+r2)*g01*g02 + (1.-r)*g12*g11 + (1.-r)*g21*g22)
        NN[2,2] = 0.5*((-r+r2)*g02*g02 + (1.-r)*g12*g12 + (1.-r)*g22*g22)
        NN[2,3] = 0.5*(( r+r2)*g00*g02 + (1.+r)*g12*g10 + (1.+r)*g20*g22)
        NN[2,4] = 0.5*(( r+r2)*g01*g02 + (1.+r)*g12*g11 + (1.+r)*g21*g22)
        NN[2,5] = 0.5*(( r+r2)*g02*g02 + (1.+r)*g12*g12 + (1.+r)*g22*g22)
        NN[2,6] = g02*(1-r2)
        for j in range(nc-3):                                          # this corresponds to the neighbor continuum elements
            NN[0,7+3*j]   = -N_[j]
            NN[1,7+3*j+1] = -N_[j]
            NN[2,7+3*j+2] = -N_[j]
        return NN

# 2D continuum elements
class ElementC2D(Element):
    def __init__(self):
        pass
    def LocatePointInElement(self, NodeList, bb, NoIndToCMInd):             # bb is target
        Finis, R, inSC, nen = False, [], [], 0
        rs = array([ 0., 0.])                                               # initial values of isoparametric coordinates
        niter = 5
        if self.Type in ["CPS4","CPS4S","CPE4","CPE4S","CPS3","CPS3S","CPE3","CPE3S","CAX4"]:
            xy = self.NodalCoordinates( NodeList, NoIndToCMInd )            # nodal coordinates
            for _ in range(niter):                                          # iteration loop to determine isoparametric coordinates / R-values of boundary control value by Newton Raphson
                r, s = rs[0], rs[1]
                XX  = self.FormX_( r, s, 0.)
                xyP = dot(XX,array(xy))                                     # undeformed global coordinates of integration point
                rxy = array([bb[0],bb[1]])-xyP                              # residuum of target coordinate to iterated coordinates
                rno = sqrt(rxy[0]**2+rxy[1]**2)
                if rno<1.e-6: 
                    Finis = True
                    break
                if self.Type in ["CPS3","CPE3"]:
                    x0, y0 = xy[0], xy[1]
                    x1, y1 = xy[2], xy[3]
                    x2, y2 = xy[4], xy[5]
                    xp = bb[0]
                    yp = bb[1]
                    deno = (-y1*x0+y1*x2+y2*x0+x1*y0-x1*y2-x2*y0)
                    rs[0] = -(x1*y2-y1*x2+y1*xp-y2*xp-x1*yp+x2*yp)/deno #(-y1*x0+y1*x2+y2*x0+x1*y0-x1*y2-x2*y0) 
                    rs[1] = (-y2*xp-x2*y0+y0*xp-yp*x0+x2*yp+y2*x0)/deno #(-y1*x0+y1*x2+y2*x0+x1*y0-x1*y2-x2*y0)                    
                else:
                    dxy = array([[  .25*(1+s)*xy[0]-.25*(1+s)*xy[2]-.25*(1-s)*xy[4]+.25*(1-s)*xy[6] ,  .25*(1+r)*xy[0]+.25*(1-r)*xy[2]-.25*(1-r)*xy[4]-.25*(1+r)*xy[6] ],
                                 [  .25*(1+s)*xy[1]-.25*(1+s)*xy[3]-.25*(1-s)*xy[5]+.25*(1-s)*xy[7] ,  .25*(1+r)*xy[1]+.25*(1-r)*xy[3]-.25*(1-r)*xy[5]-.25*(1+r)*xy[7] ]])
                    dxyI  = inv(dxy)
                    rs = rs + dot(dxyI,rxy)                                 # improved value of iso-parametric coordinates
                if rs[0]> 1.: rs[0] =  1.                                   # constraints of iso-parametric coordinates
                if rs[0]<-1.: rs[0] = -1.                                   # --> iteration will not converge if bb is outside
                if rs[1]> 1.: rs[1] =  1.
                if rs[1]<-1.: rs[1] = -1.
            if Finis:
                NN   = self.FormN( rs[0], rs[1], 0.)
                inSC = self.Inzi
                nen  = len(inSC)
                R = zeros((nen),dtype=float)
                for j_ in range(nen): R[j_] = NN[0,2*j_]
            return Finis, rs, R, inSC
        else: raise NameError("CaeFemElements::ElementC2D:LocatePointInElement: wrong element type")
    def CreateSDARankine(self, nea, MatList, NodeList, ElemList, SolSecs, NoLabToNoInd,NoIndToCMInd, ff, uu, du, TimeDelt):
        nea_= nea
        Mat = self.Material
        MatName = self.MatN 
        if Mat.RType==3:                                                            # seems to be redundant 
            if Mat.Type in ['ISODAMAGE','ELASTICLT_SDA',"MicroPl"]:
                if   Mat.Type=='ISODAMAGE':     offset = 3 + 6
                elif Mat.Type=='ELASTICLT_SDA': offset =     6
                elif Mat.Type=="MicroPl":    offset = self.Material.iS + 6          # iS = 21 + 5 = 26 -> +6  -> 32
                eps1M, nnM, ttM = self.Means( offset )                              # mean of all integration points
#                eps1M, nnM, ttM = self.Max( offset )                                # integration point with maximum     
                ln = sqrt(nnM[0]**2+nnM[1]**2)
                if ln>ZeroD: 
                    nnM[0] = nnM[0]/ln                                              # nnM should be normal --> unit normal
                    nnM[1] = nnM[1]/ln
                if Mat.PrinStrains: Flag = eps1M>Mat.eps_ct                         # eps1M written in self.Means from respective state variables
                else:               Flag = eps1M>Mat.RankineFactor*Mat.fct          # Rankine factor for, e.g., Isodam where softening occurs before uniaxial tensile strength
                if Flag:
                    xy = self.NodalCoordinates( NodeList, NoIndToCMInd )
                    if self.Type in ["CPS4","CPE4"]: 
                        XX = self.FormX_( 0., 0., 0.)
                    if self.Type in ["CPS3","CPE3","CPS6","CPE6"]:
                        XX = self.FormX_( 1./3., 1./3., 1./3.)
                        LS = [sqrt((xy[2]-xy[0])**2+(xy[3]-xy[1])**2),sqrt((xy[4]-xy[2])**2+(xy[5]-xy[3])**2),sqrt((xy[0]-xy[4])**2+(xy[1]-xy[5])**2)]                                                     # for length of edges
#                        for i in range(int(len(xy)/2)): LS += [sqrt((xy[2*i]-xy[2*i-2])**2+(xy[2*i+1]-xy[2*i-1])**2)]
                        qLS = max(LS)/min(LS)
                    else: qLS = -1                                                  # for what is this used ???  -- for output only to control proper shape
                    xyP= dot(XX,array(xy))                                          # undeformed global coordinates of crack point / element center
                    xyC= self.FindGlobalDiscontinuityEdges( xyP, nnM, NodeList, NoIndToCMInd)# global coordinates of discontinuity end point
                    dxC= xyC[0][0]-xyC[1][0]
                    dyC= xyC[0][1]-xyC[1][1]
                    dC = sqrt(dxC*dxC + dyC*dyC)                                    # looks like length of crack
                    if abs(nnM[0]*dxC + nnM[1]*dyC)>ZeroD: 
                        raise NameError("ConFem:CPS4.CreateSDARankine: cut end not perp to normal "+str(self.Label)+str(xyP)+str(nnM)+str(xyC))
                    #
                    NewLabel = Element.LastLabel+1
                    if self.Type in ["CPS4","CPE4"]:
                        l0 = NodeList[ NoIndToCMInd[self.Inzi[0]] ].Label           # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l1 = NodeList[ NoIndToCMInd[self.Inzi[1]] ].Label
                        l2 = NodeList[ NoIndToCMInd[self.Inzi[2]] ].Label
                        l3 = NodeList[ NoIndToCMInd[self.Inzi[3]] ].Label
                        RedInt = Mat.RedInt
                        ElemList += [ CPE4_SDA( NewLabel, [l0,l1,l2,l3],      self.PlSt, NodeList, [xyP[0],xyP[1],0.], [nnM[0],nnM[1],0.], SolSecs[self.Set], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, Mat, RedInt)]
                    elif self.Type in ["CPS3","CPE3"]:
                        l0 = NodeList[ NoIndToCMInd[self.Inzi[0]] ].Label           # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l1 = NodeList[ NoIndToCMInd[self.Inzi[1]] ].Label
                        l2 = NodeList[ NoIndToCMInd[self.Inzi[2]] ].Label
#                        ElemList += [CPE3_SDA( self.Label, [l0,l1,l2],         self.PlSt, NodeList, [xyP[0],xyP[1],0.], [nnM[0],nnM[1],0.], SolSecs[self.Set], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, Mat)]
                        ElemList += [CPE3_SDA( NewLabel, [l0,l1,l2],         self.PlSt, NodeList, [xyP[0],xyP[1],0.], [nnM[0],nnM[1],0.], SolSecs[self.Set], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, Mat)]
                    elif self.Type in ["CPS6","CPE6"]:
                        l0 = NodeList[ NoIndToCMInd[self.Inzi[0]] ].Label           # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l1 = NodeList[ NoIndToCMInd[self.Inzi[1]] ].Label
                        l2 = NodeList[ NoIndToCMInd[self.Inzi[2]] ].Label
                        l3 = NodeList[ NoIndToCMInd[self.Inzi[3]] ].Label           # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l4 = NodeList[ NoIndToCMInd[self.Inzi[4]] ].Label
                        l5 = NodeList[ NoIndToCMInd[self.Inzi[5]] ].Label
#                        ElemList += [CPE6_SDA( self.Label, [l0,l1,l2,l3,l4,l5],self.PlSt, NodeList, [xyP[0],xyP[1],0.], [nnM[0],nnM[1],0.], SolSecs[self.Set], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, Mat)]
                        ElemList += [CPE6_SDA( NewLabel, [l0,l1,l2,l3,l4,l5],self.PlSt, NodeList, [xyP[0],xyP[1],0.], [nnM[0],nnM[1],0.], SolSecs[self.Set], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, Mat)]
                    else: raise NameError("CaeFemElements:CPS4.CreateSDARankine: unknown element type")
                    # some discontinuity data for new element -- 1st index indicates crack
                    NewEl = ElemList[-1]
                    NewEl.xyC = xyC
                    NewEl.CrC[0,0] = xyC[0][0]
                    NewEl.CrC[0,1] = xyC[0][1]
                    NewEl.CrC[0,2] = xyC[1][0]
                    NewEl.CrC[0,3] = xyC[1][1]
                    NewEl.CrC[0,4] = dC
                    NewEl.CrTT[0,0,0] = nnM[1]                                      # transform to discontinuity orientation, discontinuity normal (alpha) and discontinuity orientation (alpha-Pi/2) are perpendicular
                    NewEl.CrTT[0,0,1] =-nnM[0]                                      # cos(alpha-Pi/2) = sin(alpha)
                    NewEl.CrTT[0,1,0] = nnM[0]                                      # sin(alpha-Pi/2) = -cos(alpha)
                    NewEl.CrTT[0,1,1] = nnM[1]                                      #
                    ttLM = dot(NewEl.CrTT[0],ttM[0:2])                              # transform global tractions to local tractions in discontinuity system - currently for control purposes only
                    ttLM_= dot(ttM,nnM)                                             # to eliminate negative values of directions / tractions - basically used for the following
        
                    NewEl.CrTL[0,0] = abs(ttLM[0])                                  # non-zero value presumably comes from viscous stress contributions
                    NewEl.CrTL[0,1] = ttLM_ # ttLM[1]
                    NewEl.wLim[0]   = 3.*NewEl.CrackEnergy/ttLM_                    # critical crack width where normal traction becomes zero derived from crack energy and actual tensile limit
                    
                    NewEl.Material = Mat
                    NewEl.MatN = MatName
                    NewEl.CrBScaleType()                                                 # only for crack band regularization but may be later needed for a combination
                    NewEl.ElemDimData( NodeList,NoIndToCMInd )
                    NewEl.InitData( NewEl.nIntL, Mat.NData, Mat.StateVar, [])
                    # complete data of new element by taking over from old element
                    for i in range(NewEl.DofI.shape[0]):
                        for j in range(NewEl.DofI.shape[1]): 
                            NewEl.DofI[i,j] = self.DofI[i,j]
                    if NewEl.Type in ["CPS3S","CPE3S"]:
                        for i in range(1):
                            for j in range(NewEl.StateVar.shape[1]):
                                NewEl.StateVar[i,j]  = self.StateVar[0,j]
                                NewEl.StateVarN[i,j] = self.StateVarN[0,j]
                    else:
                        nIP = self.StateVar.shape[0]                                # number of integration points of current element
                        nSV = self.StateVar.shape[1]                                # number of state variables per integration point of current element
#                        if NewEl.StateVar.shape[0]==nIP:                            # full integration for CPS4_SDA
                        if NewEl.nIntL==nIP:                                        # full integration for CPS4_SDA
                            for i in range(nIP):                                    # loop over integration points
                                for j in range(nSV):                                # loop over number of state variables
                                    NewEl.StateVar[i,j]  = self.StateVar[i,j]
                                    NewEl.StateVarN[i,j] = self.StateVarN[i,j]
                        else:                                                       # should be CPS4R_SDA
                            x = zeros((nSV), dtype=float)
                            for i in range(nIP):                                    # loop over integration points
                                for j in range(nSV):
                                    x[j] += self.StateVar[i,j]
                            for i in range(nSV):                                    # loop over number of state variables
                                x[i] = x[i]/nIP                                     # mean of current element
                            for j in range(nSV):
                                NewEl.StateVar[0,j] = x[j]                          # one integration point for new element
                                
#                    NewEl.ElemIntPointCoordinates(NodeList, NoIndToCMInd)
                    SolSecs[self.Set].Elems += [nea_]                                # add to SolSecs
                    nea_ += 1
                    # initial fw, fd
                    NewEl.CrGeRa[0] = 1.                                            # used for scaling between fd and fw
                    _, fd = NewEl.CrComp(TimeDelt, ff)                              # nodal forces from crack tractions
                    fw = NewEl.CrFw( uu, du, MatList, TimeDelt)                     # nodal forces from internal state of stress; uu, uo are currently identical due to position of CreateSDARankine in CaeFemMain
                    if norm(fw)>ZeroD: NewEl.CrGeRa[0] = norm(fw)/norm(fd)
                    else: raise NameError("CaeFemElements:CPS4.CreateSDARankine: no traction upon SDA initialization")
                    Echo(f'create SDA {self.Label:6d}, xy {xyP[0]:8.3f}{xyP[1]:8.3f}, LL{dC:8.3f}, nn {nnM[0]:6.3f}{nnM[1]:6.3f}, ts {NewEl.CrTL[0,0]:6.3f}, tn{NewEl.CrTL[0,1]:6.3f}, CrGeRa{NewEl.CrGeRa[0]:6.3f}, qLS{qLS:6.3f}',ff)
                    self.Active = False
#                    self.Label  = -self.Label
            else: 
                raise NameError("CaeFemElements:CPS4.CreateSDARankine: regularization type RType 3 not yet implemented for this material",Mat.Type)
        return nea_
    # means of crack tractions, crack normal, largest principal strain
    def Means( self, offset ):
        ttM, nnM, eps1M, XCounter = zeros((3), dtype= float),zeros((3), dtype= float), 0., 0  
        for i in range(self.nIntL):                                             # loop over integration sampling points
            if self.StateVar[i,offset+0]>ZeroD:
                eps1M  += self.StateVar[i,offset+0]                             # might also be largest principal stress
                nnM[0] += self.StateVar[i,offset+1]
                nnM[1] += self.StateVar[i,offset+2]
                nnM[2] += self.StateVar[i,offset+3]
                ttM[0] += self.StateVar[i,offset+4]
                ttM[1] += self.StateVar[i,offset+5]
                ttM[2] += self.StateVar[i,offset+6]
                XCounter += 1
        if XCounter>0:                                                          # means
            eps1M  = eps1M /XCounter 
            nnM[:] = nnM[:]/XCounter
            ttM[:] = ttM[:]/XCounter
#            ttM[:] = 0.9*ttM[:]/XCounter                                        # heuristic scaling down regarding fan spreading of cracked elements 
        return eps1M, nnM, ttM
    # maximum of crack tractions, crack normal, largest principal strain
    def Max( self, offset ):
        ttM, nnM, eps1M = zeros((3), dtype= float),zeros((3), dtype= float), 0.  
        for i in range(self.nIntL):                                             # loop over integration sampling points
            if self.StateVar[i,offset+0]>eps1M:
                eps1M  = self.StateVar[i,offset+0]                              # might also be largest principal stress
                nnM[0] = self.StateVar[i,offset+1]
                nnM[1] = self.StateVar[i,offset+2]
                nnM[2] = self.StateVar[i,offset+3]
                ttM[0] = self.StateVar[i,offset+4]
                ttM[1] = self.StateVar[i,offset+5]
                ttM[2] = self.StateVar[i,offset+6]
        return eps1M, nnM, ttM
    # see SDA_Quad4_DiscEdgePoints.mws        
    def FindGlobalDiscontinuityEdges( self, xy, nn, NodeList, NoIndToCMInd):
        # direction of discontinuity - determine global direction of discontinuity as normal to discontinuity normal
        if nn[0]>0.1: 
            dirDx = -nn[1]/nn[0]
            dirDy = 1.
        else: 
            dirDx = 1.
            dirDy = -nn[0]/nn[1]
        # loop over all element edges
        xyC = []
        if self.Type in ["CPS6","CPE6","CPS6S","CPE6S"]: Inzi_ = self.Inzi[0:3]
        else:                            Inzi_ = self.Inzi
        for i in range(len(Inzi_)):
            j0 = Inzi_[i-1]
            j1 = Inzi_[i] 
            node0 = NodeList[NoIndToCMInd[j0]]
            node1 = NodeList[NoIndToCMInd[j1]]
            xE0 = node0.XCo                                                 # coordinates of element edge end points
            yE0 = node0.YCo
            xE1 = node1.XCo 
            yE1 = node1.YCo
            dirEx = xE1 - xE0                                               # direction of edge
            dirEy = yE1 - yE0
            det = dirDy*dirEx-dirDx*dirEy
            if abs(det) > ZeroD:                                            # needs no else, ruled by length of xyC 
                s = (dirDy*xy[0]+dirDx*yE0-dirDy*xE0-dirDx*xy[1])/det       # s intersection coordinate along edge length
                Tol_ = 1.0e-12
                if 0.-Tol_<=s and s<=1.+Tol_:                                         # coordinate s within element
                    xC = xE0 + s*dirEx
                    yC = yE0 + s*dirEy
                    xyC += [array([xC,yC])]
        # check for duplicates
        for i in range(len(xyC)):
            for j in range(i+1,len(xyC)):
                try:
                    L_ = xyC[i]-xyC[j]
                    LL = sqrt(L_[0]**2+L_[1]**2)
                    if LL<100.*Tol_: xyC[j] = []
                except:
                    pass
#                if xyC[i]==xyC[j]: xyC[j] = []
        xyC_ = []
        for i in xyC:
            if i != []: xyC_ += [i]
        # more than two non duplicates should not make sense
        if len(xyC_)!=2:
            print(self.Label,Inzi_,'_',xy,nn,[dirDx,dirDy],'_',xyC_) 
            raise NameError("CaeFemElements:CPS4.CreateSDARankine: Error 2")
        return xyC_                                                         # returns [[x, y of 1st edge],[x, y, of 2nd edge]]

class ElementC2D_SDA():
    def __init__(self, Coor, Normal, Material):                             # performed in ElementC2D.CreateSDARankine
        self.CrXCo = Coor[0]
        self.CrYCo = Coor[1]
        self.CrZCo = Coor[2]
        lN = sqrt(Normal[0]**2+Normal[1]**2+Normal[2]**2)                               # length of crack normal        
        self.CrN = array([[Normal[0]/lN, Normal[1]/lN, Normal[2]/lN],                   # for crack normals to be completed in CreateSDA2ndCrack
                          [0.,           0.,           0.]])
        self.CrC = zeros((2,5), dtype = float)                                          # for end points of crack [x0,y0,x1,y1, length], to be filled in CreateSDARnnkine/CreateSDA2ndCrack  
        self.Mw    = self.FormMw( Coor[0], Coor[1]) 
        if self.SDANew: self.Bw = self.FormBw( 0., 0., 0.)                                  # dummy coordinates
        self.ww    = zeros((6), dtype=float)                                            # crack width dofs
        self.dw    = zeros((6), dtype=float)                                            # crack width dofs - step increment        
        self.wwOld = zeros((6), dtype=float)
        self.fwd   = zeros((6), dtype=float)                                            # storage local forces of crack width dofs
        self.Dd    = zeros((6,6), dtype=float)                                          # for discontinuity dofs from discontinuity imbalance
        self.CrTL  = zeros((2,2), dtype=float)                                          # local crack traction upon crack initiation
        self.CrTT  = zeros((2,2,2), dtype=float)                                        # coordinate transformation matrix -- for 2 cracks from global to local with local x along crack direction
        self.nIntLCr= 2                                                                 # for integration along crack contour 
        self.IntTCr= 0
        self.nIntCr= 1
        self.wwL   = zeros((2,self.nIntLCr,2), dtype=float)                             # crack width in local system -- x local direction in crack contour direction
        self.wwT   = zeros((2,self.nIntLCr,2), dtype=float)                             # crack traction in local system -- "
        self.wwLM  = zeros((2,self.nIntLCr,2), dtype=float)                             # crack width in local system -- maximum during load history
        self.wwTM  = zeros((2,self.nIntLCr,2), dtype=float)                             # crack traction in local system -- "
        self.wwV   = zeros((2,self.nIntLCr,2), dtype=float)                             # crack velocities in local system -- for viscous contributions
        self.wwVN  = zeros((2,self.nIntLCr,2), dtype=float)                             # final crack velocities in local system -- for viscous contributions
        self.CrackStatus = zeros((2,self.nIntLCr), dtype=int)                           # 0: loading; 1: unloading, reloading in tension; 2: unloading, reloading with crack closure
        self.CrackStatusOld = zeros((2,self.nIntLCr), dtype=int)
        self.wLim  = zeros((2), dtype = float)
        self.CrGeRa = zeros((2), dtype = float)                                         # used for scaling between fd and fw
        self.RotCrack = False 
        if Material.RType==3:
            self.CrackTraction = True
            if Material.Type in ['ISODAMAGE','ELASTICLT_SDA',"MicroPl" ]:
                self.CrackEnergy = Material.RegPar
                self.wLim[0], self.wLim[1]  = 0.1, 0.1                                  # default value -- may be later modified to comply with crack energy
                self.BulkEmod = Material.Emod
                self.CrackEta = Material.etaCr
                self.RotCrack = Material.RotCrack                                       # for rotating crack / discontinuity
                self.ShearRetFactor = Material.ShearRetFac                              # shear retention factor
            else:     raise NameError("CaeFemElements:CPS3/4_SDA.__init__: actual material type not yet implemented")
        else: self.CrackTraction = False
        #
        self.NoCracks = 1
    #
    def CrIni(self, eqiter, witer, xe, ff):                                                 # xe element displacement increment from iteration
        Kww,fw = zeros((6,6),dtype=float), zeros((6),dtype=float) 
        if   self.Type in ["CPS4S","CPE4S"]: Kuw,Kwu = zeros((8,6),dtype=float), zeros((6,8),dtype=float) 
        elif self.Type in ["CPS3S","CPE3S"]: Kuw,Kwu = zeros((6,6),dtype=float), zeros((6,6),dtype=float)
        elif self.Type in ["CPS6S","CPE6S"]: Kuw,Kwu = zeros((12,6),dtype=float), zeros((6,12),dtype=float)
        elif self.Type=="C3D8S":          Kuw,Kwu = zeros((24,6),dtype=float), zeros((6,24),dtype=float) 
        if eqiter>0:
            if witer==0: ddw = dot( self.Dd, (dot(self.Kwu,array(xe))+self.fwd) )# equilibrium iteration increment of discontinuity dofs -- xe = 0 for eqiter = 0
            else:        ddw = dot( self.Dd,                          self.fwd  )# " but disregard nodal displacement contribution as xe is not updated from super ordinated equilibrium iteration
        else:            ddw = zeros((6),dtype=float)   
        self.ww[:] = self.ww[:] + ddw[:] 
        self.dw[:] = self.ww[:] - self.wwOld[:]                             # crack width step increment
#        print(' dw',ddw[0:3],'\n ww',self.ww[0:3], file=ff)             # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return fw, Kww, Kuw, Kwu
    
    def CrFw(self, uu, du, MatList, TimeDelt):               # uo -> du                           # uu, uo are currently identical due to position of CreateSDARankine in CaeFemMain
        # following should be the same as in IntForces
        ue = zeros((self.DofE), dtype=float)                                            # initialization element nodal displacement vector
        de = zeros((self.DofE), dtype=float)                                            # initialization element nodal displacement load / time step increment vector
        ndof0 = 0
        for j in range(self.nNod):                                                      # loop over nodes of element
            for k in range(self.DofN[j]):                                               # loop over dofs of element node
                kk = self.DofI[j,k]                                                     # global row index of dof
                ue[ndof0+k] = uu[kk]                                                    # local displacement
                de[ndof0+k] = du[kk]                                                    # local displacement step increment
            ndof0 += self.DofN[j]
        #
        fw = zeros((6),dtype=float)                                                     
        for j in range(self.nIntL):
            rst = SamplePoints[self.IntT,self.nInt-1,j]                                 # local integration sampling point coordinates
            BB, JJ, _ = self.FormB( rst[0], rst[1], rst[2], None )
            eps = dot(BB,ue)                                                            # integration point strains
            dps = dot(BB,de)                                                            # integration point strain load / time step increments
            if not self.SDANew: BBw_ = self.FormBw( rst[0], rst[1], rst[2] )            # Mw comes here into play
            else:               BBw_ = self.Bw[j]
            eps = eps - dot(BBw_,self.ww)
            dps = dps - dot(BBw_,self.dw)
            sig, _, _ = self.Material.Sig( None, 1, TimeDelt, None, j, self, dps, eps, [], [], [])
            f   = self.Geom[1,0]*JJ*self.Geom[0,0]*SampleWeight[self.IntT,self.nInt-1,j]               # weighting factor - area contribution is in JJ (Geom[0,0]=1) or Geom[0,0] (JJ=1) in case of one poin integration
            fw  = fw  + f*dot(transpose(BBw_),sig)                                      # fw: integrated discontinuity forces from bulk stresses
#            
#            if self.Label==1: print('YYY',self.Label,j,f,sig,'\n',BBw_)
#            
#            
#            print(f'YYY,  {j:d}, sig {", ".join(z for z in [f"{x:f}" for x in sig[0:3]])}, f {f:f}, fw {", ".join(z for z in [f"{x:f}" for x in fw[0:6]])}')#
#
        return fw
    #
    def CrComp(self, TimeDelt, ff):
        Kdd,fd = zeros((6,6),dtype=float), zeros((6),dtype=float)
        for i in range(len(self.CrC)):                                                  # crack geometry data for two cracks
            CrXY = self.CrC[i]                                                          # end coordinates and length of crack [x0,y0,x1,y1, length]
            if CrXY[4]>0.:                                                              # crack is there 
                for j in range(self.nIntLCr):                                           # loop over integration points along embedded discontinuity
                    rst = SamplePoints[self.IntTCr,self.nIntCr,j]                       # local integration sampling point coordinates
                    XCut= self.FormXCut(rst[0],rst[1],0.)                               # shape function of embedded discontinuity geometry
                    xyI = dot(XCut,CrXY[0:4])                                           # global coordinates of integration points - derived from global embedded discontinuity end points
                    Nw  = self.FormNw( self.CrXCo, self.CrYCo, xyI[0], xyI[1])          # shape function for embedded discontinuity
                    ww  = dot(Nw,self.ww[i*3:(i+1)*3])
                    dw  = dot(Nw,self.dw[i*3:(i+1)*3])
                    wwL = dot(self.CrTT[i],ww) 
                    dwL = dot(self.CrTT[i],dw) 
                    self.CrStatus( i, j, wwL)
                    tr, TW = self.CrTraction( i, j, wwL, ff)
#                    
#                    print('YYY',i,j,self.CrackStatus[i,j],wwL,tr)
#                    
                    # viscous contribution to avoid roughness
                    if self.CrackEta>0.:
                        zz, Vw = self.ViscExtenCr2D( TimeDelt, self.CrackEta, dwL, 13, i, j)  # both corresponding state variable values are currently assigned to integration point 0
# shear viscosity currently disregarded 
#                        tr[0] = tr[0] + self.CrackEta*Vw[0]
#                        TW[0,0] = TW[0,0] + zz 
                        tr[1] = tr[1] + self.CrackEta*Vw[1]
                        TW[1,1] = TW[1,1] + zz 
                    self.wwL[i,j,:] = wwL[:]                                            # store for output
                    self.wwT[i,j,:] = tr[:]                                             # "
                    TW_ = dot(transpose(self.CrTT[i]),dot(TW,self.CrTT[i]))             # transformation of local tangential stiffness into global
                    f = self.Geom[1,0]*0.5*CrXY[4]*self.CrGeRa[i]*SampleWeight[self.IntTCr,self.nIntCr,j]  # weighting factor, self.dC length of embedded discontinuity in global system
                    fd[ i*3:(i+1)*3]             = fd[i*3:(i+1)*3]             +f*dot(transpose(Nw),dot(transpose(self.CrTT[i]),tr))
                    Kdd[i*3:(i+1)*3,i*3:(i+1)*3] = Kdd[i*3:(i+1)*3,i*3:(i+1)*3]+f*dot(transpose(Nw),dot(TW_,Nw))
#                    
#                    print(f'XXX, tr {j:d}, {tr[0]:f}, {tr[1]:f}, f {f:f},, fd {", ".join(z for z in [f"{x:f}" for x in fd[0:6]])} ')
#                    
            else:                                                                       # dummy values for second crack
                Kdd[3*i,3*i], Kdd[3*i+1,3*i+1], Kdd[3*i+2,3*i+2] = 1., 1., 1.
        
        return Kdd, fd
    #
    def CrStatus(self, i, j, ww): 
        # loading
        if ww[1]>=(self.wwLM[i,j,1]-1.0e-12):
            self.CrackStatus[i,j] = 0
        # unloading, reloading
        elif ww[1]>=0.:
            self.CrackStatus[i,j] = 1
        # crack closure
        else:
            self.CrackStatus[i,j] = 2
    #  local crack tractions      
    def CrTraction(self, i, j, ww, ff ):
        sf = self.ShearRetFactor                                                        # factor for local shear stiffness -- be careful with this!
        tt = zeros((2), dtype=float)
        dt = zeros((2,2), dtype=float)
        wI = 1./self.wLim[i]
        x  = wI * ww
        # normal - quadratic
        if x[1]<1.:                                                                     # critical crack width not exceeded
            tI = self.CrTL[i,1]                                                         # crack traction upon crack initiation
            x1 = self.wwLM[i,j,1]                                                       # starting point of unloading
            t1 = self.wwTM[i,j,1]                                                       # crack traction of starting point
            w  = ww[1]
            beta1= 0.1 #0.1 # 0.1                                                          # ratio of horizontal point offset (to left) to unloading starting point  
            beta = 0. # 0. # 0.1                                                             # ratio of current unloading starting point to start with crack closure transition
            dd = beta1*x1                                                               # point with horizontal course of transition
            alpha = 3.                                                                  # factor for left end (--> x1-alpha*dd) of transition
            xx = x1-alpha*dd                                                            # anchor point for anchor left of smooth unloading transition
            # loading
            if self.CrackStatus[i,j]==0:
                tt[1]   =    tI   *(x[1]*x[1] -2.*x[1] + 1.)
                dt[1,1] = 2.*tI*wI*(x[1]-1.)
                # linear
#                tt[1]   =    tI*(1.-x[1])
#                dt[1,1] =   -tI*wI
                #
            # unloading - reloading
            elif self.CrackStatus[i,j]==1:
                # smoothed transition
                if dd > 0.: # ZeroD:
                    denom = dd**2*(-6*alpha**2*dd + 2*dd-3*x1*alpha+3*alpha**2*x1+6*alpha*dd) # denom might become very small, don't use ZeroD to catch for zero division
                    x1_ = x1*wI                                                         # related starting point of unloading
                    yy1 =    tI   *(x1_*x1_ -2.*x1_ + 1.)                               # traction of unloading starting point
                    dy1 = 2.*tI*wI*(x1_-1.)                                             # derivative / tangential stiffness
                    # linear
#                    yy1 =    tI*(1.-x1_)
#                    dy1 =   -tI*wI
                    #
                    aa  = (yy1-2*dd*alpha*dy1+2*dd*dy1+x1*alpha*dy1-x1*dy1)/denom
                    bb  = -0.5*(6*x1*yy1-12*x1*alpha*dd*dy1+15*dd*x1*dy1+6*x1**2*alpha*dy1-6*x1**2*dy1-3*dd*yy1-8*dd**2*dy1+6*dd**2*dy1*alpha**2-3*dd*dy1*alpha**2*x1)/denom
                    cc  = (3*x1**2*yy1-6*x1**2*dd*alpha*dy1+9*x1**2*dy1*dd+3*x1**3*alpha*dy1-3*x1**3*dy1-3*x1*dd*yy1-8*dd**2*x1*dy1+9*x1*dd**2*dy1*alpha**2-3*dd*dy1*alpha**2*x1**2-6*dd**3*dy1*alpha**2+6*dd**3*dy1*alpha-3*dd**2*dy1*x1*alpha+2*dd**3*dy1)/denom
                    dd_ = -0.5*(2*x1**3*yy1-4*x1**3*alpha*dd*dy1+7*x1**3*dd*dy1+2*x1**4*alpha*dy1-2*x1**4*dy1-3*x1**2*dd*yy1-6*x1**2*alpha*dd**2*dy1-8*x1**2*dd**2*dy1-12*dd**3*x1*dy1*alpha**2+12*dd**3*x1*dy1*alpha+12*dd**2*x1**2*dy1*alpha**2+4*dd**3*x1*dy1+12*dd**3*yy1*alpha**2-12*dd**3*yy1*alpha-6*dd**2*yy1*alpha**2*x1+6*dd**2*yy1*x1*alpha-4*dd**3*yy1-3*x1**3*dy1*dd*alpha**2)/denom
                    # smooth transition 
                    if w > xx:                                                 # checks left end of transition from loading path  
                        tt[1]   =    aa*w**3 +    bb*w**2 + cc*w +dd_
                        dt[1,1] = 3.*aa*w**2 + 2.*bb*w    + cc
                    # intermediate linear to origin
                    elif w > beta*xx:                                        # left end of transition from loading path x beta
                        yy = aa*xx**3 + bb*xx**2 + cc*xx +dd_                           # traction at anchor point
                        CC  = yy/xx
                        tt[1]   = CC * w
                        dt[1,1] = CC
                    # transition to crack closure -- smooth for beta>0
                    else:
                        yy = aa*xx**3 + bb*xx**2 + cc*xx +dd_
                        aa_ = -1./4.*(xx*self.BulkEmod-yy)/beta/xx**2
                        bb_ =  1./2.*(xx*self.BulkEmod+yy)/xx
                        cc_ = -1./4.*self.BulkEmod*beta*xx+1./4.*yy*beta
                        tt[1]   =    aa_*w**2 + bb_*w + cc_
                        dt[1,1] = 2.*aa_*w    + bb_ 
                # sharp transition
                else:
                    C = self.wwTM[i,j,1]/self.wwLM[i,j,1]
                    tt[1]   = C * w
                    dt[1,1] = C
            # crack closure    
            elif self.CrackStatus[i,j]==2:
                # smooth transition -- beta > 0 should be fulfilled here
                if w > -beta*x1:
                    aa_ = -1./4./beta/x1**2*(x1*self.BulkEmod-t1)
                    bb_ =  1./2.*(x1*self.BulkEmod+t1)/x1
                    cc_ = -1./4.*self.BulkEmod*beta*x1+1/4*t1*beta
                    tt[1]   = aa_*w**2 + bb_*w + cc_
                    dt[1,1] = 2.*aa_*w + bb_ 
                # linear
                else:
                    tt[1]   = self.BulkEmod*w
                    dt[1,1] = self.BulkEmod
#                sf = 1.                    # experiences convergence problems uhc 190210 
        # shear -- linear
        if x[1]<1.:                                                         # critical crack width not exceeded
            if x[1]>0.: xx = 1.-x[1]
            else:       xx = 1.
            tt[0]   =  xx**2 * sf*self.BulkEmod * ww[0]
            dt[0,0] =  xx**2 * sf*self.BulkEmod
        #
        return tt, dt
    #
    def CrFin(self, fw, Kww, TimeDelt, f6):
        if self.Type in ["CPS4S","CPE4S","CPS3S","CPE3S"]:
            if Kww[2,2]==0.: Kww[2,2] = 1.                                  # ad hoc solution for reduced integration element
        if self.CrackTraction:
            Kdd, fd = self.CrComp( TimeDelt, f6 )
            self.fwd[:] = fw[:]-fd[:]                                       # residuum of local integrated discontinuity forces
#            print('fwd', self.fwd[0:3], file = f6)         # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else: 
            if self.Type in ["CPS4S","CPE4S","CPS3S","CPE3S"]:
#                if Kww[2,2]==0.: Kww[2,2] = 1.                              # ad hoc solution for reduced integration element
                if Kww[3,3]==0.: Kww[3,3], Kww[4,4], Kww[5,5] = 1., 1., 1.  # ad hoc solution for elements with one crack without traction'
                Kdd = zeros((6,6),dtype=float)
            elif self.Type in ["C3D8S"]:
                Kdd = zeros((9,9),dtype=float)
        try:        
            self.Dd[:] = inv(Kww+Kdd)[:]                                    # self.Dd carried over to next equilibrium iteration loop and is element specific
        except:
            EchoM( Kww, Label='Kww '+str(self.Label), ff=f6)
            EchoM( Kdd, Label='Kdd '+str(self.Label), ff=f6)
            EchoM(Kww+Kdd, Label='Kww+Kdd'+str(self.Label), ff=f6)
            raise NameError("ConFemElem::CrFin: Kww+Kdd nearly singular Elem",self.Label,det(Kww+Kdd))
        
    def ViscExtenCr2D(self, Dt, eta, Dww, sI, i, ipI):                      # sI base index, see also ElasticLt.__init__
        if Dt < ZeroD: return 0., [0., 0.]#raise NameError("CaeFemElements:ViscExtenCr2D: time step zero")
        # determine current scalar velocity
        def EvalV( VwOld, Dww, Dt):                                    
            if VwOld<ZeroD or Dww*VwOld<0.: VwOld = Dww/Dt
            Vw = 2.*Dww/Dt - VwOld
            return Vw
        # end def
        VwOld0 = self.wwV[i,ipI,0]
        VwOld1 = self.wwV[i,ipI,1]
        #
        Vw0 = EvalV( VwOld0, Dww[0], Dt)                                       # shear component
        Vw1 = EvalV( VwOld1, Dww[1], Dt)                                       # normal component
        zz = 2.*eta/Dt
        self.wwVN[i,ipI,0] = Vw0
        self.wwVN[i,ipI,1] = Vw1
        return zz, [ Vw0, Vw1] 
    # update discontinuity data - rotate discontinuity if scheduled
    def Update(self, Mat, NodeList, NoIndToCMInd, ff, MsgLevel):
        self.wwOld[:] = self.ww[:]                                                      # crack width in global system
        # check for discontinuity orientation and apply rotating discontinuity
        if self.RotCrack and Mat.Type in ['ISODAMAGE','ELASTICLT_SDA']:
            if   Mat.Type=='ISODAMAGE':      offset = 9
            elif Mat.Type=='ELASTICLT_SDA': offset = 6
            eps1M, nnM, ttM = self.Means( offset )
            if Mat.PrinStrains: Lim = 1.5*Mat.eps_ct
            else:               Lim = 1.5*Mat.fct
            if MsgLevel>0 and eps1M > Lim:                                                  # tensile strength strain, criterion for SDA initiation hard coded!!!
                print("-- mean bulk strain exceeds limit",self.Label,self.CrackStatus,self.CrNx,self.CrNy,'_',nnM[0],nnM[1],'__',eps1M/Mat.eps_ct)
                print("-- mean bulk strain exceeds limit",self.Label,self.CrackStatus,self.CrNx,self.CrNy,'_',nnM[0],nnM[1],'__',eps1M/Mat.eps_ct, file=ff)
            ln = sqrt(nnM[0]**2+nnM[1]**2)
            if ln>ZeroD: 
                nnM[0] = nnM[0]/ln                                                      # nnM should be normal --> unit normal
                nnM[1] = nnM[1]/ln                                                      #
                xyC= self.FindGlobalDiscontinuityEdges( [ self.CrXCo, self.CrYCo ], nnM, NodeList, NoIndToCMInd)# global coordinates of discontinuity end point
                dxC= xyC[0][0]-xyC[1][0]
                dyC= xyC[0][1]-xyC[1][1]
                dC = sqrt(dxC*dxC + dyC*dyC)
                if abs(nnM[0]*dxC + nnM[1]*dyC)>ZeroD: raise NameError("CaeFemElements:CPSx_SDA.Update: discontinuity end not perp to normal ",self.Label,nnM,xyC)
                RotDiscFlag = True
                for j in range(self.nIntLCr):
                    if self.CrackStatus[j]>1: RotDiscFlag = False                      # crack at least partially closed - strain criterion might be misleading due to Poisson effect
                if self.RotCrack and RotDiscFlag:
                    if (nnM[0]*self.CrNx+nnM[1]*self.CrNy) < 0.:
                        nnM[0] = -nnM[0]
                        nnM[1] = -nnM[1]
                    if MsgLevel>0 and (nnM[0]*self.CrNx+nnM[1]*self.CrNy) < 0.9:
                        print("-- large change of crack orientation "+str(self.Label),self.CrackStatus,self.CrNx,self.CrNy,'_',nnM[0],nnM[1])
                        print("-- large change of crack orientation "+str(self.Label),self.CrackStatus,self.CrNx,self.CrNy,'_',nnM[0],nnM[1], file=ff)
                    # transform internal state values into global system
                    LL = zeros((2,2), dtype=float)
                    TT = zeros((2,2), dtype=float)
                    LM = zeros((2,2), dtype=float)
                    TM = zeros((2,2), dtype=float)
                    for j in range(self.nIntLCr):
                        LL[j] = dot(transpose(self.CrTT),self.wwL[j])
                        TT[j] = dot(transpose(self.CrTT),self.wwT[j])
                        LM[j] = dot(transpose(self.CrTT),self.wwLM[j])
                        TM[j] = dot(transpose(self.CrTT),self.wwTM[j])
                    # quantities to change
                    if True:
                        self.CrNx  = nnM[0]          
                        self.CrNy  = nnM[1]
                        self.CrNz  = 0.
                        self.dC = dC
                        self.Mw, self.Mw_ = self.FormMw( self.CrXCo, self.CrYCo )
                        self.CrTT[0,0] =  nnM[1]
                        self.CrTT[0,1] = -nnM[0]
                        self.CrTT[1,0] =  nnM[0]
                        self.CrTT[1,1] =  nnM[1]
                    # transform internal state values into new local system
                    for j in range(self.nIntLCr):
                        self.wwL[j]  = dot(self.CrTT,LL[j])
                        self.wwT[j]  = dot(self.CrTT,TT[j])
                        self.wwLM[j] = dot(self.CrTT,LM[j])
                        self.wwTM[j] = dot(self.CrTT,TM[j])
                else:
                    pass
#                    CrNx_  = nnM[0]          
#                    CrNy_  = nnM[1]
#                    CrNz_  = nnM[2]
#                    if (self.CrNx*CrNx_+self.CrNy*CrNy_+self.CrNz*CrNz_ < 0.9) and RotDiscFlag:
#                        print("large change of disc orientation "+str(self.Label),self.CrackStatus,self.CrNx,self.CrNy,self.CrNz,'_',CrNx_,CrNy_,CrNz_)
#                        print("large change of disc orientation "+str(self.Label),self.CrackStatus,self.CrNx,self.CrNy,self.CrNz,'_',CrNx_,CrNy_,CrNz_, file=ff)
            else:                                                                       # whole element under compression
                pass
        # update max crack width
        for i in range(len(self.CrC)):
            if self.CrC[i,4]>0.:                                            # crack is there
                for j in range(self.nIntLCr):
                    if self.wwL[i,j,1]>self.wwLM[i,j,1]: 
                        self.wwLM[i,j,:] = self.wwL[i,j,:]                                          # update for largest crack width reached in loading history
                        self.wwTM[i,j,:] = self.wwT[i,j,:]
                    self.wwV[i,j,:] = self.wwVN[i,j,:]                                              # update local crack velocities
                    if MsgLevel>0 and self.CrackStatus[i,j] != self.CrackStatusOld[i,j]:
                        print("-- crack status changed ",self.Label,j,self.CrackStatusOld[i,j],"->",self.CrackStatus[i,j])
                        print("-- crack status changed ",self.Label,j,self.CrackStatusOld[i,j],"->",self.CrackStatus[i,j], file=ff)
                        self.CrackStatusOld[i,j] = self.CrackStatus[i,j]
        self.fwd[:] = 0.
    #
    def CreateSDA2ndCrack(self, MatList, NodeList, NoIndToCMInd, ff, uu, du, TimeDelt):
#        Mat = MatList[self.MatIp[0]]                                                    # material for integration point --> for whole element
        Mat = self.Material
        if Mat.S2ndCrack and Mat.RType==3: 
            if Mat.Type in ['ISODAMAGE','ELASTICLT_SDA',"MicroPl"]:
                if   Mat.Type=='ISODAMAGE':     offset = 3 + 6
                elif Mat.Type=='ELASTICLT_SDA': offset =     6
                elif Mat.Type=="MicroPl":    offset = self.Material.iS + 6                # iS = 21 + 5 = 26 -> +6  -> 32
                eps1M, nnM, ttM = self.Means( offset )                                  # means of: largest principal value (strain or stress), direction, corr. stress vector, see ElasticLT::MatC: update of state variables
                ln = sqrt(nnM[0]**2+nnM[1]**2)
                if ln>ZeroD: 
                    nnM[0] = nnM[0]/ln                                                  # nnM should be normal --> unit normal
                    nnM[1] = nnM[1]/ln
                dev = abs(nnM[0]*self.CrN[0,0]+nnM[1]*self.CrN[0,1])
#                if dev<Mat.S2ndCrackA: Flag_ = True                                     # scalar product of 1st and 2nd crack orientation (both normalized!) should be small enough that 2nd crack is more or less perpendicular to 1st
                if abs(dev)<Mat.S2ndCrackA: Flag_ = True                                     # scalar product of 1st and 2nd crack orientation (both normalized!) should be small enough that 2nd crack is more or less perpendicular to 1st
                else:                  Flag_ = False
                if Mat.PrinStrains:    Flag = eps1M>Mat.eps_ct     
                else:                  Flag = eps1M>Mat.fct
                if Flag_ and Flag:
#                if False:
                    if   self.Type in ["CPS4S","CPE4S"]: XX = self.FormX_( 0., 0., 0.)
                    elif self.Type in ["CPS3S","CPE3S"]: XX = self.FormX_( 1./3., 1./3., 1./3.)
                    xy = self.NodalCoordinates( NodeList, NoIndToCMInd )
                    if self.Type=="CPS3":
                        LS = []                                                         # for length of edges
                        for i in range(int(len(xy)/2)): LS += [sqrt((xy[2*i]-xy[2*i-2])**2+(xy[2*i+1]-xy[2*i-1])**2)]
                        qLS = max(LS)/min(LS)
                    else: qLS = -1
                    xyP= dot(XX,array(xy))                                              # undeformed global coordinates of crack point / element center
                    xyC= self.FindGlobalDiscontinuityEdges( xyP, nnM, NodeList, NoIndToCMInd)# global coordinates of discontinuity end point
                    dxC= xyC[0][0]-xyC[1][0]
                    dyC= xyC[0][1]-xyC[1][1]
                    dC = sqrt(dxC*dxC + dyC*dyC)
                    if abs(nnM[0]*dxC + nnM[1]*dyC)>ZeroD: raise NameError("CaeFemElements:CPS4.CreateSDA2ndCrack: cut end not perp to normal "+str(self.Label)+str(xyP)+str(nnM)+str(xyC))
                    # some discontinuity data
                    lN = sqrt(nnM[0]**2+nnM[1]**2)                                      # length of crack normal
                    self.CrN[ 1,0]   = nnM[0]/lN                                        # crack normal
                    self.CrN[ 1,1]   = nnM[1]/lN
                    self.Mw = self.FormMw( self.CrXCo, self.CrYCo )                     # recalculation of Mw
                    self.CrC[ 1,0]   = xyC[0][0]                                        #
                    self.CrC[ 1,1]   = xyC[0][1]
                    self.CrC[ 1,2]   = xyC[1][0]
                    self.CrC[ 1,3]   = xyC[1][1]
                    self.CrC[ 1,4]   = dC                                               # discontinuity length
                    self.CrTT[1,0,0] = nnM[1]                                           # transform to discontinuity orientation, discontinuity normal (alpha) and discontinuity orientation (alpha-Pi/2) are perpendicular
                    self.CrTT[1,0,1] =-nnM[0]                                           # cos(alpha-Pi/2) = sin(alpha)
                    self.CrTT[1,1,0] = nnM[0]                                           # sin(alpha-Pi/2) = -cos(alpha)
                    self.CrTT[1,1,1] = nnM[1]                                           #
                    ttLM = dot(self.CrTT[1],ttM[0:2])                                   # transform global tractions to local tractions in discontinuity system - currently for control purposes only
                    ttLM_= dot(ttM,nnM)                                                 # to eliminate negative values of directions / tractions - basically used for the following
                    self.CrTL[1,0]  = abs(ttLM[0])                                      # non-zero value presumably comes from viscous stress contributions
                    self.CrTL[1,1]  = ttLM_ # ttLM[1]
                    self.wLim[1] = 3.*self.CrackEnergy/ttLM_                            # critical crack width where normal traction becomes zero derived from crack energy and actual tensile limit
                    #
                    self.CrGeRa[1] = 1.                                            # used for scaling between fd and fw
                    _, fd = self.CrComp(TimeDelt, ff)
                    fw = self.CrFw( uu, du, MatList, TimeDelt)                     # uu, uo are currently identical due to position of CreateSDARankine in CaeFemMain
                    if norm(fw)>ZeroD: self.CrGeRa[1] = norm(fw)/norm(fd)
                    else: raise NameError("CaeFemElements:CPS4.CreateSDARankine: no traction upon SDA initialization")
                    print('create SDA Rankine 2nd', self.Label, self.Type, xyP, nnM, self.CrTL[1], '__', xyC, '__', dC, self.CrGeRa[1], qLS)
                    print('create SDA Rankine 2nd', self.Label, self.Type, xyP, nnM, self.CrTL[1], '__', xyC, '__', dC, self.CrGeRa[1], qLS, file=ff)
                    self.NoCracks = 2

class CPE3(ElementC2D):
    NLGeomI = False                                                         # Flag whether large deformations are already implemented for this element
    def __init__(self, Label, SetLabel, InzList, NoList, PlSt, NoLabToNoInd):
        if PlSt: XXX = "CPS3"
        else:    XXX = "CPE3"
#       Element.__init__(self,XXX,InzList, 2, 3,1,1, None, 2,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        Element.__init__(self,XXX,InzList, 2, 3,1,1, None, 2,False,Label,SetLabel,2, None,None,None,       NoList,NoLabToNoInd,[])
        self.PlSt = PlSt                                                    # flag for plane stress (True->plane stress, False->plane strain)
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        self.X0 = NoList[i0].XCo
        self.Y0 = NoList[i0].YCo
        self.X1 = NoList[i1].XCo
        self.Y1 = NoList[i1].YCo
        self.X2 = NoList[i2].XCo
        self.Y2 = NoList[i2].YCo
        #
        SolSecDic = SecDict[self.Set]
        mat = SolSecDic.Mat
        self.MatN = mat
        Material  = MaList[mat]
        self.Material = Material
        RegType  = Material.RType
        if self.RegType in [1]: self.DofT = (set([1,2,7]),set([1,2,7]),set([1,2,7]))       # tuple, type of dof for every node of this element: 1 -> u_x, 7->gradient field
        else:                   self.DofT = (set([1,2]),  set([1,2]),  set([1,2])  )
        self.nNod = len(self.DofT)                                          # number of nodes
        if RegType==1: self.DofI = zeros( (self.nNod,3), dtype=int)
        else:          self.DofI = zeros( (self.nNod,2), dtype=int)         # indices of global dofs per node
        self.RegType = RegType
        #
        self.AA = -self.Y0*self.X1+self.Y0*self.X2+self.Y2*self.X1+self.Y1*self.X0-self.Y1*self.X2-self.Y2*self.X0 # double of element area
        if self.AA<=0.: raise NameError("Something is wrong with this CPE3-element", self.Label,NoList[i0].Label,self.X0,self.Y0,'_',NoList[i1].Label,self.X1,self.Y1,'_',NoList[i2].Label,self.X2,self.Y2)
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 0.5*self.AA                                        # element area for numerical integration
        self.Geom[1,0] = SolSecDic.Val                                      # thickness
        self.Lch_ = sqrt(self.AA)                                           # characteristic length
        self.Lch = sqrt(self.AA) /self.nInt              # characteristic length - presumably obsolete -- but still needed here and there -- see Lubliner
        self.a1=(self.Y2*self.X1-self.Y1*self.X2)/self.AA
        self.a2=(self.Y0*self.X2-self.Y2*self.X0)/self.AA
        self.a3=(self.Y1*self.X0-self.Y0*self.X1)/self.AA
        self.b1=(self.Y1-self.Y2)/self.AA
        self.b2=(self.Y2-self.Y0)/self.AA
        self.b3=(self.Y0-self.Y1)/self.AA
        self.c1=(self.X2-self.X1)/self.AA
        self.c2=(self.X0-self.X2)/self.AA
        self.c3=(self.X1-self.X0)/self.AA
#        if MaList[self.MatN].Type in ['ISODAMAGE',"MicroPl"] and MaList[self.MatN].RType ==2:                             # find scaling factor for band width regularization 
#            self.CrBwS, self.CrBwS2 = MaList[self.MatN].CrackBandScale( self.Lch_)
#        Element.CharLengthCPE3 += [[self.Lch_,self.Label]]
        self.CrBScaleType()
        return [self.Lch_,self.Label]
    def FormN(self, L1, L2, L3):
        L3=1-L1-L2
        N = array([[L1,0,L2,0,L3,0],[0,L1,0,L2,0,L3]])        
        return N
    def FormB(self, L1, L2, L3, NLg):
        c1=self.c1
        c2=self.c2
        c3=self.c3
        b1=self.b1
        b2=self.b2
        b3=self.b3
        B =array([[b1,0,b2,0,b3,0],[0,c1,0,c2,0,c3],[c1,b1,c2,b2,c3,b3]])
        return B, 1, 0 
    def FormT(self, L1, L2, L3):
        L3=1-L1-L2
        T = array([L1, 0, L2, 0, L3, 0])
        return T
    def FormX(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([L1, L2, L3])
        return X
    def FormX_(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([[L1, 0, L2, 0, L3, 0],[0, L1, 0, L2, 0, L3]])
        return X
    def JacoD(self, r, s, t):
        return 1
    def NLocCoorFromGlobCoor(self, xyI, xyC):
        L1 = self.a1 + self.b1 * xyC[0] + self.c1 * xyC[1];
        L2 = self.a2 + self.b2 * xyC[0] + self.c2 * xyC[1];
        L3 = self.a3 + self.b3 * xyC[0] + self.c3 * xyC[1];
        return L1, L2, L3

class CPE3_SDA(CPE3,ElementC2D_SDA):
    def __init__(self, Label, InzList, PlSt, NoList,      Coor,               Normal,            SolSecDic, SetLabel,        MatName, NoLabToNoInd,NoIndToCMInd, Material):
        CPE3.__init__(self, Label, SetLabel, InzList, NoList, PlSt, NoLabToNoInd)
        if PlSt: self.Type = "CPS3S"
        else:    self.Type = "CPE3S"
        self.nNod = 3                                                       # number of nodes
        self.DofT = (set([1,2]),  set([1,2]),  set([1,2]))
        # replaces Ini2
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        self.X0 = NoList[i0].XCo
        self.Y0 = NoList[i0].YCo
        self.X1 = NoList[i1].XCo
        self.Y1 = NoList[i1].YCo
        self.X2 = NoList[i2].XCo
        self.Y2 = NoList[i2].YCo
        self.AA = -self.Y0*self.X1+self.Y0*self.X2+self.Y2*self.X1+self.Y1*self.X0-self.Y1*self.X2-self.Y2*self.X0 # double of element area
        if self.AA<=0.: raise NameError("Something is wrong with this CPE3-element")
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 0.5*self.AA                                        # element area for numerical integration
        self.Geom[1,0] = SolSecDic.Val                                      # thickness
        self.Lch_ = sqrt(self.AA)              # characteristic length
#        self.Lch = sqrt(self.AA) /self.nInt              # characteristic length - presumably obsolete
        self.a1=(self.Y2*self.X1-self.Y1*self.X2)/self.AA
        self.a2=(self.Y0*self.X2-self.Y2*self.X0)/self.AA
        self.a3=(self.Y1*self.X0-self.Y0*self.X1)/self.AA
        self.b1=(self.Y1-self.Y2)/self.AA
        self.b2=(self.Y2-self.Y0)/self.AA
        self.b3=(self.Y0-self.Y1)/self.AA
        self.c1=(self.X2-self.X1)/self.AA
        self.c2=(self.X0-self.X2)/self.AA
        self.c3=(self.X1-self.X0)/self.AA
#        if Material.Type in ['ISODAMAGE',"MicroPl"] and Material.RType ==2:                             # find scaling factor for band width regularization 
#            self.CrBwS, self.CrBwS2 = Material.CrackBandScale( self.Lch_)
#        self.CrBScaleType()   --> CreateSDARankine
        #
#        self.ElemDimData( NoList,NoIndToCMInd ) --> CreateSDARankine
        #
        self.SDANew = True # False #True
        ElementC2D_SDA.__init__(self, Coor, Normal, Material)
        # further SDA definitions follow in CreateSDARankine
    # Nw-matrix of embedded discontinuity
    def FormNw(self, x0, y0, x, y):
        Nw = array([[1, 0, -y+y0],[0,1,x-x0]])
#        Nw_ = array([[1, 0, -y+y0, 0, 0,     0],
#                    [0, 1,  x-x0, 0, 0,     0],
#                    [0, 0,  0,    1, 0, -y+y0],
#                    [0, 0,  0,    0, 1,  x-x0]])
        return Nw
    # Mw-matrix of embedded discontinuity
    def FormMw(self, x0, y0):
        xx, yy = [self.X0,self.X1,self.X2], [self.Y0,self.Y1,self.Y2]       # nodal coordinates
        if not self.SDANew:
            rr = 0.5
            Mw = zeros((6,6), dtype=float)                                      # 1st 6 for CPx3 nodal degrees of freedom, 2nd 6 for 2 x 3 crack degrees of freedom
            for cr in range(self.CrN.shape[0]):                                 # number of shape[0] potential cracks 
                k, HH_ = 0, []
                for i in range(len(xx)):                                        # loop over nodes  of SDA element
                    # assemble H-matrix for crack cr
                    if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                        if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.:  # determine which node is on the positive part and which on the negative for matrix H
                            for j in range(self.DofN[i]): HH_ += [rr]
                        else:
                            for j in range(self.DofN[i]): HH_ += [rr-1.]
                    else:
                        for     j in range(self.DofN[i]): HH_ += [0.]
                    # assemble Mw-matrix for crack cr
                    Nw = self.FormNw(x0,y0,xx[i],yy[i])
                    for j in range(self.DofN[i]):
                        for r in range(3): Mw[k,3*cr+r] = HH_[k]*Nw[j,r]        # 3*cr for two potential cracks, r for discontinuity degrees of freedom -- 3 in total
                        k += 1                                                                      # sums up dimension of Mw
            return Mw
        else:
            Mw = zeros((self.nIntL,2,6), dtype=float)                       # nIntL: number of integration points, (: nodal degrees of freedom, 6; discontinuity degrees of freedom
            for cr in range(self.CrN.shape[0]):                             # number of shape[0] potential cracks 
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if cr==1: raise NameError("2nd crack not yet implemented")
                    for j in range(self.nIntL):                             # integration point loop
                        r = SamplePoints[self.IntT,self.nInt-1,j][0]        # local integration point coordinates  - integration type and order of element
                        s = SamplePoints[self.IntT,self.nInt-1,j][1]
                        X = self.FormX( r, s, 0.)                           # coefficients of form function 
                        xi= dot(X,xx)                                       # integration point global coordinates
                        yi= dot(X,yy)
                        if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xi-self.CrXCo,yi-self.CrYCo]) )>0.:  # determine if ip is on positive part and which on negative - works only for centered cracks
                            H = 1.
                        else: H = 0.
                        phi = 0.
                        for i in range(len(xx)):                            # loop over nodes  of SDA element
                            if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.: # determine which node is on the positive part and which on the negative for matrix H - works only for centered cracks
                                phi += X[i]                                 # sum of shape function in positive area
                        Mw[j,:,:3] = (H - phi)*self.FormNw( x0, y0, xi, yi)
            return Mw
    # Bw-matrix of embedded discontinuity
    def FormBw(self, L1, L2, L3):
        if not self.SDANew:
            BB, _, _ = self.FormB( L1, L2, L3, None)
            return dot(BB,self.Mw)
        else:
            BB = zeros((self.nIntL,3,6), dtype=float)
            for cr in range(self.CrN.shape[0]):                                 # number of shape[0] potential cracks 
                xx, yy = [self.X0,self.X1,self.X2], [self.Y0,self.Y1,self.Y2]   # nodal coordinates
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if cr==1: raise NameError("2nd crack not yet implemented")
                    for j in range(self.nIntL):                             # integration point loop
                        r = SamplePoints[self.IntT,self.nInt-1,j][0]        # local integration point coordinates  - integration type and order of element
                        s = SamplePoints[self.IntT,self.nInt-1,j][1]
                        X = self.FormX( r, s, 0.)                           # coefficients of form function 
                        xi= dot(X,xx)                                       # integration point global coordinates
                        yi= dot(X,yy)
                        B, _, _ = self.FormB( r, s, 0, None)
                        for i in range(len(xx)):                            # loop over nodes  of SDA element
                            if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.: # determine which node is on the positive part and which on the negative for matrix H - works only for centered cracks
                                BB[j,0,0] += B[0,2*i]
                                BB[j,0,2] += B[0,2*i]  *(-yi+self.CrYCo)
                                BB[j,1,1] += B[1,2*i+1]
                                BB[j,1,2] += B[1,2*i+1]*( xi-self.CrXCo)
                                BB[j,2,0] += B[2,2*i]
                                BB[j,2,1] += B[2,2*i+1]
                                BB[j,2,2] += B[1,2*i+1]*(-yi+self.CrYCo)+B[0,2*i]*(xi-self.CrXCo)
            return BB
    # shape function of embedded discontinuity geometry
    def FormXCut(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X

class CPE6(CPE3, ElementC2D):
    NLGeomI = False                                                         # Flag whether large deformations are already implemented for this element
    def __init__(self, Label, SetLabel, InzList, MatName,Material, NoList, SolSecDic, StateV, NData, PlSt, NoLabToNoInd):
        if PlSt: XXX = "CPS6"
        else:    XXX = "CPE6"
        Element.__init__(self,XXX,InzList, 2, 3,2,3, (set([1,2]),set([1,2]),set([1,2]),set([1,2]),set([1,2]),set([1,2])), 2,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
#       Element.__init__(self,XXX,InzList, 2, 3,3,4, (set([1,2]),set([1,2]),set([1,2]),set([1,2]),set([1,2]),set([1,2])), 2,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        self.PlSt = PlSt                                                    # flag for plane stress (True->plane stress, False->plane strain)
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        self.X0 = NoList[i0].XCo
        self.Y0 = NoList[i0].YCo
        self.X1 = NoList[i1].XCo
        self.Y1 = NoList[i1].YCo
        self.X2 = NoList[i2].XCo
        self.Y2 = NoList[i2].YCo

        SolSecDic = SecDict[self.Set]
        mat = SolSecDic.Mat
        self.MatN = mat
        Material  = MaList[mat]
        self.Material = Material
        RegType  = Material.RType
        if self.RegType in [1]: self.DofT = (set([1,2,7]),set([1,2,7]),set([1,2,7]),set([1,2,7]),set([1,2,7]),set([1,2,7]))       # tuple, type of dof for every node of this element: 1 -> u_x, 7->gradient field
        else:                   self.DofT = (set([1,2]),  set([1,2]),  set([1,2]),  set([1,2]),  set([1,2]),  set([1,2])  )
        self.nNod = len(self.DofT)                                          # number of nodes
        if RegType==1: self.DofI = zeros( (self.nNod,3), dtype=int)
        else:          self.DofI = zeros( (self.nNod,2), dtype=int)         # indices of global dofs per node
        self.RegType = RegType
        
        self.AA = -self.Y0*self.X1+self.Y0*self.X2+self.Y2*self.X1+self.Y1*self.X0-self.Y1*self.X2-self.Y2*self.X0 # double of element area
        if self.AA<=0.: raise NameError("Something is wrong with this CPE3-element")
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 0.5*self.AA                                        # element area for numerical integration
        self.Geom[1,0] = SolSecDic.Val                                      # thickness
        self.Lch_ = sqrt(self.AA)                                           # characteristic length
#        self.Lch = sqrt(self.AA) /self.nInt              # characteristic length - presumably obsolete
        self.a1=(self.Y2*self.X1-self.Y1*self.X2)/self.AA
        self.a2=(self.Y0*self.X2-self.Y2*self.X0)/self.AA
        self.a3=(self.Y1*self.X0-self.Y0*self.X1)/self.AA
        self.b1=(self.Y1-self.Y2)/self.AA
        self.b2=(self.Y2-self.Y0)/self.AA
        self.b3=(self.Y0-self.Y1)/self.AA
        self.c1=(self.X2-self.X1)/self.AA
        self.c2=(self.X0-self.X2)/self.AA
        self.c3=(self.X1-self.X0)/self.AA
        self.CrBScaleType()
        return []
    def FormN(self, L1, L2, L3):
        L3=1-L1-L2
        N = array([[(2*L1-1)*L1, 0,         (2*L2-1)*L2, 0,         (2*L3-1)*L3, 0,          4*L1*L2, 0,       4*L2*L3, 0,       4*L3*L1, 0],
                   [ 0,         (2*L1-1)*L1, 0,         (2*L2-1)*L2, 0,         (2*L3-1)*L3, 0,       4*L1*L2, 0,       4*L2*L3, 0,       4*L3*L1]])        
        return N
    def FormB(self, L1, L2, L3, NLg):
        c1=self.c1
        c2=self.c2
        c3=self.c3
        b1=self.b1
        b2=self.b2
        b3=self.b3
        B = array([[b1*(4*L1-1), 0, b2*(4*L2-1), 0, b3*(4*L3-1), 0, 4*b1*L2+4*b2*L1, 0, 4*b2*L3+4*b3*L2, 0, 4*b3*L1+4*b1*L3, 0],
                   [0, c1*(4*L1-1), 0, c2*(4*L2-1), 0, c3*(4*L3-1), 0, 4*c1*L2+4*c2*L1, 0, 4*c2*L3+4*c3*L2, 0, 4*c3*L1+4*c1*L3],
                   [c1*(4*L1-1), b1*(4*L1-1), c2*(4*L2-1), b2*(4*L2-1), c3*(4*L3-1), b3*(4*L3-1), 4*c1*L2+4*c2*L1, 4*b1*L2+4*b2*L1, 4*c2*L3+4*c3*L2, 4*b2*L3+4*b3*L2, 4*c3*L1+4*c1*L3, 4*b3*L1+4*b1*L3]])
        return B, 1, 0 
    def FormT(self, L1, L2, L3):
        L3=1-L1-L2
        T = array([(2*L1-1)*L1, 0, (2*L2-1)*L2, 0, (2*L3-1)*L3, 0, 4*L1*L2, 0, 4*L2*L3, 0, 4*L3*L1, 0])
        return T
    def FormX(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([(2*L1-1)*L1, (2*L2-1)*L2, (2*L3-1)*L3, 4*L1*L2, 4*L2*L3, 4*L3*L1])
        return X
    def FormX_(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([[(2*L1-1)*L1, 0, (2*L2-1)*L2, 0, (2*L3-1)*L3, 0, 4*L1*L2, 0, 4*L2*L3, 0, 4*L3*L1, 0],
                   [0, (2*L1-1)*L1, 0, (2*L2-1)*L2, 0, (2*L3-1)*L3, 0, 4*L1*L2, 0, 4*L2*L3, 0, 4*L3*L1]])
        return X
#    def JacoD(self, r, s, t):
#        return 1
#    def NLocCoorFromGlobCoor(self, xyI, xyC):
#        L1 = self.a1 + self.b1 * xyC[0] + self.c1 * xyC[1];
#        L2 = self.a2 + self.b2 * xyC[0] + self.c2 * xyC[1];
#        L3 = self.a3 + self.b3 * xyC[0] + self.c3 * xyC[1];
#        return L1, L2, L3

class CPE6_SDA(CPE6,ElementC2D_SDA):
    def __init__(self, Label, InzList, PlSt, NoList,      Coor,               Normal,            SolSecDic, SetLabel,        MatName, NoLabToNoInd,NoIndToCMInd, Material):
        CPE6.__init__(self, Label, SetLabel, InzList, MatName,Material, NoList, SolSecDic, Material.StateVar, Material.NData, PlSt, NoLabToNoInd)
        if PlSt: self.Type = "CPS6S"
        else:    self.Type = "CPE6S"
        # replaces Ini2
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        i3 = NoIndToCMInd[self.Inzi[3]]
        i4 = NoIndToCMInd[self.Inzi[4]]
        i5 = NoIndToCMInd[self.Inzi[5]]
        self.X0 = NoList[i0].XCo
        self.Y0 = NoList[i0].YCo
        self.X1 = NoList[i1].XCo
        self.Y1 = NoList[i1].YCo
        self.X2 = NoList[i2].XCo
        self.Y2 = NoList[i2].YCo
        self.X3 = NoList[i3].XCo
        self.Y3 = NoList[i3].YCo
        self.X4 = NoList[i4].XCo
        self.Y4 = NoList[i4].YCo
        self.X5 = NoList[i5].XCo
        self.Y5 = NoList[i5].YCo
        self.AA = -self.Y0*self.X1+self.Y0*self.X2+self.Y2*self.X1+self.Y1*self.X0-self.Y1*self.X2-self.Y2*self.X0 # double of element area
        if self.AA<=0.: raise NameError("Something is wrong with this CPE3-element")
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 0.5*self.AA                                        # element area for numerical integration
        self.Geom[1,0] = SolSecDic.Val                                      # thickness
        self.Lch_ = sqrt(self.AA)                                           # characteristic length
#        self.Lch = sqrt(self.AA) /self.nInt                                 # characteristic length - presumably obsolete
        self.a1=(self.Y2*self.X1-self.Y1*self.X2)/self.AA
        self.a2=(self.Y0*self.X2-self.Y2*self.X0)/self.AA
        self.a3=(self.Y1*self.X0-self.Y0*self.X1)/self.AA
        self.b1=(self.Y1-self.Y2)/self.AA
        self.b2=(self.Y2-self.Y0)/self.AA
        self.b3=(self.Y0-self.Y1)/self.AA
        self.c1=(self.X2-self.X1)/self.AA
        self.c2=(self.X0-self.X2)/self.AA
        self.c3=(self.X1-self.X0)/self.AA
#        if Material.Type in ['ISODAMAGE',"MicroPl"] and Material.RType ==2:                                              # find scaling factor for band width regularization 
#            self.CrBwS, self.CrBwS2 = Material.CrackBandScale( self.Lch_)
#        self.CrBScaleType()
        #
        self.ElemDimData( NoList,NoIndToCMInd ) # !!!!!!! for SDANew = False 
        #
        self.SDANew = False #True
        ElementC2D_SDA.__init__(self, Coor, Normal, Material)

    # Nw-matrix of embedded discontinuity
    def FormNw(self, x0, y0, x, y):
        Nw = array([[1, 0, -y+y0],[0,1,x-x0]])
        return Nw
    # Mw-matrix of embedded discontinuity
    def FormMw(self, x0, y0):
        xx, yy = [self.X0,self.X1,self.X2,self.X3,self.X4,self.X5], [self.Y0,self.Y1,self.Y2,self.Y3,self.Y4,self.Y5] # nodal coordinates
        rr = 0.5
        Mw = zeros((12,6), dtype=float)                                     # 1st 12 for CPx6 nodal degrees of freedom, 2nd 6 for 2 x 3 crack degrees of freedom
        for cr in range(self.CrN.shape[0]):                                 # number of shape[0] potential cracks 
            k, HH_ = 0, []
            for i in range(len(xx)):                                        # loop over nodes  of SDA element
                # assemble H-matrix for crack cr
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.:  # determine which node is on the positive part and which on the negative for matrix H
                        for j in range(self.DofN[i]): HH_ += [rr]
                    else:
                        for j in range(self.DofN[i]): HH_ += [rr-1.]
                else:
                    for     j in range(self.DofN[i]): HH_ += [0.]
                # assemble Mw-matrix for crack cr
                Nw = self.FormNw(x0,y0,xx[i],yy[i])
                for j in range(self.DofN[i]):
                    for r in range(3): Mw[k,3*cr+r] = HH_[k]*Nw[j,r]        # 3*cr for two potential cracks, r for discontinuity degrees of freedom -- 3 in total
                    k += 1                                                                      # sums up dimension of Mw
        return Mw
    # Bw-matrix of embedded discontinuity
    def FormBw(self, L1, L2, L3):
        BB, _, _ = self.FormB( L1, L2, L3, None)
        return dot(BB,self.Mw)
    # shape function of embedded discontinuity geometry
    def FormXCut(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X

class CPE4(ElementC2D):
    def __init__(self, Label, SetLabel, InzList, NoList, PlSt, nI, NoLabToNoInd):
        if PlSt: XXX = "CPS4"
        else:    XXX = "CPE4"
#       if RegType==1: Element.__init__(self,XXX,InzList, 2, 1,2,4, (set([1,2,7]),set([1,2,7]),set([1,2,7]),set([1,2,7])), 2,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
#       else:          Element.__init__(self,XXX,InzList, 2, 1,2,4, (set([1,2]),  set([1,2]),  set([1,2]),  set([1,2])),   2,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
        Element.__init__(self,XXX,InzList, 2, 1,2,4, None,   2,False,Label,SetLabel,2, None,None,None, NoList,NoLabToNoInd,[])
        self.PlSt = PlSt                                            # flag for plane stress (True->plane stress, False->plane strain)
        if nI!=2:
            self.nInt = nI                                          # integration order - modification of Element.__ini__ settings
            self.nIntL= nI*nI                                       # total number of integration points
#        if RegType==1: self.DofI = zeros( (self.nNod,3), dtype=int)
#        else:          self.DofI = zeros( (self.nNod,2), dtype=int) # indices of global dofs per node
#        if Material.RType == 3: self.RegType = 3                    # regularization 3 SDA
#        else:                   self.RegType = RegType              # 0 for no regularization, 1 gradient, 2 crack band width 
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        i3 = NoIndToCMInd[self.Inzi[3]]
        self.X0 = NoList[i0].XCo
        self.Y0 = NoList[i0].YCo
        self.X1 = NoList[i1].XCo
        self.Y1 = NoList[i1].YCo
        self.X2 = NoList[i2].XCo
        self.Y2 = NoList[i2].YCo
        self.X3 = NoList[i3].XCo
        self.Y3 = NoList[i3].YCo
        if self.Type in ["CAX4"]:
            self.X = zeros((self.nIntL), dtype=float)
            for j in range(self.nIntL):  # build element stiffness with integration loop
                r = SamplePoints[self.IntT,self.nInt-1,j][0]
                s = SamplePoints[self.IntT,self.nInt-1,j][1]
                self.X[j] = 0.25 * ((1+r) * (1+s) * self.X0 + (1-r) * (1+s) * self.X1 + (1-r) * (1-s) * self.X2 + (1+r) * (1-s) * self.X3)
        SolSecDic = SecDict[self.Set]
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 1.                                                 # dummy for Area / Jacobi determinant used instead
        self.Geom[1,0]                          = SolSecDic.Val             # thickness
        if self.Type in['CAX4']: self.Geom[1,0] = 1.0                       # unit as dummy for circumference
        mat = SolSecDic.Mat
        self.MatN = mat
        Material  = MaList[mat]
        self.Material = Material
        RegType  = Material.RType
        if self.RegType in [1]: self.DofT = (set([1,2,7]),set([1,2,7]),set([1,2,7]),set([1,2,7]))       # tuple, type of dof for every node of this element: 1 -> u_x, 7->gradient field
        else:                   self.DofT = (set([1,2]),  set([1,2]),  set([1,2]),  set([1,2]))
        self.nNod = len(self.DofT)                                          # number of nodes
        if RegType==1: self.DofI = zeros( (self.nNod,3), dtype=int)
        else:          self.DofI = zeros( (self.nNod,2), dtype=int)         # indices of global dofs per node
        self.RegType = RegType
        
        self.AA = 0.5*( (self.X1-self.X3)*(self.Y2-self.Y0) + (self.X2-self.X0)*(self.Y3-self.Y1)) # element area
        if self.AA<=0.:
            raise NameError("Something is wrong with this CPE4-element",self.Label,[self.X0,self.Y0],[self.X1,self.Y1],[self.X2,self.Y2,[self.X3,self.Y3]])
        self.Lch = sqrt(self.AA)/self.nInt                                  # characteristic length --- used for ElasticLT only
        self.Lch_ = sqrt(self.AA)                                           # previous approach might be misleading
        self.CrBScaleType()
        return [self.Label,self.Lch_,self.ScaleType,self.CrBwS]
        
    def FormX(self, r, s, t):
        X = array([(1+r)*(1+s)*0.25,     (1-r)*(1+s)*0.25,    (1-r)*(1-s)*0.25,    (1+r)*(1-s)*0.25])
        return X
    def FormX_(self, r, s, t):                                                       # for geometry interpolation
        X = array([[(1+r)*(1+s)*0.25, 0, (1-r)*(1+s)*0.25,  0, (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25, 0],
                   [0, (1+r)*(1+s)*0.25,0, (1-r)*(1+s)*0.25,0, (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25]])
        return X
    def FormN(self, r, s, t):
        N = array([[(1+r)*(1+s)*0.25, 0, (1-r)*(1+s)*0.25, 0,  (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25, 0],
                   [0, (1+r)*(1+s)*0.25,0, (1-r)*(1+s)*0.25,0, (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25]])
        return N
    def FormB(self, r, s, t, NLg):
        h00= ( 1+s)*0.25
        h01= ( 1+r)*0.25
        h10=-( 1+s)*0.25
        h11= ( 1-r)*0.25
        h20= (-1+s)*0.25
        h21= (-1+r)*0.25
        h30= ( 1-s)*0.25
        h31=-( 1+r)*0.25
        JAC00 = h00*self.X0 + h10*self.X1 + h20*self.X2 + h30*self.X3
        JAC01 = h00*self.Y0 + h10*self.Y1 + h20*self.Y2 + h30*self.Y3
        JAC10 = h01*self.X0 + h11*self.X1 + h21*self.X2 + h31*self.X3
        JAC11 = h01*self.Y0 + h11*self.Y1 + h21*self.Y2 + h31*self.Y3
        det = JAC00*JAC11 - JAC01*JAC10
        deti = 1./det
        a1 = self.Y0*h01 + self.Y1*h11 + self.Y2*h21 + self.Y3*h31
        a2 = self.Y0*h00 + self.Y1*h10 + self.Y2*h20 + self.Y3*h30
        b1 = self.X0*h01 + self.X1*h11 + self.X2*h21 + self.X3*h31
        b2 = self.X0*h00 + self.X1*h10 + self.X2*h20 + self.X3*h30
        B00 = deti*( h00 * a1 - h01 * a2 )
        B10 = deti*( h10 * a1 - h11 * a2 )
        B20 = deti*( h20 * a1 - h21 * a2 )
        B30 = deti*( h30 * a1 - h31 * a2 )
        B01 = deti*(-h00 * b1 + h01 * b2 )
        B11 = deti*(-h10 * b1 + h11 * b2 )
        B21 = deti*(-h20 * b1 + h21 * b2 )
        B31 = deti*(-h30 * b1 + h31 * b2 )
        if self.RegType==1:
            N0 = (1+r)*(1+s)*0.25 
            N1 = (1-r)*(1+s)*0.25 
            N2 = (1-r)*(1-s)*0.25 
            N3 = (1+r)*(1-s)*0.25
            B = array([[ B00, 0,   0,   B10, 0,   0,   B20, 0,   0,   B30, 0   ,0  ],
                       [ 0,   B01, 0,   0,   B11, 0,   0,   B21, 0,   0,   B31, 0  ],
                       [ B01, B00, 0,   B11, B10, 0,   B21, B20, 0,   B31, B30, 0  ],
                       [ 0,   0,   B00, 0,   0,   B10, 0,   0,   B20, 0,   0,   B30],
                       [ 0,   0,   B01, 0,   0,   B11, 0,   0,   B21, 0,   0,   B31]])
            BN= array([[ B00, 0,   0,   B10, 0,   0,   B20, 0,   0,   B30, 0   ,0  ],
                       [ 0,   B01, 0,   0,   B11, 0,   0,   B21, 0,   0,   B31, 0  ],
                       [ B01, B00, 0,   B11, B10, 0,   B21, B20, 0,   B31, B30, 0  ],
                       [ 0,   0,   N0,  0,   0,   N1,  0,   0,   N2,  0,   0,   N3 ]])
            return B, BN, det, 0
        else:
            B = array([[ B00, 0,   B10, 0,   B20, 0,   B30, 0  ],
                       [ 0,   B01, 0,   B11, 0,   B21, 0,   B31],
                       [ B01, B00, B11, B10, B21, B20, B31, B30]])
            return B, det, 0
    def FormT(self, r, s, t):                                       # interpolation on temperature
        # !!! ordering might not yet be correct !!!
        if self.RegType==1: T = array([(1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25, (1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, 0, 0, 0, 0, 0, 0, 0, 0])
        else:               T = array([(1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25, (1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, 0, 0, 0, 0])
        return T
    def JacoD(self, r, s, t):
        h00= ( 1+s)*0.25
        h01= ( 1+r)*0.25
        h10=-( 1+s)*0.25
        h11= ( 1-r)*0.25
        h20= (-1+s)*0.25
        h21= (-1+r)*0.25
        h30= ( 1-s)*0.25
        h31=-( 1+r)*0.25
        JAC00 = h00*self.X0 + h10*self.X1 + h20*self.X2 + h30*self.X3
        JAC01 = h00*self.Y0 + h10*self.Y1 + h20*self.Y2 + h30*self.Y3
        JAC10 = h01*self.X0 + h11*self.X1 + h21*self.X2 + h31*self.X3
        JAC11 = h01*self.Y0 + h11*self.Y1 + h21*self.Y2 + h31*self.Y3
        return JAC00*JAC11 - JAC01*JAC10
    # determines local coordinates from global coordinates -- xyI nodal coordinates
    def NLocCoorFromGlobCoor(self, xyI, xyC):                                           # xyC should be given arbitrary global coordinate
        r, s = 0., 0.
        nI = 10
        for i in range(nI):
            xy = dot(self.FormN(r,s,0.),xyI)
            RR = xyC-xy
            rn = sqrt(RR[0]**2+RR[1]**2)
            if rn<1.e-6: break
#            N  = array([[(1+r)*(1+s)*0.25, 0, (1-r)*(1+s)*0.25, 0, (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25, 0],
#                       [0, (1+r)*(1+s)*0.25, 0, (1-r)*(1+s)*0.25, 0, (1-r)*(1-s)*0.25, 0, (1+r)*(1-s)*0.25]])
            Nr =  array([[(0+1)*(1+s)*0.25, 0, (0-1)*(1+s)*0.25, 0, (0-1)*(1-s)*0.25, 0, (0+1)*(1-s)*0.25, 0],
                       [0, (0+1)*(1+s)*0.25, 0, (0-1)*(1+s)*0.25, 0, (0-1)*(1-s)*0.25, 0, (0+1)*(1-s)*0.25]])
            Ns =  array([[(1+r)*(0+1)*0.25, 0, (1-r)*(0+1)*0.25, 0, (1-r)*(0-1)*0.25, 0, (1+r)*(0-1)*0.25, 0],
                       [0, (1+r)*(0+1)*0.25, 0, (1-r)*(0+1)*0.25, 0, (1-r)*(0-1)*0.25, 0, (1+r)*(0-1)*0.25]])
            Rr = -dot(Nr,xyI)
            Rs = -dot(Ns,xyI)
            Jaco = array([[Rr[0],Rs[0]],
                          [Rr[1],Rs[1]]])
            det_ = Jaco[0,0]*Jaco[1,1]-Jaco[0,1]*Jaco[1,0]
            if abs(det_) < ZeroD: raise NameError("CaeFemElements:CPS4.LocCoorFromGlobCoor: Something is wrong with this CP4-element")
            detI = 1./det_
            JI = detI*array([[ Jaco[1,1],-Jaco[0,1]],
                             [-Jaco[1,0], Jaco[0,0]]])
            r = r - JI[0,0]*RR[0] - JI[0,1]*RR[1]
            s = s - JI[1,0]*RR[0] - JI[1,1]*RR[1]
        if i>nI-1: raise NameError("CaeFemElements:CPS4.LocCoorFromGlobCoor: Iteration failed")
        return r, s, 0.
    def FormNI(self):                                               # N-inverse for gaussian 2D order 2x2 integration, see N_Inverse.py, to get nodal values from integration point values
        return array([[ 0.1339746, -0.5,       -0.5,        1.8660254],
                      [-0.5,        1.8660254,  0.1339746, -0.5      ],
                      [ 1.8660254, -0.5,       -0.5,        0.1339746],
                      [-0.5,        0.1339746,  1.8660254, -0.5      ]])
    
class CPE4_SDA(CPE4,ElementC2D_SDA):
    def __init__(self, Label, InzList, PlSt, NoList, Coor, Normal, SolSecDic, Set, MatName, NoLabToNoInd,NoIndToCMInd, Material, RedInt):
#       if RedInt: CPE4.__init__(self, Label, Set,      InzList, MatName,Material, NoList, SolSecDic, Material.StateVar, Material.NData, PlSt, 1, 3,       None, NoLabToNoInd) # reduced integration
#       else:      CPE4.__init__(self, Label, Set,      InzList, MatName,Material, NoList, SolSecDic, Material.StateVar, Material.NData, PlSt, 2, 3,       None, NoLabToNoInd)

        if RedInt: CPE4.__init__(self, Label, Set,      InzList, NoList, PlSt, 1, NoLabToNoInd)  # reduced integration
        else:      CPE4.__init__(self, Label, Set,      InzList, NoList, PlSt, 2, NoLabToNoInd)

        if PlSt: self.Type = "CPS4S"
        else:    self.Type = "CPE4S"
        # replaces Ini2
        if len(NoIndToCMInd)==0: Ind = [     self.Inzi[0],      self.Inzi[1],      self.Inzi[2],      self.Inzi[3]]
        else:                    Ind = [NoIndToCMInd[self.Inzi[0]],NoIndToCMInd[self.Inzi[1]],NoIndToCMInd[self.Inzi[2]],NoIndToCMInd[self.Inzi[3]]]
        self.nNod = 4                                                       # number of nodes
        self.DofT = (set([1,2]),  set([1,2]),  set([1,2]),  set([1,2]))
        self.X0 = NoList[Ind[0]].XCo
        self.Y0 = NoList[Ind[0]].YCo
        self.X1 = NoList[Ind[1]].XCo
        self.Y1 = NoList[Ind[1]].YCo
        self.X2 = NoList[Ind[2]].XCo
        self.Y2 = NoList[Ind[2]].YCo
        self.X3 = NoList[Ind[3]].XCo
        self.Y3 = NoList[Ind[3]].YCo
        self.AA = 0.5*( (self.X1-self.X3)*(self.Y2-self.Y0) + (self.X2-self.X0)*(self.Y3-self.Y1)) # element area
        if self.AA<=0.: raise NameError("CaeFemElements:CPS4.__init__: Something is wrong with this CPE4-element "+str(self.Label))
        self.Lch = sqrt(self.AA)                                            # characteristic length
        self.Lch_= self.Lch
        #
        self.Geom = zeros( (2,1), dtype=double)
        self.Geom[0,0] = 1.                                                 # dummy for Area / Jacobi determinant used instead
        self.Geom[1,0] = SolSecDic.Val                                      # thickness
#        self.Material = Material
#        self.CrBScaleType()                                                 # only for crack band regularization but may be later needed for a combination
#        self.ElemDimData( NoList,NoIndToCMInd )
        #
        self.SDANew = True # False #True
        ElementC2D_SDA.__init__(self, Coor, Normal, Material)
    # Nw-matrix of embedded discontinuity
    def FormNw(self, x0, y0, x, y):
        Nw = array([[1, 0, -y+y0],[0,1,x-x0]])
        return Nw
    # Mw-matrix of embedded discontinuity
    def FormMw(self, x0, y0):
        xx, yy = [self.X0,self.X1,self.X2,self.X3], [self.Y0,self.Y1,self.Y2,self.Y3]   # nodal coordinates
        if not self.SDANew:
            rr = 0.5
            Mw = zeros((8,6), dtype=float)
            for cr in range(self.CrN.shape[0]):                             # number of shape[0] potential cracks - crack normal CrN comes from ElementC2D_SDA,__init__  
                k, HH_ = 0, []
                for i in range(len(xx)):                                    # loop over nodes  of SDA element 
                    # assemble H-matrix for crack cr
                    if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:         # indicates whether crack cr already exists
                        if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.:  # determine which node is on the positive part and which on the negative for matrix H - works only for centered cracks
                            for j in range(self.DofN[i]): HH_ += [rr]
                        else:
                            for j in range(self.DofN[i]): HH_ += [rr-1.]
                    else:
                        for     j in range(self.DofN[i]): HH_ += [0.]
                    # assemble Mw-matrix for crack cr
                    Nw = self.FormNw(x0,y0,xx[i],yy[i])
                    for j in range(self.DofN[i]):
                        for r in range(3): Mw[k,3*cr+r] = HH_[k]*Nw[j,r]    # 3*cr for two potential cracks, r for discontinuity degrees of freedom -- 3 in total
                        k += 1                                                                      # sums up dimension of Mw
            return Mw
        else:
            Mw = zeros((self.nIntL,2,6), dtype=float)                       # nIntL: number of integration points, (: nodal degrees of freedom, 6; discontinuity degrees of freedom
            for cr in range(self.CrN.shape[0]):                             # number of shape[0] potential cracks 
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if cr==1: raise NameError("2nd crack not yet implemented")
                    for j in range(self.nIntL):                             # integration point loop
                        r = SamplePoints[self.IntT,self.nInt-1,j][0]        # local integration point coordinates  - integration type and order of element
                        s = SamplePoints[self.IntT,self.nInt-1,j][1]
                        X = self.FormX( r, s, 0.)                           # coefficients of form function 
                        xi= dot(X,xx)                                       # integration point global coordinates
                        yi= dot(X,yy)
                        if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xi-self.CrXCo,yi-self.CrYCo]) )>0.:  # determine if ip is on positive part and which on negative - works only for centered cracks
                            H = 1.
                        else: H = 0.
                        phi = 0.
                        for i in range(len(xx)):                            # loop over nodes  of SDA element
                            if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.: # determine which node is on the positive part and which on the negative for matrix H - works only for centered cracks
                                phi += X[i]                                 # sum of shape function in positive area
                        Mw[j,:,:3] = (H - phi)*self.FormNw( x0, y0, xi, yi)
            return Mw
    # Bw-matrix of embedded discontinuity
    def FormBw(self, r, s, t):
        if not self.SDANew:
            BB, _, _ = self.FormB( r, s, t, None)
            return dot(BB,self.Mw)
        else:
            BB = zeros((self.nIntL,3,6), dtype=float)
            for cr in range(self.CrN.shape[0]):                                 # number of shape[0] potential cracks 
                xx, yy = [self.X0,self.X1,self.X2,self.X3], [self.Y0,self.Y1,self.Y2,self.Y3]   # nodal coordinates
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if cr==1: raise NameError("2nd crack not yet implemented")
                    for j in range(self.nIntL):                             # integration point loop
                        r = SamplePoints[self.IntT,self.nInt-1,j][0]        # local integration point coordinates  - integration type and order of element
                        s = SamplePoints[self.IntT,self.nInt-1,j][1]
                        X = self.FormX( r, s, 0.)                           # coefficients of form function 
                        xi= dot(X,xx)                                       # integration point global coordinates
                        yi= dot(X,yy)
                        B, _, _ = self.FormB( r, s, 0, None)
                        for i in range(len(xx)):                            # loop over nodes  of SDA element
                            if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.: # determine which node is on the positive part and which on the negative for matrix H - works only for centered cracks
                                BB[j,0,0] += B[0,2*i]
                                BB[j,0,2] += B[0,2*i]  *(-yi+self.CrYCo)
                                BB[j,1,1] += B[1,2*i+1]
                                BB[j,1,2] += B[1,2*i+1]*( xi-self.CrXCo)
                                BB[j,2,0] += B[2,2*i]
                                BB[j,2,1] += B[2,2*i+1]
                                BB[j,2,2] += B[1,2*i+1]*(-yi+self.CrYCo)+B[0,2*i]*(xi-self.CrXCo)
            return BB
            
    # shape function of embedded discontinuity geometry
    def FormXCut(self, r, s, t):
        X = array([[ 0.5*(1-r), 0., 0.5*(1+r), 0.],
                   [ 0., 0.5*(1-r), 0., 0.5*(1+r)]])
        return X
    def FormH(self, x0, y0):
        xx, yy = [self.X0,self.X1,self.X2,self.X3], [self.Y0,self.Y1,self.Y2,self.Y3]   # nodal coordinates
        h = zeros((4), dtype=float)
        for cr in range(self.CrN.shape[0]):                                 # number of shape[0] potential cracks 
            for i in range(len(xx)):                                        # loop over nodes  of SDA element 
                if (self.CrN[cr,0]**2+self.CrN[cr,1]**2)>ZeroD:             # indicates whether crack cr already exists
                    if cr==1: raise NameError("2nd crack not yet implemented")
                    if dot( array([self.CrN[cr,0],self.CrN[cr,1]]), array([xx[i]-self.CrXCo,yy[i]-self.CrYCo]) )>0.:  # determine which node is on the positive part and which on the negative for matrix H
                        h[i] = 1.
    def FormNI(self):                                               # N-inverse for gaussian 2D order 1x1 integration
        return array([[ 1.0 ],
                      [ 1.0 ],
                      [ 1.0 ],
                      [ 1.0 ]])

class CAX4(CPE4): # axisymmetric 4 node bilinear
    def __init__(self, Label, SetLabel, InzList, NoList, nI, NoLabToNoInd):
        Element.__init__(self,"CAX4",InzList, 2, 1,2,4, None,   4,False,Label,SetLabel,2, None,None,None, NoList,NoLabToNoInd,[])
        self.PlSt = False                                           # not plane stress
        if nI!=2:
            self.nInt = nI                                          # integration order - modification of Element.__ini__ settings
            self.nIntL= nI*nI                                       # total number of integration points
    def FormB(self, r, s, t, NLg):
        h00 = (1 + s) * 0.25
        h01 = (1 + r) * 0.25
        h10 = -(1 + s) * 0.25
        h11 = (1 - r) * 0.25
        h20 = (-1 + s) * 0.25
        h21 = (-1 + r) * 0.25
        h30 = (1 - s) * 0.25
        h31 = -(1 + r) * 0.25
        JAC00 = h00 * self.X0 + h10 * self.X1 + h20 * self.X2 + h30 * self.X3
        JAC01 = h00 * self.Y0 + h10 * self.Y1 + h20 * self.Y2 + h30 * self.Y3
        JAC10 = h01 * self.X0 + h11 * self.X1 + h21 * self.X2 + h31 * self.X3
        JAC11 = h01 * self.Y0 + h11 * self.Y1 + h21 * self.Y2 + h31 * self.Y3
        det = JAC00 * JAC11 - JAC01 * JAC10
        deti = 1. / det
        a1 = self.Y0 * h01 + self.Y1 * h11 + self.Y2 * h21 + self.Y3 * h31
        a2 = self.Y0 * h00 + self.Y1 * h10 + self.Y2 * h20 + self.Y3 * h30
        b1 = self.X0 * h01 + self.X1 * h11 + self.X2 * h21 + self.X3 * h31
        b2 = self.X0 * h00 + self.X1 * h10 + self.X2 * h20 + self.X3 * h30
        B00 = deti * (h00 * a1 - h01 * a2)
        B10 = deti * (h10 * a1 - h11 * a2)
        B20 = deti * (h20 * a1 - h21 * a2)
        B30 = deti * (h30 * a1 - h31 * a2)
        B01 = deti * (-h00 * b1 + h01 * b2)
        B11 = deti * (-h10 * b1 + h11 * b2)
        B21 = deti * (-h20 * b1 + h21 * b2)
        B31 = deti * (-h30 * b1 + h31 * b2)
        x = 0.25*( (1+r)*(1+s)*self.X0 + (1-r)*(1+s)*self.X1 + (1-r)*(1-s)*self.X2 +(1+r)*(1-s)*self.X3 )
        x_ = 0.25/x                                                 # includes factor for B
        B = array([[B00, 0, B10, 0, B20, 0, B30, 0],
                   [0, B01, 0, B11, 0, B21, 0, B31],
                   [x_*(1+r)*(1+s),0, x_*(1-r)*(1+s),0, x_*(1-r)*(1-s),0, x_*(1+r)*(1-s),0 ],
                   [B01,         B00, B11,         B10, B21,         B20, B31,         B30]])
        return B, det, 0

# 3D continuum elements
class ElementC3D(Element):
    def __init__(self, NoIndToCMInd,NoList):
        for i in self.Inzi:
            if len(NoIndToCMInd)==0: node = NoList[i]
            else:                    node = NoList[NoIndToCMInd[i]]
            if   node.Type=="":      node.Type = "C3D"
            elif node.Type=="C3D":   pass
            else:                    raise NameError("CaeFemElements:ElementC3D: node has already different type", node.Type)

    def LocatePointInElement(self, NodeList, bb, NoIndToCMInd):
        Finis = False
        rs = array([ 0., 0., 0.])                                   # initial values of isoparametric coordinates
        R, inSC  = None, None 
        niter = 5
        if self.Type in ["C3D8","C3D8_EFG"]:
            xy = self.NodalCoordinates( NodeList, NoIndToCMInd )        # nodal coordinates
            for iter in range(niter):                               # iteration loop to determine isoparametric coordinates / R-values of boundary control value by Newton Raphson
                r, s, t  = rs[0], rs[1], rs[2]
                XX  = self.FormX_( r, s, t)
                xyP = dot(XX,array(xy))                             # undeformed global coordinates of sample point
                rxy = bb-xyP                                        # residuum of target coordinate to iterated coordinates
                rno = sqrt(rxy[0]**2+rxy[1]**2+rxy[2]**2)
                if rno<1.e-6: 
                    Finis = True
                    break
                d00 = -(1.-s)*(1.-t)*0.125*xy[0]+(1.-s)*(1.-t)*0.125*xy[3]+(1.+s)*(1.-t)*0.125*xy[6]-(1.+s)*(1.-t)*0.125*xy[9] -(1.-s)*(1.+t)*0.125*xy[12]+(1.-s)*(1.+t)*0.125*xy[15]+(1.+s)*(1.+t)*0.125*xy[18]-(1.+s)*(1.+t)*0.125*xy[21]
                d01 = -(1.-r)*(1.-t)*0.125*xy[0]-(1.+r)*(1.-t)*0.125*xy[3]+(1.+r)*(1.-t)*0.125*xy[6]+(1.-r)*(1.-t)*0.125*xy[9] -(1.-r)*(1.+t)*0.125*xy[12]-(1.+r)*(1.+t)*0.125*xy[15]+(1.+r)*(1.+t)*0.125*xy[18]+(1.-r)*(1.+t)*0.125*xy[21]
                d02 = -(1.-r)*(1.-s)*0.125*xy[0]-(1.+r)*(1.-s)*0.125*xy[3]-(1.+r)*(1.+s)*0.125*xy[6]-(1.-r)*(1.+s)*0.125*xy[9] +(1.-r)*(1.-s)*0.125*xy[12]+(1.+r)*(1.-s)*0.125*xy[15]+(1.+r)*(1.+s)*0.125*xy[18]+(1.-r)*(1.+s)*0.125*xy[21]
                d10 = -(1.-s)*(1.-t)*0.125*xy[1]+(1.-s)*(1.-t)*0.125*xy[4]+(1.+s)*(1.-t)*0.125*xy[7]-(1.+s)*(1.-t)*0.125*xy[10]-(1.-s)*(1.+t)*0.125*xy[13]+(1.-s)*(1.+t)*0.125*xy[16]+(1.+s)*(1.+t)*0.125*xy[19]-(1.+s)*(1.+t)*0.125*xy[22]
                d11 = -(1.-r)*(1.-t)*0.125*xy[1]-(1.+r)*(1.-t)*0.125*xy[4]+(1.+r)*(1.-t)*0.125*xy[7]+(1.-r)*(1.-t)*0.125*xy[10]-(1.-r)*(1.+t)*0.125*xy[13]-(1.+r)*(1.+t)*0.125*xy[16]+(1.+r)*(1.+t)*0.125*xy[19]+(1.-r)*(1.+t)*0.125*xy[22]
                d12 = -(1.-r)*(1.-s)*0.125*xy[1]-(1.+r)*(1.-s)*0.125*xy[4]-(1.+r)*(1.+s)*0.125*xy[7]-(1.-r)*(1.+s)*0.125*xy[10]+(1.-r)*(1.-s)*0.125*xy[13]+(1.+r)*(1.-s)*0.125*xy[16]+(1.+r)*(1.+s)*0.125*xy[19]+(1.-r)*(1.+s)*0.125*xy[22]
                d20 = -(1.-s)*(1.-t)*0.125*xy[2]+(1.-s)*(1.-t)*0.125*xy[5]+(1.+s)*(1.-t)*0.125*xy[8]-(1.+s)*(1.-t)*0.125*xy[11]-(1.-s)*(1.+t)*0.125*xy[14]+(1.-s)*(1.+t)*0.125*xy[17]+(1.+s)*(1.+t)*0.125*xy[20]-(1.+s)*(1.+t)*0.125*xy[23]
                d21 = -(1.-r)*(1.-t)*0.125*xy[2]-(1.+r)*(1.-t)*0.125*xy[5]+(1.+r)*(1.-t)*0.125*xy[8]+(1.-r)*(1.-t)*0.125*xy[11]-(1.-r)*(1.+t)*0.125*xy[14]-(1.+r)*(1.+t)*0.125*xy[17]+(1.+r)*(1.+t)*0.125*xy[20]+(1.-r)*(1.+t)*0.125*xy[23]
                d22 = -(1.-r)*(1.-s)*0.125*xy[2]-(1.+r)*(1.-s)*0.125*xy[5]-(1.+r)*(1.+s)*0.125*xy[8]-(1.-r)*(1.+s)*0.125*xy[11]+(1.-r)*(1.-s)*0.125*xy[14]+(1.+r)*(1.-s)*0.125*xy[17]+(1.+r)*(1.+s)*0.125*xy[20]+(1.-r)*(1.+s)*0.125*xy[23]
                dxy = array([[ d00, d01, d02 ], [ d10, d11, d12 ] , [ d20, d21, d22] ])
                dxyI  = inv(dxy)
                rs = rs + dot(dxyI,rxy)                             # improved value of iso-parametric coordinates
                if rs[0]> 1.: rs[0] =  1.                           # constraints of iso-parametric coordinates
                if rs[0]<-1.: rs[0] = -1.
                if rs[1]> 1.: rs[1] =  1.
                if rs[1]<-1.: rs[1] = -1.
                if rs[2]> 1.: rs[2] =  1.
                if rs[2]<-1.: rs[2] = -1.
            if Finis:
#                if self.Type=="C3D8_EFG":
#                    inSC, NN = self.ArbFormN( NodeList,NoIndToCMInd, CoorTree_,NoCoLitoNoLi, rs[0], rs[1], rs[2], xyP)
#                else:                                               # should be C3D8
                NN  = self.FormN( rs[0], rs[1], rs[2])
                inSC= self.Inzi                                 # node indices for rs - trivial for C3D8 but not for C3D8_EFG 
                nen = len(inSC)                                     # number of base functions for rs/xyP
                R = zeros((nen),dtype=float)
                for j_ in range(nen): R[j_] = NN[0,3*j_]
        # tetraeder
        elif self.Type in ["C3D4"]:
#            xx = [self.X0,self.X1,self.X2,self.X3]
#            yy = [self.Y0,self.Y1,self.Y2,self.Y3]
#            zz = [self.Z0,self.Z1,self.Z2,self.Z3]
            xp = bb[0]
            yp = bb[1]
            zp = bb[2]
            xx0 = self.X0 
            xx1 = self.X1 
            xx2 = self.X2 
            xx3 = self.X3 
            yy0 = self.Y0 
            yy1 = self.Y1 
            yy2 = self.Y2 
            yy3 = self.Y3 
            zz0 = self.Z0 
            zz1 = self.Z1 
            zz2 = self.Z2 
            zz3 = self.Z3
            deno = (xx1*yy2*zz3+xx1*yy0*zz2-yy1*xx3*zz0+yy0*zz1*xx3-xx1*yy3*zz2-yy0*xx3*zz2-yy3*zz0*xx2+yy3*xx0*zz2+xx1*yy3*zz0+yy3*zz1*xx2-yy3*zz1*xx0-yy0*zz1*xx2-yy1*xx0*zz2+yy1*xx3*zz2+yy2*zz1*xx0+yy0*xx2*zz3-yy2*xx0*zz3+yy2*zz0*xx3-xx1*yy2*zz0-yy2*zz1*xx3+yy1*xx0*zz3-yy1*xx2*zz3+yy1*xx2*zz0-xx1*yy0*zz3)
            rr = (zz3*xx0*yp-zz3*yy0*xp-yy0*xx3*zz2-zz0*xx3*yp+zz0*xp*yy3-zp*xx0*yy3+zp*yy0*xx3-yy3*zz0*xx2+yy3*xx0*zz2-xp*yy3*zz2+yy0*xx2*zz3+zp*yy3*xx2-zz3*yp*xx2+xx3*yp*zz2-xx3*yy2*zp+xp*yy2*zz3-yy2*xx0*zz3+yy2*zz0*xx3-yp*xx0*zz2+yp*zz0*xx2+yy0*xp*zz2-yy2*zz0*xp-yy0*xx2*zp+yy2*xx0*zp)/deno
            ss =-(xx1*yp*zz0-xx1*yy3*zz0+yy1*xx3*zz0-yy1*xp*zz0-zz0*xx3*yp+zz0*xp*yy3+xx1*yy0*zz3+xx1*zp*yy3-xx1*zz3*yp-xx1*yy0*zp-zp*yy1*xx3-zp*xx0*yy3-zz3*yy0*xp+zz3*xx0*yp+yy3*zz1*xx0-yy0*zz1*xx3-yy1*xx0*zz3+yy1*xx0*zp+zp*yy0*xx3-yp*zz1*xx0-zz1*xp*yy3+zz3*yy1*xp+zz1*xx3*yp+yy0*zz1*xp)/deno
            tt = (xx1*yy2*zp+yy0*zz1*xp-yy2*zz1*xp-xx1*yp*zz2-yy1*xx2*zp+xx1*yp*zz0+yy2*zz0*xp-yy0*zz1*xx2+yy1*xx0*zp+yp*xx0*zz2+yy1*xp*zz2+yy0*xx2*zp-xx1*yy0*zp-yy2*xx0*zp-yy1*xp*zz0+yy2*zz1*xx0+xx1*yy0*zz2-yp*zz1*xx0+yp*zz1*xx2+yy1*xx2*zz0-yy1*xx0*zz2-xx1*yy2*zz0-yy0*xp*zz2-yp*zz0*xx2)/deno
#            for i in range(len(xx)):
#                if sqrt((xx[i]-xp)**2+(yy[i]-yp)**2+(zz[i]-zp)**2) < ZeroD:
#                    if   i==0: rs = [0.,0.,0.]
#                    elif i==1: rs = [1.,0.,0.]
#                    elif i==2: rs = [0.,1.,0.]
#                    elif i==3: rs = [0.,0.,1.]
#                    print('XXX',rs[0]-rr,rs[1]-ss,rs[2]-tt)
#                    Finis = True
#                    break
            Tol = 1.0e-9
            if (0.-Tol<=rr and rr<=1.+Tol and 0.-Tol<=ss and ss<=1.+Tol and 0.-Tol<=tt and tt<=1.+Tol):
                Finis = True
                rs[0] = rr 
                rs[1] = ss
                rs[2] = tt
                NN  = self.FormN( rs[0], rs[1], rs[2], None)
                inSC= self.Inzi                                         # node indices for rs - trivial for C3D8 but not for C3D8_EFG 
                nen = len(inSC)                                         # number of base functions for rs/xyP
                R   = zeros((nen),dtype=float)
                for j_ in range(nen): R[j_] = NN[0,3*j_]
        return Finis, rs, R, inSC
    #
    def CreateSDARankine(self, nea, MatList, NodeList, ElemList, SolSecs, NoLabToNoInd,NoIndToCMInd, ff, uu, uo, TimeDelt):
        nea_= nea
        MatName = self.MatN 
        Mat = MatList[MatName]
        if Mat.RType==3: 
            if Mat.Type in ["IsoD",'ELASTICLT_SDA']:
                if   Mat.Type=="IsoD":      offset = 9
                elif Mat.Type=='ELASTICLT_SDA': offset = 6
                eps1M, nnM, ttM = self.Means( offset )                      # see below for means - mean principal strain, direction, stress vector of integration points
                ln = sqrt(nnM[0]**2+nnM[1]**2+nnM[2]**2)
                if ln>ZeroD: 
                    nnM[0] = nnM[0]/ln                                      # nnM should be normal --> unit normal
                    nnM[1] = nnM[1]/ln
                    nnM[2] = nnM[2]/ln
                if Mat.PrinStrains: Flag = eps1M>Mat.eps_ct     
                else:               Flag = eps1M>Mat.fct
                if Flag:
#                    if   self.Type=="C3D8": XX = self.FormX_( 0., 0., 0.)   # to determine the coor of the center point
                    xyz = self.NodalCoordinates( NodeList, NoIndToCMInd )   # coor of element nodes
                    xyzP= dot( self.FormX_( 0., 0., 0.), array(xyz))        # undeformed global coordinates of crack point / element center
                    xyzC, v_num = self.FindGlobalDiscontinuityEdges( xyzP, nnM, xyz)# global coordinates of discontinuity end point
                    Cx = dot( nnM, array([1., 0., 0.])) # nnM[0] ? direction cosines ?
                    Cy = dot( nnM, array([0., 1., 0.]))
                    Cz = dot( nnM, array([0., 0., 1.]))
                    T1 = array([ Cx, Cy, Cz ])            # nnM again ?
                    # create local coordinate system
                    T3 = cross( T1, array([0., 0., 1.]))
                    if norm(T3)<1.e-3: T3 = cross( T1, array([1., 0., 0.]))
                    L3 = sqrt( T3[0]**2 + T3[1]**2 + T3[2]**2 )
                    T3 = T3/L3
                    T2 = cross(T3, T1)
                    # create SDA element
                    if self.Type == "C3D8":
                        l0 = NodeList[ NoIndToCMInd[self.Inzi[0]] ].Label   # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l1 = NodeList[ NoIndToCMInd[self.Inzi[1]] ].Label
                        l2 = NodeList[ NoIndToCMInd[self.Inzi[2]] ].Label
                        l3 = NodeList[ NoIndToCMInd[self.Inzi[3]] ].Label
                        l4 = NodeList[ NoIndToCMInd[self.Inzi[4]] ].Label   # following creation of new element instance generally requires node labels for Inzi instead of indices
                        l5 = NodeList[ NoIndToCMInd[self.Inzi[5]] ].Label
                        l6 = NodeList[ NoIndToCMInd[self.Inzi[6]] ].Label
                        l7 = NodeList[ NoIndToCMInd[self.Inzi[7]] ].Label
                        ElemList += [C3D8_SDA( self.Label, [l0,l1,l2,l3,l4,l5,l6,l7],self.Geom, NodeList, v_num, [xyzP[0],xyzP[1],xyzP[2]], [nnM[0],nnM[1],nnM[2]], self.Set, MatName, NoLabToNoInd,NoIndToCMInd, MatList, SolSecs)]
                    else:
                        raise NameError("CaeFemElements:C3D8.CreateSDARankine: unknown element type")
                    # some discontinuity data
                    NewEl = ElemList[-1]
                    NewEl.xyzC = xyzC                                       # vertices coordinates of cutting plane
                    NewEl.CrTT = array([[T1[0],T1[1],T1[2]],                # rotation matrix global to local # rotation matrix global(x, y, z) to local(n, m1, m2)
                                        [T2[0],T2[1],T2[2]],
                                        [T3[0],T3[1],T3[2]]])
                    ttLM_= dot( ttM, nnM)                                   # traction perpendicular to crack only
                    NewEl.CrTL[0] = ttLM_                                   # in C3D8_SDA.__init__ -- local crack traction upon crack initiation
                    #ttLM = dot(NewEl.CrTT, ttM[0:2])          # transform global tractions to local tractions in discontinuity system
                    #NewEl.CrTL[0] = abs(ttLM[0])              # non-zero value presumably comes from viscous stress contributions
                    NewEl.wLim = 3.*NewEl.CrackEnergy/ttLM_                 # critical crack width where normal traction becomes zero derived from crack energy and actual tensile limit
                    NewEl.Material = Mat
                    NewEl.CrBScaleType()                                    # only for crack band regularization but may be later needed for a combination
                    NewEl.ElemDimData( NodeList,NoIndToCMInd )
                    NData = len(self.Data[0])                               # number of data
                    NewEl.InitData( NewEl.nIntL, NData, Mat.StateVar, [])
                    for j in range(NewEl.n_vertices):                       # loop over triangles of cutting plane
                        XCut = NewEl.FormXCut()                             # shape function of embedded discontinuity geometry
                        c = NewEl.centCoor                                  # center coordinates
                        p1 = NewEl.xyzC[j - 1]                              # nodes to form triangle -- xyzC vertices coordinates of cutting plane
                        p2 = NewEl.xyzC[j]
                        NewEl.wwArea[j] = NewEl.tria_Area(c, p1, p2)
                        tria_coor = array([c[0], c[1], c[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]])
                        NewEl.wwIntPCoor += [dot(XCut,tria_coor)]           # global coordinates of integration points - derived from global embedded discontinuity vertics points
                    # complete data of new element by taking over from old element
                    for i in range(NewEl.DofI.shape[0]):
                        for j in range(NewEl.DofI.shape[1]): 
                            NewEl.DofI[i,j] = self.DofI[i,j]
                    nIP = self.StateVar.shape[0]                            # number of integration points
                    nSV = self.StateVar.shape[1]                            # number of state variables
                    if NewEl.nIntL==nIP:                                    # full integration for C3D8_SDA
                        for i in range(nIP):                                # loop over integration points
                            for j in range(nSV):                            # loop over number of state variables
                                NewEl.StateVar[i,j]  = self.StateVar[i,j]
                                NewEl.StateVarN[i,j] = self.StateVarN[i,j]
                    else:
                        raise NameError("ConFemElements:CPS4.CreateSDARankine: some kind of inconsistency")
                    NewEl.Data[:,:] = self.Data[:,:]                        # take current data of source element uhc
                    SolSecs[self.Set].Elems += [nea_]                       # add to SolSecs
                    nea_ += 1
                    # initial fw, fd
                    NewEl.CrGeRa = 1. 
                    Kdd, fd = NewEl.CrComp(TimeDelt, ff)
                    fw = NewEl.CrFw( uu, uo, MatList, TimeDelt)                     # uu, uo are currently identical due to position of CreateSDARankine in CaeFemMain
                    if norm(fw)>ZeroD: NewEl.CrGeRa = norm(fw)/norm(fd)
                    else: raise NameError("ConFemElements:CPS4.CreateSDARankine: no traction upon SDA initialization")
                    # echo                                             element center, crack normal, crack traction, discontinuity end points, some scaling factor for crack forces
                    print('create SDA Rankine', self.Label, self.Type, xyzP, nnM, NewEl.CrTL, '__', xyzC, '__', NewEl.CrGeRa)
                    print('create SDA Rankine', self.Label, self.Type, xyzP, nnM, NewEl.CrTL, '__', xyzC, '__', NewEl.CrGeRa, file=ff)
                    self.Active = False
                else:                                                       # no crack initiated
                    pass
            else:
                raise NameError("ConFemElements:CPS4.CreateSDARankine: regularization type RType 3 not yet implemented for this material")
        return nea_
    # means of crack tractions, crack normal, largest principal strain
    def Means( self, offset ):
        ttM, nnM, eps1M, XCounter = zeros((3), dtype= float),zeros((3), dtype= float), 0., 0  
        for i in range(self.nIntL):                                         # loop over integration sampling points
            if self.StateVar[i,offset+0]>ZeroD:
                # [offset] largest principal strain, [offset 1- offset 3] largest principal strain direction, [off set4 - offset 6 stress vector in largest principal strain direction
                eps1M  += self.StateVar[i,offset+0]                         # might also be largest principal stress
                nnM[0] += self.StateVar[i,offset+1]
                nnM[1] += self.StateVar[i,offset+2]
                nnM[2] += self.StateVar[i,offset+3]
                ttM[0] += self.StateVar[i,offset+4]
                ttM[1] += self.StateVar[i,offset+5]
                ttM[2] += self.StateVar[i,offset+6]
                XCounter += 1
        if XCounter>0:                                                      # means
            eps1M  = eps1M /XCounter 
            nnM[:] = nnM[:]/XCounter
            ttM[:] = ttM[:]/XCounter
        return eps1M, nnM, ttM
    
    def FindGlobalDiscontinuityEdges( self, xyz0, nn, xyz):   # xyz0 : element center, nn strain eigen vector
        # edge lines group one:
        xg1 = [xyz[0], xyz[3], xyz[6], xyz[9] ]
        yg1 = [xyz[1], xyz[4], xyz[7], xyz[10] ]
        zg1 = [xyz[2], xyz[5], xyz[8], xyz[11] ]
        # edge lines group two:
        xg2 = [xyz[12], xyz[15], xyz[18], xyz[21] ]
        yg2 = [xyz[13], xyz[16], xyz[19], xyz[22] ]
        zg2 = [xyz[14], xyz[17], xyz[20], xyz[23] ]
        # loop over all element edges
        xyzC = []
        for i in range(4):
            xE0 = xg1[i-1]
            xE1 = xg1[i]                       # coordinates of element edge end points
            yE0 = yg1[i-1]
            yE1 = yg1[i] 
            zE0 = zg1[i-1]
            zE1 = zg1[i] 
            
            dirEx = xE1 - xE0                       # direction of edge
            dirEy = yE1 - yE0
            dirEz = zE1 - zE0
            edge = array([dirEx, dirEy, dirEz])
            
            det = dot(nn, edge)           # vector dot product to check if they are parrallel
            if abs(det) > ZeroD:          # needs no else, ruled by length of xyzC 
                w = array([xE0, yE0, zE0 ]) - xyz0  # element center to egde node
                u = edge
                
                s = -dot(nn, w ) / dot(nn, u)
                if 0.<=s and s<=1.:
                    xC = xE0 + s*dirEx
                    yC = yE0 + s*dirEy
                    zC = zE0 + s*dirEz
                    xyzC += [[xC,yC, zC]]   # coor. of plane vertices 
            
            xE0 = xg2[i-1]
            xE1 = xg2[i]                       # coordinates of element edge end points
            yE0 = yg2[i-1]
            yE1 = yg2[i] 
            zE0 = zg2[i-1]
            zE1 = zg2[i] 
            dirEx = xE1 - xE0                       # direction of edge
            dirEy = yE1 - yE0
            dirEz = zE1 - zE0
            edge = array([dirEx, dirEy, dirEz])
            
            det = dot(nn, edge)           # vector dot product to check if they are parrallel
            if abs(det) > ZeroD:          
                w = array([xE0, yE0, zE0 ]) - xyz0
                u = edge
                s = -dot(nn, w ) / dot(nn, u)
                
                if 0.<=s and s<=1.:
                    xC = xE0 + s*dirEx
                    yC = yE0 + s*dirEy
                    zC = zE0 + s*dirEz
                    xyzC += [[xC,yC, zC]]
            
        # group three:    
            xE0 = xg1[i]
            xE1 = xg2[i]                       # coordinates of element edge end points
            yE0 = yg1[i]
            yE1 = yg2[i] 
            zE0 = zg1[i]
            zE1 = zg2[i] 
            dirEx = xE1 - xE0                       # direction of edge
            dirEy = yE1 - yE0
            dirEz = zE1 - zE0
            edge = array([dirEx, dirEy, dirEz])
            
            det = dot(nn, edge)           # vector dot product to check if they are parrallel
            if abs(det) > ZeroD:          
                w = array([xE0, yE0, zE0 ]) - xyz0
                u = edge
                
                s = -dot(nn, w ) / dot(nn, u)
                if 0.<=s and s<=1.:
                    xC = xE0 + s*dirEx
                    yC = yE0 + s*dirEy
                    zC = zE0 + s*dirEz
                    xyzC += [[xC,yC, zC]]
        # Sort  of vertices coordinates  
        v_num = len(xyzC)
        a = xyzC[0]
        xyzCut = pClosest0(xyzC, a)
        xyzC.remove(a)
        xyzC.remove(xyzCut[0])
        a = xyzCut[2]
        for z in range(v_num):
            LL = pClosest1(xyzC, a)
            xyzCut += [LL]
            xyzC.remove(a)
            a = LL
            if len(xyzCut) == v_num : break
        return xyzCut, v_num

class ElementC3D_SDA( ElementC2D_SDA ):
    def __init__(self, Coor, Normal, Material, v_num):
        self.CrXCo = Coor[0]                                                # crack center coordinates
        self.CrYCo = Coor[1]
        self.CrZCo = Coor[2]
        self.centCoor = Coor
        lN = sqrt(Normal[0]**2+Normal[1]**2+Normal[2]**2)                   # length of crack normal
        self.CrNx  = Normal[0]/lN                                           # crack normal
        self.CrNy  = Normal[1]/lN
        self.CrNz  = Normal[2]/lN
        self.CrTL  = zeros((3), dtype=float)                                # local crack traction upon crack initiation
        #self.CrTT  = zeros((3,3), dtype=float)                                          # coordinate transformation matrix -- from global to local with local x along crack direction
        #self.Mw    = self.FormMw( Coor[0], Coor[1], Coor[2]) #'add Coor[2]'
        self.Dd    = zeros((6,6), dtype=float)                              # for discontinuity dofs from discontinuity imbalance               

        self.ww    = zeros((6), dtype=float)                                # crack width dofs
        self.dw    = zeros((6), dtype=float)                                # crack width dofs - step increment        
        self.fwd   = zeros((6), dtype=float)                                # storage local forces of crack width dofs
        self.wwOld = zeros((6), dtype=float)
        #self.xyC   = zeros((4), dtype=float)                            # end points of crack [x0,y0,x1,y1]
        #self.dC    = 0.                                                 # length of crack 
        #self.nIntLCr= 2                                                 # for integration along crack contour 
        #self.IntTCr= 0
        #self.nIntCr= 1
        self.n_vertices = v_num
        self.wwL   = zeros((self.n_vertices,3), dtype=float)                # crack width in local system
        self.wwLM  = zeros((self.n_vertices,3), dtype=float)                # maximum crack width in local system during load history
        self.wwTM  = zeros((self.n_vertices,3), dtype=float)                # crack traction in local system -- "
        self.wwT   = zeros((self.n_vertices,3), dtype=float)                # crack traction in local system -- "
        self.wwTn  = zeros((self.n_vertices), dtype=float)                  # normal component
        self.wwLn  = zeros((self.n_vertices), dtype=float)
        self.wwV   = zeros((self.n_vertices,3), dtype=float)                # crack velocities in local system -- for viscous contributions
        self.wwVN  = zeros((self.n_vertices,3), dtype=float)                # final crack velocities in local system -- for viscous contributions
        self.CrackStatus = zeros((self.n_vertices), dtype=int)              # 0: loading; 1: unloading, reloading in tension; 2: unloading, reloading with crack closure
        self.CrackStatusOld = zeros((self.n_vertices), dtype=int)
        self.wwArea= zeros((self.n_vertices), dtype=float)                  # area assigned to integration points
        self.wwIntPCoor = []                                                # coordinates of integration points
        if Material.RType==3:
            self.CrackTraction = True
            if Material.Type in ["IsoD",'ELASTICLT_SDA']:
                self.CrackEnergy = Material.RegPar
                self.wLim = 0.1                                             # default value -- may be later modified to comply with crack energy
                self.BulkEmod = Material.Emod
                self.CrackEta = Material.etaCr
                self.RotCrack = Material.RotCrack                           # for rotating crack / discontinuity
                self.ShearRetFactor = Material.ShearRetFac                  # shear retention factor
            else:     raise NameError("CaeFemElements:C3D8_SDA.__init__: actual material type not yet implemented")
        else: self.CrackTraction = False
    #
    def CrFw(self, uu, uo, MatList, TimeDelt):  # uu, uo are currently identical due to position of CreateSDARankine in CaeFemMain
        #
        ue = zeros((self.DofE), dtype=float)                                               # initialization element nodal displacement vector
        de = zeros((self.DofE), dtype=float)                                               # initialization element nodal displacement load / time step increment vector
        ndof0 = 0
        for j in range(self.nNod):                                                      # loop over nodes of element
            for k in range(self.DofN[j]):                                               # loop over dofs of element node               
                kk = self.DofI[j,k]                                                     # global row index of dof
                ue[ndof0+k] = uu[kk]                                                    # local displacement
                de[ndof0+k] = uu[kk]-uo[kk]                                             # local displacement step increment
            ndof0 += self.DofN[j]

        fw = zeros((6),dtype=float)
        for j in range(self.nIntL):
            rst = SamplePoints[self.IntT,self.nInt-1,j]                                   # local integration sampling point coordinates
            val = None                                                                  # some element types need additional values for B-matrix, e.g. XFEM, EFG
            BB, JJ, _ = self.FormB( rst[0], rst[1], rst[2], val )
            eps = dot(BB,ue)                                                            # integration point strains
            dps = dot(BB,de)                                                            # integration point strain load / time step increments
            BBw = self.FormBw( rst[0], rst[1], rst[2] )                                 # Mw comes here into play
            eps = eps - dot(BBw,self.ww)
            dps = dps - dot(BBw,self.dw)
#            sig, MatM, X = MatList[self.MatN].MatC( None, 1, TimeDelt, j, self, dps, eps, [] ) # material for integration point    # MatM = MatList[self.Mat].MatC( self.dim, self.PlSt)         # material for whole element
            Material = MatList[self.MatN]
            sig, _, _ = Material.Sig( None, 1, TimeDelt, None, j, self, dps, eps, [], [], [])
            f   = JJ*self.Geom[1,0]*SampleWeight[self.IntT,self.nInt-1,j]                      # weighting factor
            fw  = fw  + f*dot(transpose(BBw),sig)                                       # fw: integrated discontinuity forces from bulk stresses
        return fw
    
    def tria_Area(self, c, p1, p2):
        cp1 = array([p1[0]-c[0], p1[1]-c[1], p1[2]-c[2]])
        cp2 = array([p2[0]-c[0], p2[1]-c[1], p2[2]-c[2]])
        Area = 0.5*norm(cross(cp1, cp2))
        return Area
    #
    def CrComp(self, TimeDelt, ff):
        Kdd,fd = zeros((6,6),dtype=float), zeros((6),dtype=float) 
        for j in range(self.n_vertices):                                    # loop over triangles of cutting plane
            xyzI = self.wwIntPCoor[j]                                       # coordinates of crack integration point within crack face - comes from CreateSDARankine
            Nw  = self.FormNw( self.centCoor, xyzI)                         # shape function for embedded discontinuity
            wwL = dot(self.CrTT,dot(Nw,self.ww))                            # crack width global system =  dot(Nw,self.ww)
                                                                            # crack width local system = dot( self.CrTT,dot(Nw,self.ww))
            dwL = dot(self.CrTT,dot(Nw,self.dw))                            # crack width increment global system --> local system
            self.CrStatus( j, wwL, None)                                    # added uhc 24-01-09
            tr, TW = self.CrTraction( j, wwL, ff)
            # viscous contribution to avoid roughness
            if self.CrackEta>0.:
                zz, Vw = self.ViscExtenCr3D( TimeDelt, self.CrackEta, dwL, 13, j)        # both corresponding state variable values are currently assigned to integration point 0
                tr[0] = tr[0] + self.CrackEta*Vw[0]
                TW[0,0] = TW[0,0] + zz 
# shear viscosity currently disregarded 
#               tr[1] = tr[1] + self.CrackEta*Vw[1]
#               TW[1,1] = TW[1,1] + zz 
#               tr[2] = tr[2] + self.CrackEta*Vw[2]
#               TW[2,2] = TW[2,2] + zz 
            self.wwL[j,:] = wwL[:]                                          # store for output
            self.wwT[j,:] = tr[:]                                           # "
            self.wwLn[j]  = wwL[0]                                           #
            self.wwTn[j]  = tr[0]                                           # normal components
            TW_ = dot(transpose(self.CrTT),dot(TW,self.CrTT))               # transformation of local tangential stiffness into global
            f  = self.wwArea[j] * self.CrGeRa
            Kdd = Kdd + f*dot(transpose(Nw),dot(TW_,Nw))
            fd  = fd  + f*dot(transpose(Nw),dot(transpose(self.CrTT),tr))   # fd: integrated discontinuity forces from crack tractions
        return Kdd, fd
    #
    def CrStatus(self, j, ww, iiter ): 
        # loading
        if ww[0]>=(self.wwLM[j,1]-1.0e-12):
            self.CrackStatus[j] = 0
        # unloading, reloading
        elif ww[0]>=0.:
            self.CrackStatus[j] = 1
        # crack closure
        else:
            self.CrackStatus[j] = 2
    #  local crack tractions      
    def CrTraction(self, j, ww, ff ):
#        sf = 1.0e-6                                                                     # factor for local shear stiffness -- be careful with this!
#        sf = 1.0e-1
        sf = self.ShearRetFactor                                                                     # factor for local shear stiffness -- be careful with this!
        tt = zeros((3), dtype=float)
        dt = zeros((3,3), dtype=float)
        wI = 1./self.wLim
        x  = wI * ww
        # normal - quadratic
        if x[0]<1.:                                                                     # critical crack width not exceeded
            tI = self.CrTL[0]                                                           # crack traction upon crack initiation
            x1 = self.wwLM[j,0]                                                         # starting point of unloading
            t1 = self.wwTM[j,0]                                                         # crack traction of starting point
            w  = ww[0]
            beta1= 0.1 #0.1 # 0.1                                                       # ratio of horizontal point offset (to left) to unloading starting point  
            beta = 0. # 0. # 0.1                                                        # ratio of current unloading starting point to start with crack closure transition
            dd = beta1*x1                                                               # point with horizontal course of transition
            alpha = 3.                                                                  # factor for left end (--> x1-alpha*dd) of transition
            xx = x1-alpha*dd                                                            # anchor point for anchor left of smooth unloading transition
            # loading
            if self.CrackStatus[j]==0:
                tt[0]   =    tI   *(x[0]*x[0] -2.*x[0] + 1.)
                dt[0,0] = 2.*tI*wI*(x[0]-1.)
                # linear
#                tt[1]   =    tI*(1.-x[1])
#                dt[1,1] =   -tI*wI
                #
            # unloading - reloading
            elif self.CrackStatus[j]==1:
                # smoothed transition
                if dd > 0.: # ZeroD:
                    denom = dd**2*(-6*alpha**2*dd + 2*dd-3*x1*alpha+3*alpha**2*x1+6*alpha*dd) # denom might become very small, don't use ZeroD to catch for zero division
                    x1_ = x1*wI                                                         # related starting point of unloading
                    yy1 =    tI   *(x1_*x1_ -2.*x1_ + 1.)                               # traction of unloading starting point
                    dy1 = 2.*tI*wI*(x1_-1.)                                             # derivative / tangential stiffness
                    # linear
#                    yy1 =    tI*(1.-x1_)
#                    dy1 =   -tI*wI
                    #
                    aa  = (yy1-2*dd*alpha*dy1+2*dd*dy1+x1*alpha*dy1-x1*dy1)/denom
                    bb  = -0.5*(6*x1*yy1-12*x1*alpha*dd*dy1+15*dd*x1*dy1+6*x1**2*alpha*dy1-6*x1**2*dy1-3*dd*yy1-8*dd**2*dy1+6*dd**2*dy1*alpha**2-3*dd*dy1*alpha**2*x1)/denom
                    cc  = (3*x1**2*yy1-6*x1**2*dd*alpha*dy1+9*x1**2*dy1*dd+3*x1**3*alpha*dy1-3*x1**3*dy1-3*x1*dd*yy1-8*dd**2*x1*dy1+9*x1*dd**2*dy1*alpha**2-3*dd*dy1*alpha**2*x1**2-6*dd**3*dy1*alpha**2+6*dd**3*dy1*alpha-3*dd**2*dy1*x1*alpha+2*dd**3*dy1)/denom
                    dd_ = -0.5*(2*x1**3*yy1-4*x1**3*alpha*dd*dy1+7*x1**3*dd*dy1+2*x1**4*alpha*dy1-2*x1**4*dy1-3*x1**2*dd*yy1-6*x1**2*alpha*dd**2*dy1-8*x1**2*dd**2*dy1-12*dd**3*x1*dy1*alpha**2+12*dd**3*x1*dy1*alpha+12*dd**2*x1**2*dy1*alpha**2+4*dd**3*x1*dy1+12*dd**3*yy1*alpha**2-12*dd**3*yy1*alpha-6*dd**2*yy1*alpha**2*x1+6*dd**2*yy1*x1*alpha-4*dd**3*yy1-3*x1**3*dy1*dd*alpha**2)/denom
                    # smooth transition 
                    if w > xx:                                                 # checks left end of transition from loading path  
                        tt[0]   =    aa*w**3 +    bb*w**2 + cc*w +dd_
                        dt[0,0] = 3.*aa*w**2 + 2.*bb*w    + cc
                    # intermediate linear to origin
                    elif w > beta*xx:                                        # left end of transition from loading path x beta
                        yy = aa*xx**3 + bb*xx**2 + cc*xx +dd_                           # traction at anchor point
                        CC  = yy/xx
                        tt[0]   = CC * w
                        dt[0,0] = CC
                    # transition to crack closure -- smooth for beta>0
                    else:
                        yy = aa*xx**3 + bb*xx**2 + cc*xx +dd_
                        aa_ = -1./4.*(xx*self.BulkEmod-yy)/beta/xx**2
                        bb_ =  1./2.*(xx*self.BulkEmod+yy)/xx
                        cc_ = -1./4.*self.BulkEmod*beta*xx+1./4.*yy*beta
                        tt[0]   =    aa_*w**2 + bb_*w + cc_
                        dt[0,0] = 2.*aa_*w    + bb_ 
                # sharp transition
                else:
                    C = self.wwTM[j,1]/self.wwLM[j,1]
                    tt[0]   = C * w
                    dt[0,0] = C
            # crack closure    
            elif self.CrackStatus[j]==2:
                # smooth transition -- beta > 0 should be fulfilled here
                if w > -beta*x1:
                    aa_ = -1./4./beta/x1**2*(x1*self.BulkEmod-t1)
                    bb_ =  1./2.*(x1*self.BulkEmod+t1)/x1
                    cc_ = -1./4.*self.BulkEmod*beta*x1+1/4*t1*beta
                    tt[0]   = aa_*w**2 + bb_*w + cc_
                    dt[0,0] = 2.*aa_*w + bb_ 
                # linear
                else:
                    tt[0]   = self.BulkEmod*w
                    dt[0,0] = self.BulkEmod
#                sf = 1.                    # experiences convergence problems uhc 190210 
        # shear -- linear
        if x[1]<1.:                                                                     # critical crack width not exceeded
            tt[1]   =  sf*self.BulkEmod * ww[1]
            dt[1,1] =  sf*self.BulkEmod
        # shear -- linear
        if x[2]<1.:                                                                     # critical crack width not exceeded
            tt[2]   =  sf*self.BulkEmod * ww[2]
            dt[2,2] =  sf*self.BulkEmod
            
        return tt, dt

    def ViscExtenCr3D(self, Dt, eta, Dww, sI, ipI):                                      # sI base index, see also ElasticLt.__init__
        if Dt < ZeroD:                                                      # might also be the case with initial check of equilibrium
            Dt = 1.0
#            raise NameError("CaeFemElements:ViscExtenCr3D: time step zero")
        # determine current scalar velocity
        def EvalV( VwOld, Dww, Dt):                                    
            if VwOld<ZeroD or Dww*VwOld<0.: 
#                if Dt>ZeroD: VwOld = Dww/Dt
#                else:        VwOld = 0.
                VwOld = Dww/Dt
#            if Dt>ZeroD and Dww>ZeroD: Vw = 2.*Dww/Dt - VwOld               # actual crack width rate -- normal component
#            else:                      Vw = VwOld
            Vw = 2.*Dww/Dt - VwOld
            return Vw
        # end def
        VwOld0 = self.wwV[ipI,0]
        VwOld1 = self.wwV[ipI,1]
        VwOld2 = self.wwV[ipI,2]
        #
        Vw0 = EvalV( VwOld0, Dww[0], Dt)                 # normal component
        #print(VwOld0, VwOld1, VwOld2)
        Vw1 = EvalV( VwOld1, Dww[1], Dt)                 # shear component
        Vw2 = EvalV( VwOld2, Dww[2], Dt)                 # shear component
#        if Dt>ZeroD: zz = 2.*eta/Dt
#        else:        zz = 0.
        zz = 2.*eta/Dt
        self.wwVN[ipI,0] = Vw0
        self.wwVN[ipI,1] = Vw1
        self.wwVN[ipI,2] = Vw2
        return zz, [ Vw0, Vw1, Vw2]
    # ExtrapolateScalar3D.mws - will not work, coeff matrix is singular
#    def ExtrapolateIpToVert(self, IPCoor, IpVal, VertCoor):
#        if len(IPCoor)!= 4: raise NameError("ConFemElements::ElementC3DS.ExtrapolateIpToVert: something wrong with crack face")
#        VertVals = []
#        x1, y1, z1 = IPCoor[0][0], IPCoor[0][1], IPCoor[0][2]
#        x2, y2, z2 = IPCoor[1][0], IPCoor[1][1], IPCoor[1][2]
#        x3, y3, z3 = IPCoor[2][0], IPCoor[2][1], IPCoor[2][2]
#        x4, y4, z4 = IPCoor[3][0], IPCoor[3][1], IPCoor[3][2]
#        deno = -y3*x4*z1+y2*x3*z4-y2*x4*z3-y3*x2*z4+y3*x2*z1+y3*x4*z2-y4*x2*z1-x3*y4*z2+x3*y1*z2-x3*y1*z4-x4*y1*z2-y2*x3*z1+y2*x4*z1+y4*x3*z1+x2*y4*z3-x2*y1*z3+x2*y1*z4+x4*y1*z3+x1*y3*z4-x1*y2*z4+x1*y4*z2+x1*y2*z3-x1*y4*z3-x1*y3*z2
#        for i in IpVal:                                                     # IpVal[0]: normal crack width, [1] normal crack traction
#            v1, v2, v3, v4 = i[0], i[1], i[2], i[3]
#            a = (  -v1*y2*z4-v4*y3*z1-z3*y4*v1+z4*y3*v1+v3*y4*z1+z3*y4*v2-v3*y4*z2+y1*z3*v4+y1*z4*v2-y1*v4*z2-y1*v3*z4-y1*z3*v2+y1*v3*z2+z4*y2*v3-z4*y3*v2-z1*y2*v3-z1*y4*v2+z1*y2*v4+z1*y3*v2-v4*y2*z3+v4*y3*z2+v1*y4*z2+v1*y2*z3-v1*y3*z2) / deno
#            b = -(  x1*z3*v4+x1*z4*v2-x1*v4*z2-x1*v3*z4-x1*z3*v2+x1*v3*z2-z4*x2*v1+v4*x2*z1+z3*x4*v2-x4*z3*v1-z3*x2*v4+x4*v3*z1+v3*x2*z4-v3*x2*z1-v3*x4*z2+x3*z1*v2+x3*v4*z2-x3*z4*v2-x3*z1*v4-x4*z1*v2+x4*v1*z2+x3*v1*z4+z3*x2*v1-x3*v1*z2) / deno
#            c = (   x3*y1*v2+x4*y1*v3+y2*x3*v4+y4*x3*v1+y2*x4*v1+x1*y3*v4-x3*y4*v2+y3*x2*v1-y3*x4*v1-y2*x4*v3-x3*y1*v4-y2*x3*v1+x1*y2*v3-x1*y4*v3+x2*y4*v3-y4*x2*v1+x2*y1*v4+x1*y4*v2-x2*y1*v3-y3*x2*v4+y3*x4*v2-x1*y2*v4-x4*y1*v2-x1*y3*v2) / deno
#            d = (   v4*y3*x2*z1+v4*x3*y1*z2+z3*y4*x2*v1+z3*x4*y1*v2-z4*x3*y1*v2-z4*y3*x2*v1+x3*z1*y4*v2-x4*z1*y3*v2-v3*y4*x2*z1-v3*x4*y1*z2+x4*v1*y3*z2-x3*v1*y4*z2-y4*x1*z3*v2+y4*x1*v3*z2-y1*z3*x2*v4+y1*v3*x2*z4+z4*y2*x3*v1-z4*x1*y2*v3+z4*x1*y3*v2-z1*y2*x3*v4+z1*y2*x4*v3+v4*x1*y2*z3-v4*x1*y3*z2-v1*y2*x4*z3) / deno
#            vV0 = a*VertCoor[0][0] +b*VertCoor[0][1] + c*VertCoor[0][2] + d
#            vV1 = a*VertCoor[1][0] +b*VertCoor[1][1] + c*VertCoor[1][2] + d
#            vV2 = a*VertCoor[2][0] +b*VertCoor[2][1] + c*VertCoor[2][2] + d
#            vV3 = a*VertCoor[3][0] +b*VertCoor[3][1] + c*VertCoor[3][2] + d
#            VertVals += [[ vV0, vV1, vV2, vV3]]
#        return VertVals
    # update discontinuity data - rotate discontinuity if scheduled
    def Update(self, Mat, NodeList, NoIndToCMInd, ff, MsgLevel):
        self.wwOld[:] = self.ww[:]                                                      # crack width in global system
        # rotate discontinuity is not implemented for 3D so far
        # update max crack width
        for j in range(self.n_vertices):
            if self.wwL[j,0]>self.wwLM[j,0]: 
                self.wwLM[j,:] = self.wwL[j,:]                                          # update for largest crack width reached in loading history
                self.wwTM[j,:] = self.wwT[j,:]
            self.wwV[j,:] = self.wwVN[j,:]                                              # update local crack velocities
            if MsgLevel>0 and self.CrackStatus[j] != self.CrackStatusOld[j]:
                print("-- crack status changed ",self.Label,j,self.CrackStatusOld[j],"->",self.CrackStatus[j])
                print("-- crack status changed ",self.Label,j,self.CrackStatusOld[j],"->",self.CrackStatus[j], file=ff)
                self.CrackStatusOld[j] = self.CrackStatus[j]
        self.fwd[:] = 0.

class C3D8(ElementC3D):
    def __init__(self, Label, SetLabel, InzList, MatName,Material, NoList, StateV, NData, nI,RegType, NoLabToNoInd):
#        if RegType==1: Element.__init__(self,"C3D8",InzList, 3, 2,2,8, (set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7])), 3,True,Label,SetLabel,3, MatName,StateV,NData+2, NoList,NoLabToNoInd,[])
#        else:          Element.__init__(self,"C3D8",InzList, 3, 2,2,8, (set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3])),   3,True,Label,SetLabel,3, MatName,StateV,NData+2, NoList,NoLabToNoInd,[])
        Element.__init__(self,"C3D8",InzList, 3, 2,2,8, None,   3,True,Label,SetLabel,3, None,None,None, NoList,NoLabToNoInd,[])
        if nI!=2:
            self.nInt  = nI                                                             # integration order - modification of Element.__ini__ settings
            self.nIntL = nI*nI*nI                                                       # total number of integration points
            self.nIntLi= nI*nI*nI                                                       # total number of integration points

    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        i3 = NoIndToCMInd[self.Inzi[3]]
        i4 = NoIndToCMInd[self.Inzi[4]]
        i5 = NoIndToCMInd[self.Inzi[5]]
        i6 = NoIndToCMInd[self.Inzi[6]]
        i7 = NoIndToCMInd[self.Inzi[7]]
        n0, n1, n2, n3, n4, n5, n6, n7 = NoList[i0],NoList[i1],NoList[i2],NoList[i3],NoList[i4],NoList[i5],NoList[i6],NoList[i7]
        self.X0, self.Y0, self.Z0 = n0.XCo, n0.YCo, n0.ZCo
        self.X1, self.Y1, self.Z1 = n1.XCo, n1.YCo, n1.ZCo
        self.X2, self.Y2, self.Z2 = n2.XCo, n2.YCo, n2.ZCo
        self.X3, self.Y3, self.Z3 = n3.XCo, n3.YCo, n3.ZCo
        self.X4, self.Y4, self.Z4 = n4.XCo, n4.YCo, n4.ZCo
        self.X5, self.Y5, self.Z5 = n5.XCo, n5.YCo, n5.ZCo
        self.X6, self.Y6, self.Z6 = n6.XCo, n6.YCo, n6.ZCo
        self.X7, self.Y7, self.Z7 = n7.XCo, n7.YCo, n7.ZCo
        self.X0u = self.X0 
        self.X1u = self.X1 
        self.X2u = self.X2 
        self.X3u = self.X3 
        self.X4u = self.X4 
        self.X5u = self.X5 
        self.X6u = self.X6 
        self.X7u = self.X7 
        self.Y0u = self.Y0 
        self.Y1u = self.Y1 
        self.Y2u = self.Y2 
        self.Y3u = self.Y3 
        self.Y4u = self.Y4 
        self.Y5u = self.Y5 
        self.Y6u = self.Y6 
        self.Y7u = self.Y7 
        self.Z0u = self.Z0 
        self.Z1u = self.Z1 
        self.Z2u = self.Z2 
        self.Z3u = self.Z3 
        self.Z4u = self.Z4 
        self.Z5u = self.Z5 
        self.Z6u = self.Z6 
        self.Z7u = self.Z7 
        #
        mat = SecDict[self.Set].Mat
        self.MatN = mat
        Material  = MaList[mat]
        self.Material = Material
        self.RegType  = Material.RType
        if self.RegType in [1]: self.DofT = (set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]),set([1,2,3,7]))       # tuple, type of dof for every node of this element: 1 -> u_x, 7->gradient field
        else:                   self.DofT = (set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]),  set([1,2,3]))
        self.nNod = len(self.DofT)                                          # number of nodes
        if self.RegType in [1]: self.DofI = zeros( (self.nNod,4), dtype=int)
        else:                   self.DofI = zeros( (self.nNod,3), dtype=int)                    # indices of global dofs per node
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 1.                                                 # dummy
        self.Geom[1,0] = 1.                                                 # dummy
        #
        AA, it, io = 0., self.IntT, self.nInt-1                             # integration type, integration order
        for i in range(self.nIntLi):                                        # volume by numerical integration for characteristic length  
            r = SamplePoints[it,io,i][0]
            s = SamplePoints[it,io,i][1]
            t = SamplePoints[it,io,i][2]
            JJ = self.JacoD(r,s,t) 
            if JJ<= 0.0: raise NameError("ConFemElements::C3D84.Ini2: something wrong with geometry",self.Label,i,r,s,t,JJ)
            f = JJ*SampleWeight[it,io,i]
            AA = AA + f
        self.Lch = AA**(1./3.)/self.nInt                     # characteristic length
        self.Lch_ = AA**(1./3.)                              # previous approach might be misleading
        self.CrBScaleType()
#        if mat.RType ==2:                                    # find scaling factor for band width regularization 
#            x = mat.bw/self.Lch_
#            CrX, CrY = mat.CrX, mat.CrY                      # support points for tensile softening scaling factor 
#            i = bisect_left(CrX, x) 
#            if i>0 and i<(mat.CrBwN+1):
#                self.CrBwS = CrY[i-1] + (x-CrX[i-1])/(CrX[i]-CrX[i-1])*(CrY[i]-CrY[i-1]) # scling factor by linear interpolation
#            else:
#                print('ZZZ', self.Lch_,x,'\n', CrX, '\n', CrY, i, mat.CrBwN)  
#                raise NameError("ConFemElem:C3D8.Ini2: RType 2 - element char length exceeds scaling factor interpolation")
        return []
    def UpdateCoord(self, dis, ddis ):
        self.X0 = self.X0u + dis[0] 
        self.Y0 = self.Y0u + dis[1] 
        self.Z0 = self.Z0u + dis[2] 
        self.X1 = self.X1u + dis[3] 
        self.Y1 = self.Y1u + dis[4] 
        self.Z1 = self.Z1u + dis[5] 
        self.X2 = self.X2u + dis[6] 
        self.Y2 = self.Y2u + dis[7] 
        self.Z2 = self.Z2u + dis[8] 
        self.X3 = self.X3u + dis[9] 
        self.Y3 = self.Y3u + dis[10] 
        self.Z3 = self.Z3u + dis[11] 
        self.X4 = self.X4u + dis[12] 
        self.Y4 = self.Y4u + dis[13] 
        self.Z4 = self.Z4u + dis[14] 
        self.X5 = self.X5u + dis[15] 
        self.Y5 = self.Y5u + dis[16] 
        self.Z5 = self.Z5u + dis[17] 
        self.X6 = self.X6u + dis[18] 
        self.Y6 = self.Y6u + dis[19] 
        self.Z6 = self.Z6u + dis[20] 
        self.X7 = self.X7u + dis[21] 
        self.Y7 = self.Y7u + dis[22] 
        self.Z7 = self.Z7u + dis[23] 
        return
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((self.DofE, self.DofE), dtype=float)
        B, _ = self.Basics( r, s, t)
        for i in range(8):
            ii = 3*i
            for j in range(8):
                jj = 3*j
                HH = sig[0]*B[0,i]*B[0,j] + sig[1]*B[1,i]*B[1,j] + sig[2]*B[2,i]*B[2,j] + sig[3]*(B[1,i]*B[2,j]+B[2,i]*B[1,j]) + sig[4]*(B[0,i]*B[2,j]+B[2,i]*B[0,j]) + sig[5]*(B[0,i]*B[1,j]+B[1,i]*B[0,j])
                for k in range(3):
                    GeomK[ii+k,jj+k] = HH
        return GeomK
    def Basics(self, r, s, t):
        br0, bs0, bt0 = -0.125*(1.-s)*(1.-t), -0.125*(1.-r)*(1.-t), -0.125*(1.-r)*(1.-s)
        br1, bs1, bt1 =  0.125*(1.-s)*(1.-t), -0.125*(1.+r)*(1.-t), -0.125*(1.+r)*(1.-s)
        br2, bs2, bt2 =  0.125*(1.+s)*(1.-t),  0.125*(1.+r)*(1.-t), -0.125*(1.+r)*(1.+s)
        br3, bs3, bt3 = -0.125*(1.+s)*(1.-t),  0.125*(1.-r)*(1.-t), -0.125*(1.-r)*(1.+s)
        br4, bs4, bt4 = -0.125*(1.-s)*(1.+t), -0.125*(1.-r)*(1.+t),  0.125*(1.-r)*(1.-s)
        br5, bs5, bt5 =  0.125*(1.-s)*(1.+t), -0.125*(1.+r)*(1.+t),  0.125*(1.+r)*(1.-s)
        br6, bs6, bt6 =  0.125*(1.+s)*(1.+t),  0.125*(1.+r)*(1.+t),  0.125*(1.+r)*(1.+s)
        br7, bs7, bt7 = -0.125*(1.+s)*(1.+t),  0.125*(1.-r)*(1.+t),  0.125*(1.-r)*(1.+s)
        JJ = zeros((3,3), dtype=float)
#        JJ[0,0] = br0*self.X0 + br1*self.X1 + br2*self.X2 + br3*self.X3 + br4*self.X4 + br5*self.X5 + br6*self.X6 + br7*self.X7
#        JJ[0,1] = bs0*self.X0 + bs1*self.X1 + bs2*self.X2 + bs3*self.X3 + bs4*self.X4 + bs5*self.X5 + bs6*self.X6 + bs7*self.X7
#        JJ[0,2] = bt0*self.X0 + bt1*self.X1 + bt2*self.X2 + bt3*self.X3 + bt4*self.X4 + bt5*self.X5 + bt6*self.X6 + bt7*self.X7
#        JJ[1,0] = br0*self.Y0 + br1*self.Y1 + br2*self.Y2 + br3*self.Y3 + br4*self.Y4 + br5*self.Y5 + br6*self.Y6 + br7*self.Y7
#        JJ[1,1] = bs0*self.Y0 + bs1*self.Y1 + bs2*self.Y2 + bs3*self.Y3 + bs4*self.Y4 + bs5*self.Y5 + bs6*self.Y6 + bs7*self.Y7
#        JJ[1,2] = bt0*self.Y0 + bt1*self.Y1 + bt2*self.Y2 + bt3*self.Y3 + bt4*self.Y4 + bt5*self.Y5 + bt6*self.Y6 + bt7*self.Y7
#        JJ[2,0] = br0*self.Z0 + br1*self.Z1 + br2*self.Z2 + br3*self.Z3 + br4*self.Z4 + br5*self.Z5 + br6*self.Z6 + br7*self.Z7
#        JJ[2,1] = bs0*self.Z0 + bs1*self.Z1 + bs2*self.Z2 + bs3*self.Z3 + bs4*self.Z4 + bs5*self.Z5 + bs6*self.Z6 + bs7*self.Z7
#        JJ[2,2] = bt0*self.Z0 + bt1*self.Z1 + bt2*self.Z2 + bt3*self.Z3 + bt4*self.Z4 + bt5*self.Z5 + bt6*self.Z6 + bt7*self.Z7

        JJ[0,0] = br0*self.X0 + br1*self.X1 + br2*self.X2 + br3*self.X3 + br4*self.X4 + br5*self.X5 + br6*self.X6 + br7*self.X7
        JJ[0,1] = br0*self.Y0 + br1*self.Y1 + br2*self.Y2 + br3*self.Y3 + br4*self.Y4 + br5*self.Y5 + br6*self.Y6 + br7*self.Y7
        JJ[0,2] = br0*self.Z0 + br1*self.Z1 + br2*self.Z2 + br3*self.Z3 + br4*self.Z4 + br5*self.Z5 + br6*self.Z6 + br7*self.Z7
        JJ[1,0] = bs0*self.X0 + bs1*self.X1 + bs2*self.X2 + bs3*self.X3 + bs4*self.X4 + bs5*self.X5 + bs6*self.X6 + bs7*self.X7
        JJ[1,1] = bs0*self.Y0 + bs1*self.Y1 + bs2*self.Y2 + bs3*self.Y3 + bs4*self.Y4 + bs5*self.Y5 + bs6*self.Y6 + bs7*self.Y7
        JJ[1,2] = bs0*self.Z0 + bs1*self.Z1 + bs2*self.Z2 + bs3*self.Z3 + bs4*self.Z4 + bs5*self.Z5 + bs6*self.Z6 + bs7*self.Z7
        JJ[2,0] = bt0*self.X0 + bt1*self.X1 + bt2*self.X2 + bt3*self.X3 + bt4*self.X4 + bt5*self.X5 + bt6*self.X6 + bt7*self.X7
        JJ[2,1] = bt0*self.Y0 + bt1*self.Y1 + bt2*self.Y2 + bt3*self.Y3 + bt4*self.Y4 + bt5*self.Y5 + bt6*self.Y6 + bt7*self.Y7
        JJ[2,2] = bt0*self.Z0 + bt1*self.Z1 + bt2*self.Z2 + bt3*self.Z3 + bt4*self.Z4 + bt5*self.Z5 + bt6*self.Z6 + bt7*self.Z7
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        JI = inv(JJ)
        BB = zeros((3,8), dtype=float)
        BB[0,0] = JI[0,0]*br0+JI[0,1]*bs0+JI[0,2]*bt0
        BB[0,1] = JI[0,0]*br1+JI[0,1]*bs1+JI[0,2]*bt1
        BB[0,2] = JI[0,0]*br2+JI[0,1]*bs2+JI[0,2]*bt2
        BB[0,3] = JI[0,0]*br3+JI[0,1]*bs3+JI[0,2]*bt3
        BB[0,4] = JI[0,0]*br4+JI[0,1]*bs4+JI[0,2]*bt4
        BB[0,5] = JI[0,0]*br5+JI[0,1]*bs5+JI[0,2]*bt5
        BB[0,6] = JI[0,0]*br6+JI[0,1]*bs6+JI[0,2]*bt6
        BB[0,7] = JI[0,0]*br7+JI[0,1]*bs7+JI[0,2]*bt7
        BB[1,0] = JI[1,0]*br0+JI[1,1]*bs0+JI[1,2]*bt0
        BB[1,1] = JI[1,0]*br1+JI[1,1]*bs1+JI[1,2]*bt1
        BB[1,2] = JI[1,0]*br2+JI[1,1]*bs2+JI[1,2]*bt2
        BB[1,3] = JI[1,0]*br3+JI[1,1]*bs3+JI[1,2]*bt3
        BB[1,4] = JI[1,0]*br4+JI[1,1]*bs4+JI[1,2]*bt4
        BB[1,5] = JI[1,0]*br5+JI[1,1]*bs5+JI[1,2]*bt5
        BB[1,6] = JI[1,0]*br6+JI[1,1]*bs6+JI[1,2]*bt6
        BB[1,7] = JI[1,0]*br7+JI[1,1]*bs7+JI[1,2]*bt7
        BB[2,0] = JI[2,0]*br0+JI[2,1]*bs0+JI[2,2]*bt0
        BB[2,1] = JI[2,0]*br1+JI[2,1]*bs1+JI[2,2]*bt1
        BB[2,2] = JI[2,0]*br2+JI[2,1]*bs2+JI[2,2]*bt2
        BB[2,3] = JI[2,0]*br3+JI[2,1]*bs3+JI[2,2]*bt3
        BB[2,4] = JI[2,0]*br4+JI[2,1]*bs4+JI[2,2]*bt4
        BB[2,5] = JI[2,0]*br5+JI[2,1]*bs5+JI[2,2]*bt5
        BB[2,6] = JI[2,0]*br6+JI[2,1]*bs6+JI[2,2]*bt6
        BB[2,7] = JI[2,0]*br7+JI[2,1]*bs7+JI[2,2]*bt7
        return BB, det
    def FormN(self, r, s, t):
        N = array([[(1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125, 0, 0],
                   [0, (1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125, 0],
                   [0, 0, (1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125]])
        return N
    def FormB(self, r, s, t, NLg):
        BB, det = self.Basics( r, s, t)
        if self.RegType==1:
            B = array([[ BB[0,0], 0,   0,   0,  BB[0,1], 0,   0,   0,  BB[0,2], 0,   0,   0,  BB[0,3], 0,   0,   0,  BB[0,4], 0,   0,   0,  BB[0,5], 0,   0,   0,  BB[0,6], 0,   0,   0,  BB[0,7], 0,   0,  0],
                       [ 0,   BB[1,0], 0,   0,  0,   BB[1,1], 0,   0,  0,   BB[1,2], 0,   0,  0,   BB[1,3], 0,   0,  0,   BB[1,4], 0,   0,  0,   BB[1,5], 0,   0,  0,   BB[1,6], 0,   0,  0,   BB[1,7], 0,  0],
                       [ 0,   0,   BB[2,0], 0,  0,   0,   BB[2,1], 0,  0,   0,   BB[2,2], 0,  0,   0,   BB[2,3], 0,  0,   0,   BB[2,4], 0,  0,   0,   BB[2,5], 0,  0,   0,   BB[2,6], 0,  0,   0,   BB[2,7],0],
                       [ 0,   BB[2,0], BB[1,0], 0,  0,   BB[2,1], BB[1,1], 0,  0,   BB[2,2], BB[1,2], 0,  0,   BB[2,3], BB[1,3], 0,  0,   BB[2,4], BB[1,4], 0,  0,   BB[2,5], BB[1,5], 0,  0,   BB[2,6], BB[1,6], 0,  0,   BB[2,7], BB[1,7],0],
                       [ BB[2,0], 0,   BB[0,0], 0,  BB[2,1], 0,   BB[0,1], 0,  BB[2,2], 0,   BB[0,2], 0,  BB[2,3], 0,   BB[0,3], 0,  BB[2,4], 0,   BB[0,4], 0,  BB[2,5], 0,   BB[0,5], 0,  BB[2,6], 0,   BB[0,6], 0,  BB[2,7], 0,   BB[0,7],0],
                       [ BB[1,0], BB[0,0], 0,   0,  BB[1,1], BB[0,1], 0,   0,  BB[1,2], BB[0,2], 0,   0,  BB[1,3], BB[0,3], 0,   0,  BB[1,4], BB[0,4], 0,   0,  BB[1,5], BB[0,5], 0,   0,  BB[1,6], BB[0,6], 0,   0,  BB[1,7], BB[0,7], 0,  0],
                       [ 0,   0,   0,   BB[0,0],0,   0,   0,   BB[0,1],0,   0,   0,   BB[0,2],0,   0,   0,   BB[0,3],0,   0,   0,   BB[0,4],0,   0,   0,   BB[0,5],0,   0,   0,   BB[0,6],0,   0,   0,  BB[0,7]],
                       [ 0,   0,   0,   BB[1,0],0,   0,   0,   BB[1,1],0,   0,   0,   BB[1,2],0,   0,   0,   BB[1,3],0,   0,   0,   BB[1,4],0,   0,   0,   BB[1,5],0,   0,   0,   BB[1,6],0,   0,   0,  BB[1,7]],
                       [ 0,   0,   0,   BB[2,0],0,   0,   0,   BB[2,1],0,   0,   0,   BB[2,2],0,   0,   0,   BB[2,3],0,   0,   0,   BB[2,4],0,   0,   0,   BB[2,5],0,   0,   0,   BB[2,6],0,   0,   0,  BB[2,7]]])
            N0 = (1.-r)*(1.-s)*(1.-t)*0.125
            N1 = (1.+r)*(1.-s)*(1.-t)*0.125
            N2 = (1.+r)*(1.+s)*(1.-t)*0.125
            N3 = (1.-r)*(1.+s)*(1.-t)*0.125
            N4 = (1.-r)*(1.-s)*(1.+t)*0.125
            N5 = (1.+r)*(1.-s)*(1.+t)*0.125
            N6 = (1.+r)*(1.+s)*(1.+t)*0.125 
            N7 = (1.-r)*(1.+s)*(1.+t)*0.125
            BN =array([[BB[0,0],0,      0,   0, BB[0,1], 0,   0,   0, BB[0,2], 0,   0,   0, BB[0,3], 0,   0,   0, BB[0,4], 0,   0,   0, BB[0,5], 0,   0,   0, BB[0,6], 0,   0,   0, BB[0,7], 0,   0,  0],
                       [0,      BB[1,0],0,   0, 0,   BB[1,1], 0,   0, 0,   BB[1,2], 0,   0, 0,   BB[1,3], 0,   0, 0,   BB[1,4], 0,   0, 0,   BB[1,5], 0,   0, 0,   BB[1,6], 0,   0, 0,   BB[1,7], 0,  0],
                       [0,      0,      BB[2,0], 0, 0,   0,   BB[2,1], 0, 0,   0,   BB[2,2], 0, 0,   0,   BB[2,3], 0, 0,   0,   BB[2,4], 0, 0,   0,   BB[2,5], 0, 0,   0,   BB[2,6], 0, 0,   0,   BB[2,7],0],
                       [0,      BB[2,0],BB[1,0], 0, 0,   BB[2,1], BB[1,1], 0, 0,   BB[2,2], BB[1,2], 0, 0,   BB[2,3], BB[1,3], 0, 0,   BB[2,4], BB[1,4], 0, 0,   BB[2,5], BB[1,5], 0, 0,   BB[2,6], BB[1,6], 0, 0,   BB[2,7], BB[1,7],0],
                       [BB[2,0],0,      BB[0,0], 0, BB[2,1], 0,   BB[0,1], 0, BB[2,2], 0,   BB[0,2], 0, BB[2,3], 0,   BB[0,3], 0, BB[2,4], 0,   BB[0,4], 0, BB[2,5], 0,   BB[0,5], 0, BB[2,6], 0,   BB[0,6], 0, BB[2,7], 0,   BB[0,7],0],
                       [BB[1,0],BB[0,0], 0,   0, BB[1,1], BB[0,1], 0,   0, BB[1,2], BB[0,2], 0,   0, BB[1,3], BB[0,3], 0,   0, BB[1,4], BB[0,4], 0,   0, BB[1,5], BB[0,5], 0,   0, BB[1,6], BB[0,6], 0,   0, BB[1,7], BB[0,7], 0,  0],
                       [0,      0,   0,   N0,0,   0,   0,   N1,0,   0,   0,   N2,0,   0,   0,   N3,0,   0,   0,   N4,0,   0,   0,   N5,0,   0,   0,   N6,0,   0,   0,  N7]])
            return B, BN, det, 0
        else:
            B = array([[BB[0,0],0,      0,      BB[0,1],0,      0,      BB[0,2], 0,   0,   BB[0,3], 0,   0,   BB[0,4], 0,   0,   BB[0,5], 0,   0,   BB[0,6], 0,   0,   BB[0,7], 0,   0  ],
                       [0,      BB[1,0],0,      0,      BB[1,1],0,      0,   BB[1,2], 0,   0,   BB[1,3], 0,   0,   BB[1,4], 0,   0,   BB[1,5], 0,   0,   BB[1,6], 0,   0,   BB[1,7], 0  ],
                       [0,      0,      BB[2,0],0,      0,      BB[2,1],0,   0,   BB[2,2], 0,   0,   BB[2,3], 0,   0,   BB[2,4], 0,   0,   BB[2,5], 0,   0,   BB[2,6], 0,   0,   BB[2,7]],
                       [0,      BB[2,0],BB[1,0],0,      BB[2,1],BB[1,1],0,   BB[2,2], BB[1,2], 0,   BB[2,3], BB[1,3], 0,   BB[2,4], BB[1,4], 0,   BB[2,5], BB[1,5], 0,   BB[2,6], BB[1,6], 0,   BB[2,7], BB[1,7]],
                       [BB[2,0],0,      BB[0,0],BB[2,1],0,      BB[0,1],BB[2,2], 0,   BB[0,2], BB[2,3], 0,   BB[0,3], BB[2,4], 0,   BB[0,4], BB[2,5], 0,   BB[0,5], BB[2,6], 0,   BB[0,6], BB[2,7], 0,   BB[0,7]],
                       [BB[1,0],BB[0,0],0,      BB[1,1],BB[0,1],0,      BB[1,2], BB[0,2], 0,   BB[1,3], BB[0,3], 0,   BB[1,4], BB[0,4], 0,   BB[1,5], BB[0,5], 0,   BB[1,6], BB[0,6], 0,   BB[1,7], BB[0,7], 0  ]])
            return B, det, 0
    def FormT(self, r, s, t):                                 # interpolation on temperature
        # !!! ordering might not yet be correct !!!
        if self.RegType==1: T = array([(1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, 0,(1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, 0,(1+r)*(1+s)*(1.-t)*0.125, 0, 0, 0,(1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, 0,(1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, 0,(1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, 0,(1+r)*(1+s)*(1.+t)*0.125, 0, 0, 0,(1.-r)*(1.+s)*(1.+t)*0.125, 0, 0, 0])
        else:               T = array([(1.-r)*(1.-s)*(1.-t)*0.125, 0, 0,   (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0,   (1+r)*(1+s)*(1.-t)*0.125, 0, 0,   (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0,   (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0,   (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0,   (1+r)*(1+s)*(1.+t)*0.125, 0, 0,   (1.-r)*(1.+s)*(1.+t)*0.125, 0, 0])
        return T
    def FormX(self, r, s, t):
        X = array([(1.-r)*(1.-s)*(1.-t)*0.125, (1.+r)*(1.-s)*(1.-t)*0.125, (1+r)*(1+s)*(1.-t)*0.125, (1.-r)*(1.+s)*(1.-t)*0.125, (1.-r)*(1.-s)*(1.+t)*0.125, (1.+r)*(1.-s)*(1.+t)*0.125, (1+r)*(1+s)*(1.+t)*0.125, (1.-r)*(1.+s)*(1.+t)*0.125])
        return X
    def FormX_(self, r, s, t):
        X = array([[(1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125, 0, 0],
                   [0, (1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125, 0],
                   [0, 0, (1.-r)*(1.-s)*(1.-t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.-t)*0.125, 0, 0, (1+r)*(1+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.-t)*0.125, 0, 0, (1.-r)*(1.-s)*(1.+t)*0.125, 0, 0, (1.+r)*(1.-s)*(1.+t)*0.125, 0, 0, (1+r)*(1+s)*(1.+t)*0.125, 0, 0, (1.-r)*(1.+s)*(1.+t)*0.125]])
        return X
    def JacoD(self, r, s, t):
        br0, bs0, bt0 = -0.125*(1.-s)*(1.-t), -0.125*(1.-r)*(1.-t), -0.125*(1.-r)*(1.-s)
        br1, bs1, bt1 =  0.125*(1.-s)*(1.-t), -0.125*(1.+r)*(1.-t), -0.125*(1.+r)*(1.-s)
        br2, bs2, bt2 =  0.125*(1.+s)*(1.-t),  0.125*(1.+r)*(1.-t), -0.125*(1.+r)*(1.+s)
        br3, bs3, bt3 = -0.125*(1.+s)*(1.-t),  0.125*(1.-r)*(1.-t), -0.125*(1.-r)*(1.+s)
        br4, bs4, bt4 = -0.125*(1.-s)*(1.+t), -0.125*(1.-r)*(1.+t),  0.125*(1.-r)*(1.-s)
        br5, bs5, bt5 =  0.125*(1.-s)*(1.+t), -0.125*(1.+r)*(1.+t),  0.125*(1.+r)*(1.-s)
        br6, bs6, bt6 =  0.125*(1.+s)*(1.+t),  0.125*(1.+r)*(1.+t),  0.125*(1.+r)*(1.+s)
        br7, bs7, bt7 = -0.125*(1.+s)*(1.+t),  0.125*(1.-r)*(1.+t),  0.125*(1.-r)*(1.+s)
        JJ = zeros((3,3), dtype=float)
        JJ[0,0] = br0*self.X0 + br1*self.X1 + br2*self.X2 + br3*self.X3 + br4*self.X4 + br5*self.X5 + br6*self.X6 + br7*self.X7
        JJ[0,1] = bs0*self.X0 + bs1*self.X1 + bs2*self.X2 + bs3*self.X3 + bs4*self.X4 + bs5*self.X5 + bs6*self.X6 + bs7*self.X7
        JJ[0,2] = bt0*self.X0 + bt1*self.X1 + bt2*self.X2 + bt3*self.X3 + bt4*self.X4 + bt5*self.X5 + bt6*self.X6 + bt7*self.X7
        JJ[1,0] = br0*self.Y0 + br1*self.Y1 + br2*self.Y2 + br3*self.Y3 + br4*self.Y4 + br5*self.Y5 + br6*self.Y6 + br7*self.Y7
        JJ[1,1] = bs0*self.Y0 + bs1*self.Y1 + bs2*self.Y2 + bs3*self.Y3 + bs4*self.Y4 + bs5*self.Y5 + bs6*self.Y6 + bs7*self.Y7
        JJ[1,2] = bt0*self.Y0 + bt1*self.Y1 + bt2*self.Y2 + bt3*self.Y3 + bt4*self.Y4 + bt5*self.Y5 + bt6*self.Y6 + bt7*self.Y7
        JJ[2,0] = br0*self.Z0 + br1*self.Z1 + br2*self.Z2 + br3*self.Z3 + br4*self.Z4 + br5*self.Z5 + br6*self.Z6 + br7*self.Z7
        JJ[2,1] = bs0*self.Z0 + bs1*self.Z1 + bs2*self.Z2 + bs3*self.Z3 + bs4*self.Z4 + bs5*self.Z5 + bs6*self.Z6 + bs7*self.Z7
        JJ[2,2] = bt0*self.Z0 + bt1*self.Z1 + bt2*self.Z2 + bt3*self.Z3 + bt4*self.Z4 + bt5*self.Z5 + bt6*self.Z6 + bt7*self.Z7
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        return det
    def FormNI(self):                                         # N-inverse for gaussian 3D order 2x2x2 integration, see N_Inverse.py, to get nodal values from integration point values
        return array([[ 2.54903811, -0.6830127,  -0.6830127,   0.1830127,  -0.6830127,   0.1830127,   0.1830127,  -0.04903811],
                      [-0.6830127,   0.1830127,   0.1830127,  -0.04903811,  2.54903811, -0.6830127,  -0.6830127,   0.1830127 ],
                      [ 0.1830127,  -0.04903811, -0.6830127,   0.1830127,  -0.6830127,   0.1830127,   2.54903811, -0.6830127 ],
                      [-0.6830127,   0.1830127,   2.54903811, -0.6830127,   0.1830127,  -0.04903811, -0.6830127,   0.1830127 ],
                      [-0.6830127,   2.54903811,  0.1830127,  -0.6830127,   0.1830127,  -0.6830127,  -0.04903811,  0.1830127 ],
                      [ 0.1830127,  -0.6830127,  -0.04903811,  0.1830127,  -0.6830127,   2.54903811,  0.1830127,  -0.6830127 ],
                      [-0.04903811,  0.1830127,   0.1830127,  -0.6830127,   0.1830127,  -0.6830127,  -0.6830127,   2.54903811],
                      [ 0.1830127,  -0.6830127,  -0.6830127,   2.54903811, -0.04903811,  0.1830127,   0.1830127,  -0.6830127 ]])

class C3D8_SDA(C3D8,ElementC3D_SDA):
    def __init__(self, Label, InzList, thickness, NoList, v_num, Coor, Normal, elSet, MatName, NoLabToNoInd,NoIndToCMInd, MatList, SecDict):
        Material = MatList[MatName]
#       C3D8.__init__(self, Label, elSet,    InzList, MatName,Material, NoList, Material.StateVar, Material.NData, 3, None, NoLabToNoInd)
        C3D8.__init__(self, Label, elSet,    InzList, MatName,Material, NoList, Material.StateVar, Material.NData, 2, None, NoLabToNoInd)
        ElementC3D_SDA.__init__(self, Coor, Normal, Material, v_num)
        self.Type = "C3D8S"
        # 
        self.Ini2( NoList, NoIndToCMInd, MatList, SecDict)
        self.ElemDimData( NoList,NoIndToCMInd )
        self.SDANew = False
        ElementC3D_SDA.__init__(self, Coor, Normal, Material, v_num)

    # Nw-matrix of embedded discontinuity
    def FormNw(self, center, point):
        m1 = self.CrTT[1]
        m2 = self.CrTT[2]
        vec = array([point[0] - center[0], point[1] - center[1], point[2] - center[2]])
        s1 = dot(vec, m1)
        s2 = dot(vec, m2)
        s = array([[1.0, 0., 0., s1, s2, 0.],   # displacement discontinuity shape function in local coordinates
                   [0., 1.0, 0., 0., 0., s2],
                   [0., 0., 1.0, 0., 0., s1]])
        Nw = dot(transpose(self.CrTT), s)                   # displacement discontinuity shape function in global coordinates
        return Nw
    # Mw-matrix of embedded discontinuity
    def FormMw(self, x0, y0, z0):
        rr, HH = 0.5, []                                                                # for subdivision of jump on positive and negative part 
        Mw, k = [], 0
        xx, yy, zz = [self.X0,self.X1,self.X2,self.X3,self.X4,self.X5,self.X6,self.X7], [self.Y0,self.Y1,self.Y2,self.Y3,self.Y4,self.Y5,self.Y6,self.Y7], [self.Z0,self.Z1,self.Z2,self.Z3,self.Z4,self.Z5,self.Z6,self.Z7]   # nodal coordinates

        for i in range(len(xx)):                                                        # loop over nodes
            if dot( array([self.CrNx,self.CrNy,self.CrNz]), array([xx[i]-self.CrXCo, yy[i]-self.CrYCo, zz[i]-self.CrZCo]) )>0.:    # determine which node is on the positive part and which on the negative for matrix H
                HH += [rr, rr, rr]
            else:
                HH += [rr-1., rr-1., rr-1]
            Nw = self.FormNw([x0,y0, z0], [xx[i],yy[i],zz[i]])
            for j in range(3):
                Mw += [ HH[k]*Nw[j] ]
                k += 1
        return array(Mw)
    # Bw-matrix of embedded discontinuity
    def FormBw(self, r, s, t):
        BB, detJ, _ = self.FormB( r, s, t, None)
        return dot(BB,self.FormMw(self.centCoor[0], self.centCoor[1], self.centCoor[2]))
    # shape function of embedded discontinuity geometry
    def FormXCut(self):
        r, s = 1/3, 1/3             # local integration sampling point coordinates
        L1 = 1-r-s
        L2 = r
        L3 = s
        X = array([[L1, 0, 0, L2, 0, 0, L3, 0, 0],   # shape function of triangular element
                   [0, L1, 0, 0, L2, 0, 0, L3, 0],
                   [0, 0, L1, 0, 0, L2, 0, 0, L3]])   # geometry discontinuity shape function 
        return X

class SB3(Element):
    def __init__(self, Label, SetLabel, InzList, MatName, NoList, ShellSecDic, StateV, NData, NoLabToNoInd):
#        Element.__init__(self,"SB3",InzList,  2,                   3,2,3,                    (set([3, 4, 5]),set([3, 4, 5]),set([3, 4, 5])), 20,False, NoLabToNoInd)
        Element.__init__(self,"SB3",InzList,   2,                   3,3,4,                    (set([3, 4, 5]),set([3, 4, 5]),set([3, 4, 5])), 20,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[])
#        Element.__init__(self,"SB3",InzList,  2,                   3,4,3,                    (set([3, 4, 5]),set([3, 4, 5]),set([3, 4, 5])), 20,False, NoLabToNoInd)
#                       (self, TypeVal,nNodVal,DofEVal,nFieVal, IntTVal,nIntVal,nIntLVal, DofTVal, DofNVal, dimVal, NLGeomIVal):
        self.DofI = zeros( (self.nNod,3), dtype=int)    # indices of global dofs per node
        self.TensStiff = False                          # flag for tension stiffening
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        X0 = NoList[i0].XCo
        Y0 = NoList[i0].YCo
        X1 = NoList[i1].XCo
        Y1 = NoList[i1].YCo
        X2 = NoList[i2].XCo
        Y2 = NoList[i2].YCo
        self.AA = -Y0*X1+Y0*X2+Y2*X1+Y1*X0-Y1*X2-Y2*X0  # double of element area
        if self.AA<=0.: raise NameError("ConFemElements::C3D84.Ini2: something wrong with geometry of SB3 element",self.Label)
        elset = self.Set
        self.Geom = zeros( (2,2), dtype=double)
        self.Geom[0,0] = 0.5*self.AA                                        # element area for numerical integration
        self.Geom[1,0] = 1                                                  # dummy for thickness
        self.Geom[1,1] = SecDict[elset].Height                              # used in ConFem.Materials
        self.Lch_ = 0.3*sqrt(0.5*self.AA)                                   # characteristic length
        l1=sqrt((X2-X1)**2+(Y2-Y1)**2)
        l2=sqrt((X0-X2)**2+(Y0-Y2)**2)
        l3=sqrt((X1-X0)**2+(Y1-Y0)**2)
        self.mu1=(l3**2-l2**2)/l1**2
        self.mu2=(l1**2-l3**2)/l2**2
        self.mu3=(l2**2-l1**2)/l3**2
        self.b1=(Y1-Y2)/self.AA
        self.b2=(Y2-Y0)/self.AA
        self.b3=(Y0-Y1)/self.AA
        self.c1=(X2-X1)/self.AA
        self.c2=(X0-X2)/self.AA
        self.c3=(X1-X0)/self.AA
        self.X0 = X0
        self.Y0 = Y0
        self.X1 = X1
        self.Y1 = Y1
        self.X2 = X2
        self.Y2 = Y2
        return []
    def FormN(self, L1, L2, L3):
        L3=1-L1-L2
        c1=self.c1
        c2=self.c2
        c3=self.c3
        b1=self.b1
        b2=self.b2
        b3=self.b3
        N = array([[0,0,0,0,0,0,0,0,0],
                  [L1-L1*L2+2*L1**2*L2+L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3)+L3*L1-2*L3**2*L1-L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2),
                   c2*(L3**2*L1+0.5*L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2))+c3*(L1**2*L2+0.5*L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3))-c2*L3*L1,
                   -b2*(L3**2*L1+0.5*L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2))-b3*(L1**2*L2+0.5*L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3))+b2*L3*L1,
                   -2*L1**2*L2-L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3)-L2*L3+L1*L2+L2+2*L2**2*L3+L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1),
                   -c3*L1*L2+c3*(L1**2*L2+0.5*L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3))+c1*(L2**2*L3+0.5*L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1)),
                   b3*L1*L2-b3*(L1**2*L2+0.5*L1*L2*L3*(3*(1-self.mu3)*L1-(1+3*self.mu3)*L2+(1+3*self.mu3)*L3))-b1*(L2**2*L3+0.5*L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1)),
                   L2*L3+2*L3**2*L1+L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2)-2*L2**2*L3-L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1)+L3-L3*L1,
                   c1*(L2**2*L3+0.5*L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1))+c2*(L3**2*L1+0.5*L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2))-c1*L2*L3,
                   -b1*(L2**2*L3+0.5*L1*L2*L3*(3*(1-self.mu1)*L2-(1+3*self.mu1)*L3+(1+3*self.mu1)*L1))-b2*(L3**2*L1+0.5*L1*L2*L3*(3*(1-self.mu2)*L3-(1+3*self.mu2)*L1+(1+3*self.mu2)*L2))+b1*L2*L3]])        
        return N
    def FormB(self, L1, L2, L3, NLg):
        L3=1-L1-L2
        c1 = self.c1
        c2 = self.c2
        c3 = self.c3
        b1 = self.b1
        b2 = self.b2
        b3 = self.b3
        mu1= self.mu1
        mu2= self.mu2
        mu3= self.mu3
#        B =array([[-2*(-3*self.b1*self.b2*self.mu3+2*self.b1*self.b2-3*self.b1*self.b2*self.mu2)*L3**2+(-2*(4*self.b3*self.b2-6*self.b2*self.b3*self.mu2+2*self.b2**2-8*self.b1*self.b2-6*self.b2*self.b3*self.mu3+3*self.b2**2*self.mu3+3*self.b2**2*self.mu2-6*self.b1*self.b2*self.mu2+6*self.b1*self.b2*self.mu3)*L1-2*(-3*self.b1**2*self.mu2+4*self.b1*self.b2+4*self.b1*self.b3+6*self.b1*self.b2*self.mu2-4*self.b1**2+6*self.b1*self.b2*self.mu3-6*self.b1*self.b3*self.mu2+3*self.b1**2*self.mu3-6*self.b1*self.b3*self.mu3)*L2-8*self.b1*self.b3)*L3-2*(2*self.b1*self.b3+3*self.b1*self.b3*self.mu2+3*self.b1*self.b3*self.mu3)*L2**2+(-2*(2*self.b3**2-3*self.b3**2*self.mu2-8*self.b1*self.b3-3*self.b3**2*self.mu3+4*self.b3*self.b2+6*self.b1*self.b3*self.mu3+6*self.b2*self.b3*self.mu2+6*self.b2*self.b3*self.mu3-6*self.b1*self.b3*self.mu2)*L1+4*self.b1**2)*L2-2*self.b1*self.b2+2*self.b1*self.b3-2*(-3*self.b2*self.b3*self.mu2-4*self.b3*self.b2+3*self.b2*self.b3*self.mu3)*L1**2-2*(-4*self.b1*self.b2+2*self.b3**2)*L1,
#                   -(3*self.b1*c2*self.b2*self.mu2-self.b1*self.b2*c3-3*self.b1*c3*self.b2*self.mu3-3*self.b1*c2*self.b2)*L3**2+(-(-c2*self.b2**2+6*self.b2*c2*self.b3*self.mu2+6*self.b1*c2*self.b2*self.mu2-6*self.b1*self.b2*c3+2*self.b1*c2*self.b2-6*self.b2*self.b3*c3*self.mu3-3*c2*self.b2**2*self.mu2+6*self.b1*c3*self.b2*self.mu3+3*c3*self.b2**2*self.mu3+c3*self.b2**2-6*self.b2*c2*self.b3-2*self.b2*self.b3*c3)*L1-(3*c2*self.b1**2*self.mu2-2*self.b1*c2*self.b2+c2*self.b1**2-6*self.b1*c2*self.b2*self.mu2+6*self.b1*c2*self.b3*self.mu2-6*self.b1*c3*self.b3*self.mu3-2*self.b1*c3*self.b3+2*self.b1*self.b2*c3+6*self.b1*c3*self.b2*self.mu3-3*c3*self.b1**2+3*c3*self.b1**2*self.mu3-6*self.b1*c2*self.b3)*L2+4*self.b1*c2*self.b3)*L3-(-3*self.b1*c2*self.b3*self.mu2-self.b1*c2*self.b3+self.b1*c3*self.b3+3*self.b1*c3*self.b3*self.mu3)*L2**2+(-(2*self.b1*c2*self.b3-3*self.b3**2*c2+6*self.b2*self.b3*c3*self.mu3-2*self.b2*c2*self.b3-6*self.b2*c2*self.b3*self.mu2-c3*self.b3**2-6*self.b1*c3*self.b3+2*self.b2*self.b3*c3+6*self.b1*c2*self.b3*self.mu2+6*self.b1*c3*self.b3*self.mu3+3*c2*self.b3**2*self.mu2-3*c3*self.b3**2*self.mu3)*L1+2*c3*self.b1**2)*L2-2*self.b1*c2*self.b3-(3*self.b2*self.b3*c3*self.mu3+self.b2*c2*self.b3+3*self.b2*c2*self.b3*self.mu2-3*self.b2*self.b3*c3)*L1**2-(-4*self.b1*self.b2*c3-2*self.b3**2*c2)*L1,
#                   (-self.b1*self.b3*self.b2+3*self.b1*self.b2**2*self.mu2-3*self.b1*self.b2*self.b3*self.mu3-3*self.b1*self.b2**2)*L3**2+((6*self.b1*self.b2**2*self.mu2+6*self.b1*self.b2*self.b3*self.mu3-6*self.b2*self.b3**2*self.mu3-6*self.b1*self.b3*self.b2+2*self.b1*self.b2**2-2*self.b2*self.b3**2-3*self.b2**3*self.mu2-5*self.b3*self.b2**2-self.b2**3+3*self.b2**2*self.b3*self.mu3+6*self.b3*self.b2**2*self.mu2)*L1+(-4*self.b1*self.b3*self.b2+6*self.b1*self.b3*self.b2*self.mu2-6*self.b1*self.b2**2*self.mu2-6*self.b1*self.b3**2*self.mu3+self.b1**2*self.b2-2*self.b1*self.b3**2-2*self.b1*self.b2**2+3*self.b1**2*self.b3*self.mu3-3*self.b1**2*self.b3+3*self.b1**2*self.b2*self.mu2+6*self.b1*self.b2*self.b3*self.mu3)*L2-4*self.b1*self.b3*self.b2)*L3+(-3*self.b1*self.b3*self.b2*self.mu2-self.b1*self.b3*self.b2+3*self.b1*self.b3**2*self.mu3+self.b1*self.b3**2)*L2**2+((-self.b2*self.b3**2+6*self.b1*self.b3*self.b2*self.mu2-self.b3**3+2*self.b1*self.b3*self.b2+6*self.b2*self.b3**2*self.mu3-3*self.b3**3*self.mu3-6*self.b1*self.b3**2-2*self.b3*self.b2**2+6*self.b1*self.b3**2*self.mu3+3*self.b3**2*self.b2*self.mu2-6*self.b3*self.b2**2*self.mu2)*L1-2*self.b1**2*self.b3)*L2+2*self.b1*self.b3*self.b2+(3*self.b3*self.b2**2*self.mu2+3*self.b2*self.b3**2*self.mu3-3*self.b2*self.b3**2+self.b3*self.b2**2)*L1**2+(-4*self.b1*self.b3*self.b2-2*self.b2*self.b3**2)*L1,
#                   2*(-3*self.b1*self.b2*self.mu3-3*self.b1*self.b2*self.mu1-2*self.b1*self.b2)*L3**2+(2*(-3*self.b2**2*self.mu1+6*self.b1*self.b2*self.mu3+6*self.b1*self.b2*self.mu1-4*self.b3*self.b2+4*self.b2**2-4*self.b1*self.b2-6*self.b2*self.b3*self.mu3-6*self.b2*self.b3*self.mu1+3*self.b2**2*self.mu3)*L1+2*(3*self.b1**2*self.mu1-6*self.b1*self.b3*self.mu3-2*self.b1**2-6*self.b1*self.b2*self.mu1-6*self.b1*self.b3*self.mu1+6*self.b1*self.b2*self.mu3+3*self.b1**2*self.mu3+8*self.b1*self.b2-4*self.b1*self.b3)*L2+4*self.b2**2)*L3+2*(4*self.b1*self.b3-3*self.b1*self.b3*self.mu1+3*self.b1*self.b3*self.mu3)*L2**2+(2*(-3*self.b3**2*self.mu3-6*self.b2*self.b3*self.mu1-3*self.b3**2*self.mu1-4*self.b1*self.b3+6*self.b2*self.b3*self.mu3+6*self.b1*self.b3*self.mu3-2*self.b3**2+8*self.b3*self.b2+6*self.b1*self.b3*self.mu1)*L1+8*self.b3*self.b2-4*self.b1**2)*L2+2*self.b1*self.b2-2*self.b3*self.b2+2*(3*self.b2*self.b3*self.mu3-2*self.b3*self.b2+3*self.b2*self.b3*self.mu1)*L1**2-8*self.b1*self.b2*L1,
#                   -(-self.b1*self.b2*c3+3*self.b1*c1*self.b2*self.mu1+c1*self.b1*self.b2-3*self.b1*c3*self.b2*self.mu3)*L3**2+(-(-6*self.b1*c1*self.b2*self.mu1-6*self.b2*self.b3*c3*self.mu3-6*self.b1*self.b2*c3+c3*self.b2**2+2*c1*self.b3*self.b2+3*c3*self.b2**2*self.mu3+3*c1*self.b2**2*self.mu1-2*self.b2*self.b3*c3-3*c1*self.b2**2+6*self.b1*c3*self.b2*self.mu3+6*self.b2*self.b3*c1*self.mu1-2*c1*self.b1*self.b2)*L1-(-6*c1*self.b1*self.b2-2*self.b1*c3*self.b3-3*c3*self.b1**2+6*self.b1*c1*self.b3*self.mu1+2*self.b1*self.b2*c3-6*self.b1*c3*self.b3*self.mu3+6*self.b1*c3*self.b2*self.mu3-c1*self.b1**2+2*c1*self.b1*self.b3-3*c1*self.b1**2*self.mu1+6*self.b1*c1*self.b2*self.mu1+3*c3*self.b1**2*self.mu3)*L2+2*c1*self.b2**2)*L3-(-3*c1*self.b1*self.b3+3*self.b1*c1*self.b3*self.mu1+self.b1*c3*self.b3+3*self.b1*c3*self.b3*self.mu3)*L2**2+(-(6*self.b2*self.b3*c3*self.mu3-c3*self.b3**2-6*self.b1*c1*self.b3*self.mu1-6*c1*self.b3*self.b2-3*c3*self.b3**2*self.mu3+6*self.b2*self.b3*c1*self.mu1+3*c1*self.b3**2*self.mu1+6*self.b1*c3*self.b3*self.mu3-6*self.b1*c3*self.b3-2*c1*self.b1*self.b3+2*self.b2*self.b3*c3+c1*self.b3**2)*L1+4*c1*self.b3*self.b2+2*c3*self.b1**2)*L2-2*self.b1*self.b2*c3-(-3*self.b2*self.b3*c3-c1*self.b3*self.b2+3*self.b2*self.b3*c3*self.mu3-3*self.b2*self.b3*c1*self.mu1)*L1**2+4*self.b1*self.b2*c3*L1,
#                   -(-self.b1**2*self.b2+3*self.b1*self.b2*self.b3*self.mu3+self.b1*self.b3*self.b2-3*self.b2*self.b1**2*self.mu1)*L3**2+(-(6*self.b2*self.b1**2*self.mu1-6*self.b1*self.b2*self.b3*self.mu3-self.b3*self.b2**2+2*self.b1**2*self.b2+4*self.b1*self.b3*self.b2+2*self.b2*self.b3**2-6*self.b2*self.b3*self.b1*self.mu1-3*self.b2**2*self.b3*self.mu3-3*self.b2**2*self.b1*self.mu1+6*self.b2*self.b3**2*self.mu3+3*self.b1*self.b2**2)*L1-(-6*self.b1*self.b2*self.b3*self.mu3+self.b1**3+2*self.b1*self.b3**2-6*self.b3*self.b1**2*self.mu1-2*self.b1*self.b3*self.b2+6*self.b1*self.b3**2*self.mu3+self.b1**2*self.b3-3*self.b1**2*self.b3*self.mu3+6*self.b1**2*self.b2+3*self.b1**3*self.mu1-6*self.b2*self.b1**2*self.mu1)*L2-2*self.b1*self.b2**2)*L3-(-3*self.b3*self.b1**2*self.mu1+3*self.b1**2*self.b3-3*self.b1*self.b3**2*self.mu3-self.b1*self.b3**2)*L2**2+(-(6*self.b3*self.b1**2*self.mu1-6*self.b1*self.b3**2*self.mu3+3*self.b3**3*self.mu3+6*self.b1*self.b3*self.b2-2*self.b2*self.b3**2+self.b3**3+2*self.b1**2*self.b3+5*self.b1*self.b3**2-6*self.b2*self.b3**2*self.mu3-3*self.b3**2*self.b1*self.mu1-6*self.b2*self.b3*self.b1*self.mu1)*L1-4*self.b1*self.b3*self.b2-2*self.b1**2*self.b3)*L2+2*self.b1*self.b3*self.b2-(3*self.b2*self.b3**2+3*self.b2*self.b3*self.b1*self.mu1+self.b1*self.b3*self.b2-3*self.b2*self.b3**2*self.mu3)*L1**2-4*self.b1*self.b3*self.b2*L1,
#                   -2*(-4*self.b1*self.b2+3*self.b1*self.b2*self.mu2-3*self.b1*self.b2*self.mu1)*L3**2+(-2*(4*self.b1*self.b2+6*self.b2*self.b3*self.mu2-3*self.b2**2*self.mu1-6*self.b2*self.b3*self.mu1-8*self.b3*self.b2+2*self.b2**2+6*self.b1*self.b2*self.mu1-3*self.b2**2*self.mu2+6*self.b1*self.b2*self.mu2)*L1-2*(2*self.b1**2-8*self.b1*self.b3-6*self.b1*self.b2*self.mu1+4*self.b1*self.b2+3*self.b1**2*self.mu2+3*self.b1**2*self.mu1-6*self.b1*self.b2*self.mu2+6*self.b1*self.b3*self.mu2-6*self.b1*self.b3*self.mu1)*L2-4*self.b2**2+8*self.b1*self.b3)*L3-2*(-3*self.b1*self.b3*self.mu2+2*self.b1*self.b3-3*self.b1*self.b3*self.mu1)*L2**2+(-2*(-6*self.b2*self.b3*self.mu1-6*self.b2*self.b3*self.mu2+6*self.b1*self.b3*self.mu2+4*self.b3*self.b2+3*self.b3**2*self.mu2+4*self.b1*self.b3-3*self.b3**2*self.mu1+6*self.b1*self.b3*self.mu1-4*self.b3**2)*L1-8*self.b3*self.b2)*L2-2*self.b1*self.b3+2*self.b3*self.b2-2*(2*self.b3*self.b2+3*self.b2*self.b3*self.mu2+3*self.b2*self.b3*self.mu1)*L1**2+4*self.b3**2*L1,
#                   (-3*self.b1*c1*self.b2*self.mu1-c1*self.b1*self.b2-3*self.b1*c2*self.b2*self.mu2+3*self.b1*c2*self.b2)*L3**2+((-6*self.b2*c2*self.b3*self.mu2-6*self.b2*self.b3*c1*self.mu1+6*self.b2*c2*self.b3+6*self.b1*c1*self.b2*self.mu1-3*c1*self.b2**2*self.mu1-2*self.b1*c2*self.b2-2*c1*self.b3*self.b2-6*self.b1*c2*self.b2*self.mu2+c2*self.b2**2+3*c2*self.b2**2*self.mu2+2*c1*self.b1*self.b2+3*c1*self.b2**2)*L1+(2*self.b1*c2*self.b2-6*self.b1*c1*self.b3*self.mu1-6*self.b1*c2*self.b3*self.mu2-2*c1*self.b1*self.b3+3*c1*self.b1**2*self.mu1+6*c1*self.b1*self.b2-c2*self.b1**2+6*self.b1*c2*self.b3+6*self.b1*c2*self.b2*self.mu2+c1*self.b1**2-3*c2*self.b1**2*self.mu2-6*self.b1*c1*self.b2*self.mu1)*L2+4*self.b1*c2*self.b3+2*c1*self.b2**2)*L3+(-3*self.b1*c1*self.b3*self.mu1+3*self.b1*c2*self.b3*self.mu2+3*c1*self.b1*self.b3+self.b1*c2*self.b3)*L2**2+((-6*self.b2*self.b3*c1*self.mu1-2*self.b1*c2*self.b3-6*self.b1*c2*self.b3*self.mu2-3*c2*self.b3**2*self.mu2-3*c1*self.b3**2*self.mu1+2*c1*self.b1*self.b3+2*self.b2*c2*self.b3-c1*self.b3**2+6*c1*self.b3*self.b2+6*self.b1*c1*self.b3*self.mu1+6*self.b2*c2*self.b3*self.mu2+3*self.b3**2*c2)*L1+4*c1*self.b3*self.b2)*L2-2*c1*self.b3*self.b2+(c1*self.b3*self.b2-3*self.b2*c2*self.b3*self.mu2-self.b2*c2*self.b3+3*self.b2*self.b3*c1*self.mu1)*L1**2+2*self.b3**2*c2*L1,
#                   -(-self.b1**2*self.b2-3*self.b2*self.b1**2*self.mu1-3*self.b1*self.b2**2*self.mu2+3*self.b1*self.b2**2)*L3**2+(-(2*self.b1**2*self.b2+6*self.b2*self.b1**2*self.mu1+self.b1*self.b2**2-3*self.b2**2*self.b1*self.mu1+3*self.b2**3*self.mu2-6*self.b3*self.b2**2*self.mu2-2*self.b1*self.b3*self.b2-6*self.b2*self.b3*self.b1*self.mu1+6*self.b3*self.b2**2-6*self.b1*self.b2**2*self.mu2+self.b2**3)*L1-(-6*self.b3*self.b1**2*self.mu1+6*self.b1*self.b2**2*self.mu2-3*self.b1**2*self.b2*self.mu2-6*self.b2*self.b1**2*self.mu1+self.b1**3-2*self.b1**2*self.b3+5*self.b1**2*self.b2-6*self.b1*self.b3*self.b2*self.mu2+3*self.b1**3*self.mu1+2*self.b1*self.b2**2+6*self.b1*self.b3*self.b2)*L2-2*self.b1*self.b2**2-4*self.b1*self.b3*self.b2)*L3-(3*self.b1**2*self.b3+self.b1*self.b3*self.b2+3*self.b1*self.b3*self.b2*self.mu2-3*self.b3*self.b1**2*self.mu1)*L2**2+(-(-6*self.b1*self.b3*self.b2*self.mu2+2*self.b3*self.b2**2+2*self.b1**2*self.b3+3*self.b2*self.b3**2+6*self.b3*self.b2**2*self.mu2-self.b1*self.b3**2-3*self.b3**2*self.b1*self.mu1+6*self.b3*self.b1**2*self.mu1-3*self.b3**2*self.b2*self.mu2+4*self.b1*self.b3*self.b2-6*self.b2*self.b3*self.b1*self.mu1)*L1-4*self.b1*self.b3*self.b2)*L2+2*self.b1*self.b3*self.b2-(-self.b3*self.b2**2-3*self.b3*self.b2**2*self.mu2+3*self.b2*self.b3*self.b1*self.mu1+self.b1*self.b3*self.b2)*L1**2-2*self.b2*self.b3**2*L1],
#                  [2*(3*c1*c2*self.mu3+3*c1*c2*self.mu2-2*c2*c1)*L3**2+(2*(-3*c2**2*self.mu2-2*c2**2+8*c2*c1+6*c1*c2*self.mu2+6*c2*c3*self.mu3+6*c2*c3*self.mu2-6*c1*c2*self.mu3-3*c2**2*self.mu3-4*c3*c2)*L1+2*(-4*c3*c1-3*c1**2*self.mu3-6*c1*c2*self.mu2+3*c1**2*self.mu2+4*c1**2-4*c2*c1+6*c1*c3*self.mu3-6*c1*c2*self.mu3+6*c1*c3*self.mu2)*L2-8*c3*c1)*L3+2*(-2*c3*c1-3*c1*c3*self.mu2-3*c1*c3*self.mu3)*L2**2+(2*(-6*c1*c3*self.mu3-6*c2*c3*self.mu2+3*c3**2*self.mu2+6*c1*c3*self.mu2-6*c2*c3*self.mu3+3*c3**2*self.mu3-4*c3*c2+8*c3*c1-2*c3**2)*L1+4*c1**2)*L2+2*c3*c1-2*c2*c1+2*(-3*c2*c3*self.mu3+4*c3*c2+3*c2*c3*self.mu2)*L1**2+2*(4*c2*c1-2*c3**2)*L1,
#                   (3*c1*c2**2+c1*c3*c2-3*c1*c2**2*self.mu2+3*c1*c2*c3*self.mu3)*L3**2+((c2**3-6*c1*c2*c3*self.mu3-6*c1*c2**2*self.mu2-3*c2**2*c3*self.mu3+3*c2**3*self.mu2+6*c2*c3**2*self.mu3+5*c2**2*c3+6*c1*c3*c2-6*c2**2*c3*self.mu2-2*c1*c2**2+2*c2*c3**2)*L1+(-6*c1*c2*c3*self.mu3+6*c1*c3**2*self.mu3-3*c1**2*c3*self.mu3-3*c1**2*c2*self.mu2+6*c1*c2**2*self.mu2+4*c1*c3*c2-c1**2*c2+2*c1*c2**2+2*c1*c3**2-6*c1*c2*c3*self.mu2+3*c1**2*c3)*L2+4*c1*c3*c2)*L3+(-c1*c3**2+3*c1*c2*c3*self.mu2+c1*c3*c2-3*c1*c3**2*self.mu3)*L2**2+((-6*c2*c3**2*self.mu3+6*c2**2*c3*self.mu2+c3**3-6*c1*c3**2*self.mu3-3*c2*c3**2*self.mu2-2*c1*c3*c2+6*c1*c3**2+2*c2**2*c3+3*c3**3*self.mu3+c2*c3**2-6*c1*c2*c3*self.mu2)*L1+2*c1**2*c3)*L2-2*c1*c3*c2+(-c2**2*c3-3*c2*c3**2*self.mu3-3*c2**2*c3*self.mu2+3*c2*c3**2)*L1**2+(4*c1*c3*c2+2*c2*c3**2)*L1,
#                   (3*c1*c2*self.b2*self.mu2-3*c1*c2*self.b3*self.mu3-c1*c2*self.b3-3*c1*c2*self.b2)*L3**2+((6*c1*c2*self.b2*self.mu2+3*c2**2*self.b3*self.mu3+6*c1*c2*self.b3*self.mu3-c2**2*self.b2+c2**2*self.b3+6*c3*self.b2*c2*self.mu2-6*c1*c2*self.b3+2*c1*c2*self.b2-6*c2*self.b3*c3*self.mu3-2*c3*self.b3*c2-6*c3*self.b2*c2-3*c2**2*self.b2*self.mu2)*L1+(-3*c1**2*self.b3-6*c1*self.b2*c3+2*c1*c2*self.b3+3*c1**2*self.b2*self.mu2-6*c1*c2*self.b2*self.mu2-6*c1*c3*self.b3*self.mu3-2*c1*c3*self.b3+6*c1*c3*self.b2*self.mu2-2*c1*c2*self.b2+c1**2*self.b2+6*c1*c2*self.b3*self.mu3+3*c1**2*self.b3*self.mu3)*L2-4*c1*self.b2*c3)*L3+(-3*c1*c3*self.b2*self.mu2+c1*c3*self.b3-c1*self.b2*c3+3*c1*c3*self.b3*self.mu3)*L2**2+((-6*c3*self.b2*c2*self.mu2+6*c1*c3*self.b3*self.mu3-c3**2*self.b3-2*c3*self.b2*c2+2*c1*self.b2*c3+2*c3*self.b3*c2+6*c1*c3*self.b2*self.mu2-3*c3**2*self.b3*self.mu3-6*c1*c3*self.b3+6*c2*self.b3*c3*self.mu3+3*c3**2*self.b2*self.mu2-3*c3**2*self.b2)*L1-2*c1**2*self.b3)*L2+2*c1*self.b2*c3+(-3*c3*self.b3*c2+3*c3*self.b2*c2*self.mu2+3*c2*self.b3*c3*self.mu3+c3*self.b2*c2)*L1**2+(-4*c1*c2*self.b3-2*c3**2*self.b2)*L1,
#                   -2*(3*c1*c2*self.mu3+2*c2*c1+3*c1*c2*self.mu1)*L3**2+(-2*(6*c2*c3*self.mu3-4*c2**2-6*c1*c2*self.mu3-3*c2**2*self.mu3+4*c3*c2+3*c2**2*self.mu1-6*c1*c2*self.mu1+6*c2*c3*self.mu1+4*c2*c1)*L1-2*(6*c1*c2*self.mu1-6*c1*c2*self.mu3+6*c1*c3*self.mu3+6*c1*c3*self.mu1-3*c1**2*self.mu1+4*c3*c1-8*c2*c1-3*c1**2*self.mu3+2*c1**2)*L2+4*c2**2)*L3-2*(3*c1*c3*self.mu1-4*c3*c1-3*c1*c3*self.mu3)*L2**2+(-2*(3*c3**2*self.mu1+3*c3**2*self.mu3-6*c1*c3*self.mu1+6*c2*c3*self.mu1-6*c1*c3*self.mu3+4*c3*c1-8*c3*c2-6*c2*c3*self.mu3+2*c3**2)*L1-4*c1**2+8*c3*c2)*L2+2*c2*c1-2*c3*c2-2*(-3*c2*c3*self.mu1-3*c2*c3*self.mu3+2*c3*c2)*L1**2-8*c2*c1*L1,
#                   (c1*c3*c2-c1**2*c2-3*c1**2*c2*self.mu1+3*c1*c2*c3*self.mu3)*L3**2+((6*c1**2*c2*self.mu1-c2**2*c3-6*c2*c3*c1*self.mu1-3*c2**2*c3*self.mu3-3*c1*c2**2*self.mu1+6*c2*c3**2*self.mu3+4*c1*c3*c2+2*c2*c3**2-6*c1*c2*c3*self.mu3+3*c1*c2**2+2*c1**2*c2)*L1+(-2*c1*c3*c2+6*c1**2*c2+3*c1**3*self.mu1-6*c1*c2*c3*self.mu3+2*c1*c3**2+c1**2*c3+6*c1*c3**2*self.mu3+c1**3-6*c1**2*c3*self.mu1-6*c1**2*c2*self.mu1-3*c1**2*c3*self.mu3)*L2+2*c1*c2**2)*L3+(3*c1**2*c3-3*c1*c3**2*self.mu3-3*c1**2*c3*self.mu1-c1*c3**2)*L2**2+((5*c1*c3**2+3*c3**3*self.mu3-6*c2*c3**2*self.mu3-6*c1*c3**2*self.mu3+6*c1**2*c3*self.mu1+c3**3-6*c2*c3*c1*self.mu1+6*c1*c3*c2-2*c2*c3**2-3*c1*c3**2*self.mu1+2*c1**2*c3)*L1+4*c1*c3*c2+2*c1**2*c3)*L2-2*c1*c3*c2+(3*c2*c3**2+c1*c3*c2+3*c2*c3*c1*self.mu1-3*c2*c3**2*self.mu3)*L1**2+4*c1*c3*c2*L1,
#                   (3*c2*self.b1*c1*self.mu1-c1*c2*self.b3+c1*c2*self.b1-3*c1*c2*self.b3*self.mu3)*L3**2+((-2*c1*c2*self.b1+3*c2**2*self.b3*self.mu3+6*c1*c2*self.b3*self.mu3-6*c2*self.b1*c1*self.mu1+3*c2**2*self.b1*self.mu1+2*c3*c2*self.b1+6*c2*c3*self.b1*self.mu1-3*c2**2*self.b1+c2**2*self.b3-6*c2*self.b3*c3*self.mu3-6*c1*c2*self.b3-2*c3*self.b3*c2)*L1+(-c1**2*self.b1+2*c1*c3*self.b1-3*c1**2*self.b1*self.mu1+6*c2*self.b1*c1*self.mu1-6*c1*c2*self.b1-3*c1**2*self.b3-6*c1*c3*self.b3*self.mu3+6*c1*c2*self.b3*self.mu3+2*c1*c2*self.b3-2*c1*c3*self.b3+6*c3*self.b1*c1*self.mu1+3*c1**2*self.b3*self.mu3)*L2-2*c2**2*self.b1)*L3+(3*c3*self.b1*c1*self.mu1+c1*c3*self.b3+3*c1*c3*self.b3*self.mu3-3*c1*c3*self.b1)*L2**2+((3*c3**2*self.b1*self.mu1+c3**2*self.b1-6*c3*self.b1*c1*self.mu1+6*c1*c3*self.b3*self.mu3-2*c1*c3*self.b1+6*c2*self.b3*c3*self.mu3-6*c3*c2*self.b1+2*c3*self.b3*c2-c3**2*self.b3-3*c3**2*self.b3*self.mu3-6*c1*c3*self.b3+6*c2*c3*self.b1*self.mu1)*L1-2*c1**2*self.b3-4*c3*c2*self.b1)*L2+2*c1*c2*self.b3+(3*c2*self.b3*c3*self.mu3-c3*c2*self.b1-3*c3*self.b3*c2-3*c2*c3*self.b1*self.mu1)*L1**2-4*c1*c2*self.b3*L1,
#                   2*(-3*c1*c2*self.mu2+4*c2*c1+3*c1*c2*self.mu1)*L3**2+(2*(6*c2*c3*self.mu1-2*c2**2+3*c2**2*self.mu2-6*c2*c3*self.mu2+3*c2**2*self.mu1-6*c1*c2*self.mu1-4*c2*c1-6*c1*c2*self.mu2+8*c3*c2)*L1+2*(-2*c1**2+8*c3*c1-6*c1*c3*self.mu2-3*c1**2*self.mu1-4*c2*c1-3*c1**2*self.mu2+6*c1*c2*self.mu1+6*c1*c2*self.mu2+6*c1*c3*self.mu1)*L2-4*c2**2+8*c3*c1)*L3+2*(-2*c3*c1+3*c1*c3*self.mu2+3*c1*c3*self.mu1)*L2**2+(2*(-3*c3**2*self.mu2+3*c3**2*self.mu1-6*c1*c3*self.mu1+6*c2*c3*self.mu1-6*c1*c3*self.mu2-4*c3*c2+4*c3**2-4*c3*c1+6*c2*c3*self.mu2)*L1-8*c3*c2)*L2-2*c3*c1+2*c3*c2+2*(-3*c2*c3*self.mu2-2*c3*c2-3*c2*c3*self.mu1)*L1**2+4*c3**2*L1,
#                   -(3*c1*c2**2*self.mu2-3*c1*c2**2+3*c1**2*c2*self.mu1+c1**2*c2)*L3**2+(-(-2*c1**2*c2-6*c2**2*c3+6*c2**2*c3*self.mu2-c1*c2**2+2*c1*c3*c2-6*c1**2*c2*self.mu1-3*c2**3*self.mu2+3*c1*c2**2*self.mu1+6*c2*c3*c1*self.mu1+6*c1*c2**2*self.mu2-c2**3)*L1-(-5*c1**2*c2-3*c1**3*self.mu1+2*c1**2*c3-c1**3+6*c1*c2*c3*self.mu2-2*c1*c2**2+6*c1**2*c3*self.mu1-6*c1*c3*c2+6*c1**2*c2*self.mu1-6*c1*c2**2*self.mu2+3*c1**2*c2*self.mu2)*L2+2*c1*c2**2+4*c1*c3*c2)*L3-(-3*c1**2*c3-c1*c3*c2-3*c1*c2*c3*self.mu2+3*c1**2*c3*self.mu1)*L2**2+(-(6*c2*c3*c1*self.mu1-6*c1**2*c3*self.mu1+c1*c3**2-2*c1**2*c3-3*c2*c3**2-2*c2**2*c3+6*c1*c2*c3*self.mu2+3*c2*c3**2*self.mu2+3*c1*c3**2*self.mu1-4*c1*c3*c2-6*c2**2*c3*self.mu2)*L1+4*c1*c3*c2)*L2-2*c1*c3*c2-(-3*c2*c3*c1*self.mu1+3*c2**2*c3*self.mu2+c2**2*c3-c1*c3*c2)*L1**2+2*c2*c3**2*L1,
#                   (3*c1*c2*self.b2*self.mu2-3*c1*c2*self.b2+c1*c2*self.b1+3*c2*self.b1*c1*self.mu1)*L3**2+((-2*c1*c2*self.b1+6*c3*self.b2*c2*self.mu2-3*c2**2*self.b2*self.mu2-6*c3*self.b2*c2-6*c2*self.b1*c1*self.mu1-c2**2*self.b2+2*c1*c2*self.b2+3*c2**2*self.b1*self.mu1+2*c3*c2*self.b1+6*c2*c3*self.b1*self.mu1-3*c2**2*self.b1+6*c1*c2*self.b2*self.mu2)*L1+(6*c3*self.b1*c1*self.mu1-c1**2*self.b1+2*c1*c3*self.b1-6*c1*c2*self.b1-3*c1**2*self.b1*self.mu1-6*c1*c2*self.b2*self.mu2-6*c1*self.b2*c3-2*c1*c2*self.b2+3*c1**2*self.b2*self.mu2+6*c1*c3*self.b2*self.mu2+6*c2*self.b1*c1*self.mu1+c1**2*self.b2)*L2-2*c2**2*self.b1-4*c1*self.b2*c3)*L3+(-3*c1*c3*self.b2*self.mu2+3*c3*self.b1*c1*self.mu1-c1*self.b2*c3-3*c1*c3*self.b1)*L2**2+((3*c3**2*self.b2*self.mu2+6*c1*c3*self.b2*self.mu2+6*c2*c3*self.b1*self.mu1-2*c1*c3*self.b1+c3**2*self.b1-2*c3*self.b2*c2+2*c1*self.b2*c3-3*c3**2*self.b2-6*c3*self.b1*c1*self.mu1+3*c3**2*self.b1*self.mu1-6*c3*c2*self.b1-6*c3*self.b2*c2*self.mu2)*L1-4*c3*c2*self.b1)*L2+2*c3*c2*self.b1+(c3*self.b2*c2+3*c3*self.b2*c2*self.mu2-c3*c2*self.b1-3*c2*c3*self.b1*self.mu1)*L1**2-2*c3**2*self.b2*L1],
#                  [-(2*self.b2*c1-3*c2*self.b1*self.mu3+2*c2*self.b1-3*c1*self.b2*self.mu2-3*c2*self.b1*self.mu2-3*c1*self.b2*self.mu3)*L3**2+(-(-6*c1*self.b2*self.mu2+4*c2*self.b2-8*self.b2*c1-6*c2*self.b1*self.mu2+6*c2*self.b2*self.mu3+6*c2*self.b2*self.mu2+6*c2*self.b1*self.mu3+6*c1*self.b2*self.mu3-8*c2*self.b1-6*c3*self.b2*self.mu2+4*c2*self.b3+4*self.b2*c3-6*c2*self.b3*self.mu2-6*c2*self.b3*self.mu3-6*c3*self.b2*self.mu3)*L1-(6*c1*self.b1*self.mu3+6*c2*self.b1*self.mu3+4*self.b3*c1+6*c2*self.b1*self.mu2+6*c1*self.b2*self.mu3-6*c1*self.b3*self.mu2-6*c3*self.b1*self.mu2-6*c1*self.b1*self.mu2+4*c3*self.b1-8*c1*self.b1+6*c1*self.b2*self.mu2-6*c1*self.b3*self.mu3+4*self.b2*c1+4*c2*self.b1-6*c3*self.b1*self.mu3)*L2-4*c3*self.b1-4*self.b3*c1)*L3-(3*c3*self.b1*self.mu3+2*self.b3*c1+3*c1*self.b3*self.mu3+3*c3*self.b1*self.mu2+2*c3*self.b1+3*c1*self.b3*self.mu2)*L2**2+(-(-8*c3*self.b1-8*self.b3*c1+6*c1*self.b3*self.mu3-6*c1*self.b3*self.mu2-6*c3*self.b3*self.mu3+6*c2*self.b3*self.mu3+4*c3*self.b3+4*c2*self.b3+6*c3*self.b1*self.mu3+4*self.b2*c3+6*c3*self.b2*self.mu3-6*c3*self.b1*self.mu2+6*c3*self.b2*self.mu2+6*c2*self.b3*self.mu2-6*c3*self.b3*self.mu2)*L1+4*c1*self.b1)*L2-self.b2*c1+c3*self.b1-c2*self.b1+self.b3*c1-(-4*self.b2*c3-3*c3*self.b2*self.mu2+3*c2*self.b3*self.mu3-4*c2*self.b3+3*c3*self.b2*self.mu3-3*c2*self.b3*self.mu2)*L1**2-(4*c3*self.b3-4*self.b2*c1-4*c2*self.b1)*L1,
#                   -0.5*(-3*c2**2*self.b1+3*c2**2*self.b1*self.mu2-3*c2*self.b1*c3*self.mu3-3*c1*c3*self.b2*self.mu3-c1*self.b2*c3-c3*c2*self.b1-3*c1*c2*self.b2+3*c1*c2*self.b2*self.mu2)*L3**2+(-0.5*(6*c2*c3*self.b2*self.mu3-2*c3*self.b3*c2-6*c3*c2*self.b1-2*c2**2*self.b2+2*c2**2*self.b1-6*c2**2*self.b2*self.mu2+6*c1*c3*self.b2*self.mu3+2*c1*c2*self.b2-2*c3**2*self.b2-6*c2**2*self.b3+6*c2*self.b1*c3*self.mu3+6*c2**2*self.b1*self.mu2-4*c3*self.b2*c2-6*c3**2*self.b2*self.mu3+6*c1*c2*self.b2*self.mu2-6*c1*self.b2*c3+6*c2**2*self.b3*self.mu2-6*c2*self.b3*c3*self.mu3+6*c3*self.b2*c2*self.mu2)*L1-0.5*(6*c1*c2*self.b1*self.mu2-4*c3*c2*self.b1-2*c1*c3*self.b3+6*c3*self.b1*c2*self.mu2-6*c1*c3*self.b1-2*c2**2*self.b1+6*c2*self.b1*c3*self.mu3+6*c1*c3*self.b1*self.mu3-2*c3**2*self.b1-6*c1*c2*self.b3-6*c1*c3*self.b3*self.mu3+6*c1*c3*self.b2*self.mu3-2*c1*c2*self.b2+6*c1*c2*self.b3*self.mu2-6*c2**2*self.b1*self.mu2-6*c1*c2*self.b2*self.mu2+2*c1*c2*self.b1-6*c3**2*self.b1*self.mu3+2*c1*self.b2*c3)*L2+2*c1*c2*self.b3+2*c3*c2*self.b1)*L3-0.5*(3*c1*c3*self.b3*self.mu3+c3**2*self.b1-3*c3*self.b1*c2*self.mu2-c3*c2*self.b1-3*c1*c2*self.b3*self.mu2+c1*c3*self.b3-c1*c2*self.b3+3*c3**2*self.b1*self.mu3)*L2**2+(-0.5*(-6*c3**2*self.b1-6*c2**2*self.b3*self.mu2+6*c1*c3*self.b3*self.mu3-6*c1*c3*self.b3+2*c1*c2*self.b3-2*c3**2*self.b3-6*c3**2*self.b3*self.mu3-4*c3*self.b3*c2-6*c3*self.b2*c2*self.mu2+6*c3*self.b1*c2*self.mu2+6*c3**2*self.b1*self.mu3+6*c3**2*self.b2*self.mu3-2*c2**2*self.b3+6*c1*c2*self.b3*self.mu2+2*c3**2*self.b2+2*c3*c2*self.b1-2*c3*self.b2*c2+6*c2*self.b3*c3*self.mu3+6*c3*c2*self.b3*self.mu2)*L1+2*c1*c3*self.b1)*L2-c3*c2*self.b1-c1*c2*self.b3-0.5*(3*c3**2*self.b2*self.mu3+3*c3*self.b2*c2*self.mu2+c3*self.b2*c2+3*c2*self.b3*c3*self.mu3+3*c2**2*self.b3*self.mu2-3*c3*self.b3*c2-3*c3**2*self.b2+c2**2*self.b3)*L1**2-0.5*(-4*c3*self.b3*c2-4*c3*c2*self.b1-4*c1*self.b2*c3)*L1,
#                   0.5*(-3*self.b1*c2*self.b2-3*c1*self.b2**2-c1*self.b2*self.b3+3*c1*self.b2**2*self.mu2+3*self.b1*c2*self.b2*self.mu2-self.b1*c2*self.b3-3*c1*self.b2*self.b3*self.mu3-3*c2*self.b1*self.b3*self.mu3)*L3**2+(0.5*(2*c1*self.b2**2+2*self.b1*c2*self.b2-4*self.b2*c2*self.b3-6*c1*self.b2*self.b3-2*c2*self.b2**2+6*c1*self.b2**2*self.mu2-6*c3*self.b2**2-6*c2*self.b2**2*self.mu2+6*self.b2*c2*self.b3*self.mu2+6*c3*self.b2**2*self.mu2-6*self.b2*self.b3*c3*self.mu3-2*self.b2*self.b3*c3-6*c2*self.b3**2*self.mu3+6*c2*self.b2*self.b3*self.mu3+6*self.b1*c2*self.b2*self.mu2-2*self.b3**2*c2+6*c2*self.b1*self.b3*self.mu3+6*c1*self.b2*self.b3*self.mu3-6*self.b1*c2*self.b3)*L1+0.5*(-4*c1*self.b2*self.b3-2*c1*self.b2**2-6*c1*self.b1*self.b3+6*c1*self.b1*self.b2*self.mu2+6*c3*self.b1*self.b2*self.mu2+2*c1*self.b1*self.b2+6*c1*self.b2*self.b3*self.mu2-2*c1*self.b3**2-6*self.b1*self.b2*c3-6*self.b1*c3*self.b3*self.mu3-6*c1*self.b3**2*self.mu3-2*self.b1*c3*self.b3+6*c2*self.b1*self.b3*self.mu3-6*c1*self.b2**2*self.mu2+2*self.b1*c2*self.b3+6*c1*self.b2*self.b3*self.mu3-6*self.b1*c2*self.b2*self.mu2-2*self.b1*c2*self.b2+6*c1*self.b1*self.b3*self.mu3)*L2-2*self.b1*self.b2*c3-2*c1*self.b2*self.b3)*L3+0.5*(-3*c3*self.b1*self.b2*self.mu2+self.b1*c3*self.b3+c1*self.b3**2-self.b1*self.b2*c3+3*self.b1*c3*self.b3*self.mu3-3*c1*self.b2*self.b3*self.mu2+3*c1*self.b3**2*self.mu3-c1*self.b2*self.b3)*L2**2+(0.5*(-2*c3*self.b3**2-6*self.b1*c3*self.b3-6*self.b2*c2*self.b3*self.mu2+6*self.b1*c3*self.b3*self.mu3-6*c3*self.b3**2*self.mu3-4*self.b2*self.b3*c3+6*c3*self.b2*self.b3*self.mu2+6*c1*self.b2*self.b3*self.mu2+6*c1*self.b3**2*self.mu3-6*c3*self.b2**2*self.mu2+2*self.b1*self.b2*c3-6*c1*self.b3**2+6*c3*self.b1*self.b2*self.mu2+6*c2*self.b3**2*self.mu3+6*self.b2*self.b3*c3*self.mu3-2*c3*self.b2**2+2*self.b3**2*c2+2*c1*self.b2*self.b3-2*self.b2*c2*self.b3)*L1-2*c1*self.b1*self.b3)*L2+self.b1*self.b2*c3+c1*self.b2*self.b3+0.5*(3*self.b2*c2*self.b3*self.mu2+3*self.b2*self.b3*c3*self.mu3+3*c3*self.b2**2*self.mu2-3*self.b2*self.b3*c3-3*self.b3**2*c2+c3*self.b2**2+3*c2*self.b3**2*self.mu3+self.b2*c2*self.b3)*L1**2+0.5*(-4*self.b1*c2*self.b3-4*c1*self.b2*self.b3-4*self.b2*self.b3*c3)*L1,
#                   (-2*c2*self.b1-3*c1*self.b2*self.mu3-3*c2*self.b1*self.mu3-3*c2*self.b1*self.mu1-3*c1*self.b2*self.mu1-2*self.b2*c1)*L3**2+((6*c2*self.b2*self.mu3-6*c3*self.b2*self.mu1-4*c2*self.b3-4*c2*self.b1-6*c2*self.b2*self.mu1-4*self.b2*c1+6*c1*self.b2*self.mu3+6*c2*self.b1*self.mu3-6*c2*self.b3*self.mu1+8*c2*self.b2+6*c1*self.b2*self.mu1-6*c2*self.b3*self.mu3-4*self.b2*c3-6*c3*self.b2*self.mu3+6*c2*self.b1*self.mu1)*L1+(-4*c3*self.b1+6*c2*self.b1*self.mu3+6*c1*self.b2*self.mu3-4*self.b3*c1-6*c3*self.b1*self.mu3-6*c2*self.b1*self.mu1-6*c3*self.b1*self.mu1-6*c1*self.b3*self.mu3-6*c1*self.b2*self.mu1-6*c1*self.b3*self.mu1+6*c1*self.b1*self.mu3+8*self.b2*c1+8*c2*self.b1+6*c1*self.b1*self.mu1-4*c1*self.b1)*L2+4*c2*self.b2)*L3+(-3*c1*self.b3*self.mu1-3*c3*self.b1*self.mu1+4*self.b3*c1+3*c1*self.b3*self.mu3+3*c3*self.b1*self.mu3+4*c3*self.b1)*L2**2+((-6*c3*self.b3*self.mu3-4*self.b3*c1+8*self.b2*c3-4*c3*self.b1+6*c3*self.b2*self.mu3+6*c3*self.b1*self.mu3-6*c2*self.b3*self.mu1+6*c1*self.b3*self.mu3-4*c3*self.b3+6*c3*self.b1*self.mu1-6*c3*self.b2*self.mu1+6*c2*self.b3*self.mu3-6*c3*self.b3*self.mu1+8*c2*self.b3+6*c1*self.b3*self.mu1)*L1-4*c1*self.b1+4*self.b2*c3+4*c2*self.b3)*L2-c2*self.b3+c2*self.b1-self.b2*c3+self.b2*c1+(3*c2*self.b3*self.mu1+3*c3*self.b2*self.mu3-2*self.b2*c3+3*c2*self.b3*self.mu3-2*c2*self.b3+3*c3*self.b2*self.mu1)*L1**2+(-4*self.b2*c1-4*c2*self.b1)*L1,
#                   0.5*(-c1**2*self.b2+c3*c2*self.b1+3*c2*self.b1*c3*self.mu3-c1*c2*self.b1+c1*self.b2*c3+3*c1*c3*self.b2*self.mu3-3*self.b2*c1**2*self.mu1-3*c2*self.b1*c1*self.mu1)*L3**2+(0.5*(-6*c1*c3*self.b2*self.mu3+6*self.b2*c1**2*self.mu1+2*c1*c2*self.b1+6*c1*c2*self.b2+2*c3*self.b3*c2-2*c3*self.b2*c2-2*c1*c2*self.b3+4*c1*self.b2*c3+6*c3*c2*self.b1-6*c2*self.b2*c1*self.mu1+2*c1**2*self.b2-6*c3*self.b2*c1*self.mu1+2*c3**2*self.b2-6*c2*self.b1*c3*self.mu3-6*c2*self.b3*c1*self.mu1+6*c2*self.b3*c3*self.mu3-6*c2*c3*self.b2*self.mu3+6*c3**2*self.b2*self.mu3+6*c2*self.b1*c1*self.mu1)*L1+0.5*(-2*c3*c2*self.b1+6*c1*c2*self.b1+6*c1*c3*self.b3*self.mu3+2*c1*c3*self.b3+4*c1*c3*self.b1-2*c1*self.b2*c3+2*c3**2*self.b1-6*c3*self.b1*c1*self.mu1-6*self.b2*c1**2*self.mu1-6*c2*self.b1*c1*self.mu1-6*self.b3*c1**2*self.mu1-6*c1*c3*self.b2*self.mu3+6*c3**2*self.b1*self.mu3+6*c1**2*self.b2-6*c1*c3*self.b1*self.mu3+2*self.b1*c1**2+6*self.b1*c1**2*self.mu1-6*c2*self.b1*c3*self.mu3-2*c1**2*self.b3)*L2+2*c1*c2*self.b2)*L3+0.5*(3*c1**2*self.b3-c1*c3*self.b3+3*c1*c3*self.b1-3*c3*self.b1*c1*self.mu1-3*c3**2*self.b1*self.mu3-3*c1*c3*self.b3*self.mu3-3*self.b3*c1**2*self.mu1-c3**2*self.b1)*L2**2+(0.5*(6*c1*c2*self.b3+6*c1*self.b2*c3-6*c3*self.b3*c1*self.mu1+6*c3**2*self.b1-2*c3*self.b3*c2+6*c3*self.b1*c1*self.mu1+2*c1*c3*self.b1-6*c3**2*self.b2*self.mu3+2*c1**2*self.b3-6*c1*c3*self.b3*self.mu3-6*c3*self.b2*c1*self.mu1-2*c3**2*self.b2+2*c3**2*self.b3+6*self.b3*c1**2*self.mu1+6*c3**2*self.b3*self.mu3+4*c1*c3*self.b3-6*c2*self.b3*c1*self.mu1-6*c2*self.b3*c3*self.mu3-6*c3**2*self.b1*self.mu3)*L1+2*c1*self.b2*c3+2*c1*c3*self.b1+2*c1*c2*self.b3)*L2-c1*self.b2*c3-c3*c2*self.b1+0.5*(c1*c2*self.b3+c1*self.b2*c3+3*c3*self.b2*c1*self.mu1+3*c2*self.b3*c1*self.mu1-3*c3**2*self.b2*self.mu3+3*c3**2*self.b2+3*c3*self.b3*c2-3*c2*self.b3*c3*self.mu3)*L1**2+0.5*(4*c1*self.b2*c3+4*c3*c2*self.b1)*L1,
#                   -0.5*(c1*self.b2*self.b3-c2*self.b1**2+self.b1*c2*self.b3+3*c1*self.b2*self.b3*self.mu3-3*self.b1*c1*self.b2*self.mu1-3*c2*self.b1**2*self.mu1+3*c2*self.b1*self.b3*self.mu3-c1*self.b1*self.b2)*L3**2+(-0.5*(2*c1*self.b1*self.b2-6*c1*self.b2*self.b3*self.mu3+6*c1*self.b2*self.b3-2*self.b2*c2*self.b3+2*self.b2*self.b3*c3+4*self.b1*c2*self.b3+6*self.b1*c2*self.b2-6*c2*self.b1*self.b3*self.mu3-6*c2*self.b2*self.b3*self.mu3+6*c2*self.b3**2*self.mu3-2*self.b1*self.b2*c3+6*self.b1*c1*self.b2*self.mu1+2*self.b3**2*c2-6*c3*self.b2*self.b1*self.mu1-6*c2*self.b3*self.b1*self.mu1+2*c2*self.b1**2+6*self.b2*self.b3*c3*self.mu3-6*c2*self.b2*self.b1*self.mu1+6*c2*self.b1**2*self.mu1)*L1-0.5*(-2*self.b1*c2*self.b3+6*c2*self.b1**2-6*c1*self.b2*self.b3*self.mu3+6*c1*self.b1**2*self.mu1+6*self.b1*c3*self.b3*self.mu3-6*c2*self.b1**2*self.mu1-6*c2*self.b1*self.b3*self.mu3+2*self.b1*c3*self.b3+6*c1*self.b1*self.b2-2*c1*self.b2*self.b3-6*c3*self.b1**2*self.mu1+6*c1*self.b3**2*self.mu3+2*c1*self.b1**2+4*c1*self.b1*self.b3-6*c1*self.b1*self.b3*self.mu3+2*c1*self.b3**2-6*self.b1*c1*self.b2*self.mu1-6*self.b1*c1*self.b3*self.mu1-2*c3*self.b1**2)*L2-2*self.b1*c2*self.b2)*L3-0.5*(-3*c3*self.b1**2*self.mu1-3*c1*self.b3**2*self.mu3-c1*self.b3**2+3*c1*self.b1*self.b3-3*self.b1*c3*self.b3*self.mu3-3*self.b1*c1*self.b3*self.mu1+3*c3*self.b1**2-self.b1*c3*self.b3)*L2**2+(-0.5*(2*c3*self.b3**2+6*c3*self.b3**2*self.mu3-6*self.b2*self.b3*c3*self.mu3+6*self.b1*c2*self.b3+6*c3*self.b1**2*self.mu1-2*self.b3**2*c2+2*c1*self.b1*self.b3-6*c2*self.b3*self.b1*self.mu1-6*c3*self.b3*self.b1*self.mu1-6*c1*self.b3**2*self.mu3+6*self.b1*self.b2*c3-6*c3*self.b2*self.b1*self.mu1+4*self.b1*c3*self.b3+6*c1*self.b3**2-6*self.b1*c3*self.b3*self.mu3+2*c3*self.b1**2-2*self.b2*self.b3*c3+6*self.b1*c1*self.b3*self.mu1-6*c2*self.b3**2*self.mu3)*L1-2*c1*self.b1*self.b3-2*self.b1*self.b2*c3-2*self.b1*c2*self.b3)*L2+c1*self.b2*self.b3+self.b1*c2*self.b3-0.5*(self.b1*self.b2*c3-3*self.b2*self.b3*c3*self.mu3+3*self.b2*self.b3*c3+self.b1*c2*self.b3+3*self.b3**2*c2-3*c2*self.b3**2*self.mu3+3*c2*self.b3*self.b1*self.mu1+3*c3*self.b2*self.b1*self.mu1)*L1**2-0.5*(4*c1*self.b2*self.b3+4*self.b1*c2*self.b3)*L1,
#                   -(-3*c2*self.b1*self.mu1+3*c2*self.b1*self.mu2-4*c2*self.b1-4*self.b2*c1-3*c1*self.b2*self.mu1+3*c1*self.b2*self.mu2)*L3**2+(-(-6*c3*self.b2*self.mu1+4*self.b2*c1-8*self.b2*c3+4*c2*self.b1+6*c3*self.b2*self.mu2+6*c1*self.b2*self.mu1-6*c2*self.b2*self.mu1+6*c2*self.b3*self.mu2+6*c2*self.b1*self.mu1+6*c2*self.b1*self.mu2-6*c2*self.b3*self.mu1+6*c1*self.b2*self.mu2-8*c2*self.b3-6*c2*self.b2*self.mu2+4*c2*self.b2)*L1-(-6*c3*self.b1*self.mu1+6*c1*self.b1*self.mu2-6*c1*self.b2*self.mu1-6*c2*self.b1*self.mu2+4*c1*self.b1+6*c1*self.b3*self.mu2+6*c1*self.b1*self.mu1-6*c1*self.b3*self.mu1-8*c3*self.b1-8*self.b3*c1+4*c2*self.b1-6*c2*self.b1*self.mu1+4*self.b2*c1+6*c3*self.b1*self.mu2-6*c1*self.b2*self.mu2)*L2-4*c2*self.b2+4*c3*self.b1+4*self.b3*c1)*L3-(-3*c3*self.b1*self.mu2+2*self.b3*c1-3*c1*self.b3*self.mu1-3*c3*self.b1*self.mu1+2*c3*self.b1-3*c1*self.b3*self.mu2)*L2**2+(-(4*c3*self.b1+4*self.b3*c1+4*c2*self.b3+6*c3*self.b3*self.mu2-8*c3*self.b3+6*c3*self.b1*self.mu1-6*c2*self.b3*self.mu2-6*c2*self.b3*self.mu1+6*c1*self.b3*self.mu1-6*c3*self.b3*self.mu1+6*c1*self.b3*self.mu2-6*c3*self.b2*self.mu1+4*self.b2*c3+6*c3*self.b1*self.mu2-6*c3*self.b2*self.mu2)*L1-4*c2*self.b3-4*self.b2*c3)*L2+c2*self.b3-self.b3*c1-c3*self.b1+self.b2*c3-(2*c2*self.b3+3*c3*self.b2*self.mu1+3*c3*self.b2*self.mu2+2*self.b2*c3+3*c2*self.b3*self.mu2+3*c2*self.b3*self.mu1)*L1**2+4*c3*self.b3*L1,
#                   0.5*(-c1**2*self.b2-3*c2*self.b1*c1*self.mu1-3*c1*c2*self.b2*self.mu2+3*c1*c2*self.b2-c1*c2*self.b1-3*c2**2*self.b1*self.mu2+3*c2**2*self.b1-3*self.b2*c1**2*self.mu1)*L3**2+(0.5*(2*c1*c2*self.b1+4*c1*c2*self.b2-6*c3*self.b2*c2*self.mu2+6*self.b2*c1**2*self.mu1-2*c2**2*self.b1+2*c1**2*self.b2+6*c2**2*self.b3+6*c3*self.b2*c2-6*c2*self.b3*c1*self.mu1-2*c1*self.b2*c3-6*c2**2*self.b1*self.mu2-6*c2*self.b2*c1*self.mu1+2*c2**2*self.b2+6*c2*self.b1*c1*self.mu1-6*c3*self.b2*c1*self.mu1-2*c1*c2*self.b3-6*c2**2*self.b3*self.mu2+6*c2**2*self.b2*self.mu2-6*c1*c2*self.b2*self.mu2)*L1+0.5*(4*c1*c2*self.b1+6*c1*c2*self.b3+6*c2**2*self.b1*self.mu2-6*c3*self.b1*c2*self.mu2+6*c1**2*self.b2-6*c3*self.b1*c1*self.mu1+6*c3*c2*self.b1-2*c1*c3*self.b1-2*c1**2*self.b3-6*self.b3*c1**2*self.mu1+6*self.b1*c1**2*self.mu1-6*self.b2*c1**2*self.mu1-6*c1*c2*self.b3*self.mu2+6*c1*c2*self.b2*self.mu2-6*c2*self.b1*c1*self.mu1+2*c1*c2*self.b2-6*c1*c2*self.b1*self.mu2+2*c2**2*self.b1+2*self.b1*c1**2)*L2+2*c1*c2*self.b3+2*c1*c2*self.b2+2*c3*c2*self.b1)*L3+0.5*(-3*c3*self.b1*c1*self.mu1+3*c1*c2*self.b3*self.mu2+3*c1**2*self.b3+c3*c2*self.b1+3*c1*c3*self.b1+c1*c2*self.b3-3*self.b3*c1**2*self.mu1+3*c3*self.b1*c2*self.mu2)*L2**2+(0.5*(-2*c3*c2*self.b1+4*c1*c2*self.b3-2*c1*c3*self.b3+2*c1**2*self.b3+6*c3*self.b3*c2+6*c1*self.b2*c3-6*c2*self.b3*c1*self.mu1+2*c1*c3*self.b1+6*c3*self.b1*c1*self.mu1-6*c3*self.b3*c1*self.mu1-6*c3*self.b2*c1*self.mu1-6*c3*self.b1*c2*self.mu2+6*c3*self.b2*c2*self.mu2+2*c3*self.b2*c2+2*c2**2*self.b3+6*self.b3*c1**2*self.mu1+6*c2**2*self.b3*self.mu2-6*c1*c2*self.b3*self.mu2-6*c3*c2*self.b3*self.mu2)*L1+2*c1*self.b2*c3+2*c1*c2*self.b3)*L2-c1*self.b2*c3-c1*c2*self.b3+0.5*(c1*c2*self.b3+3*c2*self.b3*c1*self.mu1-c3*self.b2*c2+3*c3*self.b2*c1*self.mu1-3*c2**2*self.b3*self.mu2-c2**2*self.b3-3*c3*self.b2*c2*self.mu2+c1*self.b2*c3)*L1**2+2*c3*self.b3*c2*L1,
#                   -0.5*(-3*self.b1*c2*self.b2*self.mu2-c2*self.b1**2+3*c1*self.b2**2-3*c1*self.b2**2*self.mu2-c1*self.b1*self.b2-3*c2*self.b1**2*self.mu1-3*self.b1*c1*self.b2*self.mu1+3*self.b1*c2*self.b2)*L3**2+(-0.5*(-6*c1*self.b2**2*self.mu2+2*c2*self.b1**2+2*c2*self.b2**2-2*c1*self.b2**2-2*self.b1*self.b2*c3+6*c3*self.b2**2+6*c2*self.b2**2*self.mu2+4*self.b1*c2*self.b2+6*c2*self.b1**2*self.mu1-6*c2*self.b3*self.b1*self.mu1-6*c2*self.b2*self.b1*self.mu1+6*self.b1*c1*self.b2*self.mu1-6*c3*self.b2**2*self.mu2+2*c1*self.b1*self.b2-6*self.b1*c2*self.b2*self.mu2+6*self.b2*c2*self.b3-6*self.b2*c2*self.b3*self.mu2-6*c3*self.b2*self.b1*self.mu1-2*self.b1*c2*self.b3)*L1-0.5*(2*c1*self.b2**2+6*c2*self.b1**2+2*self.b1*c2*self.b2-6*c1*self.b2*self.b3*self.mu2-6*c3*self.b1**2*self.mu1+4*c1*self.b1*self.b2+6*c1*self.b1**2*self.mu1-6*c2*self.b1**2*self.mu1-6*c3*self.b1*self.b2*self.mu2+6*self.b1*self.b2*c3-6*self.b1*c1*self.b2*self.mu1+6*c1*self.b2*self.b3+2*c1*self.b1**2-6*self.b1*c1*self.b3*self.mu1-2*c1*self.b1*self.b3-2*c3*self.b1**2-6*c1*self.b1*self.b2*self.mu2+6*self.b1*c2*self.b2*self.mu2+6*c1*self.b2**2*self.mu2)*L2-2*self.b1*c2*self.b2-2*c1*self.b2*self.b3-2*self.b1*self.b2*c3)*L3-0.5*(3*c1*self.b1*self.b3+3*c3*self.b1*self.b2*self.mu2-3*self.b1*c1*self.b3*self.mu1+3*c1*self.b2*self.b3*self.mu2-3*c3*self.b1**2*self.mu1+self.b1*self.b2*c3+3*c3*self.b1**2+c1*self.b2*self.b3)*L2**2+(-0.5*(2*c3*self.b1**2+6*self.b1*c1*self.b3*self.mu1+6*c3*self.b1**2*self.mu1-6*c1*self.b2*self.b3*self.mu2-2*c1*self.b2*self.b3-6*c2*self.b3*self.b1*self.mu1+2*c3*self.b2**2+6*self.b1*c2*self.b3+6*self.b2*c2*self.b3*self.mu2-2*self.b1*c3*self.b3+2*c1*self.b1*self.b3+6*self.b2*self.b3*c3+2*self.b2*c2*self.b3-6*c3*self.b1*self.b2*self.mu2-6*c3*self.b2*self.b3*self.mu2-6*c3*self.b3*self.b1*self.mu1+4*self.b1*self.b2*c3-6*c3*self.b2*self.b1*self.mu1+6*c3*self.b2**2*self.mu2)*L1-2*self.b1*self.b2*c3-2*self.b1*c2*self.b3)*L2+self.b1*self.b2*c3+self.b1*c2*self.b3-0.5*(self.b1*c2*self.b3-c3*self.b2**2+3*c2*self.b3*self.b1*self.mu1-3*c3*self.b2**2*self.mu2-3*self.b2*c2*self.b3*self.mu2+self.b1*self.b2*c3+3*c3*self.b2*self.b1*self.mu1-self.b2*c2*self.b3)*L1**2-2*self.b2*self.b3*c3*L1]])

        B =array([[-2*(-3*b1*b2*mu3+2*b1*b2-3*b1*b2*mu2)*L3**2+(-2*(4*b3*b2-6*b2*b3*mu2+2*b2**2-8*b1*b2-6*b2*b3*mu3+3*b2**2*mu3+3*b2**2*mu2-6*b1*b2*mu2+6*b1*b2*mu3)*L1-2*(-3*b1**2*mu2+4*b1*b2+4*b1*b3+6*b1*b2*mu2-4*b1**2+6*b1*b2*mu3-6*b1*b3*mu2+3*b1**2*mu3-6*b1*b3*mu3)*L2-8*b1*b3)*L3-2*(2*b1*b3+3*b1*b3*mu2+3*b1*b3*mu3)*L2**2+(-2*(2*b3**2-3*b3**2*mu2-8*b1*b3-3*b3**2*mu3+4*b3*b2+6*b1*b3*mu3+6*b2*b3*mu2+6*b2*b3*mu3-6*b1*b3*mu2)*L1+4*b1**2)*L2-2*b1*b2+2*b1*b3-2*(-3*b2*b3*mu2-4*b3*b2+3*b2*b3*mu3)*L1**2-2*(-4*b1*b2+2*b3**2)*L1,
                   -(3*b1*c2*b2*mu2-b1*b2*c3-3*b1*c3*b2*mu3-3*b1*c2*b2)*L3**2+(-(-c2*b2**2+6*b2*c2*b3*mu2+6*b1*c2*b2*mu2-6*b1*b2*c3+2*b1*c2*b2-6*b2*b3*c3*mu3-3*c2*b2**2*mu2+6*b1*c3*b2*mu3+3*c3*b2**2*mu3+c3*b2**2-6*b2*c2*b3-2*b2*b3*c3)*L1-(3*c2*b1**2*mu2-2*b1*c2*b2+c2*b1**2-6*b1*c2*b2*mu2+6*b1*c2*b3*mu2-6*b1*c3*b3*mu3-2*b1*c3*b3+2*b1*b2*c3+6*b1*c3*b2*mu3-3*c3*b1**2+3*c3*b1**2*mu3-6*b1*c2*b3)*L2+4*b1*c2*b3)*L3-(-3*b1*c2*b3*mu2-b1*c2*b3+b1*c3*b3+3*b1*c3*b3*mu3)*L2**2+(-(2*b1*c2*b3-3*b3**2*c2+6*b2*b3*c3*mu3-2*b2*c2*b3-6*b2*c2*b3*mu2-c3*b3**2-6*b1*c3*b3+2*b2*b3*c3+6*b1*c2*b3*mu2+6*b1*c3*b3*mu3+3*c2*b3**2*mu2-3*c3*b3**2*mu3)*L1+2*c3*b1**2)*L2-2*b1*c2*b3-(3*b2*b3*c3*mu3+b2*c2*b3+3*b2*c2*b3*mu2-3*b2*b3*c3)*L1**2-(-4*b1*b2*c3-2*b3**2*c2)*L1,
                   (-b1*b3*b2+3*b1*b2**2*mu2-3*b1*b2*b3*mu3-3*b1*b2**2)*L3**2+((6*b1*b2**2*mu2+6*b1*b2*b3*mu3-6*b2*b3**2*mu3-6*b1*b3*b2+2*b1*b2**2-2*b2*b3**2-3*b2**3*mu2-5*b3*b2**2-b2**3+3*b2**2*b3*mu3+6*b3*b2**2*mu2)*L1+(-4*b1*b3*b2+6*b1*b3*b2*mu2-6*b1*b2**2*mu2-6*b1*b3**2*mu3+b1**2*b2-2*b1*b3**2-2*b1*b2**2+3*b1**2*b3*mu3-3*b1**2*b3+3*b1**2*b2*mu2+6*b1*b2*b3*mu3)*L2-4*b1*b3*b2)*L3+(-3*b1*b3*b2*mu2-b1*b3*b2+3*b1*b3**2*mu3+b1*b3**2)*L2**2+((-b2*b3**2+6*b1*b3*b2*mu2-b3**3+2*b1*b3*b2+6*b2*b3**2*mu3-3*b3**3*mu3-6*b1*b3**2-2*b3*b2**2+6*b1*b3**2*mu3+3*b3**2*b2*mu2-6*b3*b2**2*mu2)*L1-2*b1**2*b3)*L2+2*b1*b3*b2+(3*b3*b2**2*mu2+3*b2*b3**2*mu3-3*b2*b3**2+b3*b2**2)*L1**2+(-4*b1*b3*b2-2*b2*b3**2)*L1,
                   2*(-3*b1*b2*mu3-3*b1*b2*mu1-2*b1*b2)*L3**2+(2*(-3*b2**2*mu1+6*b1*b2*mu3+6*b1*b2*mu1-4*b3*b2+4*b2**2-4*b1*b2-6*b2*b3*mu3-6*b2*b3*mu1+3*b2**2*mu3)*L1+2*(3*b1**2*mu1-6*b1*b3*mu3-2*b1**2-6*b1*b2*mu1-6*b1*b3*mu1+6*b1*b2*mu3+3*b1**2*mu3+8*b1*b2-4*b1*b3)*L2+4*b2**2)*L3+2*(4*b1*b3-3*b1*b3*mu1+3*b1*b3*mu3)*L2**2+(2*(-3*b3**2*mu3-6*b2*b3*mu1-3*b3**2*mu1-4*b1*b3+6*b2*b3*mu3+6*b1*b3*mu3-2*b3**2+8*b3*b2+6*b1*b3*mu1)*L1+8*b3*b2-4*b1**2)*L2+2*b1*b2-2*b3*b2+2*(3*b2*b3*mu3-2*b3*b2+3*b2*b3*mu1)*L1**2-8*b1*b2*L1,
                   -(-b1*b2*c3+3*b1*c1*b2*mu1+c1*b1*b2-3*b1*c3*b2*mu3)*L3**2+(-(-6*b1*c1*b2*mu1-6*b2*b3*c3*mu3-6*b1*b2*c3+c3*b2**2+2*c1*b3*b2+3*c3*b2**2*mu3+3*c1*b2**2*mu1-2*b2*b3*c3-3*c1*b2**2+6*b1*c3*b2*mu3+6*b2*b3*c1*mu1-2*c1*b1*b2)*L1-(-6*c1*b1*b2-2*b1*c3*b3-3*c3*b1**2+6*b1*c1*b3*mu1+2*b1*b2*c3-6*b1*c3*b3*mu3+6*b1*c3*b2*mu3-c1*b1**2+2*c1*b1*b3-3*c1*b1**2*mu1+6*b1*c1*b2*mu1+3*c3*b1**2*mu3)*L2+2*c1*b2**2)*L3-(-3*c1*b1*b3+3*b1*c1*b3*mu1+b1*c3*b3+3*b1*c3*b3*mu3)*L2**2+(-(6*b2*b3*c3*mu3-c3*b3**2-6*b1*c1*b3*mu1-6*c1*b3*b2-3*c3*b3**2*mu3+6*b2*b3*c1*mu1+3*c1*b3**2*mu1+6*b1*c3*b3*mu3-6*b1*c3*b3-2*c1*b1*b3+2*b2*b3*c3+c1*b3**2)*L1+4*c1*b3*b2+2*c3*b1**2)*L2-2*b1*b2*c3-(-3*b2*b3*c3-c1*b3*b2+3*b2*b3*c3*mu3-3*b2*b3*c1*mu1)*L1**2+4*b1*b2*c3*L1,
                   -(-b1**2*b2+3*b1*b2*b3*mu3+b1*b3*b2-3*b2*b1**2*mu1)*L3**2+(-(6*b2*b1**2*mu1-6*b1*b2*b3*mu3-b3*b2**2+2*b1**2*b2+4*b1*b3*b2+2*b2*b3**2-6*b2*b3*b1*mu1-3*b2**2*b3*mu3-3*b2**2*b1*mu1+6*b2*b3**2*mu3+3*b1*b2**2)*L1-(-6*b1*b2*b3*mu3+b1**3+2*b1*b3**2-6*b3*b1**2*mu1-2*b1*b3*b2+6*b1*b3**2*mu3+b1**2*b3-3*b1**2*b3*mu3+6*b1**2*b2+3*b1**3*mu1-6*b2*b1**2*mu1)*L2-2*b1*b2**2)*L3-(-3*b3*b1**2*mu1+3*b1**2*b3-3*b1*b3**2*mu3-b1*b3**2)*L2**2+(-(6*b3*b1**2*mu1-6*b1*b3**2*mu3+3*b3**3*mu3+6*b1*b3*b2-2*b2*b3**2+b3**3+2*b1**2*b3+5*b1*b3**2-6*b2*b3**2*mu3-3*b3**2*b1*mu1-6*b2*b3*b1*mu1)*L1-4*b1*b3*b2-2*b1**2*b3)*L2+2*b1*b3*b2-(3*b2*b3**2+3*b2*b3*b1*mu1+b1*b3*b2-3*b2*b3**2*mu3)*L1**2-4*b1*b3*b2*L1,
                   -2*(-4*b1*b2+3*b1*b2*mu2-3*b1*b2*mu1)*L3**2+(-2*(4*b1*b2+6*b2*b3*mu2-3*b2**2*mu1-6*b2*b3*mu1-8*b3*b2+2*b2**2+6*b1*b2*mu1-3*b2**2*mu2+6*b1*b2*mu2)*L1-2*(2*b1**2-8*b1*b3-6*b1*b2*mu1+4*b1*b2+3*b1**2*mu2+3*b1**2*mu1-6*b1*b2*mu2+6*b1*b3*mu2-6*b1*b3*mu1)*L2-4*b2**2+8*b1*b3)*L3-2*(-3*b1*b3*mu2+2*b1*b3-3*b1*b3*mu1)*L2**2+(-2*(-6*b2*b3*mu1-6*b2*b3*mu2+6*b1*b3*mu2+4*b3*b2+3*b3**2*mu2+4*b1*b3-3*b3**2*mu1+6*b1*b3*mu1-4*b3**2)*L1-8*b3*b2)*L2-2*b1*b3+2*b3*b2-2*(2*b3*b2+3*b2*b3*mu2+3*b2*b3*mu1)*L1**2+4*b3**2*L1,
                   (-3*b1*c1*b2*mu1-c1*b1*b2-3*b1*c2*b2*mu2+3*b1*c2*b2)*L3**2+((-6*b2*c2*b3*mu2-6*b2*b3*c1*mu1+6*b2*c2*b3+6*b1*c1*b2*mu1-3*c1*b2**2*mu1-2*b1*c2*b2-2*c1*b3*b2-6*b1*c2*b2*mu2+c2*b2**2+3*c2*b2**2*mu2+2*c1*b1*b2+3*c1*b2**2)*L1+(2*b1*c2*b2-6*b1*c1*b3*mu1-6*b1*c2*b3*mu2-2*c1*b1*b3+3*c1*b1**2*mu1+6*c1*b1*b2-c2*b1**2+6*b1*c2*b3+6*b1*c2*b2*mu2+c1*b1**2-3*c2*b1**2*mu2-6*b1*c1*b2*mu1)*L2+4*b1*c2*b3+2*c1*b2**2)*L3+(-3*b1*c1*b3*mu1+3*b1*c2*b3*mu2+3*c1*b1*b3+b1*c2*b3)*L2**2+((-6*b2*b3*c1*mu1-2*b1*c2*b3-6*b1*c2*b3*mu2-3*c2*b3**2*mu2-3*c1*b3**2*mu1+2*c1*b1*b3+2*b2*c2*b3-c1*b3**2+6*c1*b3*b2+6*b1*c1*b3*mu1+6*b2*c2*b3*mu2+3*b3**2*c2)*L1+4*c1*b3*b2)*L2-2*c1*b3*b2+(c1*b3*b2-3*b2*c2*b3*mu2-b2*c2*b3+3*b2*b3*c1*mu1)*L1**2+2*b3**2*c2*L1,
                   -(-b1**2*b2-3*b2*b1**2*mu1-3*b1*b2**2*mu2+3*b1*b2**2)*L3**2+(-(2*b1**2*b2+6*b2*b1**2*mu1+b1*b2**2-3*b2**2*b1*mu1+3*b2**3*mu2-6*b3*b2**2*mu2-2*b1*b3*b2-6*b2*b3*b1*mu1+6*b3*b2**2-6*b1*b2**2*mu2+b2**3)*L1-(-6*b3*b1**2*mu1+6*b1*b2**2*mu2-3*b1**2*b2*mu2-6*b2*b1**2*mu1+b1**3-2*b1**2*b3+5*b1**2*b2-6*b1*b3*b2*mu2+3*b1**3*mu1+2*b1*b2**2+6*b1*b3*b2)*L2-2*b1*b2**2-4*b1*b3*b2)*L3-(3*b1**2*b3+b1*b3*b2+3*b1*b3*b2*mu2-3*b3*b1**2*mu1)*L2**2+(-(-6*b1*b3*b2*mu2+2*b3*b2**2+2*b1**2*b3+3*b2*b3**2+6*b3*b2**2*mu2-b1*b3**2-3*b3**2*b1*mu1+6*b3*b1**2*mu1-3*b3**2*b2*mu2+4*b1*b3*b2-6*b2*b3*b1*mu1)*L1-4*b1*b3*b2)*L2+2*b1*b3*b2-(-b3*b2**2-3*b3*b2**2*mu2+3*b2*b3*b1*mu1+b1*b3*b2)*L1**2-2*b2*b3**2*L1],
                  [2*(3*c1*c2*mu3+3*c1*c2*mu2-2*c2*c1)*L3**2+(2*(-3*c2**2*mu2-2*c2**2+8*c2*c1+6*c1*c2*mu2+6*c2*c3*mu3+6*c2*c3*mu2-6*c1*c2*mu3-3*c2**2*mu3-4*c3*c2)*L1+2*(-4*c3*c1-3*c1**2*mu3-6*c1*c2*mu2+3*c1**2*mu2+4*c1**2-4*c2*c1+6*c1*c3*mu3-6*c1*c2*mu3+6*c1*c3*mu2)*L2-8*c3*c1)*L3+2*(-2*c3*c1-3*c1*c3*mu2-3*c1*c3*mu3)*L2**2+(2*(-6*c1*c3*mu3-6*c2*c3*mu2+3*c3**2*mu2+6*c1*c3*mu2-6*c2*c3*mu3+3*c3**2*mu3-4*c3*c2+8*c3*c1-2*c3**2)*L1+4*c1**2)*L2+2*c3*c1-2*c2*c1+2*(-3*c2*c3*mu3+4*c3*c2+3*c2*c3*mu2)*L1**2+2*(4*c2*c1-2*c3**2)*L1,
                   (3*c1*c2**2+c1*c3*c2-3*c1*c2**2*mu2+3*c1*c2*c3*mu3)*L3**2+((c2**3-6*c1*c2*c3*mu3-6*c1*c2**2*mu2-3*c2**2*c3*mu3+3*c2**3*mu2+6*c2*c3**2*mu3+5*c2**2*c3+6*c1*c3*c2-6*c2**2*c3*mu2-2*c1*c2**2+2*c2*c3**2)*L1+(-6*c1*c2*c3*mu3+6*c1*c3**2*mu3-3*c1**2*c3*mu3-3*c1**2*c2*mu2+6*c1*c2**2*mu2+4*c1*c3*c2-c1**2*c2+2*c1*c2**2+2*c1*c3**2-6*c1*c2*c3*mu2+3*c1**2*c3)*L2+4*c1*c3*c2)*L3+(-c1*c3**2+3*c1*c2*c3*mu2+c1*c3*c2-3*c1*c3**2*mu3)*L2**2+((-6*c2*c3**2*mu3+6*c2**2*c3*mu2+c3**3-6*c1*c3**2*mu3-3*c2*c3**2*mu2-2*c1*c3*c2+6*c1*c3**2+2*c2**2*c3+3*c3**3*mu3+c2*c3**2-6*c1*c2*c3*mu2)*L1+2*c1**2*c3)*L2-2*c1*c3*c2+(-c2**2*c3-3*c2*c3**2*mu3-3*c2**2*c3*mu2+3*c2*c3**2)*L1**2+(4*c1*c3*c2+2*c2*c3**2)*L1,
                   (3*c1*c2*b2*mu2-3*c1*c2*b3*mu3-c1*c2*b3-3*c1*c2*b2)*L3**2+((6*c1*c2*b2*mu2+3*c2**2*b3*mu3+6*c1*c2*b3*mu3-c2**2*b2+c2**2*b3+6*c3*b2*c2*mu2-6*c1*c2*b3+2*c1*c2*b2-6*c2*b3*c3*mu3-2*c3*b3*c2-6*c3*b2*c2-3*c2**2*b2*mu2)*L1+(-3*c1**2*b3-6*c1*b2*c3+2*c1*c2*b3+3*c1**2*b2*mu2-6*c1*c2*b2*mu2-6*c1*c3*b3*mu3-2*c1*c3*b3+6*c1*c3*b2*mu2-2*c1*c2*b2+c1**2*b2+6*c1*c2*b3*mu3+3*c1**2*b3*mu3)*L2-4*c1*b2*c3)*L3+(-3*c1*c3*b2*mu2+c1*c3*b3-c1*b2*c3+3*c1*c3*b3*mu3)*L2**2+((-6*c3*b2*c2*mu2+6*c1*c3*b3*mu3-c3**2*b3-2*c3*b2*c2+2*c1*b2*c3+2*c3*b3*c2+6*c1*c3*b2*mu2-3*c3**2*b3*mu3-6*c1*c3*b3+6*c2*b3*c3*mu3+3*c3**2*b2*mu2-3*c3**2*b2)*L1-2*c1**2*b3)*L2+2*c1*b2*c3+(-3*c3*b3*c2+3*c3*b2*c2*mu2+3*c2*b3*c3*mu3+c3*b2*c2)*L1**2+(-4*c1*c2*b3-2*c3**2*b2)*L1,
                   -2*(3*c1*c2*mu3+2*c2*c1+3*c1*c2*mu1)*L3**2+(-2*(6*c2*c3*mu3-4*c2**2-6*c1*c2*mu3-3*c2**2*mu3+4*c3*c2+3*c2**2*mu1-6*c1*c2*mu1+6*c2*c3*mu1+4*c2*c1)*L1-2*(6*c1*c2*mu1-6*c1*c2*mu3+6*c1*c3*mu3+6*c1*c3*mu1-3*c1**2*mu1+4*c3*c1-8*c2*c1-3*c1**2*mu3+2*c1**2)*L2+4*c2**2)*L3-2*(3*c1*c3*mu1-4*c3*c1-3*c1*c3*mu3)*L2**2+(-2*(3*c3**2*mu1+3*c3**2*mu3-6*c1*c3*mu1+6*c2*c3*mu1-6*c1*c3*mu3+4*c3*c1-8*c3*c2-6*c2*c3*mu3+2*c3**2)*L1-4*c1**2+8*c3*c2)*L2+2*c2*c1-2*c3*c2-2*(-3*c2*c3*mu1-3*c2*c3*mu3+2*c3*c2)*L1**2-8*c2*c1*L1,
                   (c1*c3*c2-c1**2*c2-3*c1**2*c2*mu1+3*c1*c2*c3*mu3)*L3**2+((6*c1**2*c2*mu1-c2**2*c3-6*c2*c3*c1*mu1-3*c2**2*c3*mu3-3*c1*c2**2*mu1+6*c2*c3**2*mu3+4*c1*c3*c2+2*c2*c3**2-6*c1*c2*c3*mu3+3*c1*c2**2+2*c1**2*c2)*L1+(-2*c1*c3*c2+6*c1**2*c2+3*c1**3*mu1-6*c1*c2*c3*mu3+2*c1*c3**2+c1**2*c3+6*c1*c3**2*mu3+c1**3-6*c1**2*c3*mu1-6*c1**2*c2*mu1-3*c1**2*c3*mu3)*L2+2*c1*c2**2)*L3+(3*c1**2*c3-3*c1*c3**2*mu3-3*c1**2*c3*mu1-c1*c3**2)*L2**2+((5*c1*c3**2+3*c3**3*mu3-6*c2*c3**2*mu3-6*c1*c3**2*mu3+6*c1**2*c3*mu1+c3**3-6*c2*c3*c1*mu1+6*c1*c3*c2-2*c2*c3**2-3*c1*c3**2*mu1+2*c1**2*c3)*L1+4*c1*c3*c2+2*c1**2*c3)*L2-2*c1*c3*c2+(3*c2*c3**2+c1*c3*c2+3*c2*c3*c1*mu1-3*c2*c3**2*mu3)*L1**2+4*c1*c3*c2*L1,
                   (3*c2*b1*c1*mu1-c1*c2*b3+c1*c2*b1-3*c1*c2*b3*mu3)*L3**2+((-2*c1*c2*b1+3*c2**2*b3*mu3+6*c1*c2*b3*mu3-6*c2*b1*c1*mu1+3*c2**2*b1*mu1+2*c3*c2*b1+6*c2*c3*b1*mu1-3*c2**2*b1+c2**2*b3-6*c2*b3*c3*mu3-6*c1*c2*b3-2*c3*b3*c2)*L1+(-c1**2*b1+2*c1*c3*b1-3*c1**2*b1*mu1+6*c2*b1*c1*mu1-6*c1*c2*b1-3*c1**2*b3-6*c1*c3*b3*mu3+6*c1*c2*b3*mu3+2*c1*c2*b3-2*c1*c3*b3+6*c3*b1*c1*mu1+3*c1**2*b3*mu3)*L2-2*c2**2*b1)*L3+(3*c3*b1*c1*mu1+c1*c3*b3+3*c1*c3*b3*mu3-3*c1*c3*b1)*L2**2+((3*c3**2*b1*mu1+c3**2*b1-6*c3*b1*c1*mu1+6*c1*c3*b3*mu3-2*c1*c3*b1+6*c2*b3*c3*mu3-6*c3*c2*b1+2*c3*b3*c2-c3**2*b3-3*c3**2*b3*mu3-6*c1*c3*b3+6*c2*c3*b1*mu1)*L1-2*c1**2*b3-4*c3*c2*b1)*L2+2*c1*c2*b3+(3*c2*b3*c3*mu3-c3*c2*b1-3*c3*b3*c2-3*c2*c3*b1*mu1)*L1**2-4*c1*c2*b3*L1,
                   2*(-3*c1*c2*mu2+4*c2*c1+3*c1*c2*mu1)*L3**2+(2*(6*c2*c3*mu1-2*c2**2+3*c2**2*mu2-6*c2*c3*mu2+3*c2**2*mu1-6*c1*c2*mu1-4*c2*c1-6*c1*c2*mu2+8*c3*c2)*L1+2*(-2*c1**2+8*c3*c1-6*c1*c3*mu2-3*c1**2*mu1-4*c2*c1-3*c1**2*mu2+6*c1*c2*mu1+6*c1*c2*mu2+6*c1*c3*mu1)*L2-4*c2**2+8*c3*c1)*L3+2*(-2*c3*c1+3*c1*c3*mu2+3*c1*c3*mu1)*L2**2+(2*(-3*c3**2*mu2+3*c3**2*mu1-6*c1*c3*mu1+6*c2*c3*mu1-6*c1*c3*mu2-4*c3*c2+4*c3**2-4*c3*c1+6*c2*c3*mu2)*L1-8*c3*c2)*L2-2*c3*c1+2*c3*c2+2*(-3*c2*c3*mu2-2*c3*c2-3*c2*c3*mu1)*L1**2+4*c3**2*L1,
                   -(3*c1*c2**2*mu2-3*c1*c2**2+3*c1**2*c2*mu1+c1**2*c2)*L3**2+(-(-2*c1**2*c2-6*c2**2*c3+6*c2**2*c3*mu2-c1*c2**2+2*c1*c3*c2-6*c1**2*c2*mu1-3*c2**3*mu2+3*c1*c2**2*mu1+6*c2*c3*c1*mu1+6*c1*c2**2*mu2-c2**3)*L1-(-5*c1**2*c2-3*c1**3*mu1+2*c1**2*c3-c1**3+6*c1*c2*c3*mu2-2*c1*c2**2+6*c1**2*c3*mu1-6*c1*c3*c2+6*c1**2*c2*mu1-6*c1*c2**2*mu2+3*c1**2*c2*mu2)*L2+2*c1*c2**2+4*c1*c3*c2)*L3-(-3*c1**2*c3-c1*c3*c2-3*c1*c2*c3*mu2+3*c1**2*c3*mu1)*L2**2+(-(6*c2*c3*c1*mu1-6*c1**2*c3*mu1+c1*c3**2-2*c1**2*c3-3*c2*c3**2-2*c2**2*c3+6*c1*c2*c3*mu2+3*c2*c3**2*mu2+3*c1*c3**2*mu1-4*c1*c3*c2-6*c2**2*c3*mu2)*L1+4*c1*c3*c2)*L2-2*c1*c3*c2-(-3*c2*c3*c1*mu1+3*c2**2*c3*mu2+c2**2*c3-c1*c3*c2)*L1**2+2*c2*c3**2*L1,
                   (3*c1*c2*b2*mu2-3*c1*c2*b2+c1*c2*b1+3*c2*b1*c1*mu1)*L3**2+((-2*c1*c2*b1+6*c3*b2*c2*mu2-3*c2**2*b2*mu2-6*c3*b2*c2-6*c2*b1*c1*mu1-c2**2*b2+2*c1*c2*b2+3*c2**2*b1*mu1+2*c3*c2*b1+6*c2*c3*b1*mu1-3*c2**2*b1+6*c1*c2*b2*mu2)*L1+(6*c3*b1*c1*mu1-c1**2*b1+2*c1*c3*b1-6*c1*c2*b1-3*c1**2*b1*mu1-6*c1*c2*b2*mu2-6*c1*b2*c3-2*c1*c2*b2+3*c1**2*b2*mu2+6*c1*c3*b2*mu2+6*c2*b1*c1*mu1+c1**2*b2)*L2-2*c2**2*b1-4*c1*b2*c3)*L3+(-3*c1*c3*b2*mu2+3*c3*b1*c1*mu1-c1*b2*c3-3*c1*c3*b1)*L2**2+((3*c3**2*b2*mu2+6*c1*c3*b2*mu2+6*c2*c3*b1*mu1-2*c1*c3*b1+c3**2*b1-2*c3*b2*c2+2*c1*b2*c3-3*c3**2*b2-6*c3*b1*c1*mu1+3*c3**2*b1*mu1-6*c3*c2*b1-6*c3*b2*c2*mu2)*L1-4*c3*c2*b1)*L2+2*c3*c2*b1+(c3*b2*c2+3*c3*b2*c2*mu2-c3*c2*b1-3*c2*c3*b1*mu1)*L1**2-2*c3**2*b2*L1],
                  [2.*(-(2*b2*c1-3*c2*b1*mu3+2*c2*b1-3*c1*b2*mu2-3*c2*b1*mu2-3*c1*b2*mu3)*L3**2+(-(-6*c1*b2*mu2+4*c2*b2-8*b2*c1-6*c2*b1*mu2+6*c2*b2*mu3+6*c2*b2*mu2+6*c2*b1*mu3+6*c1*b2*mu3-8*c2*b1-6*c3*b2*mu2+4*c2*b3+4*b2*c3-6*c2*b3*mu2-6*c2*b3*mu3-6*c3*b2*mu3)*L1-(6*c1*b1*mu3+6*c2*b1*mu3+4*b3*c1+6*c2*b1*mu2+6*c1*b2*mu3-6*c1*b3*mu2-6*c3*b1*mu2-6*c1*b1*mu2+4*c3*b1-8*c1*b1+6*c1*b2*mu2-6*c1*b3*mu3+4*b2*c1+4*c2*b1-6*c3*b1*mu3)*L2-4*c3*b1-4*b3*c1)*L3-(3*c3*b1*mu3+2*b3*c1+3*c1*b3*mu3+3*c3*b1*mu2+2*c3*b1+3*c1*b3*mu2)*L2**2+(-(-8*c3*b1-8*b3*c1+6*c1*b3*mu3-6*c1*b3*mu2-6*c3*b3*mu3+6*c2*b3*mu3+4*c3*b3+4*c2*b3+6*c3*b1*mu3+4*b2*c3+6*c3*b2*mu3-6*c3*b1*mu2+6*c3*b2*mu2+6*c2*b3*mu2-6*c3*b3*mu2)*L1+4*c1*b1)*L2-b2*c1+c3*b1-c2*b1+b3*c1-(-4*b2*c3-3*c3*b2*mu2+3*c2*b3*mu3-4*c2*b3+3*c3*b2*mu3-3*c2*b3*mu2)*L1**2-(4*c3*b3-4*b2*c1-4*c2*b1)*L1),
                   2.*(-0.5*(-3*c2**2*b1+3*c2**2*b1*mu2-3*c2*b1*c3*mu3-3*c1*c3*b2*mu3-c1*b2*c3-c3*c2*b1-3*c1*c2*b2+3*c1*c2*b2*mu2)*L3**2+(-0.5*(6*c2*c3*b2*mu3-2*c3*b3*c2-6*c3*c2*b1-2*c2**2*b2+2*c2**2*b1-6*c2**2*b2*mu2+6*c1*c3*b2*mu3+2*c1*c2*b2-2*c3**2*b2-6*c2**2*b3+6*c2*b1*c3*mu3+6*c2**2*b1*mu2-4*c3*b2*c2-6*c3**2*b2*mu3+6*c1*c2*b2*mu2-6*c1*b2*c3+6*c2**2*b3*mu2-6*c2*b3*c3*mu3+6*c3*b2*c2*mu2)*L1-0.5*(6*c1*c2*b1*mu2-4*c3*c2*b1-2*c1*c3*b3+6*c3*b1*c2*mu2-6*c1*c3*b1-2*c2**2*b1+6*c2*b1*c3*mu3+6*c1*c3*b1*mu3-2*c3**2*b1-6*c1*c2*b3-6*c1*c3*b3*mu3+6*c1*c3*b2*mu3-2*c1*c2*b2+6*c1*c2*b3*mu2-6*c2**2*b1*mu2-6*c1*c2*b2*mu2+2*c1*c2*b1-6*c3**2*b1*mu3+2*c1*b2*c3)*L2+2*c1*c2*b3+2*c3*c2*b1)*L3-0.5*(3*c1*c3*b3*mu3+c3**2*b1-3*c3*b1*c2*mu2-c3*c2*b1-3*c1*c2*b3*mu2+c1*c3*b3-c1*c2*b3+3*c3**2*b1*mu3)*L2**2+(-0.5*(-6*c3**2*b1-6*c2**2*b3*mu2+6*c1*c3*b3*mu3-6*c1*c3*b3+2*c1*c2*b3-2*c3**2*b3-6*c3**2*b3*mu3-4*c3*b3*c2-6*c3*b2*c2*mu2+6*c3*b1*c2*mu2+6*c3**2*b1*mu3+6*c3**2*b2*mu3-2*c2**2*b3+6*c1*c2*b3*mu2+2*c3**2*b2+2*c3*c2*b1-2*c3*b2*c2+6*c2*b3*c3*mu3+6*c3*c2*b3*mu2)*L1+2*c1*c3*b1)*L2-c3*c2*b1-c1*c2*b3-0.5*(3*c3**2*b2*mu3+3*c3*b2*c2*mu2+c3*b2*c2+3*c2*b3*c3*mu3+3*c2**2*b3*mu2-3*c3*b3*c2-3*c3**2*b2+c2**2*b3)*L1**2-0.5*(-4*c3*b3*c2-4*c3*c2*b1-4*c1*b2*c3)*L1),
                   2.*(0.5*(-3*b1*c2*b2-3*c1*b2**2-c1*b2*b3+3*c1*b2**2*mu2+3*b1*c2*b2*mu2-b1*c2*b3-3*c1*b2*b3*mu3-3*c2*b1*b3*mu3)*L3**2+(0.5*(2*c1*b2**2+2*b1*c2*b2-4*b2*c2*b3-6*c1*b2*b3-2*c2*b2**2+6*c1*b2**2*mu2-6*c3*b2**2-6*c2*b2**2*mu2+6*b2*c2*b3*mu2+6*c3*b2**2*mu2-6*b2*b3*c3*mu3-2*b2*b3*c3-6*c2*b3**2*mu3+6*c2*b2*b3*mu3+6*b1*c2*b2*mu2-2*b3**2*c2+6*c2*b1*b3*mu3+6*c1*b2*b3*mu3-6*b1*c2*b3)*L1+0.5*(-4*c1*b2*b3-2*c1*b2**2-6*c1*b1*b3+6*c1*b1*b2*mu2+6*c3*b1*b2*mu2+2*c1*b1*b2+6*c1*b2*b3*mu2-2*c1*b3**2-6*b1*b2*c3-6*b1*c3*b3*mu3-6*c1*b3**2*mu3-2*b1*c3*b3+6*c2*b1*b3*mu3-6*c1*b2**2*mu2+2*b1*c2*b3+6*c1*b2*b3*mu3-6*b1*c2*b2*mu2-2*b1*c2*b2+6*c1*b1*b3*mu3)*L2-2*b1*b2*c3-2*c1*b2*b3)*L3+0.5*(-3*c3*b1*b2*mu2+b1*c3*b3+c1*b3**2-b1*b2*c3+3*b1*c3*b3*mu3-3*c1*b2*b3*mu2+3*c1*b3**2*mu3-c1*b2*b3)*L2**2+(0.5*(-2*c3*b3**2-6*b1*c3*b3-6*b2*c2*b3*mu2+6*b1*c3*b3*mu3-6*c3*b3**2*mu3-4*b2*b3*c3+6*c3*b2*b3*mu2+6*c1*b2*b3*mu2+6*c1*b3**2*mu3-6*c3*b2**2*mu2+2*b1*b2*c3-6*c1*b3**2+6*c3*b1*b2*mu2+6*c2*b3**2*mu3+6*b2*b3*c3*mu3-2*c3*b2**2+2*b3**2*c2+2*c1*b2*b3-2*b2*c2*b3)*L1-2*c1*b1*b3)*L2+b1*b2*c3+c1*b2*b3+0.5*(3*b2*c2*b3*mu2+3*b2*b3*c3*mu3+3*c3*b2**2*mu2-3*b2*b3*c3-3*b3**2*c2+c3*b2**2+3*c2*b3**2*mu3+b2*c2*b3)*L1**2+0.5*(-4*b1*c2*b3-4*c1*b2*b3-4*b2*b3*c3)*L1),
                   2.*((-2*c2*b1-3*c1*b2*mu3-3*c2*b1*mu3-3*c2*b1*mu1-3*c1*b2*mu1-2*b2*c1)*L3**2+((6*c2*b2*mu3-6*c3*b2*mu1-4*c2*b3-4*c2*b1-6*c2*b2*mu1-4*b2*c1+6*c1*b2*mu3+6*c2*b1*mu3-6*c2*b3*mu1+8*c2*b2+6*c1*b2*mu1-6*c2*b3*mu3-4*b2*c3-6*c3*b2*mu3+6*c2*b1*mu1)*L1+(-4*c3*b1+6*c2*b1*mu3+6*c1*b2*mu3-4*b3*c1-6*c3*b1*mu3-6*c2*b1*mu1-6*c3*b1*mu1-6*c1*b3*mu3-6*c1*b2*mu1-6*c1*b3*mu1+6*c1*b1*mu3+8*b2*c1+8*c2*b1+6*c1*b1*mu1-4*c1*b1)*L2+4*c2*b2)*L3+(-3*c1*b3*mu1-3*c3*b1*mu1+4*b3*c1+3*c1*b3*mu3+3*c3*b1*mu3+4*c3*b1)*L2**2+((-6*c3*b3*mu3-4*b3*c1+8*b2*c3-4*c3*b1+6*c3*b2*mu3+6*c3*b1*mu3-6*c2*b3*mu1+6*c1*b3*mu3-4*c3*b3+6*c3*b1*mu1-6*c3*b2*mu1+6*c2*b3*mu3-6*c3*b3*mu1+8*c2*b3+6*c1*b3*mu1)*L1-4*c1*b1+4*b2*c3+4*c2*b3)*L2-c2*b3+c2*b1-b2*c3+b2*c1+(3*c2*b3*mu1+3*c3*b2*mu3-2*b2*c3+3*c2*b3*mu3-2*c2*b3+3*c3*b2*mu1)*L1**2+(-4*b2*c1-4*c2*b1)*L1),
                   2.*(0.5*(-c1**2*b2+c3*c2*b1+3*c2*b1*c3*mu3-c1*c2*b1+c1*b2*c3+3*c1*c3*b2*mu3-3*b2*c1**2*mu1-3*c2*b1*c1*mu1)*L3**2+(0.5*(-6*c1*c3*b2*mu3+6*b2*c1**2*mu1+2*c1*c2*b1+6*c1*c2*b2+2*c3*b3*c2-2*c3*b2*c2-2*c1*c2*b3+4*c1*b2*c3+6*c3*c2*b1-6*c2*b2*c1*mu1+2*c1**2*b2-6*c3*b2*c1*mu1+2*c3**2*b2-6*c2*b1*c3*mu3-6*c2*b3*c1*mu1+6*c2*b3*c3*mu3-6*c2*c3*b2*mu3+6*c3**2*b2*mu3+6*c2*b1*c1*mu1)*L1+0.5*(-2*c3*c2*b1+6*c1*c2*b1+6*c1*c3*b3*mu3+2*c1*c3*b3+4*c1*c3*b1-2*c1*b2*c3+2*c3**2*b1-6*c3*b1*c1*mu1-6*b2*c1**2*mu1-6*c2*b1*c1*mu1-6*b3*c1**2*mu1-6*c1*c3*b2*mu3+6*c3**2*b1*mu3+6*c1**2*b2-6*c1*c3*b1*mu3+2*b1*c1**2+6*b1*c1**2*mu1-6*c2*b1*c3*mu3-2*c1**2*b3)*L2+2*c1*c2*b2)*L3+0.5*(3*c1**2*b3-c1*c3*b3+3*c1*c3*b1-3*c3*b1*c1*mu1-3*c3**2*b1*mu3-3*c1*c3*b3*mu3-3*b3*c1**2*mu1-c3**2*b1)*L2**2+(0.5*(6*c1*c2*b3+6*c1*b2*c3-6*c3*b3*c1*mu1+6*c3**2*b1-2*c3*b3*c2+6*c3*b1*c1*mu1+2*c1*c3*b1-6*c3**2*b2*mu3+2*c1**2*b3-6*c1*c3*b3*mu3-6*c3*b2*c1*mu1-2*c3**2*b2+2*c3**2*b3+6*b3*c1**2*mu1+6*c3**2*b3*mu3+4*c1*c3*b3-6*c2*b3*c1*mu1-6*c2*b3*c3*mu3-6*c3**2*b1*mu3)*L1+2*c1*b2*c3+2*c1*c3*b1+2*c1*c2*b3)*L2-c1*b2*c3-c3*c2*b1+0.5*(c1*c2*b3+c1*b2*c3+3*c3*b2*c1*mu1+3*c2*b3*c1*mu1-3*c3**2*b2*mu3+3*c3**2*b2+3*c3*b3*c2-3*c2*b3*c3*mu3)*L1**2+0.5*(4*c1*b2*c3+4*c3*c2*b1)*L1),
                   2.*(-0.5*(c1*b2*b3-c2*b1**2+b1*c2*b3+3*c1*b2*b3*mu3-3*b1*c1*b2*mu1-3*c2*b1**2*mu1+3*c2*b1*b3*mu3-c1*b1*b2)*L3**2+(-0.5*(2*c1*b1*b2-6*c1*b2*b3*mu3+6*c1*b2*b3-2*b2*c2*b3+2*b2*b3*c3+4*b1*c2*b3+6*b1*c2*b2-6*c2*b1*b3*mu3-6*c2*b2*b3*mu3+6*c2*b3**2*mu3-2*b1*b2*c3+6*b1*c1*b2*mu1+2*b3**2*c2-6*c3*b2*b1*mu1-6*c2*b3*b1*mu1+2*c2*b1**2+6*b2*b3*c3*mu3-6*c2*b2*b1*mu1+6*c2*b1**2*mu1)*L1-0.5*(-2*b1*c2*b3+6*c2*b1**2-6*c1*b2*b3*mu3+6*c1*b1**2*mu1+6*b1*c3*b3*mu3-6*c2*b1**2*mu1-6*c2*b1*b3*mu3+2*b1*c3*b3+6*c1*b1*b2-2*c1*b2*b3-6*c3*b1**2*mu1+6*c1*b3**2*mu3+2*c1*b1**2+4*c1*b1*b3-6*c1*b1*b3*mu3+2*c1*b3**2-6*b1*c1*b2*mu1-6*b1*c1*b3*mu1-2*c3*b1**2)*L2-2*b1*c2*b2)*L3-0.5*(-3*c3*b1**2*mu1-3*c1*b3**2*mu3-c1*b3**2+3*c1*b1*b3-3*b1*c3*b3*mu3-3*b1*c1*b3*mu1+3*c3*b1**2-b1*c3*b3)*L2**2+(-0.5*(2*c3*b3**2+6*c3*b3**2*mu3-6*b2*b3*c3*mu3+6*b1*c2*b3+6*c3*b1**2*mu1-2*b3**2*c2+2*c1*b1*b3-6*c2*b3*b1*mu1-6*c3*b3*b1*mu1-6*c1*b3**2*mu3+6*b1*b2*c3-6*c3*b2*b1*mu1+4*b1*c3*b3+6*c1*b3**2-6*b1*c3*b3*mu3+2*c3*b1**2-2*b2*b3*c3+6*b1*c1*b3*mu1-6*c2*b3**2*mu3)*L1-2*c1*b1*b3-2*b1*b2*c3-2*b1*c2*b3)*L2+c1*b2*b3+b1*c2*b3-0.5*(b1*b2*c3-3*b2*b3*c3*mu3+3*b2*b3*c3+b1*c2*b3+3*b3**2*c2-3*c2*b3**2*mu3+3*c2*b3*b1*mu1+3*c3*b2*b1*mu1)*L1**2-0.5*(4*c1*b2*b3+4*b1*c2*b3)*L1),
                   2.*(-(-3*c2*b1*mu1+3*c2*b1*mu2-4*c2*b1-4*b2*c1-3*c1*b2*mu1+3*c1*b2*mu2)*L3**2+(-(-6*c3*b2*mu1+4*b2*c1-8*b2*c3+4*c2*b1+6*c3*b2*mu2+6*c1*b2*mu1-6*c2*b2*mu1+6*c2*b3*mu2+6*c2*b1*mu1+6*c2*b1*mu2-6*c2*b3*mu1+6*c1*b2*mu2-8*c2*b3-6*c2*b2*mu2+4*c2*b2)*L1-(-6*c3*b1*mu1+6*c1*b1*mu2-6*c1*b2*mu1-6*c2*b1*mu2+4*c1*b1+6*c1*b3*mu2+6*c1*b1*mu1-6*c1*b3*mu1-8*c3*b1-8*b3*c1+4*c2*b1-6*c2*b1*mu1+4*b2*c1+6*c3*b1*mu2-6*c1*b2*mu2)*L2-4*c2*b2+4*c3*b1+4*b3*c1)*L3-(-3*c3*b1*mu2+2*b3*c1-3*c1*b3*mu1-3*c3*b1*mu1+2*c3*b1-3*c1*b3*mu2)*L2**2+(-(4*c3*b1+4*b3*c1+4*c2*b3+6*c3*b3*mu2-8*c3*b3+6*c3*b1*mu1-6*c2*b3*mu2-6*c2*b3*mu1+6*c1*b3*mu1-6*c3*b3*mu1+6*c1*b3*mu2-6*c3*b2*mu1+4*b2*c3+6*c3*b1*mu2-6*c3*b2*mu2)*L1-4*c2*b3-4*b2*c3)*L2+c2*b3-b3*c1-c3*b1+b2*c3-(2*c2*b3+3*c3*b2*mu1+3*c3*b2*mu2+2*b2*c3+3*c2*b3*mu2+3*c2*b3*mu1)*L1**2+4*c3*b3*L1),
                   2.*(0.5*(-c1**2*b2-3*c2*b1*c1*mu1-3*c1*c2*b2*mu2+3*c1*c2*b2-c1*c2*b1-3*c2**2*b1*mu2+3*c2**2*b1-3*b2*c1**2*mu1)*L3**2+(0.5*(2*c1*c2*b1+4*c1*c2*b2-6*c3*b2*c2*mu2+6*b2*c1**2*mu1-2*c2**2*b1+2*c1**2*b2+6*c2**2*b3+6*c3*b2*c2-6*c2*b3*c1*mu1-2*c1*b2*c3-6*c2**2*b1*mu2-6*c2*b2*c1*mu1+2*c2**2*b2+6*c2*b1*c1*mu1-6*c3*b2*c1*mu1-2*c1*c2*b3-6*c2**2*b3*mu2+6*c2**2*b2*mu2-6*c1*c2*b2*mu2)*L1+0.5*(4*c1*c2*b1+6*c1*c2*b3+6*c2**2*b1*mu2-6*c3*b1*c2*mu2+6*c1**2*b2-6*c3*b1*c1*mu1+6*c3*c2*b1-2*c1*c3*b1-2*c1**2*b3-6*b3*c1**2*mu1+6*b1*c1**2*mu1-6*b2*c1**2*mu1-6*c1*c2*b3*mu2+6*c1*c2*b2*mu2-6*c2*b1*c1*mu1+2*c1*c2*b2-6*c1*c2*b1*mu2+2*c2**2*b1+2*b1*c1**2)*L2+2*c1*c2*b3+2*c1*c2*b2+2*c3*c2*b1)*L3+0.5*(-3*c3*b1*c1*mu1+3*c1*c2*b3*mu2+3*c1**2*b3+c3*c2*b1+3*c1*c3*b1+c1*c2*b3-3*b3*c1**2*mu1+3*c3*b1*c2*mu2)*L2**2+(0.5*(-2*c3*c2*b1+4*c1*c2*b3-2*c1*c3*b3+2*c1**2*b3+6*c3*b3*c2+6*c1*b2*c3-6*c2*b3*c1*mu1+2*c1*c3*b1+6*c3*b1*c1*mu1-6*c3*b3*c1*mu1-6*c3*b2*c1*mu1-6*c3*b1*c2*mu2+6*c3*b2*c2*mu2+2*c3*b2*c2+2*c2**2*b3+6*b3*c1**2*mu1+6*c2**2*b3*mu2-6*c1*c2*b3*mu2-6*c3*c2*b3*mu2)*L1+2*c1*b2*c3+2*c1*c2*b3)*L2-c1*b2*c3-c1*c2*b3+0.5*(c1*c2*b3+3*c2*b3*c1*mu1-c3*b2*c2+3*c3*b2*c1*mu1-3*c2**2*b3*mu2-c2**2*b3-3*c3*b2*c2*mu2+c1*b2*c3)*L1**2+2*c3*b3*c2*L1),
                   2.*(-0.5*(-3*b1*c2*b2*mu2-c2*b1**2+3*c1*b2**2-3*c1*b2**2*mu2-c1*b1*b2-3*c2*b1**2*mu1-3*b1*c1*b2*mu1+3*b1*c2*b2)*L3**2+(-0.5*(-6*c1*b2**2*mu2+2*c2*b1**2+2*c2*b2**2-2*c1*b2**2-2*b1*b2*c3+6*c3*b2**2+6*c2*b2**2*mu2+4*b1*c2*b2+6*c2*b1**2*mu1-6*c2*b3*b1*mu1-6*c2*b2*b1*mu1+6*b1*c1*b2*mu1-6*c3*b2**2*mu2+2*c1*b1*b2-6*b1*c2*b2*mu2+6*b2*c2*b3-6*b2*c2*b3*mu2-6*c3*b2*b1*mu1-2*b1*c2*b3)*L1-0.5*(2*c1*b2**2+6*c2*b1**2+2*b1*c2*b2-6*c1*b2*b3*mu2-6*c3*b1**2*mu1+4*c1*b1*b2+6*c1*b1**2*mu1-6*c2*b1**2*mu1-6*c3*b1*b2*mu2+6*b1*b2*c3-6*b1*c1*b2*mu1+6*c1*b2*b3+2*c1*b1**2-6*b1*c1*b3*mu1-2*c1*b1*b3-2*c3*b1**2-6*c1*b1*b2*mu2+6*b1*c2*b2*mu2+6*c1*b2**2*mu2)*L2-2*b1*c2*b2-2*c1*b2*b3-2*b1*b2*c3)*L3-0.5*(3*c1*b1*b3+3*c3*b1*b2*mu2-3*b1*c1*b3*mu1+3*c1*b2*b3*mu2-3*c3*b1**2*mu1+b1*b2*c3+3*c3*b1**2+c1*b2*b3)*L2**2+(-0.5*(2*c3*b1**2+6*b1*c1*b3*mu1+6*c3*b1**2*mu1-6*c1*b2*b3*mu2-2*c1*b2*b3-6*c2*b3*b1*mu1+2*c3*b2**2+6*b1*c2*b3+6*b2*c2*b3*mu2-2*b1*c3*b3+2*c1*b1*b3+6*b2*b3*c3+2*b2*c2*b3-6*c3*b1*b2*mu2-6*c3*b2*b3*mu2-6*c3*b3*b1*mu1+4*b1*b2*c3-6*c3*b2*b1*mu1+6*c3*b2**2*mu2)*L1-2*b1*b2*c3-2*b1*c2*b3)*L2+b1*b2*c3+b1*c2*b3-0.5*(b1*c2*b3-c3*b2**2+3*c2*b3*b1*mu1-3*c3*b2**2*mu2-3*b2*c2*b3*mu2+b1*b2*c3+3*c3*b2*b1*mu1-b2*c2*b3)*L1**2-2*b2*b3*c3*L1)]])
        
        return B, 1, 0 
    def FormT(self, L1, L2, L3):
        L3=1-L1-L2
        T = array([L1, 0, 0, L2, 0, 0, L3, 0, 0])
        return T
    def FormX(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([L1, L2, L3])
        return X
    def FormX_(self, L1, L2, L3):
        L3=1-L1-L2
        X = array([[L1, 0, L2, 0, L3, 0],[0,L1,0,L2,0,L3]])
        return X
    def JacoD(self, r, s, t):
        return 1

    def ShearForces(self, NodeList,NoIndToCMInd, WrNoVal ):
        AI = array([[53./12.,-11./6.,-11./6.],[-11./6.,53./12.,-11./6.],[-11./6.,-11./6.,53./12.]])
        XX = array([[1./3.,1./3.,1./3.],[0.2,0.6,0.2],[0.2,0.2,0.6],[0.6,0.2,0.2]])
        for j in range(len(self.Inzi)): NodeList[NoIndToCMInd[self.Inzi[j]]].c += 1 #count assigned elements
        mx = [self.Data[0,3],self.Data[1,3],self.Data[2,3],self.Data[3,3]]
        aa = dot(AI, dot(transpose(XX),mx))                             # coefficients of linear approximation
        if WrNoVal: 
            for j in range(len(self.Inzi)): NodeList[NoIndToCMInd[self.Inzi[j]]].mx += aa[j]#extrapolation to nodes
        mx_x = aa[0]*self.b1+aa[1]*self.b2+aa[2]*self.b3                # x derivative
#            mx_y = aa[0]*self.c1+aa[1]*self.c2+aa[2]*self.c3                # y derivative
        #
        my = [self.Data[0,4],self.Data[1,4],self.Data[2,4],self.Data[3,4]]
        aa = dot(AI, dot(transpose(XX),my))
        if WrNoVal: 
            for j in range(len(self.Inzi)): NodeList[NoIndToCMInd[self.Inzi[j]]].my += aa[j]#extrapolation to nodes
#            my_x = aa[0]*self.b1+aa[1]*self.b2+aa[2]*self.b3
        my_y = aa[0]*self.c1+aa[1]*self.c2+aa[2]*self.c3
        #
        mxy = [self.Data[0,5],self.Data[1,5],self.Data[2,5],self.Data[3,5]]
        aa = dot(AI, dot(transpose(XX),mxy))
        if WrNoVal: 
            for j in range(len(self.Inzi)): NodeList[NoIndToCMInd[self.Inzi[j]]].mxy += aa[j]#extrapolation to nodes
        mxy_x = aa[0]*self.b1+aa[1]*self.b2+aa[2]*self.b3
        mxy_y = aa[0]*self.c1+aa[1]*self.c2+aa[2]*self.c3
        #
        for j in range(self.nIntL): self.Data[j,6] = -mx_x -mxy_y       # x shear force
        for j in range(self.nIntL): self.Data[j,7] = -my_y -mxy_x       # y shear force
    
class SH4(Element):
    """ Shell element 4 nodes
    """
    def LengV(self, Vec):                               # bring Vec to length 1
        LL = sqrt(Vec[0]**2+Vec[1]**2+Vec[2]**2)
        LL = 1./LL                                      # division by zero not catched !
        Vec[0]= LL*Vec[0]
        Vec[1]= LL*Vec[1]
        Vec[2]= LL*Vec[2]
    def CompNoNor(self, i0, i1, i2):                    # returns normal to i0 -> i1 and i0- -> i2
        ax = self.XX[0,i1] - self.XX[0,i0]
        ay = self.XX[1,i1] - self.XX[1,i0]
        az = self.XX[2,i1] - self.XX[2,i0]
        bx = self.XX[0,i2] - self.XX[0,i0]
        by = self.XX[1,i2] - self.XX[1,i0]
        bz = self.XX[2,i2] - self.XX[2,i0]
        nx = ay*bz - az*by
        ny = az*bx - ax*bz
        nz = ax*by - ay*bx
        ll = sqrt(nx*nx+ny*ny+nz*nz)
        if ll<ZeroD: raise NameError("ConFemElements::SHX.__ini__::VecPro: something wrong with geometry")
        return nx/ll, ny/ll, nz/ll
    def ComputeGG(self):
        for i in range(3): 
            for j in range(4): self.gg[0,i,j] = -self.a[j]*self.V2[i,j]/2.     # 1st index direction, 2nd index node 
            for j in range(4): self.gg[1,i,j] =  self.a[j]*self.V1[i,j]/2.
    def ComputeEdgeDir(self):                           # direction of 1st edge for later use
        self.EdgeDir[0] = self.XX[0,1]-self.XX[0,0]
        self.EdgeDir[1] = self.XX[1,1]-self.XX[1,0]
        self.EdgeDir[2] = self.XX[2,1]-self.XX[2,0]
        self.LengV(self.EdgeDir)
    def CompleteTriad(self, dd):                        # complete director system with V1, V2
        Tol = 0.01
        V1 = zeros((3),dtype=float)
        V2 = zeros((3),dtype=float)
        if abs(dd[2])>Tol or abs(dd[0])>Tol: V1[0]= dd[2];              V1[2]=-dd[0] # seems to regard on direction / unit vector of global coordinate system
        else:                                V1[0]=-dd[1]; V1[1]= dd[0]
        self.LengV(V1)
        V2[0] = dd[1]*V1[2]
        V2[1] = dd[2]*V1[0]-dd[0]*V1[2]
        V2[2] =            -dd[1]*V1[0]
        return V1, V2
    def ComputeTransLocal(self, i, nn, TTW):            # 
        V1_, V2_ = self.CompleteTriad( nn)              # second approach for director triad with input data for nodes
        TT = zeros((2,3),dtype=float)                   # coordinate transformation matrix for rotations regarding current node
        TT[0,0] = dot(V1_,self.V1[:,i])
        TT[0,1] = dot(V1_,self.V2[:,i])
        TT[0,2] = dot(V1_,self.Vg[:,i])
        TT[1,0] = dot(V2_,self.V1[:,i])
        TT[1,1] = dot(V2_,self.V2[:,i])
        TT[1,2] = dot(V2_,self.Vg[:,i])
        TTW[i]  = TT.copy()
        for j in range(3): self.V1[j,i] = V1_[j]       # set rest of the director triad anew
        for j in range(3): self.V2[j,i] = V2_[j]       #
    def ComputeTransLocalAll(self, TTW):
        base, base_ = 0, 0
        for i in range(4):                             # loop over nodes of element  
            for j in range(3):          self.Trans[base_+j,base+j] = 1.
            if self.SixthDoF[i]:
                for j in range(2):
                    for jj in range(3): self.Trans[base_+3+j,base+3+jj] = TTW[i][j][jj] 
#   uhc             XX = dot(transpose(TTW[i]),TTW[i])
#   uhc             YY = dot(TTW[i],transpose(TTW[i]))
#                print i, TTW[i],'\n', XX, '\n', YY
            else:
                for j in range(2):      self.Trans[base_+3+j,base+3+j] = 1.
            base = base + self.DofN[i]
            base_= base_+ self.DofNini[i]            
#        for j in xrange(self.Trans.shape[0]):
#            for jj in xrange(self.Trans.shape[1]): sys.stdout.write('%6.2f'%(self.Trans[j,jj]))
#            sys.stdout.write('\n')
#        raise NameError ("Exit")
    def __init__(self, Label, SetLabel, InzList, MatName,Material, NoList, ShellSec, StateV, NData, RCFlag, NoLabToNoInd):
#       Element.__init__(self,"SH4",InzList, 3, 4,2,16, [set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5])], 21,True,Label,SetLabel,3, MatName,StateV,NData, NoList,NoLabToNoInd,[]) # four integration points over cross section height
        Element.__init__(self,"SH4",InzList, 3, 4,5,20, [set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5])], 21,True,Label,SetLabel,3, MatName,StateV,NData, NoList,NoLabToNoInd,[]) # five integration points over cross section height
    def Ini2(self, NoList,NoIndToCMInd, MaList, SecDict):
        self.PlSt = True                                                # flag for plane stress  - should be true, otherwise inconsistencies in e.g. MISES
        self.RotM= True                                                 # Flag for elements with local coordinate system for materials
        self.TensStiff = False                                          # flag for tension stiffening
        self.ElemUpdate = True                                          # element update might be required in case of NLGEOM
        SetLabel = self.Set
        ShellSec = SecDict[SetLabel]
        RCFlag = self.ShellRCFlag                                       # set at post-DataIn
        if self.Type == 'SH3': self.a = array([ShellSec.Height,ShellSec.Height,ShellSec.Height]) # shell thickness
        else:                  self.a = array([ShellSec.Height,ShellSec.Height,ShellSec.Height,ShellSec.Height]) # shell thickness
        nRe = len(ShellSec.Reinf)                           # number of reinforcement layers
        self.Geom = zeros( (2+nRe,5), dtype=double)
        self.Geom[0,0] = 1                                  # dummy for Area / Jacobi determinant used instead
        self.Geom[1,0] = 1                                  # dummy for height / thickness
        if RCFlag and nRe>0:
            if self.Type == 'SH3':
                raise NameError("ConFemElem::Ini2: reinforcement not yet tested for SH3")
            if   self.nInt==2: i1 = 1; i2 = 16              # nInt integration order indicator # will presumably not work for sh3
            elif self.nInt==5: i1 = 4; i2 = 20
            else: raise NameError("ConFemElements::SH4.__ini__: integration order not implemented for RC",self.nInt)
            import ConFemBasics                             # to modify sampling points for numerical integration
            for j in range(nRe):
                self.Geom[2+j,0] = ShellSec.Reinf[j][0]      # reinforcement cross section
                self.Geom[2+j,1] = ShellSec.Reinf[j][1]      # " lever arm
                self.Geom[2+j,2] = ShellSec.Reinf[j][2]      # " effective reinforcement ratio for tension stiffening
                self.Geom[2+j,3] = ShellSec.Reinf[j][3]      # " parameter 0<beta<=1 for tension stiffening, beta=0 no tension stiffening
                self.Geom[2+j,4] = ShellSec.Reinf[j][4]      # " direction
                tt = 2.*ShellSec.Reinf[j][1]/ShellSec.Height # local t-coordinate for reinforcement
                # sampling points
                ConFemBasics.SamplePointsRCShell[SetLabel,4,i1,i2+4*j+0] =[-0.577350269189626,-0.577350269189626, tt] # every reinforcement layer / j gets consecutive indices in base plane
                ConFemBasics.SamplePointsRCShell[SetLabel,4,i1,i2+4*j+1] =[-0.577350269189626, 0.577350269189626, tt] #
                ConFemBasics.SamplePointsRCShell[SetLabel,4,i1,i2+4*j+2] =[ 0.577350269189626,-0.577350269189626, tt] #
                ConFemBasics.SamplePointsRCShell[SetLabel,4,i1,i2+4*j+3] =[ 0.577350269189626, 0.577350269189626, tt] #
                # sampling weights
                ConFemBasics.SampleWeightRCShell[SetLabel,4,i1,i2+4*j+0]= 2.*ShellSec.Reinf[j][0]/ShellSec.Height
                ConFemBasics.SampleWeightRCShell[SetLabel,4,i1,i2+4*j+1]= 2.*ShellSec.Reinf[j][0]/ShellSec.Height
                ConFemBasics.SampleWeightRCShell[SetLabel,4,i1,i2+4*j+2]= 2.*ShellSec.Reinf[j][0]/ShellSec.Height
                ConFemBasics.SampleWeightRCShell[SetLabel,4,i1,i2+4*j+3]= 2.*ShellSec.Reinf[j][0]/ShellSec.Height
        return []
            
    def Ini3(self, NoList, NoIndToCMInd):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        i3 = NoIndToCMInd[self.Inzi[3]]
        ni = [NoList[i0], NoList[i1],NoList[i2],NoList[i3]]
        self.XX =      array([[ni[0].XCo, ni[1].XCo, ni[2].XCo, ni[3].XCo],
                              [ni[0].YCo, ni[1].YCo, ni[2].YCo, ni[3].YCo],
                              [ni[0].ZCo, ni[1].ZCo, ni[2].ZCo, ni[3].ZCo]])# collect nodal coordinates in a compact form for later use
        self.EdgeDir = zeros((3),dtype=float)                               # initialize direction of 1st edge for later use
        self.ComputeEdgeDir()                                               # direction of 1st edge for later use
        nn = zeros((4,3), dtype = float)
        nn[0,:] = self.CompNoNor( 0, 1, 3)                                  # roughly unit normals to shell surface at nodes
        nn[1,:] = self.CompNoNor( 1, 2, 0)
        nn[2,:] = self.CompNoNor( 2, 3, 1)
        nn[3,:] = self.CompNoNor( 3, 0, 2)
        self.V1 = zeros((3,4),dtype=float)                                  # initialize director triad
        self.V2 = zeros((3,4),dtype=float)
        self.Vn = zeros((3,4),dtype=float)
        self.Vg = zeros((3,4),dtype=float)                                  # director as defined per node via input data, might not be the actually used director as is ruled in the following
        self.VD = zeros((3,4),dtype=float)                                  # director increment in time increment
        TTW = [None, None, None, None]                                      # for temporal storage of transformation coefficients
        self.SixthDoF = [False, False, False, False]
        for i in range(4):                                                  # loop over nodes of element  
            LL = sqrt(ni[i].XDi**2+ni[i].YDi**2+ni[i].ZDi**2)               # length of directors
            self.Vg[0,i] = ni[i].XDi/LL                                     # values given with input data
            self.Vg[1,i] = ni[i].YDi/LL                                     # "
            self.Vg[2,i] = ni[i].ZDi/LL                                     # "
            self.V1[:,i], self.V2[:,i] = self.CompleteTriad(self.Vg[:,i])
            if dot(nn[i],self.Vg[:,i])<0.8:
                print("ConFemElements::SH4.__ini__::ControlGeom: unfavorable shell director element %s local node index %s"%(str(self.Label),str(i)))
                self.SixthDoF[i] = True
                ni[i].DofT = ni[i].DofT.union(set([6]))                     # extend types of dof for this node
                self.DofT[i] = self.DofT[i].union(set([6]))                 # extend types of dof for this element
                self.DofE  = self.DofE + 1                                   # adapt number of dofs of whole element
                DofN_ = list(self.DofN)                                     # transform tuple into list to make it changeable
                DofN_[i] = 6                                                # one more degree of freedom for local node i
                self.DofN = tuple(DofN_)                                    # transform back into tuple
                self.Rot = True                                             # Element has formally to be rotated as a whole
                self.RotG= True                                             # geometric stiffness has also to be rotated for nonlinear geometry
                self.Trans = zeros((self.DofEini, self.DofE), dtype=float)  # initialize rotation / transformation matrix for element / coordinate transformation matrix, not quadratic anymore!
                self.ComputeTransLocal(i, nn[i], TTW)                       # transformation matrix for modified director system
                self.Vn[0,i] = nn[i,0]                                      # final director
                self.Vn[1,i] = nn[i,1]
                self.Vn[2,i] = nn[i,2]
            else:
                self.Vn[0,i] = self.Vg[0,i]                                 # final director
                self.Vn[1,i] = self.Vg[1,i]
                self.Vn[2,i] = self.Vg[2,i]
        if self.Rot: self.ComputeTransLocalAll(TTW)                         # build rotation / transformation matrix for element / coordinate transformation matrix, not quadratic anymore!
        self.XX0 = self.XX.copy() #copy(self.XX)                                            # retain initial values for NLGEOM
        self.DofI = array([[-1,-1,-1,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,-1,-1,-1,-1]],dtype=int)
        self.gg = zeros((2,3,4), dtype=float)
        self.ComputeGG()                                                    # scaled axes for rotational degrees of freedom
        #
        AA = 0.
        for i in range(16):                                                 # volume by numerical integration for characteristic length  
            r = SamplePoints[4,1,i][0]
            s = SamplePoints[4,1,i][1]
            f = self.JacoD(r,s,0)*SampleWeight[4,1,i]
            AA = AA + f
        self.Lch_ = sqrt(AA/(0.25*(self.a[0]+self.a[1]+self.a[2]+self.a[3]))) #/self.nInt   # characteristic length -> side length of square of same area of shell area
        self.Lch  = self.Lch_/self.nInt                                     # for ELASTICLT
        self.CrBScaleType()
    def Lists1(self):                                                       # indices for first integration points in base area
        if self.nInt==1:                                                    # integration order   
            Lis = [0,1, 2,3]
        elif self.nInt==2:                                                  # see ConFemBasics SamplePoints for shells + 1
            Lis = [0,4, 8,12]                                               # indices for first integration points in base area, 4 Gaussian integration points over cross section height
        elif self.nInt>=5:                                                  # # see ConFemBasics SamplePoints for shells + 1
            Lis = [0,5,10,15]                                               # "                                                  5 "
        return Lis
    def Lists2(self, nRe, j):                                               # integration point indices over height relative to base point integration index j 
        if self.nInt==1:
            Lis2 = [0]
        elif self.nInt==2: 
            Lis2, Lis3 = [0,1,2,3], []                                      # local counter for loop over cross section height
                                                                            # indices of reinforcement layers in the actual base plane integration point, needs to be filled 
            j_ = j//4                                                       # floor division -> integer value, consecutive counter for point in base plane
            Offs2 = 16
        elif self.nInt>=5:                                                  # currently only intorder 5 relevant 
            Lis2, Lis3 = [0,1,2,3,4], []                                    # "
            j_ = j//5                                                       # "
            Offs2 = 20
        for k in range(nRe): Lis3.append(Offs2+k*4+j_-j)                    # indices for reinforcement layers only, 4 should be for number of base plane integration points
                                                                            # -j compensates for adding j in next loop
        Lis2.extend(Lis3)                                                   # appends indices for reinforcement layers to local height counter 
        return Lis2, Lis3                                                   # Lis2 for bulk height integration points + reinforcement integration points, Lis3 for reinforcement integration points ONLY
    
    def FormX(self, r, s, t):                               # interpolation on geometry
        X = array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25])
        return X
#    def FormX_(self, r, s, t):                                              # for geometry interpolation --- different scheme compared to NN ??? 
#        X = array([[(1+r)*(1+s)*0.25, 0, 0, (1-r)*(1+s)*0.25, 0, 0, (1-r)*(1-s)*0.25, 0, 0, (1+r)*(1-s)*0.25, 0, 0],
#                   [0, (1+r)*(1+s)*0.25, 0, 0, (1-r)*(1+s)*0.25, 0, 0, (1-r)*(1-s)*0.25, 0, 0, (1+r)*(1-s)*0.25, 0],
#                   [0, 0, (1+r)*(1+s)*0.25, 0, 0, (1-r)*(1+s)*0.25, 0, 0, (1-r)*(1-s)*0.25, 0, 0, (1+r)*(1-s)*0.25]])
        return X
    def FormX_(self, r, s, t):                                              # for geometry interpolation 
        X = array([[(1-r)*(1-s)*0.25, 0, 0,    (1+r)*(1-s)*0.25, 0, 0,    (1+r)*(1+s)*0.25, 0, 0,    (1-r)*(1+s)*0.25, 0, 0],
                   [0, (1-r)*(1-s)*0.25, 0, 0,   (1+r)*(1-s)*0.25, 0, 0,    (1+r)*(1+s)*0.25, 0, 0,    (1-r)*(1+s)*0.25, 0],
                   [0, 0, (1-r)*(1-s)*0.25, 0, 0,  (1+r)*(1-s)*0.25, 0, 0,    (1+r)*(1+s)*0.25, 0, 0,    (1-r)*(1+s)*0.25]])
        return X
    def FormN(self, r, s, t):
        N = array([[(1-r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1+s)*0.25, 0,0,0,0, (1-r)*(1+s)*0.25, 0,0,0,0],
                   [0,(1-r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1+s)*0.25, 0,0,0,0, (1-r)*(1+s)*0.25, 0,0,0],
                   [0,0,(1-r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1-s)*0.25, 0,0,0,0, (1+r)*(1+s)*0.25, 0,0,0,0, (1-r)*(1+s)*0.25, 0,0]])
        return N
    def Basics(self, r, s, t):
        N =  array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25])
        br = array([(-1+s)*0.25,      ( 1-s)*0.25,     ( 1+s)*0.25,      -( 1+s)*0.25])
        bs = array([(-1+r)*0.25,     -( 1+r)*0.25,     ( 1+r)*0.25,       ( 1-r)*0.25])
        JJ = zeros((3,3),dtype=float)
        for k in range(3): JJ[k,0]=br[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+br[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+br[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+br[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,1]=bs[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+bs[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+bs[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+bs[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,2]=N[0]*(1/2.*self.a[0]*self.Vn[k,0])+N[1]*(1/2.*self.a[1]*self.Vn[k,1])+N[2]*(1/2.*self.a[2]*self.Vn[k,2])+N[3]*(1/2.*self.a[3]*self.Vn[k,3])
        JI = inv(JJ)
        ll=sqrt(JJ[0,2]**2+JJ[1,2]**2+JJ[2,2]**2)   
        vv = array([[0.,0.,JJ[0,2]/ll],[0.,0.,JJ[1,2]/ll],[0.,0.,JJ[2,2]/ll]]) # unit normal of local coordinate system, 3rd column
        LoC = False                                                        
        if LoC:                                             # local right handed coordinate system with 1st direction / column aligned to element edge
            x0 = self.EdgeDir[1]*vv[2,2]-self.EdgeDir[2]*vv[1,2]
            x1 = self.EdgeDir[2]*vv[0,2]-self.EdgeDir[0]*vv[2,2]
            x2 = self.EdgeDir[0]*vv[1,2]-self.EdgeDir[1]*vv[0,2]
            xx = sqrt(x0**2+x1**2+x2**2)
            vv[0,1] = -x0/xx                                # 2nd column, approx perp. to element edge, reversed in sign to preserve right handedness
            vv[1,1] = -x1/xx
            vv[2,1] = -x2/xx 
            x0 = vv[1,2]*vv[2,1]-vv[2,2]*vv[1,1]
            x1 = vv[2,2]*vv[0,1]-vv[0,2]*vv[2,1]
            x2 = vv[0,2]*vv[1,1]-vv[1,2]*vv[0,1]
            xx = sqrt(x0**2+x1**2+x2**2)
            vv[0,0] = -x0/xx                                # 1st column, approx aligned to element edge, sign reversal of 2nd column is implicitely corrected
            vv[1,0] = -x1/xx
            vv[2,0] = -x2/xx 
        else:                                               # local coordinate system with one axis aligned to global axis
            if abs(vv[1,2])<0.99:                           # local coordinate system V_1 from cross product of V_n and e_y ( V_1 in e_x - e_z plane) if V_n is not to much aligned to e_y  
                ll = sqrt(vv[2,2]**2+vv[0,2]**2)            # length of V_1
                vv[0,0] = vv[2,2]/ll                        # V_1[0] normalized;  V1[1] = 0
                vv[2,0] =-vv[0,2]/ll                        # V_1[2] normalized

                vv[0,1] = vv[1,2]*vv[2,0]                   # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
                vv[1,1] = vv[2,2]*vv[0,0]-vv[0,2]*vv[2,0]
                vv[2,1] =-vv[1,2]*vv[0,0]
            else:                                           # local coordinate system V_1 from cross product of V_n and e_x ( V_1 in e_y - e_z plane)
                ll = sqrt(vv[2,2]**2+vv[1,2]**2)            # length of V_1
                vv[1,0] =-vv[2,2]/ll                        # V_1[0] normalized;  V1[0] = 0
                vv[2,0] = vv[1,2]/ll                        # V_1[2] normalized

                vv[0,1] = vv[1,2]*vv[2,0]-vv[2,2]*vv[1,0]   # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
                vv[1,1] =-vv[0,2]*vv[2,0]
                vv[2,1] = vv[0,2]*vv[1,0]
        return N, br, bs, JJ, JI, vv
    def FormB(self, r, s, t_, NLg):
        t = t_
        N, br, bs, JJ, JI, vv = self.Basics( r, s, t)
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        HH = zeros((3,2,4),dtype=float)
        for k in range(4):
            HH[0,0,k]=JJ[0,0]*self.gg[0,0,k]+JJ[1,0]*self.gg[0,1,k]+JJ[2,0]*self.gg[0,2,k]
            HH[0,1,k]=JJ[0,0]*self.gg[1,0,k]+JJ[1,0]*self.gg[1,1,k]+JJ[2,0]*self.gg[1,2,k]
            HH[1,0,k]=JJ[0,1]*self.gg[0,0,k]+JJ[1,1]*self.gg[0,1,k]+JJ[2,1]*self.gg[0,2,k]
            HH[1,1,k]=JJ[0,1]*self.gg[1,0,k]+JJ[1,1]*self.gg[1,1,k]+JJ[2,1]*self.gg[1,2,k]
            HH[2,0,k]=JJ[0,2]*self.gg[0,0,k]+JJ[1,2]*self.gg[0,1,k]+JJ[2,2]*self.gg[0,2,k]
            HH[2,1,k]=JJ[0,2]*self.gg[1,0,k]+JJ[1,2]*self.gg[1,1,k]+JJ[2,2]*self.gg[1,2,k]
        BB = zeros((6,20),dtype=float)
        for k in range(4):
            BB[0,k*5+0]=JJ[0,0]*br[k]
            BB[0,k*5+1]=JJ[1,0]*br[k]
            BB[0,k*5+2]=JJ[2,0]*br[k]
            BB[0,k*5+3]=               t*br[k]*HH[0,0,k]
            BB[0,k*5+4]=               t*br[k]*HH[0,1,k]
            BB[1,k*5+0]=JJ[0,1]*bs[k]
            BB[1,k*5+1]=JJ[1,1]*bs[k]
            BB[1,k*5+2]=JJ[2,1]*bs[k]
            BB[1,k*5+3]=               t*bs[k]*HH[1,0,k]
            BB[1,k*5+4]=               t*bs[k]*HH[1,1,k]
            BB[2,k*5+3]=N[k]*HH[2,0,k]
            BB[2,k*5+4]=N[k]*HH[2,1,k]
            BB[5,k*5+0]=JJ[0,0]*bs[k]  +JJ[0,1]*br[k]
            BB[5,k*5+1]=JJ[1,0]*bs[k]  +JJ[1,1]*br[k]
            BB[5,k*5+2]=JJ[2,0]*bs[k]  +JJ[2,1]*br[k]
            BB[5,k*5+3]=     t*(bs[k]*HH[0,0,k]+br[k]*HH[1,0,k])
            BB[5,k*5+4]=     t*(bs[k]*HH[0,1,k]+br[k]*HH[1,1,k])
        flag = True    #  True --> Assumed-Natural-Strain-Method to avoid transverse shear locking
        if not flag:
            for k in range(4):
                BB[3,(k-0)*5+0]=JJ[0,2]*bs[k]
                BB[3,(k-0)*5+1]=JJ[1,2]*bs[k]
                BB[3,(k-0)*5+2]=JJ[2,2]*bs[k]
                BB[3,(k-0)*5+3]=N[k]*HH[1,0,k]+t*bs[k]*HH[2,0,k]
                BB[3,(k-0)*5+4]=N[k]*HH[1,1,k]+t*bs[k]*HH[2,1,k]
                BB[4,(k-0)*5+0]=JJ[0,2]*br[k]
                BB[4,(k-0)*5+1]=JJ[1,2]*br[k]
                BB[4,(k-0)*5+2]=JJ[2,2]*br[k]
                BB[4,(k-0)*5+3]=N[k]*HH[0,0,k]+t*br[k]*HH[2,0,k]
                BB[4,(k-0)*5+4]=N[k]*HH[0,1,k]+t*br[k]*HH[2,1,k]
        else:
#                t = 0
                J01D=-0.5*(self.XX[0,1]-self.XX[0,2])-0.25*t*(self.a[1]*self.Vn[0,1]-self.a[2]*self.Vn[0,2])
                J11D=-0.5*(self.XX[1,1]-self.XX[1,2])-0.25*t*(self.a[1]*self.Vn[1,1]-self.a[2]*self.Vn[1,2])
                J21D=-0.5*(self.XX[2,1]-self.XX[2,2])-0.25*t*(self.a[1]*self.Vn[2,1]-self.a[2]*self.Vn[2,2])
                J01B=-0.5*(self.XX[0,0]-self.XX[0,3])-0.25*t*(self.a[0]*self.Vn[0,0]-self.a[3]*self.Vn[0,3])
                J11B=-0.5*(self.XX[1,0]-self.XX[1,3])-0.25*t*(self.a[0]*self.Vn[1,0]-self.a[3]*self.Vn[1,3])
                J21B=-0.5*(self.XX[2,0]-self.XX[2,3])-0.25*t*(self.a[0]*self.Vn[2,0]-self.a[3]*self.Vn[2,3])
                J00A= 0.5*(self.XX[0,2]-self.XX[0,3])+0.25*t*(self.a[2]*self.Vn[0,2]-self.a[3]*self.Vn[0,3])
                J10A= 0.5*(self.XX[1,2]-self.XX[1,3])+0.25*t*(self.a[2]*self.Vn[1,2]-self.a[3]*self.Vn[1,3])
                J20A= 0.5*(self.XX[2,2]-self.XX[2,3])+0.25*t*(self.a[2]*self.Vn[2,2]-self.a[3]*self.Vn[2,3])
                J00C=-0.5*(self.XX[0,0]-self.XX[0,1])-0.25*t*(self.a[0]*self.Vn[0,0]-self.a[1]*self.Vn[0,1])
                J10C=-0.5*(self.XX[1,0]-self.XX[1,1])-0.25*t*(self.a[0]*self.Vn[1,0]-self.a[1]*self.Vn[1,1])
                J20C=-0.5*(self.XX[2,0]-self.XX[2,1])-0.25*t*(self.a[0]*self.Vn[2,0]-self.a[1]*self.Vn[2,1])
                J02D=0.25*(self.a[1]*self.Vn[0,1]+self.a[2]*self.Vn[0,2])
                J12D=0.25*(self.a[1]*self.Vn[1,1]+self.a[2]*self.Vn[1,2])
                J22D=0.25*(self.a[1]*self.Vn[2,1]+self.a[2]*self.Vn[2,2])
                J02B=0.25*(self.a[0]*self.Vn[0,0]+self.a[3]*self.Vn[0,3])
                J12B=0.25*(self.a[0]*self.Vn[1,0]+self.a[3]*self.Vn[1,3])
                J22B=0.25*(self.a[0]*self.Vn[2,0]+self.a[3]*self.Vn[2,3])
                J02A=0.25*(self.a[2]*self.Vn[0,2]+self.a[3]*self.Vn[0,3])
                J12A=0.25*(self.a[2]*self.Vn[1,2]+self.a[3]*self.Vn[1,3])
                J22A=0.25*(self.a[2]*self.Vn[2,2]+self.a[3]*self.Vn[2,3])
                J02C=0.25*(self.a[0]*self.Vn[0,0]+self.a[1]*self.Vn[0,1])
                J12C=0.25*(self.a[0]*self.Vn[1,0]+self.a[1]*self.Vn[1,1])
                J22C=0.25*(self.a[0]*self.Vn[2,0]+self.a[1]*self.Vn[2,1])
                BB[3,0] =-0.25*(1-r)*J02B
                BB[3,1] =-0.25*(1-r)*J12B
                BB[3,2] =-0.25*(1-r)*J22B
                BB[3,3] = 0.25*(1-r)*((J01B*self.gg[0,0,0]+J11B*self.gg[0,1,0]+J21B*self.gg[0,2,0])-t*(J02B*self.gg[0,0,0]+J12B*self.gg[0,1,0]+J22B*self.gg[0,2,0]))
                BB[3,4] = 0.25*(1-r)*((J01B*self.gg[1,0,0]+J11B*self.gg[1,1,0]+J21B*self.gg[1,2,0])-t*(J02B*self.gg[1,0,0]+J12B*self.gg[1,1,0]+J22B*self.gg[1,2,0]))
                BB[3,5] =-0.25*(1+r)*J02D 
                BB[3,6] =-0.25*(1+r)*J12D
                BB[3,7] =-0.25*(1+r)*J22D
                BB[3,8] = 0.25*(1+r)*((J01D*self.gg[0,0,1]+J11D*self.gg[0,1,1]+J21D*self.gg[0,2,1])-t*(J02D*self.gg[0,0,1]+J12D*self.gg[0,1,1]+J22D*self.gg[0,2,1]))
                BB[3,9]= 0.25*(1+r)*((J01D*self.gg[1,0,1]+J11D*self.gg[1,1,1]+J21D*self.gg[1,2,1])-t*(J02D*self.gg[1,0,1]+J12D*self.gg[1,1,1]+J22D*self.gg[1,2,1]))
                BB[3,10]= 0.25*(1+r)*J02D 
                BB[3,11]= 0.25*(1+r)*J12D
                BB[3,12]= 0.25*(1+r)*J22D
                BB[3,13]= 0.25*(1+r)*((J01D*self.gg[0,0,2]+J11D*self.gg[0,1,2]+J21D*self.gg[0,2,2])+t*(J02D*self.gg[0,0,2]+J12D*self.gg[0,1,2]+J22D*self.gg[0,2,2]))
                BB[3,14]= 0.25*(1+r)*((J01D*self.gg[1,0,2]+J11D*self.gg[1,1,2]+J21D*self.gg[1,2,2])+t*(J02D*self.gg[1,0,2]+J12D*self.gg[1,1,2]+J22D*self.gg[1,2,2]))
                BB[3,15]= 0.25*(1-r)*J02B
                BB[3,16]= 0.25*(1-r)*J12B
                BB[3,17]= 0.25*(1-r)*J22B
                BB[3,18]= 0.25*(1-r)*((J01B*self.gg[0,0,3]+J11B*self.gg[0,1,3]+J21B*self.gg[0,2,3])+t*(J02B*self.gg[0,0,3]+J12B*self.gg[0,1,3]+J22B*self.gg[0,2,3]))
                BB[3,19]= 0.25*(1-r)*((J01B*self.gg[1,0,3]+J11B*self.gg[1,1,3]+J21B*self.gg[1,2,3])+t*(J02B*self.gg[1,0,3]+J12B*self.gg[1,1,3]+J22B*self.gg[1,2,3]))
                BB[4,0] =-0.25*(1-s)*J02C
                BB[4,1] =-0.25*(1-s)*J12C
                BB[4,2] =-0.25*(1-s)*J22C
                BB[4,3] = 0.25*(1-s)*((J00C*self.gg[0,0,0]+J10C*self.gg[0,1,0]+J20C*self.gg[0,2,0])-t*(J02C*self.gg[0,0,0]+J12C*self.gg[0,1,0]+J22C*self.gg[0,2,0]))
                BB[4,4] = 0.25*(1-s)*((J00C*self.gg[1,0,0]+J10C*self.gg[1,1,0]+J20C*self.gg[1,2,0])-t*(J02C*self.gg[1,0,0]+J12C*self.gg[1,1,0]+J22C*self.gg[1,2,0]))
                BB[4,5] = 0.25*(1-s)*J02C 
                BB[4,6] = 0.25*(1-s)*J12C
                BB[4,7] = 0.25*(1-s)*J22C
                BB[4,8] = 0.25*(1-s)*((J00C*self.gg[0,0,1]+J10C*self.gg[0,1,1]+J20C*self.gg[0,2,1])+t*(J02C*self.gg[0,0,1]+J12C*self.gg[0,1,1]+J22C*self.gg[0,2,1]))
                BB[4,9] = 0.25*(1-s)*((J00C*self.gg[1,0,1]+J10C*self.gg[1,1,1]+J20C*self.gg[1,2,1])+t*(J02C*self.gg[1,0,1]+J12C*self.gg[1,1,1]+J22C*self.gg[1,2,1]))
                BB[4,10]= 0.25*(1+s)*J02A 
                BB[4,11]= 0.25*(1+s)*J12A 
                BB[4,12]= 0.25*(1+s)*J22A
                BB[4,13]= 0.25*(1+s)*((J00A*self.gg[0,0,2]+J10A*self.gg[0,1,2]+J20A*self.gg[0,2,2])+t*(J02A*self.gg[0,0,2]+J12A*self.gg[0,1,2]+J22A*self.gg[0,2,2]))
                BB[4,14]= 0.25*(1+s)*((J00A*self.gg[1,0,2]+J10A*self.gg[1,1,2]+J20A*self.gg[1,2,2])+t*(J02A*self.gg[1,0,2]+J12A*self.gg[1,1,2]+J22A*self.gg[1,2,2]))
                BB[4,15]=-0.25*(1+s)*J02A 
                BB[4,16]=-0.25*(1+s)*J12A
                BB[4,17]=-0.25*(1+s)*J22A
                BB[4,18]= 0.25*(1+s)*((J00A*self.gg[0,0,3]+J10A*self.gg[0,1,3]+J20A*self.gg[0,2,3])-t*(J02A*self.gg[0,0,3]+J12A*self.gg[0,1,3]+J22A*self.gg[0,2,3]))
                BB[4,19]= 0.25*(1+s)*((J00A*self.gg[1,0,3]+J10A*self.gg[1,1,3]+J20A*self.gg[1,2,3])-t*(J02A*self.gg[1,0,3]+J12A*self.gg[1,1,3]+J22A*self.gg[1,2,3]))

#       td=array([[JI[0,0]*vv[0,0]+JI[0,1]*vv[0,1]+JI[0,2]*vv[0,2],JI[0,0]*vv[1,0]+JI[0,1]*vv[1,1]+JI[0,2]*vv[1,2],JI[0,0]*vv[2,0]+JI[0,1]*vv[2,1]+JI[0,2]*vv[2,2]],
#                 [JI[1,0]*vv[0,0]+JI[1,1]*vv[0,1]+JI[1,2]*vv[0,2],JI[1,0]*vv[1,0]+JI[1,1]*vv[1,1]+JI[1,2]*vv[1,2],JI[1,0]*vv[2,0]+JI[1,1]*vv[2,1]+JI[1,2]*vv[2,2]],
#                 [JI[2,0]*vv[0,0]+JI[2,1]*vv[0,1]+JI[2,2]*vv[0,2],JI[2,0]*vv[1,0]+JI[2,1]*vv[1,1]+JI[2,2]*vv[1,2],JI[2,0]*vv[2,0]+JI[2,1]*vv[2,1]+JI[2,2]*vv[2,2]]])
        td=array([[JI[0,0]*vv[0,0]+JI[0,1]*vv[1,0]+JI[0,2]*vv[2,0],JI[0,0]*vv[0,1]+JI[0,1]*vv[1,1]+JI[0,2]*vv[2,1],JI[0,0]*vv[0,2]+JI[0,1]*vv[1,2]+JI[0,2]*vv[2,2]],
                  [JI[1,0]*vv[0,0]+JI[1,1]*vv[1,0]+JI[1,2]*vv[2,0],JI[1,0]*vv[0,1]+JI[1,1]*vv[1,1]+JI[1,2]*vv[2,1],JI[1,0]*vv[0,2]+JI[1,1]*vv[1,2]+JI[1,2]*vv[2,2]],
                  [JI[2,0]*vv[0,0]+JI[2,1]*vv[1,0]+JI[2,2]*vv[2,0],JI[2,0]*vv[0,1]+JI[2,1]*vv[1,1]+JI[2,2]*vv[2,1],JI[2,0]*vv[0,2]+JI[2,1]*vv[1,2]+JI[2,2]*vv[2,2]]])
        TD=array([[td[0,0]**2,         td[1,0]**2,       td[2,0]**2,     td[1,0]*td[2,0],                td[0,0]*td[2,0],                td[0,0]*td[1,0]],
                  [td[0,1]**2,         td[1,1]**2,       td[2,1]**2,     td[1,1]*td[2,1],                td[0,1]*td[2,1],                td[0,1]*td[1,1]],
                  [td[0,2]**2,         td[1,2]**2,       td[2,2]**2,     td[1,2]*td[2,2],                td[0,2]*td[2,2],                td[0,2]*td[1,2]],
                  [2*td[0,1]*td[0,2],2*td[1,1]*td[1,2],2*td[2,1]*td[2,2],td[1,1]*td[2,2]+td[2,1]*td[1,2],td[0,1]*td[2,2]+td[0,2]*td[2,1],td[0,1]*td[1,2]+td[1,1]*td[0,2]],
                  [2*td[0,0]*td[0,2],2*td[1,0]*td[1,2],2*td[2,0]*td[2,2],td[1,0]*td[2,2]+td[2,0]*td[1,2],td[0,0]*td[2,2]+td[0,2]*td[2,0],td[0,0]*td[1,2]+td[1,0]*td[0,2]],
                  [2*td[0,0]*td[0,1],2*td[1,0]*td[1,1],2*td[2,0]*td[2,1],td[1,0]*td[2,1]+td[2,0]*td[1,1],td[0,0]*td[2,1]+td[0,1]*td[2,0],td[0,0]*td[1,1]+td[1,0]*td[0,1]]])
#            for kk in xrange(3): sys.stdout.write('%10.4f'%(JJ[k,kk]))
#            sys.stdout.write('\n')  
#        for k in xrange(3): 
#            for kk in xrange(3): sys.stdout.write('%10.4f'%(vv[k,kk]))
#            sys.stdout.write('\n')  
#        for k in xrange(3): 
#            for kk in xrange(3): sys.stdout.write('%10.4f'%(td[k,kk]))
#            sys.stdout.write('\n')  
        return BB, det, TD
    
    def FormT(self, r, s, t):                               # interpolation on temperature - currently not used
        T = array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        return T
    def UpdateElemData(self):
        for i in range(4):                                 # loop over nodes
            self.Vn[:,i] = self.Vn[:,i] + self.VD[:,i] 
            self.Vg[:,i] = self.Vg[:,i] + self.VD[:,i]
            self.V1[:,i], self.V2[:,i] = self.CompleteTriad(self.Vg[:,i]) # may be overwritten in the following in case of a local system
        if self.Rot:
            TTW = [None, None, None, None]                  # for temporal storage of transformation coefficients
            for i in range(4):                             # loop over nodes of element  
                if self.SixthDoF[i]: self.ComputeTransLocal(i,self.Vn[:,i],TTW)  # transformation matrix for modified director system, completed director triad for local system if required
            self.ComputeTransLocalAll(TTW)                  # build rotation / transformation matrix for element / coordinate transformation matrix, not quadratic anymore!
        self.ComputeGG()                                    # scaled axes for rotational degrees of freedom
    def UpdateCoord(self, dis, ddis ):
        # length of dis, ddis must be 20 here!
        self.XX[0,0] = self.XX0[0,0] + dis[0]               # x-displacement 1st node
        self.XX[1,0] = self.XX0[1,0] + dis[1]               # y-displacement
        self.XX[2,0] = self.XX0[2,0] + dis[2]               # z-displacement
        self.XX[0,1] = self.XX0[0,1] + dis[5]               # x-displacement 2nd node
        self.XX[1,1] = self.XX0[1,1] + dis[6]               # y-displacement
        self.XX[2,1] = self.XX0[2,1] + dis[7]
        self.XX[0,2] = self.XX0[0,2] + dis[10]              # x-displacement 3rd node
        self.XX[1,2] = self.XX0[1,2] + dis[11]              # y-displacement
        self.XX[2,2] = self.XX0[2,2] + dis[12]
        self.XX[0,3] = self.XX0[0,3] + dis[15]              # x-displacement 4th node
        self.XX[1,3] = self.XX0[1,3] + dis[16]              # y-displacement
        self.XX[2,3] = self.XX0[2,3] + dis[17]
        self.ComputeEdgeDir()
        self.VD[0,0] =  - self.V2[0,0]*ddis[3]  + self.V1[0,0]*ddis[4] # ddis is increment of rotation angle in time increment; V1, V2 are displaced director triad at beginning of time increment
        self.VD[1,0] =  - self.V2[1,0]*ddis[3]  + self.V1[1,0]*ddis[4]
        self.VD[2,0] =  - self.V2[2,0]*ddis[3]  + self.V1[2,0]*ddis[4]
        self.VD[0,1] =  - self.V2[0,1]*ddis[8]  + self.V1[0,1]*ddis[9]
        self.VD[1,1] =  - self.V2[1,1]*ddis[8]  + self.V1[1,1]*ddis[9]
        self.VD[2,1] =  - self.V2[2,1]*ddis[8]  + self.V1[2,1]*ddis[9]
        self.VD[0,2] =  - self.V2[0,2]*ddis[13] + self.V1[0,2]*ddis[14]
        self.VD[1,2] =  - self.V2[1,2]*ddis[13] + self.V1[1,2]*ddis[14]
        self.VD[2,2] =  - self.V2[2,2]*ddis[13] + self.V1[2,2]*ddis[14]
        self.VD[0,3] =  - self.V2[0,3]*ddis[18] + self.V1[0,3]*ddis[19]
        self.VD[1,3] =  - self.V2[1,3]*ddis[18] + self.V1[1,3]*ddis[19]
        self.VD[2,3] =  - self.V2[2,3]*ddis[18] + self.V1[2,3]*ddis[19]
        for i in range(4): self.LengV(self.Vn[:,i])
    def GeomStiff(self, r, s, t, sig):
        GeomK = zeros((20, 20), dtype=float)
        N =  array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25])
        br = array([(-1+s)*0.25,      ( 1-s)*0.25,     ( 1+s)*0.25,      -( 1+s)*0.25])
        bs = array([(-1+r)*0.25,     -( 1+r)*0.25,     ( 1+r)*0.25,       ( 1-r)*0.25])
        t2 = t*t
#        bi = 0                                              # base for dof index
        for i in range(4):                                 # loop over nodes
            bi = i*5                                        # index for dofs
#            bj = 0                                          # base for dof index
            for j in range(4):                             # loop over nodes
                bj = j*5                                    # index for dofs
                S11   = sig[0] * br[i]*br[j]                #  s_11 sequence according to voigt notation
                S22   = sig[1] * bs[i]*bs[j]                #  s_22
                S33   = sig[2] * N[i] *N[j]                 #  s_33
                S23si = sig[3] * bs[i]*N[j]                 #  s_23
                S23sj = sig[3] * N[i] *bs[j]                #  s_23
                S13ri = sig[4] * br[i]*N[j]                 #  s_13
                S13rj = sig[4] * N[i] *br[j]                #  s_13
                S12   = sig[5] *(br[i]*bs[j]+bs[i]*br[j])   #  s_12     
                ggg00 = (self.gg[0,0,i]*self.gg[0,0,j] + self.gg[0,1,i]*self.gg[0,1,j] + self.gg[0,2,i]*self.gg[0,2,j])
                ggg11 = (self.gg[1,0,i]*self.gg[1,0,j] + self.gg[1,1,i]*self.gg[1,1,j] + self.gg[1,2,i]*self.gg[1,2,j])
                ggg01 = (self.gg[0,0,i]*self.gg[1,0,j] + self.gg[0,1,i]*self.gg[1,1,j] + self.gg[0,2,i]*self.gg[1,2,j])
                ggg10 = (self.gg[1,0,i]*self.gg[0,0,j] + self.gg[1,1,i]*self.gg[0,1,j] + self.gg[1,2,i]*self.gg[0,2,j])
                
                if False:
                    GeomK[bi,  bj]   = S11+S22+S12
                    GeomK[bi+1,bj+1] = S11+S22+S12
                    GeomK[bi+2,bj+2] = S11+S22+S12
                    GeomK[bi+3,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33 + S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg00
                    GeomK[bi+4,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33 + S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg11
                    GeomK[bi,  bj+1] = 0. 
                    GeomK[bi,  bj+2] = 0.
                    GeomK[bi,  bj+3] = (S11*t + S22*t + S13ri + S23si + S12*t)*self.gg[0,0,j]
                    GeomK[bi,  bj+4] = (S11*t + S22*t + S13ri + S23si + S12*t)*self.gg[1,0,j]
                    GeomK[bi+1,bj]   = 0.
                    GeomK[bi+1,bj+2] = 0.
                    GeomK[bi+1,bj+3] = (S11*t + S22*t + S13ri + S23si + S12*t)*self.gg[0,1,j]
                    GeomK[bi+1,bj+4] = (S11*t + S22*t + S13ri + S23si + S12*t)*self.gg[1,1,j]
                    GeomK[bi+2,bj]   = 0.
                    GeomK[bi+2,bj+1] = 0.
                    GeomK[bi+2,bj+3] = (S11*t + S22*t + S12*t + S13ri + S23si)*self.gg[0,2,j]
                    GeomK[bi+2,bj+4] = (S11*t + S22*t + S12*t + S13ri + S23si)*self.gg[1,2,j]
                    GeomK[bi+3,bj]   = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[0,0,i]
                    GeomK[bi+3,bj+1] = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[0,1,i]
                    GeomK[bi+3,bj+2] = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[0,2,i]
                    GeomK[bi+3,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33 + S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg01
                    GeomK[bi+4,bj]   = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[1,0,i]
                    GeomK[bi+4,bj+1] = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[1,1,i]
                    GeomK[bi+4,bj+2] = (S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[1,2,i]
                    GeomK[bi+4,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33 + S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg10
                else:
                    GeomK[bi,  bj]   = S11+S22+S12
                    GeomK[bi+1,bj+1] = S11+S22+S12
                    GeomK[bi+2,bj+2] = S11+S22+S12
                    GeomK[bi+3,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg00
                    GeomK[bi+4,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg11
                    GeomK[bi,  bj+1] = 0. 
                    GeomK[bi,  bj+2] = 0.
                    GeomK[bi,  bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,0,j]
                    GeomK[bi,  bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,0,j]
                    GeomK[bi+1,bj]   = 0.
                    GeomK[bi+1,bj+2] = 0.
                    GeomK[bi+1,bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,1,j]
                    GeomK[bi+1,bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,1,j]
                    GeomK[bi+2,bj]   = 0.
                    GeomK[bi+2,bj+1] = 0.
                    GeomK[bi+2,bj+3] = (S11*t + S22*t + S12*t)*self.gg[0,2,j]
                    GeomK[bi+2,bj+4] = (S11*t + S22*t + S12*t)*self.gg[1,2,j]
                    GeomK[bi+3,bj]   = (S11*t + S22*t + S12*t)*self.gg[0,0,i]
                    GeomK[bi+3,bj+1] = (S11*t + S22*t + S12*t)*self.gg[0,1,i]
                    GeomK[bi+3,bj+2] = (S11*t + S22*t + S12*t)*self.gg[0,2,i]
                    GeomK[bi+3,bj+4] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg01
                    GeomK[bi+4,bj]   = (S11*t + S22*t + S12*t)*self.gg[1,0,i]
                    GeomK[bi+4,bj+1] = (S11*t + S22*t + S12*t)*self.gg[1,1,i]
                    GeomK[bi+4,bj+2] = (S11*t + S22*t + S12*t)*self.gg[1,2,i]
                    GeomK[bi+4,bj+3] = (S11*t2 + S22*t2 + S12*t2 + S33)*ggg10
                    if True:
                        GeomK[bi+3,bj+3] += (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg00
                        GeomK[bi+4,bj+4] += (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg11
                        GeomK[bi,  bj+3] += (S13ri + S23si)*self.gg[0,0,j]
                        GeomK[bi,  bj+4] += (S13ri + S23si)*self.gg[1,0,j]
                        GeomK[bi+1,bj+3] += (S13ri + S23si)*self.gg[0,1,j]
                        GeomK[bi+1,bj+4] += (S13ri + S23si)*self.gg[1,1,j]
                        GeomK[bi+2,bj+3] += (S13ri + S23si)*self.gg[0,2,j]
                        GeomK[bi+2,bj+4] += (S13ri + S23si)*self.gg[1,2,j]
                        GeomK[bi+3,bj]   += (S13rj + S23sj)*self.gg[0,0,i]
                        GeomK[bi+3,bj+1] += (S13rj + S23sj)*self.gg[0,1,i]
                        GeomK[bi+3,bj+2] += (S13rj + S23sj)*self.gg[0,2,i]
                        GeomK[bi+3,bj+4] += (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg01
                        GeomK[bi+4,bj]   += (S13rj + S23sj)*self.gg[1,0,i]
                        GeomK[bi+4,bj+1] += (S13rj + S23sj)*self.gg[1,1,i]
                        GeomK[bi+4,bj+2] += (S13rj + S23sj)*self.gg[1,2,i]
                        GeomK[bi+4,bj+3] += (S13ri*t + S13rj*t + S23si*t + S23sj*t)*ggg10
                        
#                if self.Label==4030:
#                    if bi==0: # or bi==5:
#                        if bj==0 or bj==5 or bj==10 or bj==15:
#                            sys.stdout.write('%4i %4i%4i %4i%4i\n'%(self.Label,bi,bj,i,j))
#                            sys.stdout.write('%10.4f%10.4f%10.4f  %10.4f\n'%(S11,S22,S12,GeomK[bi,bj]))
                
#        GeomK = zeros((self.DofE, self.DofE), dtype=float)
#        if self.Label==6:
#            for bi in [0,5]:
#                for bj in [10]:
#                    sys.stdout.write('%4i%4i%4i\n'%(self.Label,bi,bj))
#                    sys.stdout.write('%10.4f%10.4f%10.4f%10.4f\n'%(GeomK[bi+1,bj+3],(S11*t + S22*t + S13ri + S23si + S12*t)*self.gg[0,1,j],\
#                                                                   GeomK[bi+3,bj+1],(S11*t + S22*t + S12*t + S13rj + S23sj)*self.gg[0,1,i]))
#                    for kk in xrange(5):
#                        for ll in xrange(5): sys.stdout.write('%10.4f'%(GeomK[bi+kk,bj+ll])) 
#                        sys.stdout.write('\n')

        return GeomK
    def JacoD(self, r, s, t):
        N = array([(1-r)*(1-s)*0.25, (1+r)*(1-s)*0.25, (1+r)*(1+s)*0.25, (1-r)*(1+s)*0.25])
        br = array([(-1+s)*0.25,( 1-s)*0.25,( 1+s)*0.25,-( 1+s)*0.25])
        bs = array([(-1+r)*0.25,-( 1+r)*0.25,( 1+r)*0.25,( 1-r)*0.25])
        JJ = zeros((3,3),dtype=float)
        for k in range(3): JJ[k,0]=br[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+br[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+br[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+br[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,1]=bs[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+bs[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+bs[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])+bs[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,2]=N[0]*(1/2.*self.a[0]*self.Vn[k,0])+N[1]*(1/2.*self.a[1]*self.Vn[k,1])+N[2]*(1/2.*self.a[2]*self.Vn[k,2])+N[3]*(1/2.*self.a[3]*self.Vn[k,3])
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        return det
    def StressIntegration(self, j, offset, nRe, ElemResults, FF):           # j is index of base integration point
        r = SamplePoints[self.IntT,self.nInt-1,j][0]                        # local integration point coordinates in reference surface
        s = SamplePoints[self.IntT,self.nInt-1,j][1]                        # "
        aa = dot( self.FormX(r,s,0), self.a)                                # interpolated shell thickness from node thicknesses
        # see WriteElemData
        if   self.Type=='SH4': Corr, CorM =0.5, 0.5             # to compensate for local iso-parametric coordinate t in range [-1..1] --> 2
        elif self.Type=='SH3': Corr, CorM =3.0, 0.5             # " + Corr -->  0.5 * 8 *3/4 = 3, see ConFemBasic::SampleWeights SH3
        #  
        if ElemResults == None:
            r = SamplePoints[self.IntT,self.nInt-1,j][0]                    # local integration point coordinates in reference surface
            s = SamplePoints[self.IntT,self.nInt-1,j][1]                    # "
            aa = dot( self.FormX(r,s,0), self.a)                            # interpolated shell thickness from node thicknesses
            nx, ny, nxy, qx, qy, mx, my, mxy = 0., 0., 0., 0., 0., 0., 0., 0.
#            nxR,nyR,nxyR,qxR,qyR,mxR,myR,mxyR= 0., 0., 0., 0., 0., 0., 0., 0.
#            nxC,nyC,nxyC,qxC,qyC,mxC,myC,mxyC= 0., 0., 0., 0., 0., 0., 0., 0.
            Lis2, Lis3 = self.Lists2( nRe, j)                               # integration point indices specific for base point
            for k in Lis2:                                                  # loop over height integration points - Lis2 should also hold reinforcement layers
                jj = j+k
                ppp = self.Data[jj]
                #
                if (k in Lis3):                                             # reinforcement layers only 
                    t  =         SamplePointsRCShell[self.Set,self.IntT,self.nInt-1,jj][2]
                    ff = Corr*aa*SampleWeightRCShell[self.Set,self.IntT,self.nInt-1,jj]  # 0.5 seems to compensate for SampleWeight local coordinates
#                    nxR = nxR + ff*ppp[0+offset]
#                    nyR = nyR + ff*ppp[1+offset]
#                    qyR = qyR + ff*ppp[3+offset]
#                    qxR = qxR + ff*ppp[4+offset]
#                    nxyR= nxyR+ ff*ppp[5+offset]
#                    mxR = mxR + 0.5*aa*t*ff*ppp[0+offset]
#                    myR = myR + 0.5*aa*t*ff*ppp[1+offset]
#                    mxyR= mxyR+ 0.5*aa*t*ff*ppp[5+offset]
                else:                                                       # bulk only        
                    t  =         SamplePoints[self.IntT,self.nInt-1,jj][2]
                    ff = Corr*aa*SampleWeight[self.IntT,self.nInt-1,jj]
#                    nxC = nxC + ff*ppp[0+offset]
#                    nyC = nyC + ff*ppp[1+offset]
#                    qyC = qyC + ff*ppp[3+offset]
#                    qxC = qxC + ff*ppp[4+offset]
#                    nxyC= nxyC+ ff*ppp[5+offset]
#                    mxC = mxC + 0.5*aa*t*ff*ppp[0+offset]
#                    myC = myC + 0.5*aa*t*ff*ppp[1+offset]
#                    mxyC= mxyC+ 0.5*aa*t*ff*ppp[5+offset]
                #
                nx = nx + ff*ppp[0+offset]                                  # do currently not know for what offset is uhc 200630
                ny = ny + ff*ppp[1+offset]
                qy = qy + ff*ppp[3+offset]
                qx = qx + ff*ppp[4+offset]
                nxy= nxy+ ff*ppp[5+offset]
                mx = mx + CorM*aa*t*ff*ppp[0+offset]
                my = my + CorM*aa*t*ff*ppp[1+offset]
                mxy= mxy+ CorM*aa*t*ff*ppp[5+offset]
                
#            nx_,ny_,nxy_ = nx-nxR-nxC, ny-nyR-nyC, nxy-nxyR-nxyC 
#            mx_,my_,mxy_ = mx-mxR-mxC, my-myR-myC, mxy-mxyR-mxyC
#            qx_, qy_ = qx-qxR-qxC, qy-qyR-qyC
                
#            Echo(f"X {self.Label:d},{j:d}\n{nxR:9.4e},{nyR:9.4e},{nxyR:9.4e},{qxR:9.4e},{qyR:9.4e},{mxR:9.4e},{myR:9.4e},{mxyR:9.4e}", FF)
#            Echo(                        f"{nxC:9.4e},{nyC:9.4e},{nxyC:9.4e},{qxC:9.4e},{qyC:9.4e},{mxC:9.4e},{myC:9.4e},{mxyC:9.4e}", FF)
#            Echo(                        f"{nx:9.4e},{ny:9.4e},{nxy:9.4e},{qx:9.4e},{qy:9.4e},{mx:9.4e},{my:9.4e},{mxy:9.4e}", FF)
#            Echo(                        f"{nx_:9.4e},{ny_:9.4e},{nxy_:9.4e},{qx_:9.4e},{qy_:9.4e},{mx_:9.4e},{my_:9.4e},{mxy_:9.4e}", FF)
#            val = 1.0e-15
#            if abs(nx_)>val or abs(ny_)>val or abs(nxy_)>val or abs(mx_)>val or abs(my_)>val or abs(mxy_)>val or abs(qx_)>val or abs(qy_)>val:
#                raise NameError("XXX")   
        else:
            key = str( self.Label)+self.Set+str(j)
            if key in ElemResults:
                nx, ny, nxy, mx, my, mxy, qx, qy = ElemResults[key][5][3],ElemResults[key][5][4],ElemResults[key][5][5],\
                                                   ElemResults[key][5][6],ElemResults[key][5][7],ElemResults[key][5][8],\
                                                   ElemResults[key][5][9],ElemResults[key][5][10] # Data
            else: 
                raise NameError("ConFemElem:StressIntegration: unknown key for data 2",self.Label)
        return nx, ny, nxy, qx, qy, mx, my, mxy, aa 
    
class SH3( SH4 ):
    """ Shell element 3 nodes
    """
    def ComputeGG(self):
        for i in range(3): 
            for j in range(3): self.gg[0,i,j] = -self.a[j]*self.V2[i,j]/2.     # 1st index direction, 2nd index node 
            for j in range(3): self.gg[1,i,j] =  self.a[j]*self.V1[i,j]/2.
    def ComputeTransLocalAll(self, TTW):
        base, base_ = 0, 0
        for i in range(3):                                 # loop over nodes of element  
            for j in range(3):          self.Trans[base_+j,base+j] = 1.
            if False:
#            if self.SixthDoF[i]:
                for j in range(2):
                    for jj in range(3): self.Trans[base_+3+j,base+3+jj] = TTW[i][j][jj] 
#   uhc             XX = dot(transpose(TTW[i]),TTW[i])
#   uhc             YY = dot(TTW[i],transpose(TTW[i]))
#                print i, TTW[i],'\n', XX, '\n', YY
            else:
                for j in range(2):      self.Trans[base_+3+j,base+3+j] = 1.
            base = base + self.DofN[i]
            base_= base_+ self.DofNini[i]            
#        for j in xrange(self.Trans.shape[0]):
#            for jj in xrange(self.Trans.shape[1]): sys.stdout.write('%6.2f'%(self.Trans[j,jj]))
#            sys.stdout.write('\n')
#        raise NameError ("Exit")
    def __init__(self, Label, SetLabel, InzList, MatName,Material, NoList, ShellSec, StateV, NData, RCFlag, NoLabToNoInd):
#       Element.__init__(self,"SH3",InzList, 3, 5,2,15, [set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5])], 21,False,Label,SetLabel,2, MatName,StateV,NData, NoList,NoLabToNoInd,[]) # five integration points over cross section height
        Element.__init__(self,"SH3",InzList, 3, 5,2,15, [set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5]),set([1, 2, 3, 4, 5])], 21,False,Label,SetLabel,3, MatName,StateV,NData, NoList,NoLabToNoInd,[]) # five integration points over cross section height
    def Lists1(self):                                       # indices for first integration points in base area
        if self.nInt==2:                                    # integration order 
            Lis = [0,5,10]                                  # indices for first integration points in base area, 4 Gaussian integration points over cross section height
        return Lis
    def Lists2(self, nRe, j):                               # integration point indices specific for base point
        if self.nInt==2: 
            Lis2, Lis3 = [0,1,2,3,4], []            
        return Lis2, Lis3                                   # RC not yet implemented
    def Ini2_(self, NoList,NoIndToCMInd, MaList, SecDict):  # Ini2 from SH4
        return []
    def Ini3(self, NoList,NoIndToCMInd ):
        i0 = NoIndToCMInd[self.Inzi[0]]
        i1 = NoIndToCMInd[self.Inzi[1]]
        i2 = NoIndToCMInd[self.Inzi[2]]
        ni = [NoList[i0], NoList[i1], NoList[i2]]
        self.XX =      array([[ni[0].XCo, ni[1].XCo, ni[2].XCo],
                              [ni[0].YCo, ni[1].YCo, ni[2].YCo],
                              [ni[0].ZCo, ni[1].ZCo, ni[2].ZCo]])           # collect nodal coordinates in a compact form for later use
        self.EdgeDir = zeros((3),dtype=float)                               # initialize direction of 1st edge for later use
        self.ComputeEdgeDir()                                               # direction of 1st edge for later use
        nn = zeros((3,3), dtype = float)                                    # 1st index for nodes, 2nd for directions
        nn[0,:] = self.CompNoNor( 0, 1, 2)                                  # roughly unit normals to shell surface at nodes
        nn[1,:] = self.CompNoNor( 1, 2, 0)
        nn[2,:] = self.CompNoNor( 2, 0, 1)
        self.V1 = zeros((3,3),dtype=float)                                  # initialize director triad, 1st index for directions, 2nd index for nodes
        self.V2 = zeros((3,3),dtype=float)
        self.Vn = zeros((3,3),dtype=float)
        self.Vg = zeros((3,3),dtype=float)                                  # director as defined per node via input data, might not be the actually used director as is ruled in the following
#        self.VD = zeros((3,4),dtype=float)                  # director increment in time increment
        TTW = [None, None, None, None]                                      # for temporal storage of transformation coefficients
        self.SixthDoF = [False, False, False]
        for i in range(3):                                                  # loop over nodes of element  
            LL = sqrt(ni[i].XDi**2+ni[i].YDi**2+ni[i].ZDi**2)               # length of directors
            self.Vg[0,i] = ni[i].XDi/LL                                     # values given with input data
            self.Vg[1,i] = ni[i].YDi/LL                                     # "
            self.Vg[2,i] = ni[i].ZDi/LL                                     # "
            self.V1[:,i], self.V2[:,i] = self.CompleteTriad(self.Vg[:,i])
            self.Vn[0,i] = self.Vg[0,i]                                 # final director
            self.Vn[1,i] = self.Vg[1,i]
            self.Vn[2,i] = self.Vg[2,i]
        if self.Rot: self.ComputeTransLocalAll(TTW)                         # ??? build rotation / transformation matrix for element / coordinate transformation matrix, not quadratic anymore!
        self.XX0 = self.XX.copy() #copy(self.XX)                                            # retain initial values for NLGEOM
        self.DofI = array([[-1,-1,-1,-1,-1,-1],[-1,-1,-1,-1,-1,-1],[-1,-1,-1,-1,-1,-1]],dtype=int)
        self.gg = zeros((2,3,3), dtype=float)                               # initialization of scaled rotation axes; 1st index for rotation dofs, 2nd index for directions of rotation axes, 3rd for nodes 
        self.ComputeGG()                                                    # scaled axes for rotational degrees of freedom
        AA = 0.
        for i in range(15):                                                 # volume by numerical integration for characteristic length  
            r = SamplePoints[5,1,i][0]
            s = SamplePoints[5,1,i][1]
            f = self.JacoD(r,s,0)*SampleWeight[5,1,i]
            AA = AA + f
        self.Lch_ = sqrt(AA/(0.25*(self.a[0]+self.a[1]+self.a[2])))         #/self.nInt   # characteristic length -> side length of square of same area of shell area
        # --> common method with sh4
        self.CrBScaleType()
#        if MaList[self.MatN].RType ==2:                 # find scaling factor for band width regularization
#            x = MaList[self.MatN].bw/self.Lch_          # ratio of crack band width to characteristic length
#            CrX, CrY = MaList[self.MatN].CrX, MaList[self.MatN].CrY # support points for tensile softening scaling factor 
#            i = bisect_left(CrX, x) 
#            if i>0 and i<(MaList[self.MatN].CrBwN+1): 
#                self.CrBwS = CrY[i-1] + (x-CrX[i-1])/(CrX[i]-CrX[i-1])*(CrY[i]-CrY[i-1]) # scaling factor by linear interpolation
#            else:
#                print('ZZZ', AA, self.Lch_,x,'\n', CrX, '\n', CrY, i, MaList[self.MatN].CrBwN)  
#                raise NameError("ConFemElem:SH3.Ini2: RType 2 - element char length exceeds scaling factor interpolation")

    def Basics(self, r, s, t):
        N =  array([ 1. -r-s, r, s])
        br = array([-1., 1., 0.])
        bs = array([-1., 0., 1.])
        JJ = zeros((3,3),dtype=float)
        # following may be simplified due to br[2]=bs[1]=0 and so on
        for k in range(3): JJ[k,0]=br[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+br[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+br[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2]) #+br[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,1]=bs[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+bs[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+bs[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2]) #+bs[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,2]=N[0]*(1/2.*self.a[0]*self.Vn[k,0])+N[1]*(1/2.*self.a[1]*self.Vn[k,1])+N[2]*(1/2.*self.a[2]*self.Vn[k,2]) #+N[3]*(1/2.*self.a[3]*self.Vn[k,3])
        JI = inv(JJ)
        ll=sqrt(JJ[0,2]**2+JJ[1,2]**2+JJ[2,2]**2)   
        vv = array([[0.,0.,JJ[0,2]/ll],[0.,0.,JJ[1,2]/ll],[0.,0.,JJ[2,2]/ll]]) # normal of local coordinate system, 3RD COLUMN
        Loc = False
        if Loc:                                         # local right handed coordinate system with 1st direction / column aligned to element edge
            x0 = self.EdgeDir[1]*vv[2,2]-self.EdgeDir[2]*vv[1,2]
            x1 = self.EdgeDir[2]*vv[0,2]-self.EdgeDir[0]*vv[2,2]
            x2 = self.EdgeDir[0]*vv[1,2]-self.EdgeDir[1]*vv[0,2]
            xx = sqrt(x0**2+x1**2+x2**2)
            vv[0,1] = -x0/xx                            # 2ND COLUMN, approx perp. to element edge, reversed in sign to preserve right handedness
            vv[1,1] = -x1/xx
            vv[2,1] = -x2/xx 
            x0 = vv[1,2]*vv[2,1]-vv[2,2]*vv[1,1]
            x1 = vv[2,2]*vv[0,1]-vv[0,2]*vv[2,1]
            x2 = vv[0,2]*vv[1,1]-vv[1,2]*vv[0,1]
            xx = sqrt(x0**2+x1**2+x2**2)
            vv[0,0] = -x0/xx                            # 1ST COLUMN, approx aligned to element edge, sign reversal of 2nd column is implicitly corrected
            vv[1,0] = -x1/xx
            vv[2,0] = -x2/xx 
        else:
            if abs(vv[1,2])<0.99:                       # local coordinate system V_1 from cross product of V_n and e_y ( V_1 in e_x - e_z plane) if V_n is not to much aligned to e_y  
                ll = sqrt(vv[2,2]**2+vv[0,2]**2)        # length of V_1
                vv[0,0] = vv[2,2]/ll                    # V_1[0] normalized;  V1[1] = 0
                vv[2,0] =-vv[0,2]/ll                    # V_1[2] normalized

                vv[0,1] = vv[1,2]*vv[2,0]               # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
                vv[1,1] = vv[2,2]*vv[0,0]-vv[0,2]*vv[2,0]
                vv[2,1] =-vv[1,2]*vv[0,0]
            else:                                       # local coordinate system V_1 from cross product of V_n and e_x ( V_1 in e_y - e_z plane)
                ll = sqrt(vv[2,2]**2+vv[1,2]**2)        # length of V_1
                vv[1,0] =-vv[2,2]/ll                    # V_1[0] normalized;  V1[0] = 0
                vv[2,0] = vv[1,2]/ll                    # V_1[2] normalized

                vv[0,1] = vv[1,2]*vv[2,0]-vv[2,2]*vv[1,0] # as V_n and V_1 are orthogonal and both have unit length V_2 also should have unit length
                vv[1,1] =-vv[0,2]*vv[2,0]
                vv[2,1] = vv[0,2]*vv[1,0]
        return N, br, bs, JJ, JI, vv

    def FormB(self, r, s, t, NLg):
        N, _, _, JJ, JI, vv = self.Basics( r, s, t)
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        HH = zeros((3,2,3),dtype=float)
        for k in range(3):
            HH[0,0,k]=JJ[0,0]*self.gg[0,0,k]+JJ[1,0]*self.gg[0,1,k]+JJ[2,0]*self.gg[0,2,k]
            HH[0,1,k]=JJ[0,0]*self.gg[1,0,k]+JJ[1,0]*self.gg[1,1,k]+JJ[2,0]*self.gg[1,2,k]
            HH[1,0,k]=JJ[0,1]*self.gg[0,0,k]+JJ[1,1]*self.gg[0,1,k]+JJ[2,1]*self.gg[0,2,k]
            HH[1,1,k]=JJ[0,1]*self.gg[1,0,k]+JJ[1,1]*self.gg[1,1,k]+JJ[2,1]*self.gg[1,2,k]
            HH[2,0,k]=JJ[0,2]*self.gg[0,0,k]+JJ[1,2]*self.gg[0,1,k]+JJ[2,2]*self.gg[0,2,k]
            HH[2,1,k]=JJ[0,2]*self.gg[1,0,k]+JJ[1,2]*self.gg[1,1,k]+JJ[2,2]*self.gg[1,2,k]
        BB = zeros((6,15),dtype=float)
        BB[0,0] =-JJ[0,0]
        BB[0,1] =-JJ[1,0]
        BB[0,2] =-JJ[2,0]
        BB[0,3] =-t*HH[0,0,0]
        BB[0,4] =-t*HH[0,1,0]
        BB[0,5] = JJ[0,0]
        BB[0,6] = JJ[1,0]
        BB[0,7] = JJ[2,0]
        BB[0,8] = t*HH[0,0,1]
        BB[0,9] = t*HH[0,1,1]
        BB[0,10]= 0
        BB[0,11]= 0
        BB[0,12]= 0
        BB[0,13]= 0
        BB[0,14]= 0
        BB[1,0] =-JJ[0,1]
        BB[1,1] =-JJ[1,1]
        BB[1,2] =-JJ[2,1]
        BB[1,3] =-t*HH[1,0,0]
        BB[1,4] =-t*HH[1,1,0]
        BB[1,5] = 0
        BB[1,6] = 0
        BB[1,7] = 0
        BB[1,8] = 0
        BB[1,9] = 0
        BB[1,10]= JJ[0,1]
        BB[1,11]= JJ[1,1]
        BB[1,12]= JJ[2,1]
        BB[1,13]= t*HH[1,0,2]
        BB[1,14]= t*HH[1,1,2]
        BB[2,0] = 0
        BB[2,1] = 0
        BB[2,2] = 0
        BB[2,3] = N[0]*HH[2,0,0]
        BB[2,4] = N[0]*HH[2,1,0]
        BB[2,5] = 0
        BB[2,6] = 0
        BB[2,7] = 0
        BB[2,8] = N[1]*HH[2,0,1]
        BB[2,9] = N[1]*HH[2,1,1]
        BB[2,10]= 0
        BB[2,11]= 0
        BB[2,12]= 0
        BB[2,13]= N[2]*HH[2,0,2]
        BB[2,14]= N[2]*HH[2,1,2]
        BB[5,0] =-JJ[0,0]-JJ[0,1]
        BB[5,1] =-JJ[1,0]-JJ[1,1]
        BB[5,2] =-JJ[2,0]-JJ[2,1]
        BB[5,3] = t*(-HH[0,0,0]-HH[1,0,0])
        BB[5,4] = t*(-HH[0,1,0]-HH[1,1,0])
        BB[5,5] = JJ[0,1]
        BB[5,6] = JJ[1,1]
        BB[5,7] = JJ[2,1]
        BB[5,8] = t*HH[1,0,1]
        BB[5,9] = t*HH[1,1,1]
        BB[5,10]= JJ[0,0]
        BB[5,11]= JJ[1,0]
        BB[5,12]= JJ[2,0]
        BB[5,13]= t*HH[0,0,2]
        BB[5,14]= t*HH[0,1,2]
        BB[3,0] =-(1.-r)*JJ[0,2]-r*JJ[0,2]+s*(JJ[0,2]-JJ[0,2])
        BB[3,1] =-(1.-r)*JJ[1,2]-r*JJ[1,2]+s*(JJ[1,2]-JJ[1,2])
        BB[3,2] =-(1.-r)*JJ[2,2]-r*JJ[2,2]+s*(JJ[2,2]-JJ[2,2])
        BB[3,3] = (1.-r)*(0.5*HH[1,0,0]-t*HH[2,0,0])-r*t*HH[2,0,0]+s*(t*HH[2,0,0]+0.5*HH[1,0,0]-t*HH[2,0,0])
        BB[3,4] = (1.-r)*(0.5*HH[1,1,0]-t*HH[2,1,0])-r*t*HH[2,1,0]+s*(t*HH[2,1,0]+0.5*HH[1,1,0]-t*HH[2,1,0])
        BB[3,5] = 0
        BB[3,6] = 0
        BB[3,7] = 0
        BB[3,8] = 0.5*r*HH[1,0,1]-0.5*s*HH[1,0,1]
        BB[3,9] = 0.5*r*HH[1,1,1]-0.5*s*HH[1,1,1]
        BB[3,10]= (1-r)*JJ[0,2]+r*JJ[0,2]+s*(-JJ[0,2]+JJ[0,2])
        BB[3,11]= (1-r)*JJ[1,2]+r*JJ[1,2]+s*(-JJ[1,2]+JJ[1,2])
        BB[3,12]= (1-r)*JJ[2,2]+r*JJ[2,2]+s*(-JJ[2,2]+JJ[2,2])
        BB[3,13]= (1-r)*(0.5*HH[1,0,2]+t*HH[2,0,2])-r*(-0.5*HH[1,0,2]-t*HH[2,0,2])+s*(-0.5*HH[1,0,2]-t*HH[2,0,2]+0.5*HH[1,0,2]+t*HH[2,0,2])
        BB[3,14]= (1-r)*(0.5*HH[1,1,2]+t*HH[2,1,2])-r*(-0.5*HH[1,1,2]-t*HH[2,1,2])+s*(-0.5*HH[1,1,2]-t*HH[2,1,2]+0.5*HH[1,1,2]+t*HH[2,1,2])
        BB[4,0] =-r*(-JJ[0,2]+JJ[0,2])-(1-s)*JJ[0,2]-s*JJ[0,2]
        BB[4,1] =-r*(-JJ[1,2]+JJ[1,2])-(1-s)*JJ[1,2]-s*JJ[1,2]
        BB[4,2] =-r*(-JJ[2,2]+JJ[2,2])-(1-s)*JJ[2,2]-s*JJ[2,2]
        BB[4,3] =-r*(-t*HH[2,0,0]-0.5*HH[0,0,0]+t*HH[2,0,0])+(1-s)*(0.5*HH[0,0,0]-t*HH[2,0,0])-s*t*HH[2,0,0]
        BB[4,4] =-r*(-t*HH[2,1,0]-0.5*HH[0,1,0]+t*HH[2,1,0])+(1-s)*(0.5*HH[0,1,0]-t*HH[2,1,0])-s*t*HH[2,1,0]
        BB[4,5] =-r*(JJ[0,2]-JJ[0,2])+(1-s)*JJ[0,2]+s*JJ[0,2]
        BB[4,6] =-r*(JJ[1,2]-JJ[1,2])+(1-s)*JJ[1,2]+s*JJ[1,2]
        BB[4,7] =-r*(JJ[2,2]-JJ[2,2])+(1-s)*JJ[2,2]+s*JJ[2,2]
        BB[4,8] =-r*(0.5*HH[0,0,1]+t*HH[2,0,1]-0.5*HH[0,0,1]-t*HH[2,0,1])+(1-s)*(0.5*HH[0,0,1]+t*HH[2,0,1])+s*(0.5*HH[0,0,1]+t*HH[2,0,1])
        BB[4,9] =-r*(0.5*HH[0,1,1]+t*HH[2,1,1]-0.5*HH[0,1,1]-t*HH[2,1,1])+(1-s)*(0.5*HH[0,1,1]+t*HH[2,1,1])+s*(0.5*HH[0,1,1]+t*HH[2,1,1])
        BB[4,10]= 0
        BB[4,11]= 0
        BB[4,12]= 0
        BB[4,13]=-0.5*r*HH[0,0,2]+0.5*s*HH[0,0,2]
        BB[4,14]=-0.5*r*HH[0,1,2]+0.5*s*HH[0,1,2]

        td=array([[JI[0,0]*vv[0,0]+JI[0,1]*vv[1,0]+JI[0,2]*vv[2,0],JI[0,0]*vv[0,1]+JI[0,1]*vv[1,1]+JI[0,2]*vv[2,1],JI[0,0]*vv[0,2]+JI[0,1]*vv[1,2]+JI[0,2]*vv[2,2]],
                  [JI[1,0]*vv[0,0]+JI[1,1]*vv[1,0]+JI[1,2]*vv[2,0],JI[1,0]*vv[0,1]+JI[1,1]*vv[1,1]+JI[1,2]*vv[2,1],JI[1,0]*vv[0,2]+JI[1,1]*vv[1,2]+JI[1,2]*vv[2,2]],
                  [JI[2,0]*vv[0,0]+JI[2,1]*vv[1,0]+JI[2,2]*vv[2,0],JI[2,0]*vv[0,1]+JI[2,1]*vv[1,1]+JI[2,2]*vv[2,1],JI[2,0]*vv[0,2]+JI[2,1]*vv[1,2]+JI[2,2]*vv[2,2]]])
        TD=array([[td[0,0]**2,         td[1,0]**2,       td[2,0]**2,     td[1,0]*td[2,0],                td[0,0]*td[2,0],                td[0,0]*td[1,0]],
                  [td[0,1]**2,         td[1,1]**2,       td[2,1]**2,     td[1,1]*td[2,1],                td[0,1]*td[2,1],                td[0,1]*td[1,1]],
                  [td[0,2]**2,         td[1,2]**2,       td[2,2]**2,     td[1,2]*td[2,2],                td[0,2]*td[2,2],                td[0,2]*td[1,2]],
                  [2*td[0,1]*td[0,2],2*td[1,1]*td[1,2],2*td[2,1]*td[2,2],td[1,1]*td[2,2]+td[2,1]*td[1,2],td[0,1]*td[2,2]+td[0,2]*td[2,1],td[0,1]*td[1,2]+td[1,1]*td[0,2]],
                  [2*td[0,0]*td[0,2],2*td[1,0]*td[1,2],2*td[2,0]*td[2,2],td[1,0]*td[2,2]+td[2,0]*td[1,2],td[0,0]*td[2,2]+td[0,2]*td[2,0],td[0,0]*td[1,2]+td[1,0]*td[0,2]],
                  [2*td[0,0]*td[0,1],2*td[1,0]*td[1,1],2*td[2,0]*td[2,1],td[1,0]*td[2,1]+td[2,0]*td[1,1],td[0,0]*td[2,1]+td[0,1]*td[2,0],td[0,0]*td[1,1]+td[1,0]*td[0,1]]])

        return BB, det, TD
    def FormT(self, r, s, t):                                               # interpolation on temperature - currently not used
        T = array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        return T
    def FormX(self, r, s, t):                                               # interpolation on geometry ??? see below
        X = array([ 1.-r-s, r, s])
        return X
    def FormX_(self, L1, L2, L3):
        L3=1-L1-L2
#        X = array([[L1, 0, L2, 0, L3, 0],
#                   [0, L1, 0, L2, 0, L3]])
        X = array([[L1, 0, 0,  L2, 0, 0,  L3, 0, 0],
                   [0, L1, 0,  0, L2, 0,  0, L3, 0],
                   [0,  0, L1, 0,  0, L2, 0, 0, L3]])
        return X
    def FormN(self, L1, L2, L3):                                            # supresses rotational degrees of freedom
        L3 = 1 - L1 - L2
        N = array([[L1, 0,  0,  0, 0, L2, 0,  0,  0, 0,  L3, 0,  0,  0, 0 ],
                   [0,  L1, 0,  0, 0, 0,  L2, 0,  0, 0,  0,  L3, 0,  0, 0 ],
                   [0,  0,  L1, 0, 0, 0,  0,  L2, 0, 0,  0,  0,  L3, 0, 0 ]])
        return N
    def JacoD(self, r, s, t):
        N =  array([ 1. -r-s, r, s])
        br = array([-1., 1., 0.])
        bs = array([-1., 0., 1.])
        JJ = zeros((3,3),dtype=float)
        # following may be simplified due to br[2]=bs[1]=0 and so on
        for k in range(3): JJ[k,0]=br[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+br[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+br[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])#+br[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,1]=bs[0]*(self.XX[k,0]+t/2.*self.a[0]*self.Vn[k,0])+bs[1]*(self.XX[k,1]+t/2.*self.a[1]*self.Vn[k,1])+bs[2]*(self.XX[k,2]+t/2.*self.a[2]*self.Vn[k,2])#+bs[3]*(self.XX[k,3]+t/2.*self.a[3]*self.Vn[k,3])
        for k in range(3): JJ[k,2]=N[0]*(1/2.*self.a[0]*self.Vn[k,0])+N[1]*(1/2.*self.a[1]*self.Vn[k,1])+N[2]*(1/2.*self.a[2]*self.Vn[k,2])#+N[3]*(1/2.*self.a[3]*self.Vn[k,3])
        det = JJ[0,0]*JJ[1,1]*JJ[2,2]-JJ[0,0]*JJ[1,2]*JJ[2,1]-JJ[1,0]*JJ[0,1]*JJ[2,2]+JJ[1,0]*JJ[0,2]*JJ[2,1]+JJ[2,0]*JJ[0,1]*JJ[1,2]-JJ[2,0]*JJ[0,2]*JJ[1,1]
        return det
        
   
