# CaeFemMain -- 2022-09-27
# Copyright (C) [2022] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from numpy import dot, array, mean, sqrt, arcsin, zeros, std, linspace, amin, amax
from scipy import stats as st
from math import cos, sin, pi, fabs
#from matplotlib.pyplot import figure, plot, grid, title, text, contour, clabel, show, axis, xticks, yticks, ylim, annotate
#from numpy.linalg import norm, det
from ConFemElem import ZeroD, CPE4, T3D3I
from ConFemInOut import ReadOptionsFile #, PlotNodes
from ConFemBasics import FindIndexByLabel, DirFilFromName
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from os import path as pth
import numpy as np
import pandas as pd
import matplotlib.colors as colors
import mpl_toolkits.mplot3d as a3
import pylab as pl
import time as ti
import sys
import os
sys.path.insert(1, '../src')
#import imp
from ConFemBasics import SamplePoints, _Truss2DAll_, _TriNodes_, _Bond2D_, _Length2D_,  _FourNodes2D_,\
                        _Shells3D_, _BeamsAll_, __Length1D__
from ConFemBasics import __Length2D__, __SpringsAll__, __BeamsBern__, __BeamsTimo__, __TriNodes__

Colors = ['blue','gray','black','darkviolet','magenta','indianred','darkred','cyan','yellow','blue']
#Colors = ['tab:green','tab:brown','magenta','darkviolet','blue','green','yellow']
#Colors = ['royalblue','darkviolet','green','magenta','blue','darkcyan','yellow']
Colors2 = ['red','green','blue','magenta','cyan','darkred']
ColorTable = {}
#ColorTable['IsoD'] = 'white'
ColorTable['ELASTIC'] = 'black'
ColorTable['BOND'] = 'green'
ColorTable['ELASTICLT'] = 'gray'
ColorTable['ElasticSDA'] = 'black'
ColorTable['MISES'] = 'red'
ColorTable['MICRODAMAGE'] = 'magenta'
ColorTable['Lubliner'] = 'cyan'
ColorTable['ISODAMAGE'] = 'cyan'
ColorTable['ELASTICLT_SDA'] = 'magenta'
ColorTable['RESHEET'] = 'red'
ColorTable['NLSLAB'] = 'red'

#Colors = ['tab:green','tab:brown','magenta','darkviolet','blue','green','yellow']
#Colors = ['royalblue','magenta','darkviolet','darkcyan','blue','green','yellow']
FonSizTi='x-large'   # fontsize title x-large
FonSizAx='x-large'     # fontsize axis
#FonSizAx='medium'     # fontsize axis
LiWiCrv=3            # line width curve in pt
LiWiCrvMid=2            # line width curve in pt
LiWiCrvThin=1.0        # line width curve in pt

class ValMinMax():
    def __init__(self, ni, nf):                                             # ni: dim of int values; nf: dim of float values
        self.Max = -1.0e+9
        self.Min =  1.0e+9
        self.datafMax = zeros((nf), dtype=float)
        self.datafMin = zeros((nf), dtype=float)
        self.dataiMax = zeros((ni), dtype=int)
        self.dataiMin = zeros((ni), dtype=int)
        self.used = False                                                   # might not be appropriate for all plot types
    def CheckMax(self, val, datai, dataf):
        if val>self.Max:
            self.Max = val
            for i in range(len(self.dataiMax)): self.dataiMax[i] = datai[i]
            for i in range(len(self.datafMax)): self.datafMax[i] = dataf[i]
    def CheckMin(self, val, datai, dataf):
        if val < self.Min:
            self.Min = val
            for i in range(len(self.dataiMin)): self.dataiMin[i] = datai[i]
            for i in range(len(self.datafMax)): self.datafMin[i] = dataf[i]
    def Annotate(self, pp):
        if len(self.datafMax)>0:
            pp.annotate('%6i'%self.dataiMax[0]+' max '+"%6.1e"%self.Max+' / '+"%6.1e"%self.datafMax[2],xy=(self.datafMax[0],self.datafMax[1]),
                        xycoords='data', xytext=(0.20, 0.87), textcoords='figure fraction',
                        arrowprops=dict(arrowstyle="->"), fontsize='medium')
        if len(self.datafMin)>0:
            pp.annotate('%6i'%self.dataiMin[0]+' min '+"%6.1e"%self.Min+' / '+"%6.1e"%self.datafMin[2],xy=(self.datafMin[0],self.datafMin[1]),
                        xycoords='data', xytext=(0.60, 0.87), textcoords='figure fraction',
                        arrowprops=dict(arrowstyle="->"), fontsize='medium')

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='large') # x-large
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 

# from Ahmad 2020
#def vtkPrincipal( NoList, ElemList, Name, timeStr, IntPointCoor, IntPointDispl, IntPointPrinStress, IntPointPrinVctr):
#    # collect data before writing to vtu
#    class Inptdata():
#        def __init__(self, coord, intpdisp, sigma, vctr):
#            self.coord    = coord
#            self.intpdisp = intpdisp
#            self.sigma    = sigma
#            self.vctr     = vctr
#    def WriteprincFile( fileOut_princ, intP):
#        f1 = open(fileOut_princ, 'w')
#        f1.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
#        f1.write('  <UnstructuredGrid>\n')
#        f1.write('    <Piece NumberOfPoints="' + str(len(intP)) + '" NumberOfCells="' + str(len(intP)) + '">\n')
#        #
#        f1.write('      <Points>\n')
#        f1.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.coord[0], i.coord[1], i.coord[2]))
#        f1.write('      </Points>\n')
#        f1.write('      <Cells>\n')
#        f1.write('        <DataArray type="Int32" Name="connectivity" format="ascii"/>\n')
#        for i in range(len(intP)): f1.write('     %3i\n' % (i))
#        f1.write('        <DataArray type="Int32" Name="offsets" format="ascii"/>\n')
#        for i in range(len(intP)): f1.write('     %3i\n' % (i))
#        f1.write('        <DataArray type="Int32" Name="types" format="ascii"/>\n')
#        for i in range(len(intP)): f1.write('       1\n')
#        f1.write('      </Cells>\n')
#        #
#        f1.write('      <PointData>\n')
#        f1.write('        <DataArray type="Float64" Name="Princpal_Stresses" NumberOfComponents=" 3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.sigma[0], i.sigma[1], i.sigma[2]))
#        f1.write('        <DataArray type="Float64" Name="IntP_displ" NumberOfComponents=" 3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.intpdisp[0], i.intpdisp[1], i.intpdisp[2]))
#        f1.write('        <DataArray type="Float64" Name="Vector_1" NumberOfComponents=" 3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[0][0], i.vctr[0][1], i.vctr[0][2]))
#        f1.write('        <DataArray type="Float64" Name="Vector_2" NumberOfComponents=" 3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[1][0], i.vctr[1][1], i.vctr[1][2]))
#        f1.write('        <DataArray type="Float64" Name="Vector_3" NumberOfComponents=" 3" format="ascii"/>\n')
#        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[2][0], i.vctr[2][1], i.vctr[2][2]))
##        f1.write('        <DataArray type="Float64" Name="Princpal_Stresses" NumberOfComponents=" 3" format="ascii"/>\n')
##        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.sigma[0], i.sigma[1], i.sigma[2]))
#        f1.write('      </PointData>\n')
#        f1.write('    </Piece>\n')
#        f1.write('  </UnstructuredGrid>\n')
#        f1.write('</VTKFile>\n')
#        f1.close()
#    # main part
#    intP = []
#    for a, b, c, d in zip(IntPointCoor, IntPointDispl, IntPointPrinStress, IntPointPrinVctr):
#        intP += [ Inptdata( a, b, c, d) ]
#    FFF = Name+'_Principal_ParaView'
#    if pth.exists(FFF): pass
#    else: os.mkdir(FFF)
#    fileOut_princ = FFF + '/Principal_'+timeStr+'.vtu'
#    WriteprincFile(fileOut_princ, intP)                                     # write data

# from Ahmad 2020
def vtkCrack(ElemList, fileName, timeStr):
        CracNodes, Tr, WL = [], [], []
        ne = 0
        Flag3D = False
        for i in ElemList:
                tr, wl = 0, 0
                if i.Type in ["C3D8S"]:
                    ne += 1
                    Flag3D = True
                    for a in range(len(i.xyzC)):
                        tr += i.wwT[a][0]   # position of normal traction in 3D = [0]
                        wl += i.wwL[a][0]   
                    tr = tr/len(i.xyzC)
                    wl = wl/len(i.xyzC)
                    for b in i.xyzC:
                        CracNodes +=  [[b[0], b[1], b[2] ]]
                        Tr += [tr]              # same traction for all nodes of the crack is assumed - Ahmad
                        WL += [wl]
                elif i.Type in ["CPS4S"]:
                    ne += 1
                    for a in range(len(i.xyC)):                             # xyC coordinates of crack endpoints
                        tr += i.wwT[0][a][1]                                # traction from CrComp; [0] 1st crack, [a] int point, [1] normal conponent
                        wl += i.wwL[0][a][1]                                # crack width
                    tr = tr/len(i.xyC)
                    wl = wl/len(i.xyC)
                    for b in i.xyC:
                        CracNodes +=  [[b[0], b[1], 0. ]]
                        Tr += [tr]              # same traction for all nodes of the crack is assumed - Ahmad
                        WL += [wl]
                        
        FFF = fileName+'_Cracks_ParaView'
        if pth.exists(FFF): pass
        else: os.mkdir(FFF)
        fN = FFF+ "/" + 'Cracks_'+timeStr+'.vtu'
        fp = open(fN,'w')
        fp.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
        fp.write('  <UnstructuredGrid>\n')
        fp.write('    <Piece NumberOfPoints='+'"'+str(len(CracNodes))+'"'+' NumberOfCells='+'"'+str(ne)+'">\n')
            #
        fp.write('      <Points>\n')
        fp.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
        for i in CracNodes:
                fp.write("%12.4e%12.4e%12.4e\n"%(i[0],i[1],i[2]))
        fp.write('      </Points>\n')
            # cell connectivity, offset, types
        fp.write('      <Cells>\n')
        fp.write('        <DataArray type="Int32" Name="connectivity" format="ascii"/>\n')
        c = 0
        for i in range(ne):
            if Flag3D:
                fp.write('%8i%8i%8i%8i'%(0+c,1+c,2+c,3+c))         # only VTK_QUAD is implimented for crack surface !!! - Ahmad
                c += 4 
                fp.write('\n') 
            else:
                fp.write('%8i%8i'%(0+c,1+c))
                c += 2
                fp.write('\n')
        fp.write('        <DataArray type="Int32" Name="offsets" format="ascii"/>\n') 
        m = 1
        for i in range(ne):
            if Flag3D:
                fp.write('%8i'%(4*m))
                m += 1
                fp.write('\n')
            else:
                fp.write('%8i'%(2*m))
                m += 1
                fp.write('\n')
        fp.write('        <DataArray type="Int32" Name="types" format="ascii"/>\n')
        for i in range(ne):
            if Flag3D: fp.write('%8i\n'%(9))
            else:      fp.write('%8i\n'%(3))
        fp.write('      </Cells>\n')
            #
        fp.write('      <PointData>\n')
        fp.write('        <DataArray type="Float64" Name="Cra_traction" NumberOfComponents=" 1" format="ascii"/>\n')
        for i in range(len(CracNodes)):
                fp.write('%12.4e'%(Tr[i]))
                fp.write('\n')
        fp.write('        <DataArray type="Float64" Name="Cra_opening" NumberOfComponents=" 1" format="ascii"/>\n')
        for i in range(len(CracNodes)):
                fp.write('%12.4e'%(WL[i]))
                fp.write('\n')
        fp.write('      </PointData>\n')
            # 
        fp.write('    </Piece>\n')
        fp.write('  </UnstructuredGrid>\n')
        fp.write('</VTKFile>\n')
        fp.close()

def vtkFile(fileName, ElemList,NodeList, timeStr, elTypesExc,noTypesInc, offset, numEl):
    # compute number of respective elements
    ne = 0
    for i in ElemList:
        if i.Type not in elTypesExc: ne += 1
    nn = 0
    for i in NodeList:
        if i.Type in noTypesInc: nn += 1
    # write data
    marker = noTypesInc[0]
    FFF = fileName+'_'+marker+'_ParaView'
    if pth.exists(FFF): pass
    else: os.mkdir(FFF)
    fN = FFF+ "/" + marker+'_'+timeStr+'.vtu'
    fp = open(fN,'w')
    fp.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
    fp.write('  <UnstructuredGrid>\n')
#    fp.write('    <Piece NumberOfPoints='+'"'+str(nn)+'"'+' NumberOfCells='+'"'+str(ne-numSDA)+'">\n')
    fp.write('    <Piece NumberOfPoints='+'"'+str(nn)+'"'+' NumberOfCells='+'"'+str(numEl)+'">\n')
    #
    # write node coordinates to vtk
    #
    fp.write('      <Points>\n')
    fp.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
    for i in NodeList:
        if i.Type in noTypesInc: fp.write("%12.4e%12.4e%12.4e\n"%(i.XCo,i.YCo,i.ZCo))
    fp.write('      </Points>\n')
    #
    # cell connectivity, offset, types to vtk
    #
    fp.write('      <Cells>\n')
    fp.write('        <DataArray type="Int32" Name="connectivity" format="ascii"/>\n')
    for i in ElemList:
        if i.Active and not i.Type in elTypesExc:
            for j in i.Inzi: fp.write('%8i'%(j-offset))    # this requires collection of nodes belonging to a common type of elements in input file !!! 
            fp.write('\n')
    fp.write('        <DataArray type="Int32" Name="offsets" format="ascii"/>\n') 
    off = 0
    for i in ElemList:
        if i.Active and not i.Type in elTypesExc:
            fp.write('%8i\n'%(len(i.Inzi)+off))
            off += len(i.Inzi)
    fp.write('        <DataArray type="Int32" Name="types" format="ascii"/>\n')
    for i in ElemList:
        if i.Active and not i.Type in elTypesExc:
            if   i.Type in ["C3D8","C3D8_EFG", "C3D8S"]: type_ = 12
            elif i.Type in ["C3D8_IGA"]:                    type_ = 25 
            elif i.Type in ["C3D4"]:                        type_ = 10
            elif i.Type in ["T3D3I"]:                       type_ = 21
            elif i.Type in ["T3D2I", "T2D2I", "T2D2", "B23I", "B23"]:    type_ = 3
            elif i.Type in ["CPS4", "CPS4S"]:            type_ = 9
            elif i.Type in ["CPS3", "CPS3S"]:            type_ = 5
            else:  raise NameError("CaeFemInOut::DataOutvtk: unknown element type 2",i.Type)
            fp.write('%8i\n'%(type_))
    fp.write('      </Cells>\n')
    #
    # point/node meta data
    #
    fp.write('      <PointData>\n')
    stressName = noTypesInc[0]+'_stress'
    if len(noTypesInc)>1: raise NameError("ConFemInOut::DataOutVTK: 1")
    if   noTypesInc[0] in ["T2DE", "T3D"]: num=1
    elif noTypesInc[0] in ["T2DE", "T3D"]: num=2  # ???????????????????????????????????????
    else:                                   num=6
    #
    # write node stresses to vtk, BufSig comes from first part of PostElemVTK
    #
    fp.write('        <DataArray type="Float64" Name="'+stressName+'" NumberOfComponents="'+str(num)+'" format="ascii"/>\n')
    BondFlag = False
    mattype = False
    for i in NodeList:
        if i.Type in noTypesInc:
            if i.Type in ["T2DE", "T3D"]:                                  # Type determined in EvNodeEl
                fp.write('%12.4e'%(i.BufSig[1]))
                BondFlag = True
            elif i.Type in ["B23I"]:
                fp.write('%12.4e%12.4e'%(i.BufSig[0], i.BufSig[1]))
                BondFlag = True
            elif i.Type in ["C2D"]:
                for j in i.BufSig: fp.write('%12.4e'%(j))
                mattype = True
            else:
                for j in i.BufSig: fp.write('%12.4e'%(j))
            fp.write('\n')
    #
    # write material to vtk
    #        
    if mattype:                                                             # for flag mattype see above
        MatNamee = noTypesInc[0]+'_Mat'
        fp.write('        <DataArray type="Float64" Name="'+MatNamee+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in NodeList:
            if i.Type in noTypesInc:
#                if i.material == "IsoDam": fp.write('%12.4e'%(i.BufSig[3])) # currently for C2D4/CPE4
#                else:                      fp.write('%12.4e'%(1))
                fp.write('%12.4e'%(i.BufSig[3]))                            # currently for C2D4/CPE4
                fp.write('\n')
    #
    #
    #
    if BondFlag:
        bondName = noTypesInc[0]+'_bond'
        fp.write('        <DataArray type="Float64" Name="'+bondName+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in NodeList:
            if i.Type in noTypesInc:
#                fp.write('%12.4e'%(i.BufBond))
                fp.write('%12.4e'%(i.BufSig[5]))
                fp.write('\n')
        
        slipName= noTypesInc[0]+'_slip'
        fp.write('        <DataArray type="Float64" Name="'+slipName+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in NodeList:
            if i.Type in noTypesInc:
#                fp.write('%12.4e'%(i.Bufslip))
                fp.write('%12.4e'%(i.BufSig[2]))
#                    if i.Label == 641:  fslip.write('%12.4e\n'%(i.Bufslip))
                fp.write('\n')
#            fslip.flush()
    #
    # write displacements to vtk
    #    
    displName = noTypesInc[0]+'_displ'
    fp.write('        <DataArray type="Float64" Name="'+displName+'" NumberOfComponents=" 3" format="ascii"/>\n')
    for i in NodeList:
        if i.Type in noTypesInc:
            if i.Type in ["B23I"]: fp.write("%12.4e%12.4e%12.4e"%(i.BufDis[0],i.BufDis[1],0))
            else:
                for j in i.BufDis: fp.write('%12.4e'%(j))
                if len(i.BufDis)==2: fp.write('%12.4e'%(0.))
            fp.write('\n')
    fp.write('      </PointData>\n')
    # 
    fp.write('    </Piece>\n')
    fp.write('  </UnstructuredGrid>\n')
    fp.write('</VTKFile>\n')
    fp.close()
    return nn

def DataFromSlices(ElData, ElDaSl):
    Data = []
    for sl in ElDaSl:                                                       # loop over sclices - maybe shis may be moved one level obove for si, sj transforming
        sl_ = sl.split(":")
        if sl_ == ['']: continue
        Data += [ElData[int(sl_[0]): int(sl_[1]) + 1]]                      # slicing range
    return Data

def WriteVTKFile( fileName, noType,timeStr, VTKNoData,VTKElData, numStress,numStressAll):
    numNo = len(VTKNoData)
    numEl = len(VTKElData)
    dir = fileName+'_'+noType+'_ParaView_'
    if pth.exists(dir): pass
    else:               os.mkdir(dir)
    fp = open( dir+ "/" + noType+'_'+timeStr+'.vtu','w')
    fp.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
    fp.write('  <UnstructuredGrid>\n')
    fp.write('    <Piece NumberOfPoints='+'"'+str(numNo)+'"'+' NumberOfCells='+'"'+str(numEl)+'">\n')
    fp.write('      <Points>\n')
    fp.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
    for i in VTKNoData:
#        coor = VTKNoData[i]["coord"]
        for j, j_ in enumerate(VTKNoData[i]["coord"]): fp.write("%12.4e"%( j_ ))
        if j==1:                                       fp.write('   0.0')   # in case of 2D
        fp.write('\n')
#        fp.write("%12.4e%12.4e%12.4e\n"%( coor[0], coor[1], coor[2]))
    fp.write('      </Points>\n')
    #
    fp.write('      <Cells>\n')
    fp.write('        <DataArray type="Int32" Name="connectivity" format="ascii"/>\n')
    for i in VTKElData:
        for j in VTKElData[i]["inzi"]: fp.write('%8i'%(j))
        fp.write('\n')
    fp.write('        <DataArray type="Int32" Name="offsets" format="ascii"/>\n')
    off = 0
    for i in VTKElData:
        delt = len(VTKElData[i]["inzi"])
        fp.write('%8i\n' % (delt + off))
        off += delt
    fp.write('        <DataArray type="Int32" Name="types" format="ascii"/>\n')
    for i in VTKElData:
        Type = VTKElData[i]["type"]
        if Type in ["C3D8", "C3D8_EFG", "C3D8S"]:               type_ = 12
        elif Type in ["C3D8_IGA"]:                              type_ = 25
        elif Type in ["C3D4"]:                                  type_ = 10
        elif Type in ["T3D3I"]:                                 type_ = 21
        elif Type in ["T3D2I", "T2D2I", "T2D2", "B23I", "B23"]: type_ = 3
        elif Type in ["CPS4", "CPS4S"]:                         type_ = 9
        elif Type in ["CPS3", "CPS3S"]:                         type_ = 5
        elif Type in ["IP"]:                                    type_ = 1
        elif Type in ["C3DCR"]:                                 type_ = 9
        else: raise NameError("ConFemInOut::WriteVTKFile: unknown element type 2", Type)
        fp.write('%8i\n' % (type_))
    fp.write('      </Cells>\n')
    #
    fp.write('      <PointData>\n')
    VTKname = noType + "_stress"
    fp.write('        <DataArray type="Float64" Name="'+VTKname+'" NumberOfComponents="'+str(numStress)+'" format="ascii"/>\n')
    for i in VTKNoData:
        if   noType in ["C3D","C3DIP"]:                                     # for element nodes and integration points
            for j in VTKNoData[i]["stress"]:
                fp.write('%12.4e'%(j))
        elif noType in ["C2D","C2DIP"]:
            for j in VTKNoData[i]["stress"]:
                fp.write('%12.4e'%(j))
        elif noType in ["T3D"]:
            fp.write('%12.4e'%(VTKNoData[i]["stress"][1]))
        elif noType in ["C3DCR"]:
            fp.write('%12.4e'%(VTKNoData[i]["stress"][0]))
        fp.write('\n')
    if noType in ["T3D"] and numStressAll==6:                               # embedded element
        VTKname = noType + "_bond"
        fp.write('        <DataArray type="Float64" Name="'+VTKname+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in VTKNoData: fp.write('%12.4e\n'% (VTKNoData[i]["stress"][5]))
        VTKname = noType + "_slip"
        fp.write('        <DataArray type="Float64" Name="' + VTKname + '" NumberOfComponents="1" format="ascii"/>\n')
        for i in VTKNoData:
            fp.write('%12.4e\n' % (VTKNoData[i]["stress"][2]))
    #
    for i in VTKNoData:                                                     # don't know valid index - returns 1st entry in dic
        vlen = len(VTKNoData[i]["vector"])                                  # list of lists
        dlen = len(VTKNoData[i]["displ"])
        break
    VTKname = noType + "_displ"
#    fp.write('        <DataArray type="Float64" Name="'+VTKname+'" NumberOfComponents=" 3" format="ascii"/>\n')
    fp.write('        <DataArray type="Float64" Name="'+VTKname+'" NumberOfComponents=" '+str(dlen)+'" format="ascii"/>\n')
    for i in VTKNoData:
        for j, j_ in enumerate(VTKNoData[i]["displ"]):
            fp.write('%12.4e'%(j_))
        if j==1: fp.write('   0.0')                                         # in case of 2D
        fp.write('\n')
    for j in range(vlen):
        VTKname = "Vector_" + str(j + 1)
        fp.write('        <DataArray type="Float64" Name="'+VTKname+'" NumberOfComponents=" '+str(vlen)+'" format="ascii"/>\n')
        for i in VTKNoData:
            for vec in VTKNoData[i]["vector"][j]:                           # j adresses a list
                fp.write('%12.4e'%(vec))
            fp.write('\n')
    fp.write('      </PointData>\n')
    fp.write('    </Piece>\n')
    fp.write('  </UnstructuredGrid>\n')
    fp.write('</VTKFile>\n')
    fp.close()

def PostElemVTK( fileName, timeStr, ElList, NodeList,NoIndToCMInd, ElResults,NodeResults):
    def Princstresses3D(sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy): # from Ahmad
        AA = np.array([[sigma_x, tau_xy, tau_xz],[tau_xy, sigma_y, tau_yz],[tau_xz, tau_yz, sigma_z]])
        eigenvalue, eigenvector = np.linalg.eig(AA)
        eigenvalue_list = [eigenvalue[0], eigenvalue[1], eigenvalue[2]]
        ma, mi = eigenvalue_list.index(max(eigenvalue)), eigenvalue_list.index(min(eigenvalue))
        v1 = np.array(np.multiply(eigenvalue[ma], eigenvector[:, ma]))
        v3 = np.array(np.multiply(eigenvalue[mi], eigenvector[:, mi]))
        for i in range(3):
            if i != ma and i != mi:
                v2 = np.array(np.multiply(eigenvalue[i], eigenvector[:, i]))
                sig_2 = eigenvalue[i]
        sig_1, sig_3 = eigenvalue[ma], eigenvalue[mi]
        return sig_1, sig_2, sig_3, v1, v2, v3
    def Princstresses2D(xx, yy, xy):  # calculation of principal stresses and vectors in 2D
        delta = sqrt((-xx - yy) ** 2 - 4 * (-xy ** 2 + xx * yy))
        b = -xx - yy
        s1 = (-b + delta) / 2
        s2 = (-b - delta) / 2
        if abs(s1 - xx) < ZeroD:  # if first principle equals to sig xx
            vctr_1 = [1., 0., 0.]
            vctr_2 = [0., 1., 0.]
        elif abs(s1 - yy) < ZeroD:  # if first principle equals to sig yy
            vctr_1 = [0., 1., 0.]
            vctr_2 = [1., 0., 0.]
        else:
            a1 = (s1 - xx) / xy
            l1 = sqrt(1 ** 2 + a1 ** 2)
            vctr_1 = [1 / l1, a1 / l1, 0.]
            a2 = (s2 - xx) / xy
            l2 = sqrt(1 ** 2 + a2 ** 2)
            vctr_2 = [1 / l2, a2 / l2, 0.]
        # if s1 < 0: s1 = 0.                  # only tension
        # if s2 > 0: s2 = 0.                  # only compression
        return s1, [s1 * vctr_1[0], s1 * vctr_1[1], s1 * vctr_1[2]], s2, [s2 * vctr_2[0], s2 * vctr_2[1], s2 * vctr_2[2]]
    def DataProcessor( Res, Slice, DatLen):
        VTKNoData,VTKElData, VTKIpData,VTKElIpInz, VTKCrData,VTKCrInzi = {},{}, {},{}, {},{}        # VTKElIpInz is dummy for "inzidence" of integration points
        numEl, numNo, numIp, numVert, numCrInz = 0, 0, 0, 0, 0              # numCrInz initialize connectivity for crack vertices - VTK label for vertices
        Labl = Res["Label"]
        Data = Res["Data"]
        Slic = Res["Slices"]
        IpIn = Res["IntPointIndex"]
        IpCo = Res["IntPointCoordinates"]
        for Labl, Data, Slices, IpIn, IpCo in zip(Labl, Data, Slic, IpIn, IpCo):
            el  = ElList[FindIndexByLabel(ElList, Labl)]
            elT = el.Type
            if IpIn.strip() in ["Cra","Cra2"]:
                if numVert not in VTKCrData:
                    VTKCrData[numVert] = {}
                if Labl not in VTKCrInzi:
                    VTKCrInzi[Labl] = {}
                    VTKCrInzi[Labl]["inzi"]   = []
                if elT in ["C3D8", "C3D8S"]:
                    VTKCrData[numVert]["coord"]  = [ IpCo[0], IpCo[1], IpCo[2]]
                    VTKCrData[numVert]["displ"]  = [Data[0]]
                    VTKCrData[numVert]["stress"] = [Data[1]]
                    VTKCrData[numVert]["vector"] = []                       # dummy
                    VTKCrInzi[Labl]["inzi"]  += [numCrInz]
                    VTKCrInzi[Labl]["type"]   = "C3DCR"
                else:
                    pass
                numCrInz += 1
                numVert +=1
            else:
                IpI = int(IpIn)
                # for whole element
                if IpI==0:
                    elNoCoor = el.NodalCoordinates( NodeList, NoIndToCMInd) # nodal coordinates
                    IntLen = el.nIntL                                       # number of integration points of element
                    if elT in ["C3D8","C3D8S"]:
                        if IntLen != 8: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 3")
                    elif elT in ["CPS4","CPE4"]:
                        if IntLen != 4: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 4")
                    elif elT in ["CPS4S", "CPE4S"]:
                        if IntLen != 1: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 4")
                    elif elT in ["T3D3I","T3D2I"]:
                        if IntLen != 2: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 2")
                        if elT=="T3D3I":
                            elNoCoor = elNoCoor[0:7]                        # discard dummy coordinates of enhancing node as enhancing node has only one dof
                            elNoCoor[6] = 0.5*(elNoCoor[0]+elNoCoor[3])     # enhanced node may have dummy coordinates
                            no0 = NodeList[NoIndToCMInd[el.Inzi[0]]]
                            no1 = NodeList[NoIndToCMInd[el.Inzi[1]]]
                            no2 = NodeList[NoIndToCMInd[el.Inzi[2]]]
                            no2.XCo = 0.5*(no0.XCo+no1.XCo)
                            no2.YCo = 0.5*(no0.YCo+no1.YCo)
                            no2.ZCo = 0.5*(no0.ZCo+no1.ZCo)
                    # element node data part 1
                    elNoDisp = []
                    for ni in el.Inzi:
                        no = NodeList[NoIndToCMInd[ni]]
                        for u in no.BufDis: elNoDisp += [u]                 # must be here - otherwise displacements of element nodes might not be complete
                        if ni not in VTKNoData:
                            VTKNoData[ni] = {}
                            VTKNoData[ni]["VTKid"] = numNo                  # reaasign node-id for VTK as only subsets of nodes may be indluded in actual VTK data
                            VTKNoData[ni]["coord"] = [ no.XCo, no.YCo, no.ZCo]
                            VTKNoData[ni]["stress"] = zeros((DatLen),dtype=float)    # to collect sigma-data
                            VTKNoData[ni]["bufCo"]  = 0                     # counter for buffer
                            VTKNoData[ni]["vector"] = []                    # dummy
                            numNo +=1
                    Buf, X = zeros((DatLen,IntLen), dtype=float), []
                # element integration point data
                r = SamplePoints[el.IntT,el.nInt-1, IpI][0]
                s = SamplePoints[el.IntT,el.nInt-1, IpI][1]
                t = SamplePoints[el.IntT,el.nInt-1, IpI][2]
                NN = el.FormN( r, s, t)
                VTKIpData[numIp], VTKElIpInz[numIp] = {}, {}
                VTKIpData[numIp]["coord"] = dot( NN, elNoCoor)
                VTKIpData[numIp]["displ"] = dot( NN, elNoDisp)
                Data_ = DataFromSlices(Data, Slice)[0]                      # integration point all stress data
                if DatLen != len(Data_): raise NameError("ConFemPostProcStandAalone::PostElemVTK: 4", DatLen, len(Data_), Data_)
                if el.Type in ["C3D8", "C3D8S"]:
                    sig_1, sig_2, sig_3, v1, v2, v3 = Princstresses3D( Data_[0], Data_[1], Data_[2], Data_[3], Data_[4], Data_[5])
                    VTKIpData[numIp]["stress"] = [ sig_1, sig_2, sig_3]
                    VTKIpData[numIp]["vector"] = [ v1, v2, v3]
                elif el.Type in ["CPS4","CPE4","CPS4S","CPE4S"]:
                    sig_1, v1, sig_2, v2 =Princstresses2D( Data_[0], Data_[1], Data_[3])
                    VTKIpData[numIp]["stress"] = [ sig_1, sig_2, 1.]
                    VTKIpData[numIp]["vector"] = [ v1, v2]
                else:
                    pass
                VTKElIpInz[numIp]["inzi"] = [ numIp ]
                VTKElIpInz[numIp]["type"] = "IP"
                numIp += 1
                for i in range(DatLen): Buf[ i, IpI] += Data_[i]            # reorder rows of Data in columns of Buf so that rows of Buf hold values of same type for all integration points
                # for whole element
                if IpI==IntLen-1:
                    VTKElData[numEl] = {}
                    VTKElData[numEl]["type"] = el.Type
                    VTKElData[numEl]["inzi"] = []
                    for i in el.Inzi:                                       # i is global index
                        VTKElData[numEl]["inzi"] += [ VTKNoData[i]["VTKid"] ]   # numEl, VTKNoData[i]["VTKid"] are local - should not be confused with VTKElIpInz[numIp]["inzi"] above
                    # element node data part 2
                    for k, i in enumerate(el.Inzi):                         # loop over element nodes
                        for j in range(DatLen):                             # loop over value types
                            Z = dot( el.FormNI(), Buf[j])                   # Buf[j] has values of same type for all ip --> extrapolate to nodes with 'inverse' of shape function
                            VTKNoData[i]["stress"][j] += Z[k]               # [j] is index for value type; [k], [i] are indices for nodes
                        VTKNoData[i]["bufCo"] += 1
                    if el.Type in ["T3D3I"]:                                # special treatment for T3D3I
                        no2.BufDis = dot( el.FormN(0.,0.,0.), elNoDisp)     # no2 from above, redefines BufDis
                        VTKElData[numEl]["inzi"] = [ VTKNoData[el.Inzi[0]]["VTKid"], VTKNoData[el.Inzi[2]]["VTKid"], VTKNoData[el.Inzi[1]]["VTKid"] ]
                    numEl += 1
                # end of integration point
            # end of integration point / crack loop
        # for whole system
        for ni in VTKNoData:                                                # ni is internal global node index
            no = NodeList[NoIndToCMInd[ni]]
            VTKNoData[ni]["displ"] = no.BufDis
            VTKNoData[ni]["stress"] = VTKNoData[ni]["stress"] / VTKNoData[ni]["bufCo"]
        #
        return VTKNoData,VTKElData, VTKIpData,VTKElIpInz, VTKCrData,VTKCrInzi
    # end def
    for i in NodeList:                                                      # initialize and fill displacement buffer
#        i.BufSig = zeros((6),dtype=float)                                   # buffer for nodal stress values used for VTK - be cautious with its length, must include all components
#        i.BufBond, i.Bufslip = 0., 0.
#        i.BufBondCount = 0
#        i.BufCount = 0
#        i.material = ""
#        i.BufDis = zeros((3), dtype = float)
        try:
            i.BufDis = NodeResults[i.Label][2]                          # displacement buffer
        except: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 5")
    #
    Res = ElResults[ElResults["Type"].isin(["C3D8", "C3D8S"])]
    if len(Res) > 0:
        VTKNoData,VTKElData,VTKIpData,VTKElIpInz, VTKCrData,VTKCrInzi = DataProcessor( Res, ["6:11"], 6 ) # starts with 0 index after integration point
        WriteVTKFile(fileName, "C3D",   timeStr, VTKNoData,VTKElData, 6,6)
        WriteVTKFile(fileName, "C3DIP", timeStr, VTKIpData,VTKElIpInz, 3,3)
        if len(VTKCrData)>0: WriteVTKFile(fileName, "C3DCR", timeStr, VTKCrData,VTKCrInzi, 1, 1)
#        vtkCrack(ElList, fileName, timeStr)
    #
    Res = ElResults[ElResults["Type"].isin(["CPS4","CPE4","CPS4S","CPE4S"])]
    if len(Res)>0:
        VTKNoData,VTKElData,VTKIpData,VTKElIpInz, _, _ = DataProcessor( Res, ["4:7"], 4 )
        WriteVTKFile(fileName, "C2D",   timeStr, VTKNoData,VTKElData, 4, 4)
        WriteVTKFile(fileName, "C2DIP", timeStr, VTKIpData,VTKElIpInz, 3, 3)
#        vtkCrack( ElList, fileName, timeStr)
    #
    Res = ElResults[ElResults["Type"].isin(["T3D3I","T3D2I"])]
    if len(Res)>0:
        DatLen = 6  # corresponds to value types - actually strain, stress for T3D3I, 6 bond items: bond item 1 for long slip, item 4 for long stress, other for lateral which are not considered
        VTKNoData,VTKElData, _, _, _, _ = DataProcessor( Res, ["0:5"], DatLen )
        WriteVTKFile(fileName, "T3D", timeStr, VTKNoData,VTKElData, 1,6)

def PostElem1D_( f2, ResultTypes, PlotTimes ):  # ResultTypes defined in ConFemInOut considering element and material type
    DataDict, elsetType = {}, {}
    # fill dictionaries
    with ( f2 ) as ff:
        for z1 in ff:
            z2 = z1.split(',')
            # new time
            if z2[0].strip()=="Time":
                Time = float(z2[1])
                if Time in PlotTimes or len(PlotTimes)==0: PlotTFlag = True
                else:                                      PlotTFlag = False
                continue                                            # next line
            elif PlotTFlag:
                ElSet  = z2[1].strip()
#                ElMat  = z2[2].strip()
                ElType = z2[3].split(":")[0].strip()
                elsetType[ElSet] = ElType
                if ElType in ['CPS4','CPS4R','CPE4','CPE4R', 'CPS3','CPE3','CPS4S','CPE4S','CPS3S','CPE3S','CAX4',
                              'SB3','SH4','SH3',
                              'T2D2','T2D2I','T2D3I','T3D2','T3D3','T3D2I','T3D3I',
                              "Bond2D2","Bond2D3","Bond3D2","Bond3D3","BondAX2","BondAX3",
                              'C3D8','C3D8S']:
                    continue
                else:
                    key = ElSet #+ "*" + ElMat + "*" + ElType
                    if key not in DataDict:
                        DataDict[key] = {}
                    if Time not in DataDict[key]:
                        DataDict[key][Time] = []
                    data = [ int(z2[4]), float(z2[5]), float(z2[6])]
                    dataLen = len(ResultTypes[ElSet])                       # number of result items without coordinates
                    for i in range(7,7+dataLen): data += [ float(z2[i])]    # first 7 items are elnumber, elset, mat, eltype, ip, x, y
                    DataDict[key][Time] += [ data ]
    # plot
    # no rotation in this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for elset in DataDict:
        ElType = elsetType[elset]
        pp, YminAll, YmaxAll = [], [], []
        for r in ResultTypes[elset]:
            P0 = plt.figure()
            p0 = P0.add_subplot(111)
            p0.set_title( elset + ': ' + r, fontsize=FonSizTi)
            p0.tick_params(axis='x', labelsize=FonSizAx)
            p0.tick_params(axis='y', labelsize=FonSizAx)
            p0.grid()
            pp += [ p0 ]
            YminAll += [ 0 ]
            YmaxAll += [ 0 ]
        for j, time in enumerate(DataDict[elset]):
            AllData = DataDict[elset][time]
            offset = 3                                                      # ip, x, y
            Colr = Colors[j % len(Colors)]
            for i, r in enumerate(ResultTypes[elset]):
#                Ymin, Ymax = 0,0
                X, Y, p0 = [], [], pp[i]
                for d in AllData:
                    factor = 1.0
                    if False:
                        if ElType in ["BAX21E"]:
                            if r=="normal force":
                                if d[0]!= 1:                                # integration point
                                    continue
                                else:
                                    factor = 1. / 2.46                              # quick and dirty for special case 1-M0-25-1.23.elemout to have stresses
                    X += [ d[1] ]
                    Y += [ factor*d[offset+i] ]
                p0.plot(X, Y, '-', color=Colr, linewidth=LiWiCrv)
                p0.plot(X, Y, 'o', color='tab:red', ms=5.)
                if min(Y) < YminAll[i]: YminAll[i] = min(Y)                 # i is for result types
                if max(Y) > YmaxAll[i]: YmaxAll[i] = max(Y)
#                Ymin = min(Ymin, min(Y))
#                Ymax = max(Ymax, max(Y))
#                if abs(Ymax-Ymin)>ZeroD: p0.set_ylim(1.1*Ymin,1.1*Ymax)
        for i, r in enumerate(ResultTypes[elset]):
            p0 = pp[i]
            Ymin = YminAll[i]
            Ymax = YmaxAll[i]
            if abs(Ymax-Ymin)>ZeroD: p0.set_ylim(1.1*Ymin,1.1*Ymax)

# uses data from *.elemout.txt
def PostElem1D( f2, ResultTypes, PlotTimes ):                       # ResultTypes defined in ConFemInOut considering element and material type
    def Rotate( dataXe, dataYe):
        n = len(dataXe)
        if n>1:
            LL = ((dataXe[n-1]-dataXe[0])**2+(dataYe[n-1]-dataYe[0])**2)**0.5
            cosA = (dataXe[n-1]-dataXe[0])/LL
            sinA = (dataYe[n-1]-dataYe[0])/LL
            return cosA, sinA
        else: return 1, 0
    dataXe = []                                                     # list for ip x per element
    dataYe = []                                                     # list for ip y per element
    dataAe = []                                                     # list for ip data per element
    dataXt = []                                                     # list for x per time step
    dataYt = []                                                     # list for y per time step
    dataAt = []                                                     # list for A per time step
    dataSt = []                                                     # list for sets per time step
    dataX  = []                                                     # list for all x
    dataY  = []                                                     # list for all y
    dataA  = []                                                     # list for all A
    dataS  = []                                                     # list for all sets
    Times  = []                                                     # list for time values
#    FlagData = False
    FlagElAll = False                                                       # to indicate whether at least one 1D element was hit
    # read all data from file
    z1 = f2.readline()
    z2 = z1.split()
    while z1!="":
        # new time
        if z2[0]=="Time":
            Time = float(z2[1])
            if Time in PlotTimes or len(PlotTimes)==0: PlotTFlag = True
            else:                                      PlotTFlag = False
            if PlotTFlag:
                ElSetI = 0
                ElSets = {}                                                 # dictionary: element set label -> index
                Times += [Time]
                if len(dataXe)>0:
                    dataXt += [dataXe]                                      # ip x coordinates: append last element of previous time
                    dataYt += [dataYe]                                      # ip y coordinates: append last element
                    dataAt += [dataAe]                                      # result data:   append last element of previous time
                    dataX  += [dataXt]                                      # append to list of all x
                    dataY  += [dataYt]                                      # append to list of all y
                    dataA  += [dataAt]                                      # append to list of all data
                    dataS  += [dataSt]                                      # 
                    dataXt, dataYt, dataAt, dataXe, dataYe, dataAe = [], [], [], [], [], [] # reset lists for elements
        # new element
        elif PlotTFlag and z2[0]=="El":
            ElType = z2[2]
            if ElType in ['CPS4','CPS4R','CPE4','CPE4R', 'CPS3','CPE3','CPS4S','CPE4S','CPS3S','CPE3S','CAX4',
                          'SB3','SH4','SH3',
#                          'T2D2','T2D3','T2D2I','T2D3I','T3D2','T3D3','T3D2I','T3D3I',
                          'T2D2','T2D2I','T2D3I','T3D2','T3D3','T3D2I','T3D3I',
                          "Bond2D2","Bond2D3","Bond3D2","Bond3D3","BondAX2","BondAX3",
                          'C3D8','C3D8S']:
                FlagEl = False 
            else:
                FlagElAll = True
                FlagEl = True
                ElSet  = z2[3]
                if ElSet not in ElSets:                                     # new elset
                    ElSetI = ElSetI+1                                       # count elsets for later plot loop
                    ElSets[ElSet] = ElSetI                                  # index for new elset
                dataSt += [ElSetI]                                          # assign elset index to current element data
                if len(dataXe)>0:
                    dataXt += [dataXe]                                      # append to list of elements per time step
                    dataYt += [dataYe]                                      # append to list of elements per time step
                    dataAt += [dataAe]                                      # append to list of data
                    dataXe = []                                             # reset list for element
                    dataYe = []                                             # reset list for element
                    dataAe = []                                             # reset list for element
        # integration point data
        elif PlotTFlag and FlagEl:                                          # ip data line
            dataXe += [float(z2[0])]                                        # append ip x-coordinate to element list
            dataYe += [float(z2[1])]                                        # append ip y-coordinate to element list
            dataAe += [[float(z2[pos]) for pos in range(2,len(z2))]]        # append ip data to element list
        z1 = f2.readline()
        z2 = z1.split()
        # data reading finished
    #
    if FlagElAll: 
        dataXt += [dataXe]                                                  # append last element
        dataYt += [dataYe]                                                  # append last element
        dataAt += [dataAe]                                                  # append last element
        dataX  += [dataXt]                                                  # append to list of all
        dataY  += [dataYt]                                                  # append to list of all
        dataA  += [dataAt]                                                  # append to list of all
        dataS  += [dataSt]
#        if len(dataX)==0 or not FlagData: return 0
        if len(dataX)==0: return 0
    
        # plot whole stuff
        scal = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
    #    if raw_input('Separate element sets for plot? (y/n):')=='y':
        if len(ElSets)>0: SFlag=True
        else:             SFlag=False                                           # does not separate for elsets
        for i in range(len(Times)): print("PostElem1D process plot for time ", Times[i],Colors[i % len(Colors)]) # pre user info
        #
        for sI, ElS in zip( list(ElSets.values()) , list(ElSets.keys()) ):      # loop over element sets; value is elset index/counter, key is elset label
            titleList = ResultTypes[ElS]                                        # ResultTypes set per elset in ConFemInOut::DataInput post
            for I, tL in enumerate(titleList):                                  # loop over result types
                Ymin, Ymax = 0,0
                if tL == None: continue
                P0 = plt.figure()
                p0 = P0.add_subplot(111)
                if SFlag: p0.set_title(ElS+': '+tL, fontsize=FonSizTi)          # ElS is elset label
                else:     p0.title(tL, fontsize=FonSizTi)
                p0.tick_params(axis='x', labelsize=FonSizAx)
                p0.tick_params(axis='y', labelsize=FonSizAx)
                p0.grid()
                #
                for i in range(len(Times)):                                     # loop over time steps
                    Colr = Colors[i % len(Colors)]
                    X_, Y_ = [], []
                    dA, dX, dY, dS = dataA[i], dataX[i], dataY[i], dataS[i]     # element data, p x-coordinate, ip y-coordinate, elset-label for time i
                    ne = len(dA)                                                # number of elements
                    PlotThrough = False
                    #
                    for j in range(ne):                                         # loop over elements
                        if dS[j] == sI or not SFlag:                            # found element in current elset to plot  or   do not separate for elsets
                            neD = len(dA[j])                                    # number of data items per element
                            if neD==1: PlotThrough = True                       # only one data item per element
                            cosA, sinA = Rotate( dX[j], dY[j])                  # ip coordinates
                            X                   = [ dX[j][k] + sinA*(scal[I]*dA[j][k][I]) for k in range(neD)] # shifted origin and and rotation applied to x = 0
                            if PlotThrough: X_ += [ dX[j][0] + sinA*(scal[I]*dA[j][0][I]) ]
                            if fabs(sinA)>ZeroD:
                                p0.plot( dX[j], dY[j], '-k')                    # integration point coordinates  undeformed structure connected by line
                                Y =                  [ dY[j][k] + cosA*(scal[I]*dA[j][k][I]) for k in range(neD)]
                                if PlotThrough: Y_+= [ dY[j][0] + cosA*(scal[I]*dA[j][0][I])]
                            else: 
                                Y =  [                       (scal[I]*dataA[i][j][k][I]) for k in range(neD)]
                                if PlotThrough: Y_+= [       (scal[I]*dataA[i][j][0][I])]
                            Ymin=min(Ymin,min(Y))
                            Ymax=max(Ymax,max(Y))
                            p0.plot( X, Y, '-', color=Colr,linewidth=LiWiCrv)
                            p0.plot( X, Y, 'o', color='tab:red',ms=5.) 
                    if PlotThrough: p0.plot( X_, Y_, '-', color=Colr,linewidth=LiWiCrvThin)
                if abs(Ymax-Ymin)>ZeroD: p0.set_ylim(1.1*Ymin,1.1*Ymax)
                else:                    P0.set_visible(False)
                P0.autofmt_xdate()
    return 0

def PostScales( SecDic, ElemList, ElResults, ScaleStress):
    Scale = {}
    # loop over element sets
    for S in SecDic:
        ScaleElSet = 1.
        # take scaling factors from options, if there are any
        for sc in ScaleStress:                                              # from options
            if sc[0].upper()==S:
                ScaleElSet = sc[1]                                          # superposed scaling factor for element set
                break
        # loop over elements of Solid Section for mean characteristic element length
        nEl, meanCharL = 0, 0.
        for e in SecDic[S].Elems:
            i = ElemList[e]
#            if i.Type in ["SH3","SH4"]: Lis = i.Lists1()
#            else:                       Lis = range(i.nIntL)
            meanCharL += (i.Lch_/sqrt(float(i.nIntL)))
            nEl       += 1
        meanCharL = meanCharL/nEl
        # extract element set data from dataframe
        SetResults = ElResults[ElResults["ElSet"] == S]
        SetResultsData   = SetResults["Data"]
        SetResultsSlices = SetResults["Slices"]
        SetResultsNoEl   = SetResults["NoEl"]
        if SetResultsNoEl.iloc[0]: sce = nEl                                # number of elements in set
        else:                      sce = 1
        nn = SetResultsData.shape[0]                                        # number of rows
        Sc_ = []
        # loop over slices
        for s_ in SetResultsSlices.iloc[0]:                                 # extracts first row of sub-dataframe to determine slices for the following
            sl = s_.split(":")
            if sl == ['']:
                continue
            si, sj = int(sl[0]), int(sl[1])                                 # slicing range
            if len(sl)!= 2 or sj<si: raise NameError("ConFemPost::PostScales: definition of slices uncorrect",sl)
            Data = zeros((nn,sj-si+1), dtype=float)
            for i in range(nn):                                             # shuffle data to numpy array
                Data[i] = SetResultsData.iloc[i][si:sj+1]
            # determine min, max, scaling
            sMin = amin(Data)
            sMax = amax(Data)
            sAbs = max(sMax, -sMin)
            if sAbs>ZeroD: sf = sce*ScaleElSet*meanCharL/sAbs
            else:          sf = 1.
            if sf > 1.e6: sf = 1.
            Sc_ += [ sf ]
        Scale[S] = Sc_
    return Scale

def PloPrin(pp, scaleP, XX, YY, p0, MaPlLib,ffred,ffgre ):
    SS = pp[0]                                                              # larger principal stress - signed
    if SS >= 0:
        if MaPlLib:
            p0.plot([XX - SS * scaleP * pp[1], XX + SS * scaleP * pp[1]],
                    [YY - SS * scaleP * pp[2], YY + SS * scaleP * pp[2]], 'r-')  # plot red from point to point
        else:
            ffred.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n" % (
            XX - SS * scaleP * pp[1], YY - SS * scaleP * pp[2], XX + SS * scaleP * pp[1], YY + SS * scaleP * pp[2]))
    else:
        if MaPlLib:
            p0.plot([XX - SS * scaleP * pp[1], XX + SS * scaleP * pp[1]],
                    [YY - SS * scaleP * pp[2], YY + SS * scaleP * pp[2]], 'g-')  # same as above but green
        else:
            ffgre.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n" % (
            XX - SS * scaleP * pp[1], YY - SS * scaleP * pp[2], XX + SS * scaleP * pp[1], YY + SS * scaleP * pp[2]))
    S2 = pp[3]                                                              # smaller principal stress - signed
    if S2 >= 0:
        if MaPlLib:
            p0.plot([XX - S2 * scaleP * pp[2], XX + S2 * scaleP * pp[2]],
                    [YY + S2 * scaleP * pp[1], YY - S2 * scaleP * pp[1]], 'r-')  # plot red from point to point
        else:
            ffred.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n" % (
            XX - S2 * scaleP * pp[2], YY + S2 * scaleP * pp[1], XX + S2 * scaleP * pp[2], YY - S2 * scaleP * pp[1]))
    else:
        if MaPlLib:
            p0.plot([XX - S2 * scaleP * pp[2], XX + S2 * scaleP * pp[2]],
                    [YY + S2 * scaleP * pp[1], YY - S2 * scaleP * pp[1]], 'g-')  # same as above but green
        else:
            ffgre.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n" % (
            XX - S2 * scaleP * pp[2], YY + S2 * scaleP * pp[1], XX + S2 * scaleP * pp[2], YY - S2 * scaleP * pp[1]))
    return SS, S2

def PlotRotLinesText_( S, el,j,xyP,angle, sign_,sig,ScaleSec, p1,pDaX,pDaY, cosA,sinA):
    lsc, lsi = len(ScaleSec), len(sig)
    if lsc != lsi: raise NameError("ConFemPy3::PostElem2D: number of scaling factors does not match for this section",S,ScaleSec,sig)
    phiC, phiS = cosA, sinA
    pX, pY, ci = [], [], 0
    for jj in range(lsc):
        for s_ in sig[jj]:
            dX = xyP[0] - sign_[jj]*ScaleSec[jj]*s_ * phiS
            dY = xyP[1] + sign_[jj]*ScaleSec[jj]*s_ * phiC
            pX += [ dX ]
            pY += [ dY ]
            p1.text( dX,dY,format(s_,".2f"),ha='left',va='top',rotation=angle,color=Colors2[ci],fontsize='medium')
            ci += 1
    pDaX[j] = pX                                                            # j is integration point index
    pDaY[j] = pY

def PostElem2D( Header,ElemList,NodeList, ScaleU,Contour2D, SolSecs, MaPlLib, NoIndToCMInd, NodeResults,ElResults,SetSlices, TimeEl, FilDir, SCales ):
    # start def
    def Plot2DCrackContour(CrXY0, CrXY1):
        r0, s0, t0 = el.NLocCoorFromGlobCoor(xy, CrXY0)  # local coordinate 'left' edge of crack
        r1, s1, t1 = el.NLocCoorFromGlobCoor(xy, CrXY1)  # local coordinate 'right' edge of crack
        CrdXY0 = dot(el.FormN(r0, s0, t0), array(uxy))  # displacement 'left' edge of crack
        CrdXY1 = dot(el.FormN(r1, s1, t1), array(uxy))  # displacement 'right' edge of crack
        xx = [CrXY0[0] + CrdXY0[0], CrXY1[0] + CrdXY1[0]]
        yy = [CrXY0[1] + CrdXY0[1], CrXY1[1] + CrdXY1[1]]
        if MaPlLib:
            p1.plot([xx[0], xx[1]], [yy[0], yy[1]], color='gray')
        else:
            ffc.write("%8.4f  %8.4f\n%8.4f  %8.4f\n\n" % (xx[0], yy[0], xx[1], yy[1]))
    # end of def
    if not MaPlLib: import PyGnuplot as gp                                       # use gnuplot instead of matplotlib
    # store values of ndof per node -- if NodeResults exists use data from nodeout files, otherwise from actual displacement vector
    DisplList = NodalDisplList( NodeList, NodeResults)

#    # echo bond elements
#    if False:
#        FF = 0.
#        for i in ElemList:
#            if i.Type in ["Bond2D2","Bond2D3"]:
#                uElC, uElT, uElA = [], [], []
#                eC   = ElemList[i.iC]
#                if eC.Type=="CPS4_IGA": xyC = eC.ElCPLi                           # geometry control point coordinates 
#                else:                    xyC = eC.NodalCoordinates( NodeList, NoIndToCMInd )
#                XXC  = eC.FormX( i.rsC[0], i.rsC[1], 0.)                               
#                xyPC = dot(XXC,array(xyC))                                          # undeformed global coordinates of integration point
#                for j in eC.Inzi:                                                   # loop over nodes affected by element -- equivalent to lines before
#                    j_ = NoIndToCMInd[j]                                                   
#                    for k in range(len(DisplList[j_])):
#                        uElC += [DisplList[j_][k]]                                   # collects nodal displacements per element
#                val = None
#                if eC.Type=="CPS4_EFG":
#                    NNC = eC.ArbFormN( NodeList,NoIndToCMInd, CoorTree_,NoCoLitoNoLi, LevSetList, i.rsC[0], i.rsC[1], xyPC)
#                else:
#                    if eC.Type=="CPS4_XFEM" and i.DofE>8: val, valR, valT = eC.LevSetVal( i.rsC[0], i.rsC[1], LevSetList)
#                    NNC = eC.FormN( i.rsC[0], i.rsC[1], 0., val)
#                uC  = dot(NNC,uElC)
#                uC_ = dot(i.Rot,uC)
#                #
#                eT   = ElemList[i.iT]
#                xyT  = eT.NodalCoordinates( NodeList, NoIndToCMInd )                    # undeformed nodal coordinates
#                XXT  = eT.FormX( i.rsT[0], 0., 0.)                               
#                xyPT = dot(XXT,array(xyT))                                          # undeformed global coordinates of integration point
#                for j in eT.Inzi:                                                   # loop over nodes affected by element -- equivalent to lines before
#                    j_ = NoIndToCMInd[j]                                                   
#                    for k in range(len(DisplList[j_])):
#                        uElT += [DisplList[j_][k]]                                  # collects nodal displacements per element
#                NNT = eT.FormN( i.rsT[0], 0., 0., val)
#                uT  = dot(NNT,uElT)
#                uT_ = dot(i.Rot,uT)
#                for j in i.Inzi:
#                    j_ = NoIndToCMInd[j]                                                   
#                    for k in range(len(DisplList[j_])):
#                        uElA += [DisplList[j_][k]]                                  # collects nodal displacements per element
#                rL = zeros((i.DofE), dtype=float)
#                for j in range(i.nIntL):
#                    rst = SamplePoints[i.IntT,i.nInt,j]                             # local integration sampling point coordinates
#                    BB, JJ = i.FormB( rst[0], rst[1], rst[2], val )
#                    sl = dot(BB,uElA)
#                    tt, MatM, X = MatList[i.MatIp[j]].MatC( None, 1, 1., j, i, None, sl, None )                 # material for integration point
#                    tt = dot(MatM,sl)
#                    f = JJ*i.Geom*SampleWeight[i.IntT,i.nInt,j]
#                    rL = rL + f*dot(transpose(BB),tt)
#                    FF += MatM[0,0]*(uT_[0]-uC_[0])*f
#                    print('XXX', j, xyPC, xyPC-xyPT, '__', norm(uT_-uC_- sl), sl)

    ContourSets, ContourX, ContourY, ContourF, ContourM, ContourP = [], {}, {}, {}, {}, {}
#    pDaX, pDaY = {}, {}                                 # for permanent storage of line plot data

    # loop over all element sections
    for s_, S in enumerate(SolSecs):
        sec = SolSecs[S]
        if sec.ElemTypes in __BeamsBern__ or sec.ElemTypes in __BeamsTimo__ \
                or sec.ElemTypes in __SpringsAll__ or sec.ElemTypes in __Length1D__ or sec.ElemTypes in [["C3D8"]]:
#        if sec.ElemTypes in [["S1D2"],["T1D2"],["B23"],["B23E"],["B21"],["B21E"],["TAX2"],['TAX3'],["BAX21E"],["BAX21"]]:
            continue                                                        # ??? this routine is not appropriate for this element types
        print('PostElem2D process plot for section', S, sec.ElemTypes, "time",TimeEl)
        try:
            ScaleSec =SCales[S]
        except:
            ScaleSec = [0]
            print("PostElem2D scaling factors for this section not provided, used 0",S)
        # matplotlib
        if MaPlLib:                                                         # init for principal value plots
            if SolSecs[S].ElemTypes in [["SH4"],["SH3"]]:
                _, p1 = DefPlot("normal forces "  + SolSecs[S].Set+ " time "+"%6.3f"%TimeEl) #+ " scale "+"%7.4f"%ScaleSec[0])
                _, p2 = DefPlot("bending moments "+ SolSecs[S].Set+ " time "+"%6.3f"%TimeEl) #+ " scale "+"%7.4f"%ScaleSec[1])
                p2.axis('equal')
            elif SolSecs[S].ElemTypes in [["SB3"]]:
                _, p1 = DefPlot("bending moments "+ SolSecs[S].Set+ " time "+"%6.3f"%TimeEl) #+ " scale "+"%7.4f"%ScaleSec[0])
                scaleQ = 2.5*ScaleSec[1]
                _, p2 = DefPlot("shear forces "   + SolSecs[S].Set+ " time "+"%6.3f"%TimeEl) #+ " scale "+"%7.4f"%scaleQ)
                p2.axis('equal')
            else:
                try:    head = Header[0]
                except: head = " "
                _, p1 = DefPlot(head+" - "+SolSecs[S].Set+ " - time "+"%6.3f"%TimeEl) #+ " scale "+"%7.4f"%ScaleSec[0])
            p1.axis('equal')
            ffred = None
            ffgre = None
        # gnuplot
        else:
            p1 = None
            dfg   = FilDir+"tmp"+str(S)+".dat"                              # temporary file for gnuplot elements
            dfc   = FilDir+"tmp"+str(S)+"c.dat"                             # temporary file for gnuplot cracks
            dfred = FilDir+"tmp"+str(S)+"red.dat"                           # temporary file for gnuplot principal stress tensile
            dfgre = FilDir+"tmp"+str(S)+"gre.dat"                           # temporary file for gnuplot principal stress compressive
            dfIP  = FilDir+"tmp"+str(S)+"IP.dat"                            # temporary file for gnuplot integration points + color for stress level
            dfIP2 = FilDir+"tmp"+str(S)+"IP2.dat"                           # temporary file for gnuplot integration points 2 + color for stress level
            ffg=open(dfg,'w'); ffc=open(dfc,'w'); ffred=open(dfred,'w'); ffgre=open( dfgre,'w'); ffIP=open(dfIP,'w'); ffIP2=open(dfIP2,'w')
            Truss2DCoor, Truss2DStr, Truss2DStrBond, = [], [], []           # for bar/beam elements
            gi = 0
        xxMin, xxMax, yyMin, yyMax = 999., -999., 999., -999.
        #
        NorMinMax = ValMinMax(1, 3)                                         # to indicate extremal values in plot -- 1 int value, 3 float values
        MomMinMax = ValMinMax(1, 3)                                         #
        # initializations for Contour plots and "D line plots
        if MaPlLib:
            for c in Contour2D:                                             # init for 2D contour plots -- if there are any defined: Contour2D comes from plot options input
                S_ = c.split()[0]
                if S_.upper() == S.upper():                                 # S is elset label, c.split()[0] is 1st part of key
                    ContourSets += [S]
                    ContourX[c], ContourY[c], ContourF[c] = [], [], []
                    try:    head = Header[0]
                    except: head = " "
                    PC, pc = DefPlot("Contour " + head + ' ' + c.split()[1] + " " + SolSecs[S].Set + " time " + "%6.3f" % TimeEl)
                    pc.axis('equal')
                    ContourP[c] = [PC, pc]
                    ContourM[c] = ValMinMax(1, 3)                           # for min / max values -- 1 int value, 3 float values
            #
            if SolSecs[S].ElemTypes in __Length2D__ or SolSecs[S].ElemTypes in [["T3D2"],["T3D3I"]]:
                LinePlotData = {}
                for i in Colors2:
                    LinePlotData[i] = {}
                    LinePlotData[i]['xData'] = []                           # collects all x-values of this set
                    LinePlotData[i]['yData'] = []
        # collect data from dataframe
        SetResults = ElResults[ElResults["ElSet"] == S]                     # subset of dataframe
        SetResultsLabl = SetResults["Label"]                                # series of subset of dataframe
        SetResultsIPCo = SetResults["IntPointCoordinates"]
        SetResultsData = SetResults["Data"]
        SetResultsIpIn = SetResults["IntPointIndex"]
        SetResultsSlic = SetResults["Slices"]
        #
        pDaX, pDaY = {}, {}  # for permanent storage of line plot data
        # loop over all elements with each integration point -- basically a loop over integration points
        for ElLabel,ElIpCo,ElData,ElIpIn,ElDaSl in zip(SetResultsLabl, SetResultsIPCo, SetResultsData, SetResultsIpIn, SetResultsSlic):
            #
            # starts with deformed mesh plot + cracks - evaluation using one ip is sufficient
            if ElIpIn.strip() in ["0","Cra","Cra2"]:
                el  = ElemList[ FindIndexByLabel( ElemList, ElLabel) ]
                if el.Type in _Bond2D_: continue
                xy = el.NodalCoordinates(NodeList, NoIndToCMInd)
                _, uxy = ElementDisplList(el, DisplList, NoIndToCMInd, ScaleU)  # def ElementDisplList see below
                linew=LiWiCrvThin
                if el.Type in _Length2D_: #_Truss2DAll_ or el.Type in _BeamsAll_:
                    samp, xx, yy = [[-1.,0.,0.],[1.,0.,0.]], [], []          # samp: local sampling points for nodes, ZeroD for isogeom
                    linew=LiWiCrvMid
                elif el.Type in _TriNodes_:
                    samp, xx, yy = [[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]], [], [] # samp: local sampling points for nodes CPS3
                else:
                    samp, xx, yy = [[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.]], [], []
                for j in range(len(samp)):
                    r, s, t = samp[j][0], samp[j][1], samp[j][2]
                    XX  = el.FormX_( r, s, t)
                    xyP = dot(XX,array(xy))                                     # undeformed coordinates
                    if el.Type in ["SH4","SH3"]: NN  = el.FormX_( r, s, 0.)
                    else:                        NN  = el.FormN( r, s, 0.)
                    try:                                                        # e.g. quad element with overlayed beam element will not work here
                        xyP = xyP + dot(NN,array(uxy))                          # displaced coordinates
                    except:
                        pass
                    xx += [xyP[0]]                                              # collects displaced nodal position of element
                    yy += [xyP[1]]
                xc, yc = mean(xx), mean(yy)                                     # coordinates of element center
                xx+= [xx[0]]                                                    # to close the sequence of edges
                yy+= [yy[0]]                                                    # to close the sequence of edges
                if min(xx) < xxMin: xxMin=min(xx)
                if max(xx) > xxMax: xxMax=max(xx)
                if min(yy) < yyMin: yyMin=min(yy)
                if max(yy) > yyMax: yyMax=max(yy)
                if MaPlLib:
                    p1.plot(                                   xx,yy,'--',color=Colors[s_],linewidth=linew) # matplotlib plot deformed element geometry
                    if el.Type in ["SH4","SH3","SB3"]: p2.plot(xx,yy,'--',color=Colors[s_],linewidth=linew)
                else:
                    for j in range(len(xx)): ffg.write("%8.4f  %8.4f\n"%(xx[j],yy[j]))  # gunplot write temp deformed mesh data
                    ffg.write("\n")
                # crack plotting
                if el.Type in ['CPS4S','CPE4S','CPS3S','CPE3S','CPS6S','CPE6S']:
                    sns = ElIpIn.strip()
                    if sns == 'Cra':
                        crX, crY, crNx, crNy = ElIpCo[0], ElIpCo[1], ElData[0], ElData[1]       # center x, center y, nx, ny
                        if MaPlLib: p1.plot([xc], [yc], 'o', color='gray')      # matplotlib mark cracked elements
                        xyC = el.FindGlobalDiscontinuityEdges([crX, crY], [crNx, crNy], NodeList, NoIndToCMInd)  # xyC: [[x, y of 1st edge],[x, y, of 2nd edge]]
                        Plot2DCrackContour(xyC[0], xyC[1])
                    elif sns == 'Cra2':
                        crX, crY, crNx, crNy = ElIpCo[0], ElIpCo[1], ElData[0], ElData[1]       # center x, center y, nx, ny
                        if sqrt(crNx ** 2 + crNy ** 2) > ZeroD:  # crack is there
                            if MaPlLib: p1.plot([xc], [yc], 'o', color='red')
                            xyC = el.FindGlobalDiscontinuityEdges([crX, crY], [crNx,crNy],NodeList,NoIndToCMInd)  # xyC: [[x, y of 1st edge],[x, y, of 2nd edge]]
                            Plot2DCrackContour(xyC[0], xyC[1])
            # that's it for deformed mesh plot
            #
            # continues with plotting of data within deformed mesh
            try:
                j = int(ElIpIn)                                             # j becomes current integration point
                r = SamplePoints[el.IntT,el.nInt-1,j][0]
                s = SamplePoints[el.IntT,el.nInt-1,j][1]
                t = SamplePoints[el.IntT,el.nInt-1,j][2]
            except:
                continue
            #
            if j==0:                                                        # for whole element
                # different plot approaches for continuum elements and line elements
                AreaF, LineF = False, False
                if   el.Type in _FourNodes2D_ or el.Type in _Shells3D_ or el.Type in _TriNodes_ :
                    AreaF = True
                elif el.Type in _Length2D_ or el.Type in ["T3D2","T3D3I"]:
                    LineF = True
                    x0, x1, y0, y1 = xx[0], xx[1], yy[0], yy[1]             # should be displaced coordinates of end nodes in base / x-y -plane --hopefully
                    Lxy = sqrt((x1-x0)**2+(y1-y0)**2)
                    sinA = ( y1 - y0 )/Lxy                                  # in base x-y plane
                    cosA = ( x1 - x0 )/Lxy
                    if cosA < 0:  angle =  90 - arcsin(sinA) * 180 / 3.141593
                    else:         angle = -90 + arcsin(sinA) * 180 / 3.141593
                    if angle < -91.: angle = angle + 180.
                else:
                    raise NameError("ConFemPost::PostElem2D: element type not yet implemented",el.Type)
            #
            sig_ = []
            for sl in ElDaSl:                                               # loop over sclices - maybe shis may be moved one level obove for si, sj transforming
                sl_ = sl.split(":")
                if sl_ == ['']: continue
                sig_ += [ ElData[int(sl_[0]) : int(sl_[1])+1] ]             # slicing range
            # following for global coordinates of integration points
            if el.Type in ["SH4","SH3"]: NN  = el.FormX_( r, s, 0.)
            else:                        NN  = el.FormN( r, s, 0.)
            # deformed global ip coordinates
            try:                                                            # e.g. quad element with overlayed beam element will not work here
                ElIpCo = ElIpCo + dot(NN,array(uxy))[0:2]                   # displaced global coordinates of integration point
            except:
                pass
            # plot / assemble data in ip
            if AreaF:                                                       # plot for continuum/area elements
#                if el.Type == 'SB3':           pp = PrinC(-sig_[0][0],-sig_[0][1],-sig_[0][2]) # to make positive values -- compression upper side -- 'green' --> confusing plot
                if el.Type == 'SB3':           pp = PrinC( sig_[0][0], sig_[0][1], sig_[0][2]) #
                elif el.Type in ['SH4','SH3']: pp = PrinC( sig_[0][0], sig_[0][1], sig_[0][2]) # slice for normal forces
                else:                          pp = PrinC( sig_[0][0], sig_[0][1], sig_[0][3]) # slice for 2D stresses
                SS, S2 = PloPrin(pp, ScaleSec[0], ElIpCo[0],ElIpCo[1], p1, MaPlLib,ffred,ffgre) # SS larger principal value, S2 smaller principal value
                NorMinMax.used = True                                       # also used for SB3 moments
                NorMinMax.CheckMax(SS, [el.Label], [ElIpCo[0], ElIpCo[1], S2])
                NorMinMax.CheckMin(S2, [el.Label], [ElIpCo[0], ElIpCo[1], SS])
                if MaPlLib:
                    if el.Type in ['SH4','SH3']:
                        pp = PrinC(sig_[1][0], sig_[1][1], sig_[1][2])      # slice for bending moments
                        SS, S2 = PloPrin(pp,ScaleSec[1], ElIpCo[0],ElIpCo[1],p2,MaPlLib, ffred,ffgre)  # SS larger principal value, S2 smaller principal value
                        MomMinMax.used = True
                        MomMinMax.CheckMax(SS, [el.Label], [ElIpCo[0], ElIpCo[1], S2])
                        MomMinMax.CheckMin(S2, [el.Label], [ElIpCo[0], ElIpCo[1], SS])
                    if el.Type == 'SB3' and ElIpIn.strip() in ["0"]:        # for shear forces -- determined for center of elmement only
                        xS = (NodeList[NoIndToCMInd[el.Inzi[0]]].XCo + NodeList[NoIndToCMInd[el.Inzi[1]]].XCo + NodeList[NoIndToCMInd[el.Inzi[2]]].XCo) / 3.
                        yS = (NodeList[NoIndToCMInd[el.Inzi[0]]].YCo + NodeList[NoIndToCMInd[el.Inzi[1]]].YCo + NodeList[NoIndToCMInd[el.Inzi[2]]].YCo) / 3.
                        qx, qy = sig_[1][0], sig_[1][1]                     # second slice
                        if qx>0: p2.plot([xS, xS + scaleQ * qx], [yS, yS], '-', color='red')
                        else:    p2.plot([xS, xS + scaleQ * qx], [yS, yS], '-', color='green')
                        if qy>0: p2.plot([xS, xS], [yS, yS + scaleQ * qy], '-', color='red')
                        else:    p2.plot([xS, xS], [yS, yS + scaleQ * qy], '-', color='green')
#                       p2.plot([xS], [yS], 'o', color='gray')
                    if S in ContourSets:                                    # ContourSets collects labels of element sets
                        for c in Contour2D:                                 # Contour2D collects input from .plt
                            ContourX[c] += [ ElIpCo[0] ]                    # x, y, and value data for contour plots
                            ContourY[c] += [ ElIpCo[1] ]
                            ir = Contour2D[c][0]  - len(ElIpCo)                 # Contour2D[c] --> index in ElemResults - 0 is starter for 1st float
                            ContourF[c] += [ ElData[ir] ]
                            ContourM[c].used = True
                            ContourM[c].CheckMax( ElData[ir], [el.Label], [ ElIpCo[0], ElIpCo[1], 0.])
                            ContourM[c].CheckMin( ElData[ir], [el.Label], [ ElIpCo[0], ElIpCo[1], 0.])
            elif LineF:                                                     # plotting for line elements
                if MaPlLib:
                    if el.Type in ["BAX23I"]:
                        sign_ = [1, 1, 1, 1, 1]
                    else:
                        sign_ = [1, 1, -1, 1, 1]
                    PlotRotLinesText_( S,el,j,ElIpCo,angle, sign_,sig_,ScaleSec, p1,pDaX,pDaY, cosA,sinA ) # j is integration point index
                    if j == el.nIntL-1:                                         # element completed - rearrange data for plotting
                        for j_ in range(len(pDaX[j])):                          # loop over colors
                            col = Colors2[j_]
                            X_, Y_ = [], []
                            for k in pDaX:                                      # loop over integration points
                                X_ += [ pDaX[k][j_] ]
                                Y_ += [ pDaY[k][j_] ]
                            LinePlotData[col]['xData'] += [ X_ ]
                            LinePlotData[col]['yData'] += [ Y_ ]
                else:                                                           # for gnuplot
                    if el.Type in _BeamsAll_:
                        Truss2DCoor    += [[xyP[0], xyP[1]]]
                        Truss2DStr     += [sig_[0][0]]                          # normal force
                        Truss2DStrBond += [sig_[1][0]]                          # bending moment
                    else:
                        pass
#                p1.plot([ElIpCo[0]], [ElIpCo[1]], 'x', color='darkred')         # cross in integration point

            # that's it for plotting of data within deformed mesh
            # end of loop over elements with each integration point
        # plot lines, contour and gnuplot
        if MaPlLib:
            # plot all 2D line data
            if LineF:
                for pp in LinePlotData:                                         # plot line data of whole set - indexed by color
                    r = len(LinePlotData[pp]['xData'])
                    for j_ in range(r):
                        xp = LinePlotData[pp]['xData'][j_]
                        yp = LinePlotData[pp]['yData'][j_]
                        p1.plot( xp, yp, color=pp)
                if S in SetSlices:                                              # work around for labels
                    for j_, sl in enumerate(SetSlices[S]):
                        p1.plot([0.,1.],[0.,0.], color=Colors2[j_], label=sl)
                    print(SetSlices[S])
                p1.legend()
            #
            if NorMinMax.used: NorMinMax.Annotate(p1)                           # method from class ValMinMax for 2D elements
            if MomMinMax.used: MomMinMax.Annotate(p2)
            # plot contour plots
            if S in ContourSets:
                for c in Contour2D:
                    x = np.array(ContourX[c][:])
                    y = np.array(ContourY[c][:])
                    z = np.array(ContourF[c][:])
                    triang = tri.Triangulation( x, y)
                    def apply_mask(triang, alpha=1.0):
                        # https://stackoverflow.com/questions/42426095/matplotlib-contour-contourf-of-concave-non-gridded-data
                        # Mask triangles with side length bigger some alpha
                        triangles = triang.triangles
                        # Mask off unwanted triangles.
                        xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
                        ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
                        maxi = np.max(np.sqrt(xtri ** 2 + ytri ** 2), axis=1)
                        # apply masking
                        triang.set_mask(maxi > alpha)
                        def CheckMask():                    # mask inidicated as True
                            for m in triang.mask:
                                if not m: return m          # one not aasked found
                            return True
                        if CheckMask():
                            raise NameError(f"PostElem2D: length {alpha:f} masks whole triangulization for {c:s}")
                    apply_mask(triang, alpha=float(Contour2D[c][1]))
                    P_, p_ = ContourP[c][0], ContourP[c][1]
                    p_.set_aspect('equal')
                    tri_ = p_.tricontourf(triang, z, cmap='coolwarm')
                    P_.colorbar(tri_, ax=p_)
                    p_.scatter( x, y, s=2, color='black', marker="+")
                    if Contour2D[c][2].upper() == 'H':  # for histogram
                        elset = c.split()[0]
                        resTy = c.split()[1]
                        RF = ContourF[c]
                        mRF, sRF = mean(RF), std(RF)
                        Ph, ph = DefPlot(elset + ' ' + resTy + " mean " + "%.3f" % mRF + ' std ' + "%.3f" % sRF)
                        ph.hist(RF, bins=20, facecolor='lightblue', alpha=0.2, histtype='stepfilled', density=True)
                        x = linspace(mRF - 4 * sRF, mRF + 4 * sRF, num=100)
                        yy = st.norm.pdf(x, mRF, sRF)
                        Ph.autofmt_xdate()
                        ph.plot(x, yy)
                    if ContourM[c].used:
                        ContourM[c].Annotate( ContourP[c][1])
            # end of MaPlLib ifs
        else:                                                               # gnuplot
            if len(Truss2DCoor)>0:                                          # gnuplot truss element data
                Truss2DStrDel = max(Truss2DStr) - min(Truss2DStr)
                if (Truss2DStrDel)>ZeroD:
                    for j in range(len(Truss2DCoor)):
                        ffIP.write("%12.4e  %12.4e  %12.4e\n"%(Truss2DCoor[j][0],Truss2DCoor[j][1],Truss2DStr[j]))
                Truss2DStrDel = max(Truss2DStrBond) - min(Truss2DStrBond)
                if (Truss2DStrDel)>ZeroD:
                    for j in range(len(Truss2DCoor)):
                        ffIP2.write("%12.4e  %12.4e  %12.4e\n"%(Truss2DCoor[j][0],Truss2DCoor[j][1],Truss2DStrBond[j]))
            ffg.close()
            ffc.close()
            ffred.close()
            ffgre.close()
            ffIP.close()
            ffIP2.close()
            figure1 = gp.gp()
            figure1.c('set term wxt ' + str(gi))
            figure1.c('set size square')
            figure1.c('set grid')
            figure1.c('unset key')
            figure1.c('set title "Set '+el.Set+' -- window '+str(gi)+' -- Time '+str(TimeEl)+' -- Scale '+str(ScaleU)+' " ')
            dxx, dyy = max(1., xxMax - xxMin), max(1., yyMax - yyMin)
            figure1.c('set xrange[' + str(xxMin - 0.08 * dxx) + ':' + str(xxMax + 0.08 * dxx) + ']')
            figure1.c('set yrange[' + str(yyMin - 0.08 * dyy) + ':' + str(yyMax + 0.08 * dyy) + ']')
            figure1.c('plot "' + dfg + '" with lines linestyle 1 linecolor rgb "#0060ad" linetype 1 linewidth 0.5 ')  # deformed elements - dfg data file
            if os.path.getsize(dfc) > 0:   figure1.c('replot "' + dfc + '" with lines linestyle 1 linecolor rgb "gray30" ')  # 2D crack contours
            if os.path.getsize(dfred) > 0: figure1.c('replot "' + dfred + '" with lines linestyle 1 linecolor rgb "red" ')  # tensile principal stresses
            if os.path.getsize(dfgre) > 0: figure1.c('replot "' + dfgre + '" with lines linestyle 1 linecolor rgb "dark-green" ')  # commpressive principal stresses
            if os.path.getsize(dfIP) > 0:  figure1.c('replot "' + dfIP + '" using 1:2:3 w p pointtype 7 pointsize 0.3 lc palette')
            #                figure1.p(filename=dfg+'.eps', width=20, height=12, fontsize=12, term='wxt')
            #                if os.path.getsize('tmp'+str(S)+'IP2.dat')>0:
            if os.path.getsize(dfIP2) > 0:
                gi += 1
                figure1.c('set term wxt ' + str(gi))
#                figure1.c('set size square')
                figure1.c('set grid')
                figure1.c('unset key')
                figure1.c('set title "Set '+el.Set+' -- window '+str(gi)+' -- Time '+str(TimeEl)+' -- Scale '+str(ScaleU)+' " ')
                figure1.c('plot "'+dfg+'" with lines linestyle 1 linecolor rgb "#0060ad" linetype 1 linewidth 0.5 ')  # deformed elements
                figure1.c('replot "' + dfIP2 + '" using 1:2:3 w p pointtype 7 pointsize 0.3 lc palette')
            #                    figure1.p(filename=FilDir+'tmp'+str(S)+'_.eps', width=20, height=12, fontsize=12, term='wxt')
            #    #            figure1.c('pause -1')
            #    #            figure1.c('a = system("read a; echo $a")')
            ti.sleep(1)
            _ = str(input('continue: '))
            os.remove(dfg); os.remove(dfc); os.remove(dfred); os.remove(dfgre); os.remove(dfIP); os.remove(dfIP2);
    #
    return 0
 
# for 3D plor of continuum shell elements with material *RCSHELL
def PostElemSH4( SolSec, Time, ScaleShellL, ElResults ):
    def IniP(Text, Time, ElSet, proj):
        if proj=='3d':
            P0 = plt.figure().add_subplot(111,title=Text+'Time '+str(Time)+', Element Set '+ElSet, projection=proj)
        else:
            P0 = plt.figure().add_subplot(111,title=Text+'Time '+Time+', Element Set '+ElSet)
            P0.axis('equal')
            P0.grid()
        P0.tick_params(axis='x', labelsize=FonSizAx) #yticks(fontsize=FonSizAx)
        P0.tick_params(axis='y', labelsize=FonSizAx) #xticks(fontsize=FonSizAx)
        P0.set_xlabel('x')
        P0.set_ylabel('y')
        return P0
    for s_, S in enumerate(SolSec):
        SetResults = ElResults[ElResults["ElSet"] == S]                     # subset of dataframe
        ShRes = SetResults[SetResults["Type"].isin(["SH4", "SH3"])]
        if len(ShRes) == 0: return 0
        ShResLabl = ShRes["Label"]
        ShResData = ShRes["Data"]
        ShResIpIn = ShRes["IntPointIndex"]
        ShResIpCo = ShRes["IntPointCoordinates"]
        ShResType = ShRes["Type"]
        ShResMark = ShRes["Marker"]
        P2 = IniP('strains: ', Time, S, '3d')
        for Labl,Data,IpIn,IpCo,Type,Marker in zip( ShResLabl,ShResData,ShResIpIn,ShResIpCo,ShResType,ShResMark):  # loop over integration points
            XX = IpCo[0]
            YY = IpCo[1]
            ZZ = IpCo[2]                                                    # local local isoparametric t-coordinate
            if Marker=='R':                                                 # reinforcement
                col='blue'
#                    P2.plot([XX],[YY],[ZZ],'o',color=col)                  # single ip point plot
            elif ZZ<-0.9 or ZZ>0.9:
#                else:
                scaleX = ScaleShellL # SP*Scales[ScaleKey][2]
                pp = PrinC( Data[0], Data[1], 0.5*Data[2] )                 # principal strains concrete
                linesty = 'dotted'
                SS = pp[0]                                                  # larger principal strain
                if SS > 0.1e-3: linesty = 'solid'
                if SS>=0: P2.plot([XX-SS*scaleX*pp[1],XX+SS*scaleX*pp[1]],[YY-SS*scaleX*pp[2],YY+SS*scaleX*pp[2]],[ZZ,ZZ],'r',linestyle=linesty)
                else:     P2.plot([XX-SS*scaleX*pp[1],XX+SS*scaleX*pp[1]],[YY-SS*scaleX*pp[2],YY+SS*scaleX*pp[2]],[ZZ,ZZ],'g',linestyle=linesty)
                linesty = 'dotted'
                SS = pp[3]
                if SS > 0.1e-3: linesty = 'solid'
                if SS>=0: P2.plot([XX-SS*scaleX*pp[2],XX+SS*scaleX*pp[2]],[YY+SS*scaleX*pp[1],YY-SS*scaleX*pp[1]],[ZZ,ZZ],'r',linestyle=linesty)
                else:     P2.plot([XX-SS*scaleX*pp[2],XX+SS*scaleX*pp[2]],[YY+SS*scaleX*pp[1],YY-SS*scaleX*pp[1]],[ZZ,ZZ],'g',linestyle=linesty)
                col='blue' # 'red'
                P2.plot([XX],[YY],[ZZ],'o',color=col,markersize=3)           # single ip point plot

        print('PostElemSH process plot for time ', Time)

def PostNode( ElList, NodeList,NICM, VecU, SN, ft):
    LX, LY, LZ, LU = [], [], [], []
    for i in NodeList:
        LX.append(i.XCo)
        LY.append(i.YCo)
        LZ.append(i.ZCo)
    for i in VecU: 
        LU.append(fabs(i))
    Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = min(LX), max(LX), min(LY), max(LY), min(LZ), max(LZ)
    D = max(Xmax-Xmin,Ymax-Ymin,Zmax-Zmin)
    U = max(LU)
    if U > ZeroD: 
        scale = SN*0.5*D/U
    else:
        print("PostNode   Zero Displacements")
        scale = 0.                                              # estimation for scaling factor
    P0 = plt.figure().add_subplot(111,title='Deformed mesh final step (scale '+'%7.4f'%(scale)+')')
    P0.axis('equal')
    P0.grid()
    P0.tick_params(axis='x', labelsize=FonSizAx) #yticks(fontsize=FonSizAx) 
    P0.tick_params(axis='y', labelsize=FonSizAx) #xticks(fontsize=FonSizAx)
    P0.set_xlabel('x')
    P0.set_ylabel('y')
    # uncomment to have node labels ---------------------------------
#    for i in NodeList: P0.text(i.XCo, i.YCo, "%i" % i.Label, ha='left', va='bottom', rotation=0., color='black', fontsize=10)
    # ---------------------------------------------------------------
    #
    ElSets, colI = {}, 0                                            # for different colors for different element sets
    for Elem in ElList:
        xS, yS, xP, yP, xN, yN = [], [], [], [], [], []
        Set = Elem.Set
        if Set not in ElSets:
            ElSets[Set] = colI
            colI = colI+1
        if Elem.Type=='SB3':
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo,NodeList[NICM[Elem.Inzi[0]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo,NodeList[NICM[Elem.Inzi[0]]].YCo]
            # plot labels
            lN = [NodeList[NICM[Elem.Inzi[0]]].Label,NodeList[NICM[Elem.Inzi[1]]].Label,NodeList[NICM[Elem.Inzi[2]]].Label]
            xS_ = (NodeList[NICM[Elem.Inzi[0]]].XCo+NodeList[NICM[Elem.Inzi[1]]].XCo+NodeList[NICM[Elem.Inzi[2]]].XCo)/3.
            yS_ = (NodeList[NICM[Elem.Inzi[0]]].YCo+NodeList[NICM[Elem.Inzi[1]]].YCo+NodeList[NICM[Elem.Inzi[2]]].YCo)/3.
            eL = "%i/%i"%(i,Elem.Label)
            plt.text(xS_,yS_,eL,ha='center',va='center')                # elements
            for xL, yL, lL in zip(xN,yN,lN): plt.text(xL, yL, lL,ha='right',va='top',color='red') # nodes
        elif Elem.Type in ['CPE4','CPS4','CAX4']:
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo,NodeList[NICM[Elem.Inzi[3]]].XCo,NodeList[NICM[Elem.Inzi[0]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo,NodeList[NICM[Elem.Inzi[3]]].YCo,NodeList[NICM[Elem.Inzi[0]]].YCo]
            xP = [NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[1]]].XCo+scale*VecU[Elem.DofI[1,0]],NodeList[NICM[Elem.Inzi[2]]].XCo+scale*VecU[Elem.DofI[2,0]],NodeList[NICM[Elem.Inzi[3]]].XCo+scale*VecU[Elem.DofI[3,0]],NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]]]
            yP = [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]],NodeList[NICM[Elem.Inzi[1]]].YCo+scale*VecU[Elem.DofI[1,1]],NodeList[NICM[Elem.Inzi[2]]].YCo+scale*VecU[Elem.DofI[2,1]],NodeList[NICM[Elem.Inzi[3]]].YCo+scale*VecU[Elem.DofI[3,1]],NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]]]
        elif Elem.Type=='CPE3' or Elem.Type=='CPS3':
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo,NodeList[NICM[Elem.Inzi[0]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo,NodeList[NICM[Elem.Inzi[0]]].YCo]
            xP = [NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[1]]].XCo+scale*VecU[Elem.DofI[1,0]],NodeList[NICM[Elem.Inzi[2]]].XCo+scale*VecU[Elem.DofI[2,0]],NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]]]
            yP = [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]],NodeList[NICM[Elem.Inzi[1]]].YCo+scale*VecU[Elem.DofI[1,1]],NodeList[NICM[Elem.Inzi[2]]].YCo+scale*VecU[Elem.DofI[2,1]],NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]]]
        elif Elem.Type=='SH4':
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo,NodeList[NICM[Elem.Inzi[3]]].XCo,NodeList[NICM[Elem.Inzi[0]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo,NodeList[NICM[Elem.Inzi[3]]].YCo,NodeList[NICM[Elem.Inzi[0]]].YCo]
            xP = [NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[1]]].XCo+scale*VecU[Elem.DofI[1,0]],NodeList[NICM[Elem.Inzi[2]]].XCo+scale*VecU[Elem.DofI[2,0]],NodeList[NICM[Elem.Inzi[3]]].XCo+scale*VecU[Elem.DofI[3,0]],NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]]]
            yP = [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]],NodeList[NICM[Elem.Inzi[1]]].YCo+scale*VecU[Elem.DofI[1,1]],NodeList[NICM[Elem.Inzi[2]]].YCo+scale*VecU[Elem.DofI[2,1]],NodeList[NICM[Elem.Inzi[3]]].YCo+scale*VecU[Elem.DofI[3,1]],NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]]]
        elif Elem.Type=='SH3':
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo,NodeList[NICM[Elem.Inzi[0]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo,NodeList[NICM[Elem.Inzi[0]]].YCo]
            xP = [NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[1]]].XCo+scale*VecU[Elem.DofI[1,0]],NodeList[NICM[Elem.Inzi[2]]].XCo+scale*VecU[Elem.DofI[2,0]],NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]]]
            yP = [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]],NodeList[NICM[Elem.Inzi[1]]].YCo+scale*VecU[Elem.DofI[1,1]],NodeList[NICM[Elem.Inzi[2]]].YCo+scale*VecU[Elem.DofI[2,1]],NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]]]
        elif Elem.Type in ['T1D2']:
            xN = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo]
            yN = [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[1]]].YCo]
            xP = [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[1]]].XCo]
            yP = [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[1]]].YCo+scale*VecU[Elem.DofI[1,0]]]
        elif Elem.Type in ['T2D2','T2D3','T2D2I','T2D3I','TAX2','TAX2I','TAX3','TAX3I'] or Elem.Type in ['T3D2','T3D3','T3D2I','T3D3I']:
#        elif Elem.Type in _Truss2DAll_ or Elem.Type in ['T3D2', 'T3D3', 'T3D2I', 'T3D3I']:
            xN = [ NodeList[NICM[Elem.Inzi[0]]].XCo, NodeList[NICM[Elem.Inzi[1]]].XCo ]             # initial x-position start node, end node
            yN = [ NodeList[NICM[Elem.Inzi[0]]].YCo, NodeList[NICM[Elem.Inzi[1]]].YCo ]             # y
            xP = [ xN[0] + scale*VecU[Elem.DofI[0,0]], xN[1] + scale*VecU[Elem.DofI[1,0]] ]
            yP = [ yN[0] + scale*VecU[Elem.DofI[0,1]], yN[1] + scale*VecU[Elem.DofI[1,1]] ]
        elif Elem.Type=='S1D2':
            j=0
            k0 = Elem.DofI[j,0]
            k1 = Elem.DofI[j+1,0]
            xS += [NodeList[NICM[Elem.Inzi[j]]].XCo]
            yS += [NodeList[NICM[Elem.Inzi[j]]].YCo+scale*(VecU[k1]-VecU[k0])]
        elif Elem.Type in ['B23E','B23EI','BAX23EI']:
            xN += [NodeList[NICM[Elem.Inzi[0]]].XCo,NodeList[NICM[Elem.Inzi[2]]].XCo]
            yN += [NodeList[NICM[Elem.Inzi[0]]].YCo,NodeList[NICM[Elem.Inzi[2]]].YCo]
            xP =  [NodeList[NICM[Elem.Inzi[0]]].XCo+scale*VecU[Elem.DofI[0,0]],NodeList[NICM[Elem.Inzi[2]]].XCo+scale*VecU[Elem.DofI[2,0]]]
            yP =  [NodeList[NICM[Elem.Inzi[0]]].YCo+scale*VecU[Elem.DofI[0,1]],NodeList[NICM[Elem.Inzi[2]]].YCo+scale*VecU[Elem.DofI[2,1]]]
        elif Elem.Type not in ['SB3','CPE4','CPS4','SH4', 'SH3', 'T1D2','S1D2','CPE3','CPS3','B23E',
                               "Bond2D2","Bond2D3","Bond3D2","Bond3D3","BondAX2","BondAX3"]:
            for j in range(Elem.nNod):
                k0 = Elem.DofI[j,0]
                k1 = Elem.DofI[j,1]
                xN += [NodeList[NICM[Elem.Inzi[j]]].XCo]                    # undeformed 2D geometry
                yN += [NodeList[NICM[Elem.Inzi[j]]].YCo]                    # undeformed 2D geometry
                xP += [NodeList[NICM[Elem.Inzi[j]]].XCo+scale*VecU[Elem.DofI[j,0]]]
                yP += [NodeList[NICM[Elem.Inzi[j]]].YCo+scale*VecU[Elem.DofI[j,1]]]
        P0.plot(xN,yN, 'b--')                                               # undeformed geometry for current element
        if (len(xP)+len(yP))>0:                                             # deformed geometry for current element
            P0.plot(xP,yP, 'b-',linewidth=LiWiCrv)
#            P0.plot(xP,yP, 'b-',linewidth=0.5) #LiWiCrvThin)
#            P0.plot(xP, yP, 'b-', linewidth=0.5)
            P0.plot(xP,yP, 'o', color=Colors[ElSets[Set]])
#            P0.scatter(xP,yP, s=10, c=Colors[ElSets[Set]]) #, alpha=0.5)
#            P0.scatter(xP,yP, s=100, c='red', alpha=0.5)
            if ft != None:
                print(Elem.Type,',',scale,',',xN,',',yN,',',xP,',',yP, file=ft )
        if len(xS)>0:                                                       # for S1D2 only
            P0.plot(xS,yS, 'b-',linewidth=LiWiCrv)
            P0.plot(xS,yS, 'o', color=Colors[ElSets[Set]])
    return 0

def PostElem3D( ElList, NodeList,NoIndToCMInd, ElResults, Time):
    XList, YList, ZList, LisLen = [], [], [], 0
    for i in NodeList:
        XList += [i.XCo]
        YList += [i.YCo]
        ZList += [i.ZCo]
    maxX, minX, maxY, minY, maxZ, minZ = max(XList), min(XList), max(YList), min(YList), max(ZList), min(ZList)
    delX, delY, delZ = maxX-minX, maxY-minY, maxZ-minZ
    delR2  = 0.5*max(delX, delY, delZ)
    meanX, meanY, meanZ = 0.5*(minX+maxX), 0.5*(minY+maxY), 0.5*(minZ+maxZ)
    #
    nxLi_,nyLi_,nxyLi_,mxLi_,myLi_,mxyLi_,qxLi_,qyLi_, aaLi_ = [], [], [], [], [], [], [], [], [] # collect internal force / thickness data for all elements
    ShRes = ElResults[ElResults["Type"].isin(["SH4","SH3"])]
    if len(ShRes) == 0: return 0
    ShResLabl = ShRes["Label"]
    ShResData = ShRes["Data"]
    ShResIpIn = ShRes["IntPointIndex"]
    ic = None
    for ElLabl, ElData, ElIpIn in zip( ShResLabl, ShResData, ShResIpIn):
        if ElIpIn.strip() == '0':
            try:
                nxLi_.append( nx/ic);nyLi_.append(ny/ic);nxyLi_.append(nxy/ic);mxLi_.append(mx/ic);myLi_.append(my/ic);mxyLi_.append(mxy/ic);qxLi_.append(qx/ic);qyLi_.append(qy/ic);aaLi_.append(aa_/ic)
            except:
                pass
            nx, ny, nxy, qy, qx, mx, my, mxy, aa_ = 0., 0., 0., 0., 0., 0., 0., 0.,  0.
            ic = 0
        nx  += ElData[0]
        ny  += ElData[1]
        nxy += ElData[2]
        mx  += ElData[3]
        my  += ElData[4]
        mxy += ElData[5]
        qx  += ElData[6]
        qy  += ElData[7]
        aa_ += ElData[8]
        ic  += 1
    nxLi_.append(nx/ic);nyLi_.append(ny/ic);nxyLi_.append(nxy/ic);mxLi_.append(mx/ic);myLi_.append(my/ic);mxyLi_.append(mxy/ic);qxLi_.append(qx/ic);qyLi_.append(qy/ic);aaLi_.append(aa_/ic)
    # end of element loop
    # plot whole stuff
    if ic!=None:
        TolP = 1.0e-9
        plotVal3D(aaLi_, 'thickness_', minX,maxX,minY,maxY,minZ,maxZ,ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        if max(map(abs,nxLi_)) > TolP: plotVal3D(nxLi_, 'n_x_', minX, maxX, minY, maxY, minZ, maxZ, ElList, NodeList, NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        if max(map(abs,nyLi_)) > TolP: plotVal3D(nyLi_, 'n_y_', minX,maxX,minY,maxY,minZ,maxZ,ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        if max(map(abs,nxyLi_))> TolP: plotVal3D(nxyLi_,'n_xy_',minX,maxX,minY,maxY,minZ,maxZ,ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        plotVal3D(mxLi_, 'm_x_', minX,maxX,minY,maxY,minZ,maxZ, ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        plotVal3D(myLi_, 'm_y_', minX,maxX,minY,maxY,minZ,maxZ, ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        plotVal3D(mxyLi_,'m_xy_',minX,maxX,minY,maxY,minZ,maxZ, ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        plotVal3D(qxLi_, 'q_x_', minX,maxX,minY,maxY,minZ,maxZ, ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)
        plotVal3D(qyLi_, 'q_y_', minX,maxX,minY,maxY,minZ,maxZ, ElList,NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2)

    return ic

def plotVal3D( valList,Label, minX,maxX,minY,maxY,minZ,maxZ, ElList, NodeList,NoIndToCMInd, Time, meanX,meanY,meanZ,delR2):
    fig = pl.figure() #figsize=(15.,10.))
    fig.text( 0.05, 0.93, Label+' time '+str(Time), size='x-large')
    maxVal = max(valList)
    minVal = min(valList)
    if abs(maxVal-minVal)<1.e-6: 
        aVal = abs(maxVal)
        maxVal = maxVal + 0.2*aVal
        minVal = minVal - 0.2*aVal
    norm_ = colors.Normalize(vmin=minVal, vmax=maxVal)
    ax1, orie = fig.add_axes([0.90, 0.05, 0.02, 0.8]), 'vertical' # left, bottom, width, height] where all quantities are in fractions of figure width and height
    ColM = plt.get_cmap('brg')
    mpl.colorbar.ColorbarBase(ax1, cmap=ColM, norm=norm_, orientation=orie)
    ax = plt.subplot(projection='3d')
    for i in range(len(ElList)):
        elem, CoordList = ElList[i], []
        for j in range(len(elem.Inzi)):
            ni = NoIndToCMInd[elem.Inzi[j]]
            node = NodeList[ni]
            CoordList += [[node.XCo, node.YCo, node.ZCo]]
        tri = a3.art3d.Poly3DCollection([CoordList])
        tri.set_color(ColM((valList[i] - minVal)/(maxVal-minVal)))                      # --> valList
        tri.set_edgecolor('k')
        ax.add_collection3d(tri)
    ax.set_xlim3d( meanX-delR2, meanX+delR2)
    ax.set_ylim3d( meanY-delR2, meanY+delR2)
    ax.set_zlim3d(meanZ-delR2, meanZ+delR2)
    ax.set_xlabel('x',fontsize='x-large')
    ax.set_ylabel('y',fontsize='x-large')
    ax.set_zlabel('z',fontsize='x-large')
    ax.tick_params(axis='x', labelsize='x-large')
    ax.tick_params(axis='y', labelsize='x-large')# ,labelcolor='white')
    ax.tick_params(axis='z', labelsize='x-large')#small',labelcolor='white')
    ax.set_xticks([minX,meanX,maxX])
    ax.set_xticklabels( ["%.2f"%(minX),"%.2f"%(meanX),"%.2f"%(maxX)] )
    ax.set_yticks([minY,meanY,maxY])
    ax.set_yticklabels( ["%.2f"%(minY),"%.2f"%(meanY),"%.2f"%(maxY)] )
    ax.set_zticks([minZ,meanZ,maxZ])
    ax.set_zticklabels( ["%.2f"%(minZ),"%.2f"%(meanZ),"%.2f"%(maxZ)] )
    ax.set_aspect('auto') #,'box', 'equal')
    ax.grid()
    return 0

# deformed shell structure 3D
def PostNode3D( ElList,NodeList,NoIndToCMInd, NodeResults,VecU, Time, ScaleDis):
    print('PostNode3D process plot for time ', Time)
    XList, YList, ZList, DisElemList = [], [], [], []
    for i in NodeList:
        XList += [i.XCo]
        YList += [i.YCo]
        ZList += [i.ZCo]
    xMin, xMax, XM = np.min(XList), np.max(XList), np.mean(XList)
    yMin, yMax, YM = np.min(YList), np.max(YList), np.mean(YList)
    zMin, zMax, ZM = np.min(ZList), np.max(ZList), np.mean(ZList)
    delX, delY, delZ = xMax-xMin, yMax-yMin, zMax-zMin
    delRef = max(delX, delY, delZ)
    delR2  = 0.5*delRef
    xMean, yMean, zMean = 0.5*(xMin+xMax), 0.5*(yMin+yMax), 0.5*(zMin+zMax) 
    # for measure of absolute displacement
    for elem in ElList:                                                     # loop over elements
        if elem.Type in ["TAX2","BAX21E","CAX4","BAX23EI"]:
            continue
        DisList = []
        for o in elem.Inzi:
            ni = NoIndToCMInd[o]
            node = NodeList[ni]
            if len(VecU)==0:
                d = NodeResults[node.Label][2]
                if   elem.Type in ['SH4','SH3']: DisList += [ sqrt(d[0]**2 + d[1]**2 + d[2]**2) ]
                elif elem.Type=='SB3':           DisList += [ sqrt(d[0]**2) ]
                elif elem.Type in _Length2D_:    DisList += [ sqrt(d[0]**2 + d[1]**2) ]
            else:
                DofInd = node.GlobDofStart
                if   elem.Type in ['SH4','SH3']: DisList += [sqrt(VecU[DofInd]**2+VecU[DofInd+1]**2+VecU[DofInd+2]**2)]
                elif elem.Type=='SB3':           DisList += [VecU[DofInd]]
                elif elem.Type in _Length2D_:                               # uhc quick and dirty
                    if elem.Type in ['B23E']:
                        if len(node.DofT)==3:   DisList += [sqrt(VecU[DofInd]**2+VecU[DofInd+1]**2)] # to exclude enhance node
                    else:
                        raise NameError("ConFemPost::PostNode3D: element ype not yet implemented",elem.Type)
        DisElemList += [np.mean(DisList)]                                   # this is always positive

    if len(DisElemList)>0:
        dMax = max(DisElemList)
        #
        scale_ = 10.                                                            # 20.  # displayed deformation should be roughly 1/scale_ of parts dimensions
        scale = ScaleDis*max(XM,YM,ZM)/(dMax*scale_)
        ColM = plt.get_cmap('summer')
        fig = pl.figure()
        fig.text( 0.05, 0.95, f"deformed mesh time {Time:7.3f} - scale {scale:7.3f}", size='x-large') # , ha='right')
        norm_ = colors.Normalize(vmin=0., vmax=dMax)
        ax1, orie = fig.add_axes([0.90, 0.05, 0.02, 0.8]), 'vertical'           # left, bottom, width, height] where all quantities are in fractions of figure width and height
    #    ax1, orie = fig.add_axes([0.05, 0.05, 0.9, 0.02]), 'horizontal' # left, bottom, width, height] where all quantities are in fractions of figure width and height
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=ColM, norm=norm_, orientation=orie)
        ax = plt.subplot(projection='3d')
        #
        for i, elem in enumerate( ElList ):
            if elem.Type in ["TAX2","BAX21E","CAX4"]:
                continue
            CoordList  = []
            for j in elem.Inzi:
                node = NodeList[ NoIndToCMInd[j] ]
                if len(VecU) == 0:
                    nL = node.Label
                    if nL in NodeResults: noRe = NodeResults[nL][2]
                    else:                 raise NameError("ConFemPost::PostNode3D: node results not available",nL)
                    if   elem.Type in ['SH4','SH3']:
                        CoordList += [[node.XCo+scale*noRe[0],  node.YCo+scale*noRe[1],  node.ZCo+scale*noRe[2]]]
                    elif elem.Type=='SB3':
                        CoordList += [[node.XCo-xMin, node.YCo-yMin, node.ZCo-zMin+scale*noRe[0]]]
                    else:
                        raise NameError ("ConFemPost::PostNode3D: element type not yet implemented",elem.Label,elem.Type)
                else:
                    DofInd = node.GlobDofStart
                    if   elem.Type in ['SH4','SH3']:
                        CoordList += [[node.XCo+scale*VecU[DofInd], node.YCo+scale*VecU[DofInd+1], node.ZCo+scale*VecU[DofInd+2]]]
                    elif elem.Type=='SB3':
                        CoordList += [[node.XCo-xMin, node.YCo-yMin, node.ZCo-zMin+scale*VecU[DofInd]]]
                    elif elem.Type in _Length2D_:                               # uhc quick and dirty
                        if elem.Type in ['B23E']:
                            if len(node.DofT)==3:   CoordList += [[node.XCo+scale*VecU[DofInd], node.YCo+scale*VecU[DofInd+1], node.ZCo]] # to exclude enhance node
                        else:                       raise NameError("ConFemPost::PostNode3D: element ype not yet implemented",elem.Type)

            tri = a3.art3d.Poly3DCollection([CoordList])
            tri.set_color(ColM( DisElemList[i]/dMax ))
            tri.set_edgecolor('k')
            ax.add_collection3d(tri)
        #
        ax.set_xlim3d( xMean-delR2, xMean+delR2 )
        ax.set_ylim3d( yMean-delR2, yMean+delR2 )
        ax.set_zlim3d( zMean-delR2, zMean+delR2 )
        ax.set_xlabel('x',fontsize='x-large')
        ax.set_ylabel('y',fontsize='x-large')
        ax.set_zlabel('z',fontsize='x-large')
        ax.tick_params(axis='x', labelsize='x-large')
        ax.tick_params(axis='y', labelsize='x-large')
        ax.tick_params(axis='z', labelsize='x-large')
        ax.set_xticks([xMin,xMean,xMax])
        ax.set_xticklabels( ["%.2f"%(xMin),"%.2f"%(xMean),"%.2f"%(xMax)] )
        ax.set_yticks([yMin,yMean,yMax])
        ax.set_yticklabels( ["%.2f"%(yMin),"%.2f"%(yMean),"%.2f"%(yMax)] )
        ax.set_zticks([zMin,zMean,zMax])
        ax.set_zticklabels( ["%.2f"%(zMin),"%.2f"%(zMean),"%.2f"%(zMax)] )
        ax.set_aspect('auto')
        ax.grid()

    return 0

def PlotHist( f5, WrNodes ):
    if WrNodes[0][3]=='p-u':   label='load-displacement node '+str(WrNodes[0][0])
    elif WrNodes[0][3]=='p-t': label='load-time node '+str(WrNodes[0][0])
    elif WrNodes[0][3]=='u-t': label='displacement-time (node '+str(WrNodes[0][0])+')'
    elif WrNodes[0][3]=='b-t': label='reaction force-time (node '+str(WrNodes[0][0])+')'
    elif WrNodes[0][3]=='b-u': label='reaction force-displacement (node '+str(WrNodes[0][0])+')'
    else: raise NameError("ConFemInOut:PlotNodes: unknown key")
    sc1, sc2 = 1., 1.
#    sc1, sc2 = -1., -1.
    L1, L2 = [0.], [0.]
    z1 = f5.readline()                                                      # 1st input line
    z2 = z1.split()
    while z1!="":
#        zLen = len(z2)                                                    
        zLen = len(z2)-2                                                    # for adding step, counter to end of line
        if (zLen % 2) != 0: raise NameError("ConFemInOut::PlotNodes: input line inconsistent",zLen,z2)
        sum_ = 0.
#        for i in range(int(zLen/2)): sum_ += float(z2[2*i-1])
        for i in range(int(zLen/2)): sum_ += float(z2[2*i+1])
        L1 += [sc1*float(z2[0])]
        L2 += [sc2*sum_]
        z1 = f5.readline()                                                  # read next line
        z2 = z1.split()
    PP = plt.figure()
    P0 = PP.add_subplot(111,title=label) 
    P0.set_title(label,fontsize=FonSizTi)
    P0.plot(L1,L2,linewidth=LiWiCrv)
    P0.set_ylim(1.05*min(L2),1.05*max(L2))
    P0.tick_params(axis='x', labelsize=FonSizAx)
    P0.tick_params(axis='y', labelsize=FonSizAx) #'x-large')
    P0.grid()
    PP.autofmt_xdate()
    return 0

# store values of ndof per node 
def NodalDisplList( NodeList, NodeResults):
    DisplList = [] 
    for i in NodeList:
        key = i.Label
        nDof = NodeResults[key][1]
        Data = NodeResults[key][2]
        dL = []
#        for j in range(nDof): dL += [Data[j]]
        for j in range(min(nDof,3)): dL += [Data[j]]                        # use coordinates only
        DisplList += [dL]
        i.dL = dL                                                           # for label text on deformed structure
    return DisplList
# store values of ndof per element 
def ElementDisplList( i, DisplList, NoIndToCMInd, ScaleU):
    uEl, uxy = [], []
    for j in i.Inzi:                                                        # loop over nodes affected by element
        j_ = NoIndToCMInd[j]
        for k in range(len(DisplList[j_])):
            uEl += [DisplList[j_][k]]                                       # collects nodal displacements per element
            uxy += [ScaleU*DisplList[j_][k]]                                # collects scaled nodal displacements per element
    return uEl, uxy

def PrinC( xx, yy, xy ):                                                    # calculation of principal values
    if ZeroD < abs(xy):
        h1 = xx+yy
        h2 = sqrt((xx-yy)**2+4*xy**2);
        s1 = .5*h1+.5*h2
        s2 = .5*h1-.5*h2
        h = (s1-xx)/xy
        L = sqrt(h**2+1)
        n11 = 1./L
        n12 = h/L
        h = (s2-xx)/xy
        L = sqrt(h**2+1)
        n21 = 1./L
        n22 = h/L 
    else:
        s1 = xx
        n11 = 1.
        n12 = 0
        s2 = yy
        n21 = 0
        n22 = 1.
    if s1>=s2: return( [s1, n11, n12, s2, n21, n22] )
    else:      return( [s2, n21, n22, s1, n11, n12 ])

def ReadLine(ff):
    z1 = ff.readline()
    z2 = z1.strip()
    return z2.split(',')

def ReadResultsFromFile( MatList, ff, ffN):
    NodeResults, EndFlag, EndFlagN  = {}, False, False
    # specifications for plotting
    ElResultsL = { # no of integration points, slices for data items to plot - each slice with own scaling, not consider elements for scaling
                                                                            # slices refer to elemout-data, index counts zero from 1st data - not from integration point coordinates
        "ELASTIC" :  {  "CPS4":    [2, ["4:7"], False],                     # 4 stress components
                        "CPS3":    [2, ["4:7"], False],                     # "
                        "CAX4":    [2, ["4:7"], False],                     # "
                        "C3D8":    [3, ["6:11"], False],                    # 6 stress compmonents
                        "BAX21EI": [2, ["3:3","4:4","5:5","8:8","9:9"], True],
                        "B23E":    [2,[], True],                            # currently dummy
                        "SB3":     [2, ["3:5","6:7"], False],               # bending moments, shear forces
                        "SH4":     [3, ["0:2","3:5","6:7"], False],         # normal forces,bending moments, shear forces  ? tlhickness nx, ny, nxy, mx, my, mxy, qx, qy, aa
                        "SH3":     [3, ["0:2","3:5","6:7"], False],         # "
                        "T2D2I":   [2, ["1:1","4:4","5:5"], True],          # uniaxial stress, bond stress, penalty stress
                        "T2D3I":   [2, ["1:1","4:4","5:5"], True],          # uniaxial stress, bond stress, penalty stress
                        "T3D3I":   [3, ["1:1", "5:5"], True],               # uniaxial stress, bond stress
                        "T3D2I":   [3, ["1:1", "5:5"], True],               # "
                        "T2D2":    [2, ["1:1"], True],
                        "TAX2":    [2, ["1:1"], True],
                        "TAX2I":   [2, ["2:2", "6:6", "7:7"], True],
                        "TAX3":    [2, ["1:1"], True],
                        "T3D2":    [3, ["1:1"], True],
                        "B23I":    [2, ["4:5", "8:8", "9:9"], True],
                        "B23EI":   [2, ["4:5", "8:8", "9:9"], True],
                        "BAX23EI": [2, ["4:5", "8:8", "9:9"], True],
                        "BAX21E":  [2, ["3:3","4:4","5:5"], True],          # normal force, bending moment, shear force
                        "BAX21":   [2, ["3:3", "4:4", "5:5"], True],
                        "B21":     [2, ["3:3", "4:4", "5:5"], True],
                        "BAX23":   [2, [], True],                # dummy
                        "B23":     [2, [], True]
#                        "BondAX3": [2, ["2:2", "3:3"], True]                # [] lateral, longitudinal bond stress
                     },
        "ELASTIC1DR":{  "TAX2I": [2, ["2:2", "6:6", "7:7"], True],
                        "TAX3I": [2, ["2:2", "6:6", "7:7"], True]
                     },
        "ELASTICC1D":{  "TAX2":  [2, ["3:3"], True],                        # circumferential uniaxial stress
                     },
        "ELASTICLT": {  "T1D2":  [2, ["1:1"], True]                         # uniaxial stress

                     },
        "ELASTICLT_SDA": { "CPS4":    [2, ["4:7"], False],
                           "CPS4S":   [2, ["4:7"], False],
                           "CPS3":    [2, ["4:7"], False],
                           "CPS3S":   [2, ["4:7"], False],
                           "C3D8":    [3, ["6:11"],False],
                           "C3D8S":   [3, ["6:11"],False]
                    },
        "ELASTIC_VISCO": { "T1D2":  [2, ["1:1"], True]                    # uniaxial stress

                     },
        "ELASTIC_PHASEFIELD": { "T1D2":  [2, ["1:1"], True]                    # uniaxial stress

                    },
        "ElasticOrtho": {"S2D6": [2, [], True]},
        "ISODAMAGE": {  "CPS4":  [2, ["4:7"],  False],
                        "CPS3":  [2, ["4:7"],  False],
                        "CPS4S": [2, ["4:7"],  False],
                        "CPS3S": [2, ["4:7"],  False],
                        "CAX4":  [2, ["4:7"],  False],
                        "C3D8":  [3, ["6:11"], False],
                        "T1D2":  [2, ["1:1"],  False]                     # uniaxial stress
                     },
        "MISES":     {  "B23EI": [2,["4:5","8:8","9:9"], True],          # edge stresses, bond stress, lateral constraint stress
                        "B23I":  [2, ["4:5", "8:8", "9:9"], True],
                        "T1D2":  [2,["1:1"], True],                      # uniaxial stress
                        "B23E":  [2, [], True],                          # dummy
                        "T2D2":  [2, ["1:1"], True],                      # uniaxial stress
                        "T2D2I": [2, ["1:1"], True],
                        "SH4":   [ 3, ["0:2","3:5","6:7"], False]           # internal forces, see above
                        },
        "MISESBEAM":  { "B23E":    [2, [], True],                          # dummy
                        "B23I":    [2, ["2:2","3:3","12:12","13:13"], True],
                        "B23EI":   [2, ["2:2","3:3", "5:5","12:12"], True],     # normal force, bending moment, bond stress
                        "BAX23I":  [2, ["2:2","3:3","12:12","13:13"], True],
                        "BAX23EI": [2, ["2:2","3:3","12:12","13:13"], True],    # normal force, bending moment
                        "BAX23EI": [2, ["2:2","3:3","12:12","13:13"], True],    # normal force, bending moment
#                        "BAX23EI": [2, ["2:2","3:3"], True],    # normal force, bending moment
                        "BAX21EI": [2, ["3:3","4:4", "5:5","12:12"], True],        # normal force, bending moment, shear force, bond stress
                        "BAX21E":  [2, [], True]
                        #                        "BAX23EI": [2, ["2:2", "3:3"], True]
                        },
        "SPRING":    {  "S1D2":  [2, ["1:1"], True]                      # force

                     },
        "RCBEAM":    { "B23E":  [2, [], True],                           # currently dummy
                       "B23":   [2, [], True],
                       "B21":   [2, [], True],
                       "B21E":  [2, [], True]
                     },
        "NLSLAB":    { "SB3":   [2, ["3:5","6:7"], False]               # bending moments, shear forces

                    },
        "RCSHELL":  {  "SH4":   [ 3, ["0:2","3:5","6:7"], False]        # marker Z: nx, ny, nxy, mx, my, mxy, qx, qy
                                                                        # marker C: inplane strains, inplane stresses, principal stress values & directions
                                                                        # marker R: local uniaxial strain stress yield limmit, principal stress values & directions
                       },
        "MICRODAMAGE": { "C3D8":  [3, ["6:11"], False ],
                         "CPS4":  [2, ["4:7"], False ],
                         "CAX4":  [2, ["4:7"], False ]
        },
        "RESHEET":     { "CPS4":  [2, ["4:7"], False ]

        }
#        "BOND":      { "Bond2D3": [2,4, ["",""]] }
    }
    # data frame for plotting
    ElemResults = pd.DataFrame({
        "Label" : [ -1 ] ,
        "ElSet" : [""],
        "Material": [""],
        "Type":  [""],
        "IntPointIndex": [-1] ,
        "IntPointCoordinates": [[]],
        "Data": [[]],
        "Slices": [[]],
        "NoEl": [False],
        "Marker": [""]                                                      # relevant for SH3, SH4 only
#        "Active":[False]
    })
    ElemResultsRCShell = ElemResults.copy()

    SetSlices = {}                                                          # to hold slices of a set
    # read from elemout -- current line should be line with time value
    z1 = ff.readline()
    z2 = z1.strip()
    z3 = z2.split(',')
    if z3[0]!="Time": raise NameError("ConFemPostProcStandAlone: something wrong with elemout")
    else:             Time = float(z3[1])
    while z1!="":                                                           # current line is dataline
        pos= ff.tell()
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split(',')
        if z1   =="":                   break
        if z3[0]=='Time': ff.seek(pos); break

        # assign data to data frame
        mat = z3[2].strip()                                                 # material name
        mty = MatList[mat].Type                                             # material type
        eltype = z3[3].split(":")[0].strip()                                # element type - to get rid of markers in case of SH3, SH4
        try:    marker = z3[3].split(":")[1].strip()                        # Z, R, C in case of  SH3, SH4 only
        except: marker = ""
        if mty in ElResultsL:                                               # material type from dictionary above
            if eltype in ["Bond2D2","Bond2D3","Bond3D2","Bond3D3","BondAX2","BondAX3"]:
                pass
            elif eltype in ElResultsL[mty]:                                 # refers to element type in dictionary above
                l1 = ElResultsL[mty][eltype][0]                             # items for integration point
                if z3[-1]=="": l2 = len(z3)-l1-6
                else:          l2 = len(z3)-l1-5
                sl = ElResultsL[mty][eltype][1]                             # slices for data
                ef = ElResultsL[mty][eltype][2]                             # consider no of elements for plot scaling
                se = z3[1].strip()                                          # set label
                # sequence of data corresponds to definition of data frame
                if marker not in ["R","C"]:                                 # excludes SH4, SH3 with material *RCSHELL
                    #                                           Label,     ElSet,                   IntPointIndex, IntPointCoordinates,         Data - includes all slices
                    ElemResults.loc[len(ElemResults.index)] = [ int(z3[0]),z3[1].strip(),mat,eltype,z3[4],[float(z3[i]) for i in range(5,5+l1)],
                                                                [float(z3[i]) for i in range(5+l1,5+l1+l2)],sl,ef,marker]
                else:
                    ElemResultsRCShell.loc[len(ElemResultsRCShell.index)] = [ int(z3[0]),z3[1].strip(),mat,eltype,z3[4],[float(z3[i]) for i in range(5,5+l1)],
                                                                [float(z3[i]) for i in range(5+l1,5+l1+l2)],sl,ef,marker]
                if se not in SetSlices:
                    SetSlices[se] = sl
            else:
                raise NameError("ConFemPost::ReadResultsFromFile: not defined in dictionary for element results",mty, eltype)
        elif mty in ["BOND"]:
            pass
        else:
            raise NameError("ConFemPost::ReadResultsFromFile: not defined in dictionary for element results", mty)

    if z1=="": EndFlag = True
    # read from nodeout
    z1 = ffN.readline()
    z2 = z1.strip()
    z3N= z2.split(',')
    if z3N[0]!="Time": raise NameError("ConFemPostProcStandAlone: something wrong with nodeout ",z3N)
    else:             TimeN = float(z3N[1])
    while z1!="":
        pos= ffN.tell()
        z1 = ffN.readline()
        z2 = z1.strip()
        z3N= z2.split(',')
        if z1    =="":                    break
        if z3N[0]=='Time': ffN.seek(pos); break
        key  = int(z3N[0])                                                  # int node label, should be unique
        nDof = int(z3N[4])                                                  # number of degrees of freedom
        Data = []
        for i in z3N[5:5+nDof]: Data += [float(i)]
        NodeResults[key] = [ key, nDof, Data]
    if z1=="": EndFlagN = True
    #
    return EndFlag, EndFlagN, Time,TimeN, NodeResults,ElemResults,SetSlices,ElemResultsRCShell

def PickleLoad( FilDir, FilName, StepCounter):
    import pickle
#    from scipy.spatial import KDTree
    import scipy.spatial as spatial
#    from os import path as pth
    if pth.isfile(FilDir+FilName+".pkl"):
        fd = open(FilDir+FilName+'.pkl', 'rb')                              # has to be in sync with pickle.dumo in ConFem, ConSimFem
        NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);FlElasticLT=pickle.load(fd);\
            VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
            VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
            TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);Step=pickle.load(fd);                Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);\
            NoLabToNoInd=pickle.load(fd);NoIndToCMInd=pickle.load(fd);ContinuumNodes=pickle.load(fd);CoNoToNoLi=pickle.load(fd);SecDic=pickle.load(fd);LinAlgFlag=pickle.load(fd);ResultTypes=pickle.load(fd);Header=pickle.load(fd); \
            StepRestart = pickle.load(fd);MaxType = pickle.load(fd);MaxEquiIter=pickle.load(fd);StabSys=pickle.load(fd);StabTolF=pickle.load(fd);SoftSys=pickle.load(fd);SoftRed=pickle.load(fd);
        fd.close()
#        s = StepList[StepCounter]
        if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes ) # for search purposes, e.g. for EFG or aggregates or embedded truss elements
        else:                     CoorTree = None
    else: 
        raise NameError ("ConFemPostProcStandAlone: .pkl does not exit - correct name ?",FilDir+FilName)
    return NodeList, ElList, MatList, CoorTree, CoNoToNoLi, VecU, NoIndToCMInd, NoLabToNoInd, SecDic, LinAlgFlag, ResultTypes, Time, Header

def PlotContour( Label, x, y, z):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    p0.axis('equal')
    p0.tricontour(x, y, z, linewidths=0.5, colors='k')
    P0.colorbar(p0.tricontourf(x, y, z, cmap="bone"), ax=p0)
    p0.plot(x, y, 'ko', ms=1)

def PlotResiduals( LogName, incr, iter, scale):
    def Stuff( Pl, Name, scale, scale2):
        ff, xx, yy, rr = open( Name, 'r'), [], [], 0.
        z1 = ff.readline()                                  # 1st input line
        z2 = z1.split()
        while z1!="":
            label, x, y, rx, ry = int(z2[0]), float(z2[1]),float(z2[2]), float(z2[5]), float(z2[6])
            rr = rr + rx**2 + ry**2
            xx += [x]
            yy += [y]
            Pl.plot([x,x+scale*rx],[y,y+scale*ry],color='red')
            Pl.text(x,y,"%i"%label,ha='left',va='center',rotation=0.,color='black',fontsize=6) #7)
#            Pl.plot([x,x+scale*ry],[y,y+scale*rx],color='red')
    #        Pl.plot([x,x+scale2*ax],[y,y+scale2*ay],color='blue')
            z1 = ff.readline()                                  # 1st input line
            z2 = z1.split()
        ff.close()
        print(sqrt(rr))
        Pl.plot(xx,yy,'x')
        Pl.set_aspect('equal')
        Pl.grid()
        return
    Name = LogName+"/log-"+str(incr)+"-"+str(iter)+".txt"
    Pl1 = plt.figure().add_subplot(111,title='residuum increment '+str(incr)+' iteraton '+str(iter))
    Stuff( Pl1, Name, scale, 1.0e1) # 1e1
#    Stuff( Pl1, Name, 2.0e+3, 1.0e3) # 1e1
#    Name = LogName+"/logD-"+str(counter)+"-"+str(i)+".txt"
#    Pl2 = plt.figure().add_subplot(111,title='Displ '+str(counter)+'/'+str(i))
#    Stuff( Pl2, Name, 1.0e2, 1.0e1)
    plt.show()
    return

def WriteResiduals( LogName, incr, iter, N=100):
    class Residual():                                  # 1D Truss 2 nodes
        def __init__(self, Label, X, Y, Z, Res):
            self.Label = Label
            self.X   = X
            self.Y   = Y
            self.Z   = Z
            self.Res = Res
    Name = LogName+"/log-"+str(incr)+"-"+str(iter)+".txt"
    Nam1 = LogName+"/logWr-"+str(incr)+"-"+str(iter)+".txt"
    ResList = []
    ff  = open( Name, 'r')
    f1  = open( Nam1, 'w')
    z1 = ff.readline()
    z2 = z1.split()
    while z1!="":
        nL, x, y, z, _, rr  = int(z2[0]), float(z2[1]), float(z2[2]), float(z2[3]), int(z2[4]), float(z2[-1])
        ResList += [Residual( nL, x, y, z, rr)]
#        f1.write("%5i ; %8.2f ; %8.2f ; %8.2f ; %16.6e\n"%(nL,x,y,z,rr))
        z1 = ff.readline()
        z2 = z1.split()
    ff.close()
    n = len(ResList)
    ResList.sort(key=lambda t: t.Res)
    print(Name)
    for i in range(N):
        r = ResList[n-1-i]
        f1.write("%5i ; %8.2f ; %8.2f ; %8.2f ; %16.6e\n"%(r.Label,r.X,r.Y,r.Z,r.Res))
#        print(i+1, r.Label, r.Res)
    f1.close()

class ConFemPost:
    def __init__(self):
        pass
    def Run(self, FilDir,FilName, Name,StepCounter, MaPlLib, VTK, Version, ResiData):
    # load system data
#        NodeList,ElemList,MatList,CoorTree,CoNoToNoLi,VecU,NoIndToCMInd,NoLabToNoInd, SecDic, LinAlgFlag, ResultTypes, LastTime, Header = PickleLoad( FilDir, Name, StepCounter)
        NodeList,ElemList,MatList,CoorTree,CoNoToNoLi,VecU,NoIndToCMInd,NoLabToNoInd, SecDic, LinAlgFlag, ResultTypes, LastTime, Header = PickleLoad( "", Name, StepCounter)
        # process plot options
        ScaleDis, PE2DFlag, PE3DFlag = 1.0, True, False
        Post1DFlag, PlotTimes, ScaleStress, PostNodeFlag, ScaleDis2D, ShellL, Contour2D = True, [], [], True, 0., False, {}
        #
        if pth.isfile(Name+".plt.txt"):                                     # read plot options file if there is any
            f1=open( Name+".plt.txt", 'r')
            from ConFemInOut import ReadPlotOptionsFile
            ScaleDis, ScaleDis2D, PE2DFlag, PE3DFlag, PlotTimes, Post1DFlag, ScaleStress, PostNodeFlag, ShellL, ScaleShellL, Contour2D = ReadPlotOptionsFile(f1, SecDic) # to modify plot scaling factors
            f1.close()
        # 
        if Version == 1:
            if Post1DFlag:
#                f2= open( Name+".elemout.txt", "r")
#                PostElem1D( f2,  ResultTypes, PlotTimes)
#                f2.close()
                f2= open( Name+".elemout.txt", "r")
                PostElem1D_(f2, ResultTypes, PlotTimes)
                f2.close()
            #
            if PostNodeFlag:
                PostNode( ElemList, NodeList,NoIndToCMInd, VecU, ScaleDis, None)
                print('PostNode   process plot for last time', "%6.4f"%LastTime)
            #
#            try:
            if True:
                ff  = open(Name+".elemout"+".txt",'r')
                ffN = open(Name+".nodeout"+".txt",'r')
                EndFlag, EndFlagN = False, False
                # loop over times in result data sets
                while not EndFlag and not EndFlagN:
                    EndFlag, EndFlagN, Time, TimeN, NodeResults, ElResults, SetSlices, ElResultsRCShell = ReadResultsFromFile( MatList, ff, ffN)
                    if Time != TimeN: raise NameError("ConFemPostProcStandAlone: async time for elements and nodes",Time, TimeN)
                    print("ElemData found time",Time)
                    if Time in PlotTimes or len(PlotTimes)==0:
                        if PE2DFlag:
                            SCa_= PostScales( SecDic, ElemList, ElResults, ScaleStress )
                            _ = PostElem2D( Header,ElemList,NodeList, ScaleDis2D,Contour2D, SecDic, MaPlLib, NoIndToCMInd, NodeResults,ElResults,SetSlices, Time, FilDir+"/", SCa_)
                        if PE3DFlag:
                            PostNode3D( ElemList,NodeList,NoIndToCMInd, NodeResults,[], Time, ScaleDis)
                            l = PostElem3D( ElemList, NodeList,NoIndToCMInd, ElResults, Time)
                            if l>0: print('PostElem3D process plot for time ', Time)
                        if VTK:
                            PostElemVTK( Name, str(int(1000*round(Time,5))), ElemList, NodeList,NoIndToCMInd, ElResults,NodeResults)
                        if ShellL:
                            PostElemSH4(SecDic, Time, ScaleShellL, ElResultsRCShell)
                ff.close()
                ffN.close()
            #
            if pth.isfile(Name+".opt.txt"):
                f4=open( Name+".opt.txt", 'r')
                WrNodes,_,_,_,_,_,_,_ = ReadOptionsFile(f4, NodeList, NoLabToNoInd, NoIndToCMInd)
                f4.close()
            else: WrNodes = []
            if len(WrNodes)>0:
                f5=open( Name+".timeout.txt", 'r')
                print('PostHist process plot for history')
                PlotHist( f5, WrNodes )
                f5.close()
            #
            return True
        elif Version == 2:
            LogName="../LogFiles"                                               # to log temporary data
            incr  = ResiData[0]
            iter  = ResiData[1]
            scale = ResiData[2]
            PlotResiduals( LogName, incr, iter, scale)
            WriteResiduals( LogName, incr, iter)
            return True
        else:
            print("nothing to do")

    ############################################################################################################################### # ts
if __name__ == "__main__":
    StepCounter = 0
    VTK,VTK3D = False, False
    # ScP scaling factor principal stresses, ScU scaling factor displacements
    #                                                                  [output times], scaling factor displacements, [[Section 1 scaling factors], [Section 2 scaling factors], .... ]
    #                                                                                                                each Section may display more than one result depending on element type
    # DataExamples
    FilDir = ""
#    Name = "../DataExamples/E03/E3-01"  # 01, 04
#    Name = "../DataExamples/E04/E4-06".
#    Name="../DataExamples/E04/E3-02_CircMises"                              # 2_CircMises"
#    Name ="../DataExamples/E05/E5-02"
#    Name = "../DataExamples/E07/E7-01b"
#    Name = "../DataExamples/E07/E7-05"                                   # 1D phase field
#    Name, VTK = "../DataExamples/E07-06/E7-06a", True                               # SDA
#    Name = "../DataExamples/E08-02/ElasticLT/E8-02"                      # deep beam nonlinear
#    Name, VTK = "../DataExamples/E08/E8-03B23e", False
#    Name, VTK = "../DataExamples/E08/E8-04", False
#    Name, VTK = "../DataExamples/E08-02/E8-02", False
#    Name = "../DataExamples/E09/E9-01"
#    Name = "../DataExamples/E09/E9-05"
#    Name = "../DataExamples/E10/E10-02a"
#    Name = "../DataExamples/E10-02/E10-02"
#    Name = "../_DataTrusses/staebe3D2"    # staebe_4, staebe3D2
#    Name = "../_DataBenchmarks/Arrea/arrea2D"
#    Name = "../_DataBenchmarks/Nooru/Nooru-1550-5"
#    Name = "../_DataTmp/Deep_beam_AacNL"
#    Name = "../_DataBenchmarks/Aachentests/Deep_beam"
#    Name = "../_DataShellsSlabs/IuriisExamples/slab_SH4"
#    Name ="../DataExamples/E10-02/IsoDam/E10-02"

#    Name, VTK,VTK3D = "../_DataC3D/Cube8", True, True
#    Name, VTK, VTK3D = "../_DataC3D/Deep_beam_SDA", True, True
    #
#    Name = "../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load" # Shell, bridge_el05m, c_1461(0.08)_2.1e5_0.3_segment load
    #
#    Name = "C:/Users/uhc/Documents/Work/FoilPap/2023/Note_ShearPlateRandom/ConFem/ShearPanel/Restart/ShearPanelR.89"
#    Name = "C:/Users/uhc/Desktop/Note_FlatSlab_Comp/ExpDataLandler_25-02-07/1-M0-25-1.23/1-M0-25-1.23"
#    Name = "C:/Users/uhc/Desktop/Note_FlatSlab_Comp/ExpDataLandler_25-01-30/1-M0-25-1.23/1-M0-25-1.23"
#    Name = "../_DataBond/PulloutAxiSym"
    Name = "../_DataTmp/E10-02_20/E10-02"
#    Name = "../_DataTmp/1-M0-25-1.23"

    DirName, FilName = DirFilFromName( Name )
    #
    print('ConFemPost for ',Name)
    Version = 1                                                             # 1: plot results, 2: plot residuals
    ResiData = [ 30, 18, 5.0e+01]                                            # incr, iter, scale
    ComFemPost_ = ConFemPost()
    MaPlLib = True                                                          # flag for using matplotlib
    Flag = ComFemPost_.Run( DirName,FilName, Name,StepCounter, MaPlLib, VTK, Version, ResiData)
    print('finish processing')
    if MaPlLib or Flag: plt.show() 
    print('finish all')
