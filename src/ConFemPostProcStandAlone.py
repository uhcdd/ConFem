# CaeFemMain -- 2017-04-19
# Copyright (C) [2017] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
from numpy import dot, array, mean, sqrt, arcsin, zeros
#from numpy.linalg import norm, det
from ConFemElem import ZeroD, C3D8
from ConFemInOut import ReadOptionsFile, PlotNodes
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path as pth
import numpy as np
import matplotlib.colors as colors
import mpl_toolkits.mplot3d as a3
import pylab as pl
import time as ti
import sys
import os
sys.path.insert(1, '../src')
#import imp
from ConFemBasics import SamplePoints

Colors = ['blue','black','darkviolet','magenta','indianred','darkred','cyan','yellow','blue']
Colors2 = ['red','green','magenta','cyan']
ColorTable = {}
ColorTable['IsoD'] = 'white'
ColorTable['Elastic'] = 'black'
ColorTable['Bond'] = 'green'
ColorTable['ElasticLT'] = 'gray'
ColorTable['Mises'] = 'red'
ColorTable['MicroPl'] = 'magenta'
ColorTable['Lubliner'] = 'cyan'
ColorTable['IsoDam'] = 'cyan'
ColorTable['ElasticSDA'] = 'magenta'

def DefPlot(Label):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    return P0, p0 

def Post2D( ElemList,NodeList,MatList, ScaleU,ScaleS, SolSecs, PostFlagStress,PostFlagMaPlLi, NoIndToCMInd, ElemResults,NodeResults, TimeEl, pdir ):
    def Annot( ):
        if Xmax!=None:
            p1.annotate(    'max s1 '+"%6.1e"%Smax+' s2 '+"%6.1e"%S2max+' '+'%6i'%(ELaMax),xy=(Xmax, Ymax),xycoords='data',xytext=(0.40,0.87),textcoords='figure fraction',arrowprops=dict(arrowstyle="->"),fontsize='medium')
#            if ElType=='SH4': 
#                P1.annotate('max s1 '+'%7.4f'%(Smax_)+' s2 '+'%7.4f'%(S2max_),xy=(Xmax_, Ymax_),xycoords='data',xytext=(0.50,0.87),textcoords='figure fraction',arrowprops=dict(arrowstyle="->"),fontsize='large')
        if Xmin!=None:
            p1.annotate(    'min s1 '+"%6.1e"%Smin+' s2 '+"%6.1e"%S2min+' '+'%6i'%(ELaMin),xy=(Xmin, Ymin),xycoords='data',xytext=(0.40,0.84),textcoords='figure fraction',arrowprops=dict(arrowstyle="->"),fontsize='medium')
#            if ElType=='SH4': 
#                P1.annotate('min s1 '+'%7.4f'%(Smin_)+' s2 '+'%7.4f'%(S2min_),xy=(Xmin_, Ymin_),xycoords='data',xytext=(0.50,0.84),textcoords='figure fraction',arrowprops=dict(arrowstyle="->"),fontsize='large')
    def PloPrin( pp, scaleP, XX, YY, p0, Smax, Smin, ElLabel, ELaMax_, ELaMin_):
        ELaMax, ELaMin, Xmax, Ymax, Xmin, Ymin = ELaMax_, ELaMin_, 0, 0, 0, 0
        Smax_ , Smin_, S2min, S2max = Smax, Smin, 0, 0
        SS = pp[0]
        if SS>=0: 
            if PostFlagMaPlLi: p0.plot([                                          XX-SS*scaleP*pp[1],XX+SS*scaleP*pp[1]],[YY-SS*scaleP*pp[2],YY+SS*scaleP*pp[2]],'r-')
            else:              ffred.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n"%( XX-SS*scaleP*pp[1],YY-SS*scaleP*pp[2],  XX+SS*scaleP*pp[1],YY+SS*scaleP*pp[2] ))
        else:     
            if PostFlagMaPlLi: p0.plot([                                          XX-SS*scaleP*pp[1],XX+SS*scaleP*pp[1]],[YY-SS*scaleP*pp[2],YY+SS*scaleP*pp[2]],'g-')
            else:              ffgre.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n"%( XX-SS*scaleP*pp[1],YY-SS*scaleP*pp[2],  XX+SS*scaleP*pp[1],YY+SS*scaleP*pp[2] ))
        if SS>Smax: ELaMax, Xmax, Ymax, Smax_, S2max = ElLabel, XX, YY, SS, pp[3]
        if SS<Smin: ELaMin, Xmin, Ymin, Smin_, S2min = ElLabel, XX, YY, SS, pp[3]
        SS = pp[3]
        if SS>=0: 
            if PostFlagMaPlLi: p0.plot(                                          [XX-SS*scaleP*pp[2],XX+SS*scaleP*pp[2]],[YY+SS*scaleP*pp[1],YY-SS*scaleP*pp[1]],'r-')
            else:              ffred.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n"%( XX-SS*scaleP*pp[2],YY+SS*scaleP*pp[1],  XX+SS*scaleP*pp[2],YY-SS*scaleP*pp[1] ))
        else:     
            if PostFlagMaPlLi: p0.plot(                                          [XX-SS*scaleP*pp[2],XX+SS*scaleP*pp[2]],[YY+SS*scaleP*pp[1],YY-SS*scaleP*pp[1]],'g-')
            else:              ffgre.write("%12.4e  %12.4e\n%12.4e  %12.4e\n\n"%( XX-SS*scaleP*pp[2],YY+SS*scaleP*pp[1],  XX+SS*scaleP*pp[2],YY-SS*scaleP*pp[1] ))
        if SS>Smax: ELaMax, Xmax, Ymax, Smax_, S2max = ElLabel, XX, YY, SS, pp[0]
        if SS<Smin: ELaMin, Xmin, Ymin, Smin_, S2min = ElLabel, XX, YY, SS, pp[0]
        return ELaMax, ELaMin, Xmax, Ymax, Smax_, S2max, Xmin, Ymin, Smin_, S2min
    
    if not PostFlagMaPlLi: import PyGnuplot as gp                                       # use gnuplot instead of matplotlib
    # store values of ndof per node -- if NodeResults exists use data from nodeout files, otherwise from actual displacement vector
    DisplList = NodalDisplList( NodeList, NodeResults)

#    # echo bond elements
#    if False:
#        FF = 0.
#        for i in ElemList:
#            if i.Type in ["B2D2E","B2D3E"]:
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
                    
    P1_, p1_, gi = [], [], 0                                                # gi counter for gnu plot windows
    # loop over all sections
    for s_, S in enumerate(SolSecs):
        print('post process solid section',s_,SolSecs[S].Set,"Time",TimeEl)
        try:    ScaleSec =ScaleS[s_]
        except: raise NameError("ConFemPy3::Post2D: should provide scaling factors for this section",S)
        SigMax, SigMin, angle, maxX, minX = -9999, 9999, 0., -9999., 9999.
        ELaMax, ELaMin, Xmax, Ymax, Xmin, Ymin = 0, 0, 0, 0, 0, 0
        if PostFlagMaPlLi:
            P1, p1 = DefPlot("forces / stresses "+ "%2i "%s_ +SolSecs[S].Set+ " Time "+"%7.4f"%TimeEl )
            p1.axis('equal')
            P1_ += [P1]
            p1_ += [p1]
        else: 
            p1 = None
            ffg = open("tmp"+str(S)+".dat", 'w')                            # temporary file for gnuplot elements
            ffc = open("tmp"+str(S)+"c.dat", 'w')                           # temporary file for gnuplot cracks
            ffred = open("tmp"+str(S)+"red.dat", 'w')                       # temporary file for gnuplot principal stress tensile
            ffgre = open("tmp"+str(S)+"gre.dat", 'w')                       # temporary file for gnuplot principal stress compressive
            ffIP = open("tmp"+str(S)+"IP.dat", 'w')                         # temporary file for gnuplot integration points + color for stress level
            ffIP2= open("tmp"+str(S)+"IP2.dat", 'w')                        # temporary file for gnuplot integration points 2 + color for stress level
        xxMin, xxMax, yyMin, yyMax = 999., -999., 999., -999.
        # deformed mesh + cracks
        for jj in SolSecs[S].Elems:                                         # loop over element indices of Solid Section
            i = ElemList[jj]
            if i.Active and not i.Type in ["B2D2E","B2D3E"]:
                key = str(i.Label)+i.Set+"0"                                # refer to 1st integration point
                xy = i.NodalCoordinates( NodeList, NoIndToCMInd )
                _, uxy =  ElementDisplList( i, DisplList, NoIndToCMInd, ScaleU)
                #
                if i.Type in ["T2D2","T2D2E","T2D3","T2D3E"]:
                    samp, xx, yy = [[-1.,0.,0.],[1.,0.,0.]], [], []         # samp: local sampling points for nodes, ZeroD for isogeom
                elif i.Type in ["CPS3",'CPS3S',"CPE3",'CPE3S','CPS6','CPE6','CPS6S','CPE6S']:
                    samp, xx, yy = [[1., 0., 0.],[0., 1., 0.],[0., 0., 1.]], [], [] # samp: local sampling points for nodes CPS3
                else:    
                    samp, xx, yy = [[-1.+ZeroD,-1.+ZeroD,0.],[1.-ZeroD,-1.+ZeroD,0.],[1.-ZeroD,1.-ZeroD,0.],[-1.+ZeroD,1.-ZeroD,0.]], [], []    # samp: local sampling points for nodes, ZeroD for isogeom
                for j in range(len(samp)):
                    r, s, t = samp[j][0], samp[j][1], samp[j][2] 
                    XX  = i.FormX_( r, s, t)
                    xyP = dot(XX,array(xy))                                 # undeformed coordinates
                    NN  = i.FormN( r, s, 0.)
                    xyP = xyP + dot(NN,array(uxy))                          # displaced coordinates
                    xx += [xyP[0]]
                    yy += [xyP[1]]
                    if xyP[0]<minX: minX=xyP[0]
                    if xyP[0]>maxX: maxX=xyP[0]
                xc, yc = mean(xx), mean(yy)                                 # coordinates of element center
                xx+= [xx[0]]                                                # to close the sequence of edges
                yy+= [yy[0]]                                                # to close the sequence of edges
                if min(xx) < xxMin: xxMin=min(xx)
                if max(xx) > xxMax: xxMax=max(xx)
                if min(yy) < yyMin: yyMin=min(yy)
                if max(yy) > yyMax: yyMax=max(yy)
                # start def
                def Plot2DCrackContour(CrXY0, CrXY1):
                    r0, s0, t0 = i.NInverse( xy, CrXY0)                     # local coordinate 'left' edge of crack
                    r1, s1, t1 = i.NInverse( xy, CrXY1)                     # local coordinate 'right' edge of crack
                    CrdXY0 = dot(i.FormN( r0, s0, t0),array(uxy))           # displacement 'left' edge of crack
                    CrdXY1 = dot(i.FormN( r1, s1, t1),array(uxy))           # displacement 'right' edge of crack
                    xx = [ CrXY0[0]+CrdXY0[0], CrXY1[0]+CrdXY1[0] ]
                    yy = [ CrXY0[1]+CrdXY0[1], CrXY1[1]+CrdXY1[1] ]
                    if PostFlagMaPlLi: p1.plot(                                    [ xx[0],xx[1] ], [ yy[0],yy[1] ],color='gray')
                    else:              ffc.write("%8.4f  %8.4f\n%8.4f  %8.4f\n\n"%(  xx[0],yy[0],     xx[1],yy[1] ))
                # end of def
                if key in ElemResults:
                    if PostFlagMaPlLi: 
                        p1.plot(xx,yy,'--',color=Colors[s_])                # matplotlib plot deformed element geometry
#                        p1_[0].plot(xx,yy,'--',color=Colors[s_])           # 1st element set !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                        try:    p1_[1].plot(xx,yy,'--',color=Colors[s_])           # 1st element set !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                        except: pass
                    else:
                        for j in range(len(xx)): ffg.write("%8.4f  %8.4f\n"%(xx[j],yy[j]))  # gunplot write temp deformed mesh data 
                        ffg.write("\n")
                    # 2D SDA elements
                    if ElemResults[key][3] in ['CPS4S','CPE4S','CPS3S','CPE3S','CPS6S','CPE6S']: 
                        key = str(i.Label)+i.Set+"Cra"                      # key first crack
                        if key in ElemResults:
                            Data = ElemResults[key][5]
                            crX, crY, crNx, crNy = Data[0], Data[1], Data[2], Data[3]  # center x, center y, nx, ny
                            if PostFlagMaPlLi: p1.plot([xc], [yc],'o',color='gray') # matplotlib mark cracked elements
                            xyC = i.FindGlobalDiscontinuityEdges( [crX,crY], [crNx,crNy], NodeList, NoIndToCMInd) # xyC: [[x, y of 1st edge],[x, y, of 2nd edge]]
                            Plot2DCrackContour( xyC[0], xyC[1] )
                        key = str(i.Label)+i.Set+"Cra2"                     # key second crack
                        if key in ElemResults:
                            Data = ElemResults[key][5]
                            crX, crY, crNx, crNy = Data[0], Data[1], Data[2], Data[3]  
                            if sqrt(crNx**2+crNy**2)>ZeroD:                 # crack is there
                                if PostFlagMaPlLi: p1.plot([xc], [yc],'o',color='red')
                                xyC = i.FindGlobalDiscontinuityEdges( [crX,crY], [crNx,crNy], NodeList, NoIndToCMInd) # xyC: [[x, y of 1st edge],[x, y, of 2nd edge]]
                                Plot2DCrackContour( xyC[0], xyC[1] )
        # end of loop over deformed elements and cracks
#    label for elements, generally commented out
#                p1.text(xc,yc,"%i"%i.Label,ha='center',va='center',rotation=0.,color='black',fontsize=7)
        # principal stresses on deformed geometry
        if PostFlagStress:
            pDataX, pDataY = {}, {}
            Truss2DCoor, Truss2DStr, Truss2DStrBond,  = [], [], []
            for jj in SolSecs[S].Elems:                                     # loop over element indices of Solid Section
                i = ElemList[jj]
                if i.Active and not i.Type in ["B2D2E","B2D3E"]:
                    for j in range(len(ScaleSec)): pDataX[j], pDataY[j] = [], []  # storage for x,y of 2D plot per each scaling factor
                    xy = i.NodalCoordinates( NodeList, NoIndToCMInd )
                    _, uxy =  ElementDisplList( i, DisplList, NoIndToCMInd, ScaleU)
                    # principal stresses in integration points
                    if i.Type in ['T2D2','T2D2E','T2D3E','B23I']: 
                        if i.cosA<0: angle =  90 - arcsin(i.sinA)*180/3.141593
                        else:        angle = -90 + arcsin(i.sinA)*180/3.141593
#                    if i.Type in ['CPS4S','CPE4S','CPS3S','CPE3S']: ww = i.ww
                    nInt, nIntL = i.nInt-1, i.nIntL 
                    for j in range(nIntL):                                  # loop over integration points of element
                        r = SamplePoints[i.IntT,nInt,j][0]
                        s = SamplePoints[i.IntT,nInt,j][1]
                        t = SamplePoints[i.IntT,nInt,j][2]
#                        val = None                                          # some element types need additional values for B-matrix, e.g. XFEM
                        key = str(i.Label)+i.Set+str(j)
                        Data= ElemResults[key][5]
                        # data specific for element type
                        if   ElemResults[key][3] in ["CPS4","CPE4","CPS4S","CPE4S","CPS3","CPE3","CPS3S","CPE3S","CPS6","CPE6","CPS6S","CPE6S"]: 
                            sig = [Data[6],Data[7],Data[9]]
                        elif ElemResults[key][3] in ["T2D2","T2D3"]: 
                            sig = [Data[3]]                                 # Data should hold all floating point from elemout, [3] long. stress, [12] bond stress 
                        elif ElemResults[key][3] in ["T2D2E","T2D3E"]: 
                            sig = [Data[3],Data[12],Data[13]]               # Data should hold all floating point from elemout, [3] long. stress, [12] bond stress, [13] penalty stress 
                        elif ElemResults[key][3] in ["B23I"]: 
                            sig = [Data[6],Data[7],Data[12],Data[13]]  # low and up stress, long. and lat. bond stress
                        else: raise NameError("CaeFem3::Post2D: unknown element type for reading from files",ElemResults[key][3])
                        #
                        XX  = i.FormX_( r, s, t)
                        xyP = dot(XX,array(xy))                             # undeformed global coordinates of integration point
                        NN  = i.FormN( r, s, t)
                        xyP = xyP + dot(NN,array(uxy))                      # displaced global coordinates of integration point
                        CoLor = ColorTable[MatList[i.MatN].Type]
                        if PostFlagMaPlLi:                                         
#                            p1.plot([xyP[0]], [xyP[1]],'x',color = 'darkred')       # matplotlib cross in integration point
#                            p1_[0].plot([xyP[0]], [xyP[1]],'x',color = 'darkred')  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            pass
                        elif ElemResults[key][3] in ["T2D2","T2D2E","T2D3","T2D3E"]:
                            Truss2DCoor += [[xyP[0], xyP[1]]] 
                            Truss2DStr  += [sig[0]]                         # longitudinal stress
                            Truss2DStrBond  += [abs(sig[1])]                # bond stress
                        def PlotRotLinesText( ind, sign_, pDataX, pDataY):
                            if len(ScaleSec) != ind: raise NameError("ConFemPy3::Post2D: number of scaling factors does not match for this section",S)
                            phiC, phiS = i.cosA, i.sinA
                            sign_ = [1,-1,-1]
                            for jj in range(ind):
                                pDataX[jj] += [ xyP[0] - sign_[jj]*ScaleSec[jj]    *phiS*sig[jj] ]
                                pDataY[jj] += [ xyP[1] + sign_[jj]*ScaleSec[jj]    *phiC*sig[jj] ]
                                p1.text(pDataX[jj][-1],pDataY[jj][-1],format(sig[jj],".2f"),ha='left',va='top',rotation=angle,color=Colors2[jj],fontsize='medium')
#                            return PDataX, PDataY
                        if i.Type in ['T2D2','T2D3']:
                            if PostFlagMaPlLi: PlotRotLinesText( 1, [1], pDataX, pDataY)
                        elif i.Type in ['T2D2E','T2D3E']:
                            if PostFlagMaPlLi:                              # matplotlib truss elements
                                PlotRotLinesText( 3, [1, -1, -1], pDataX, pDataY)
#                                if len(ScaleSec) != 3: raise NameError("ConFemPy3::Post2D: number of scaling factors does not match for this section",S)
#                                phiC, phiS = i.cosA, i.sinA
#                                sign_ = [1,-1,-1]
#                                for jj in range(3):
#                                    pDataX[jj] += [ xyP[0] - sign_[jj]*ScaleSec[jj]    *phiS*abs(sig[jj]) ]
#                                    pDataY[jj] += [ xyP[1] + sign_[jj]*ScaleSec[jj]    *phiC*abs(sig[jj]) ]
#                                    p1.text(pDataX[jj][-1],pDataY[jj][-1],format(sig[jj],".2f"),ha='left',va='top',rotation=angle,color=Colors2[jj],fontsize='medium')
                        elif i.dim==2:                                              # 2D continuum elements
                            pp = PrinC( sig[0], sig[1], sig[2] )
                            ELaMax_, ELaMin_, Xmax_, Ymax_, Smax, S2max_, Xmin_, Ymin_, Smin, S2min_ = PloPrin( pp, ScaleSec[0], xyP[0], xyP[1], p1, SigMax, SigMin, i.Label, ELaMax, ELaMin)
                            if ELaMax != ELaMax_: SigMax, ELaMax, Xmax, Ymax, S2max = Smax, ELaMax_, Xmax_, Ymax_, S2max_ 
                            if ELaMin != ELaMin_: SigMin, ELaMin, Xmin, Ymin, S2min = Smin, ELaMin_, Xmin_, Ymin_, S2min_
                        elif i.Type in ['B23I']:                                             # bernoulli beam
                            if PostFlagMaPlLi:
                                if len(ScaleSec) != 4: raise NameError("ConFemPy3::Post2D: number of scaling factors does not match for this section",S.Set)
                                phiC, phiS = i.cosA, i.sinA
                                sign_ = [1,1,-1,-1]
                                for jj in range(4):
                                    pDataX[jj] += [ xyP[0] - sign_[jj]*ScaleSec[jj]    *phiS*abs(sig[jj]) ]
                                    pDataY[jj] += [ xyP[1] + sign_[jj]*ScaleSec[jj]    *phiC*abs(sig[jj]) ]
                                    p1.text(pDataX[jj][-1],pDataY[jj][-1],format(sig[jj],".2f"),ha='left',va='top',rotation=angle,color=Colors2[jj],fontsize='medium')
                    if PostFlagMaPlLi:
                        for a in pDataX: p1.plot(pDataX[a],pDataY[a], color = Colors2[a])
                    
            # gnuplot truss element data
            if len(Truss2DCoor)>0:                                          
                Truss2DStrDel = max(Truss2DStr) - min(Truss2DStr) 
                if (Truss2DStrDel)>ZeroD:
                    for j in range(len(Truss2DCoor)): 
                        ffIP.write("%12.4e  %12.4e  %12.4e\n"%(Truss2DCoor[j][0],Truss2DCoor[j][1],Truss2DStr[j]))
                Truss2DStrDel = max(Truss2DStrBond) - min(Truss2DStrBond) 
                if (Truss2DStrDel)>ZeroD:
                    for j in range(len(Truss2DCoor)): 
                        ffIP2.write("%12.4e  %12.4e  %12.4e\n"%(Truss2DCoor[j][0],Truss2DCoor[j][1],Truss2DStrBond[j]))

        # gnuplot
        if not PostFlagMaPlLi:
            ffg.close()
            ffc.close()
            ffred.close()
            ffgre.close()
            ffIP.close()
            ffIP2.close()
            gp.c('set term wxt '+str(gi))
            gp.c('set size square')
            gp.c('set grid')
            gp.c('unset key')
            gp.c('set title "Set '+i.Set+' -- Time '+str(TimeEl)+' -- Scale '+str(ScaleU)+' " ')
            import os
            dxx, dyy = max(1.,xxMax-xxMin), max(1.,yyMax-yyMin)
            gp.c('set xrange['+str(xxMin-0.08*dxx)+':'+str(xxMax+0.08*dxx)+']')
            gp.c('set yrange['+str(yyMin-0.08*dyy)+':'+str(yyMax+0.08*dyy)+']')
            gp.c('plot "tmp'+str(S)+'.dat" with lines linestyle 1 linecolor rgb "#0060ad" linetype 1 linewidth 0.5 ')              # deformed elements
            if os.path.getsize('tmp'+str(S)+'c.dat')>0:   gp.c('replot "tmp'+str(S)+'c.dat" with lines linestyle 1 linecolor rgb "gray30" ') # 2D crack contours 
            if os.path.getsize('tmp'+str(S)+'red.dat')>0: gp.c('replot "tmp'+str(S)+'red.dat" with lines linestyle 1 linecolor rgb "red" ') # tensile principal stresses
            if os.path.getsize('tmp'+str(S)+'gre.dat')>0: gp.c('replot "tmp'+str(S)+'gre.dat" with lines linestyle 1 linecolor rgb "dark-green" ') # commpressive principal stresses
            if os.path.getsize('tmp'+str(S)+'IP.dat')>0:  gp.c('replot "tmp'+str(S)+'IP.dat" using 1:2:3 w p pointtype 7 pointsize 0.3 lc palette')
            gp.p(filename=pdir+'tmp'+str(S)+'.eps', width=20, height=12, fontsize=12, term='wxt')
            if os.path.getsize('tmp'+str(S)+'IP2.dat')>0:
                gi += 1
                gp.c('set term wxt '+str(gi))
                gp.c('set size square')
                gp.c('set grid')
                gp.c('unset key')
                gp.c('set title "Set '+i.Set+' -- Time ' +str(TimeEl)+' -- Scale '+str(ScaleU)+' " ')
                gp.c('plot "tmp'+str(S)+'.dat" with lines linestyle 1 linecolor rgb "#0060ad" linetype 1 linewidth 0.5 ')              # deformed elements
                gp.c('replot "tmp'+str(S)+'IP2.dat" using 1:2:3 w p pointtype 7 pointsize 0.3 lc palette')
                gp.p(filename=pdir+'tmp'+str(S)+'_.eps', width=20, height=12, fontsize=12, term='wxt')
#            gp.c('pause -1')
#            gp.c('a = system("read a; echo $a")')
            ti.sleep(2)
            _ = str(input('\ncontinue: '))
            os.remove('tmp'+str(S)+'.dat');os.remove('tmp'+str(S)+'c.dat');os.remove('tmp'+str(S)+'red.dat');os.remove('tmp'+str(S)+'gre.dat');os.remove('tmp'+str(S)+'IP.dat');os.remove('tmp'+str(S)+'IP2.dat'); 
        else:
#    label for nodes, generally commented out
#            for i in NodeList: 
#                if len(i.dL)==2: p1.text(i.XCo+ScaleU*i.dL[0],i.YCo+ScaleU*i.dL[1],"%i"%(i.Label-0),ha='left',va='bottom',rotation=0.,color='red',fontsize=7)
            for p in p1_: p.plot([minX-0.05*abs(maxX-minX),maxX+0.05*abs(maxX-minX)],[0,0],color='lightgrey')
            if i.dim==2 and PostFlagStress: Annot()
        gi += 1
    return 0
 
def Principal( NoList, ElemList, Name, timeStr):   
                        
    def Princstresses3D(sigma_x, sigma_y, sigma_z, tau_yz, tau_xz, tau_xy):
        
        AA = np.array([[sigma_x, tau_xy, tau_xz],
                       [tau_xy, sigma_y, tau_yz],
                       [tau_xz, tau_yz, sigma_z]])
            
        eigenvalue, eigenvector = np.linalg.eig(AA)
        eigenvalue_list = [eigenvalue[0], eigenvalue[1], eigenvalue[2]] 
        ma, mi = eigenvalue_list.index(max(eigenvalue)), eigenvalue_list.index(min(eigenvalue)) 
        v1= np.array(np.multiply(eigenvalue[ma], eigenvector[:,ma]))
        v3= np.array(np.multiply(eigenvalue[mi], eigenvector[:,mi]))
        for i in range(3):
            if i != ma and i != mi:
                v2= np.array(np.multiply(eigenvalue[i], eigenvector[:,i]))
                sig_2 = eigenvalue[i]
        sig_1 , sig_3 = eigenvalue[ma], eigenvalue[mi] 
        return sig_1, sig_2, sig_3, v1, v2, v3
    
    def Princstresses2D( xx, yy, xy ):                # calculation of principal stresses and vectors in 2D 
        delta = sqrt((-xx-yy)**2 - 4*(-xy**2 + xx*yy))
        b = -xx-yy
        s1 = (-b + delta)/2
        s2 = (-b - delta)/2
        
        if abs(s1 - xx) < ZeroD:                                # if first principle equals to sig xx
            vctr_1 = [1. , 0., 0.]
            vctr_2 = [0. , 1., 0.]
        elif abs(s1 - yy) < ZeroD:                              # if first principle equals to sig yy
            vctr_1 = [0. , 1., 0.]
            vctr_2 = [1. , 0., 0.]
        else:
            a1 = (s1-xx)/xy
            l1 = sqrt(1**2 + a1**2)
            vctr_1 = [1/l1 , a1/l1, 0.]
                
            a2 = (s2-xx)/xy
            l2 = sqrt(1**2 + a2**2)
            vctr_2 = [1/l2 , a2/l2, 0.]
        #if s1 < 0: s1 = 0.                  # only tension 
        #if s2 > 0: s2 = 0.                  # only compression
        return s1, [s1*vctr_1[0],s1*vctr_1[1],s1*vctr_1[2]], s2, [s2*vctr_2[0],s2*vctr_2[1],s2*vctr_2[2]]
    
    class Inptdata():
        def __init__(self, coord, intpdisp, sigma, vctr):
            self.coord    = coord
            self.intpdisp = intpdisp
            self.sigma    = sigma
            self.vctr     = vctr
                  
    def ElemIntPointCoordinates( Elem, NoList):                             # integration point coordinates of element
        xyP, z, xy = [], [], []
        for i in Elem.Inzi:                                                 # ordered (without CM etc) nodal coordinates 
            node = NoList[i]
            if   Elem.SDim==1: xy += [node.XCo]
            elif Elem.SDim==2: xy += [node.XCo,node.YCo]
            elif Elem.SDim==3: xy += [node.XCo,node.YCo,node.ZCo]
        for j in range(Elem.nIntL):                                         # loop over integration points of element
            rst = SamplePoints[Elem.IntT,Elem.nInt-1,j]                     # local integration point coordinates - SamplePoints different to SamplePoints_
            XX  = Elem.FormX_( rst[0], rst[1], rst[2] )
            xyP = dot(XX,xy)                                                # global coordinates of integration point XXX
            z.append(xyP.tolist())
        return z                  

    def intpoint( NoList, ElemList):
        intP = []
        for i in ElemList:
            if i.Type in ["CPS3", "CPS3_SDA", "CPS4", "CPS4_SDA", "C3D8", "C3D8_SDA"]:
                noDisp = []
                for j in i.Inzi:
                    for k in range(len(NoList[j].BufDis)): 
                        noDisp.append(NoList[j].BufDis[k])                  # collect dofs of element from data read in
                IntpGlobal = ElemIntPointCoordinates( i, NoList)            # integration point coordinates
                for j in range(i.nIntL): 
                    Sig  = i.Data[j].tolist()    
                    Intp = SamplePoints[i.IntT,i.nInt-1,j]
                    intpdisp = dot(i.FormN(Intp[0], Intp[1], Intp[2]), noDisp)
                    if i.Type in ["CPS4", "CPS4_SDA","CPS3", "CPS3_SDA"]:   
                        if i.ElSet == "Supp": sig_1, v1, sig_2, v2 = 0., [0., 0., 0.], 0., [0., 0., 0.]
                        else : sig_1, v1, sig_2, v2 = Princstresses2D(float(Sig[4]), float(Sig[5]), float(Sig[7]))
                        intP += [Inptdata([float(IntpGlobal[0]), float(IntpGlobal[1]), 0.], [ intpdisp[0], intpdisp[1], 0. ], [sig_1, sig_2, 0],[v1,v2, [0., 0., 0.]])] 
                    elif i.Type in ["C3D8", "C3D8_SDA"]: 
                        if i.Set in ["Supp", "Supp1"]: sig_1, v1, sig_2, v2, sig_3, v3 = 0., [0., 0., 0.], 0., [0., 0., 0.], 0., [0., 0., 0.] 
                        else:                          sig_1,sig_2,sig_3,v1,v2,v3 = Princstresses3D(float(Sig[0]), float(Sig[1]), float(Sig[2]), float(Sig[3]), float(Sig[4]), float(Sig[5]))
                        intP += [Inptdata([float(IntpGlobal[j][0]), float(IntpGlobal[j][1]), float(IntpGlobal[j][2])], [ intpdisp[0], intpdisp[1], intpdisp[2] ], [sig_1, sig_2, sig_3],[v1,v2,v3])] 
        return intP  
        
    def WriteprincFile(fileOut_princ):
        f1 = open(fileOut_princ, 'w')
        f1.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
        f1.write('  <UnstructuredGrid>\n')
        f1.write('    <Piece NumberOfPoints="' + str(len(intP)) + '" NumberOfCells="' + str(len(intP)) + '">\n')
        #
        f1.write('      <Points>\n')
        f1.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.coord[0], i.coord[1], i.coord[2]))
        f1.write('      </Points>\n')
        f1.write('      <Cells>\n')
        f1.write('        <DataArray type="Int32" Name="connectivity" format="ascii"/>\n')
        for i in range(len(intP)): f1.write('       %3i\n' % (i))
        f1.write('        <DataArray type="Int32" Name="offsets" format="ascii"/>\n')
        for i in range(len(intP)): f1.write('       %3i\n' % (i))
        f1.write('        <DataArray type="Int32" Name="types" format="ascii"/>\n')
        for i in range(len(intP)): f1.write('       1\n')
        f1.write('      </Cells>\n')
        #
        f1.write('      <PointData>\n')
        f1.write('        <DataArray type="Float64" Name="IntP_displ" NumberOfComponents=" 3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.intpdisp[0], i.intpdisp[1], i.intpdisp[2])) 
        f1.write('        <DataArray type="Float64" Name="Vector_1" NumberOfComponents=" 3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[0][0], i.vctr[0][1], i.vctr[0][2]))
        f1.write('        <DataArray type="Float64" Name="Vector_2" NumberOfComponents=" 3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[1][0], i.vctr[1][1], i.vctr[1][2]))
        f1.write('        <DataArray type="Float64" Name="Vector_3" NumberOfComponents=" 3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.vctr[2][0], i.vctr[2][1], i.vctr[2][2]))  
        f1.write('        <DataArray type="Float64" Name="Princpal_Stresses" NumberOfComponents=" 3" format="ascii"/>\n')
        for i in intP:             f1.write('%12.4e%12.4e%12.4e\n' % (i.sigma[0], i.sigma[1], i.sigma[2])) 
        f1.write('      </PointData>\n')
        f1.write('    </Piece>\n')
        f1.write('  </UnstructuredGrid>\n')
        f1.write('</VTKFile>\n')
        f1.close()          
    
    intP = intpoint( NoList, ElemList)
    FFF = Name+'_Principal_ParaView'
    if pth.exists(FFF): pass
    else: os.mkdir(FFF)
    fileOut_princ = FFF + '/Principal_'+timeStr+'.vtu'
    WriteprincFile(fileOut_princ)     

def Crack(ElemList, fileName, timeStr):
        CracNodes, Tr, WL = [], [], []
        ne = 0
        Flag3D = False
        for i in ElemList:
                tr, wl = 0, 0
                if i.Type in ["C3D8_SDA"]:
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
                elif i.Type in ["CPS4_SDA"]:
                    ne += 1
                    for a in range(len(i.xyC)):
                        tr += i.wwT[0][a][1]   # position of normal traction in 2D = [1]
                        wl += i.wwL[0][a][1]   
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

def VtkFile(fileName, ElemList, VTKNodeList, marker, timeStr, elTypesExc, noTypesInc, offset, numSDA):
    # compute number of respective elements
    ne = 0
    for i in ElemList:
        if i.Type not in elTypesExc: ne += 1
    nn, nn_ = 0, 0
    for i in VTKNodeList:
        if i.Type in noTypesInc: nn += 1
    # write data
    FFF = fileName+'_'+marker+'_ParaView'
    if pth.exists(FFF): pass
    else: os.mkdir(FFF)
    fN = FFF+ "/" + marker+'_'+timeStr+'.vtu'
    fp = open(fN,'w')
    fp.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
    fp.write('  <UnstructuredGrid>\n')
    fp.write('    <Piece NumberOfPoints='+'"'+str(nn)+'"'+' NumberOfCells='+'"'+str(ne-numSDA)+'">\n')
    # points
    fp.write('      <Points>\n')
    fp.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>\n')
    for i in VTKNodeList:
        if i.Type in noTypesInc: fp.write("%12.4e%12.4e%12.4e\n"%(i.XCo,i.YCo,i.ZCo))
    fp.write('      </Points>\n')
    # cell connectivity, offset, types
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
            if   i.Type in ["C3D8","C3D8_EFG", "C3D8_SDA"]: type_ = 12
            elif i.Type in ["C3D8_IGA"]:                    type_ = 25 
            elif i.Type in ["C3D4"]:                        type_ = 10
            elif i.Type in ["T3D3E"]:                       type_ = 21
            elif i.Type in ["T3D2E", "T2D2E", "T2D2", "BM2D2E", "BM2D2"]:    type_ = 3
            elif i.Type in ["CPS4", "CPS4_SDA"]:            type_ = 9
            elif i.Type in ["CPS3", "CPS3_SDA"]:            type_ = 5
            else:  raise NameError("CaeFemInOut::DataOutVtk: unknown element type 2")
            fp.write('%8i\n'%(type_))
    fp.write('      </Cells>\n')
    # point data
    fp.write('      <PointData>\n')
    stressName = noTypesInc[0]+'_stress'
    if len(noTypesInc)>1: raise NameError("ConFemInOut::DataOutVTK: 1")
    if   noTypesInc[0] in ["T2DE", "T3DE"]: num=1
    elif noTypesInc[0] in ["T2DE", "T3DE"]: num=2  # ???????????????????????????????????????
    else:                                   num=6
    fp.write('        <DataArray type="Float64" Name="'+stressName+'" NumberOfComponents="'+str(num)+'" format="ascii"/>\n')
    BondFlag = False
    mattype = False
    for i in VTKNodeList:
        if i.Type in noTypesInc:
            if i.Type in ["T2DE", "T3DE"]:                                  # Type determined in EvNodeEl
                fp.write('%12.4e'%(i.BufSig[1]))
                BondFlag = True
            elif i.Type in ["BM2DE"]:
                fp.write('%12.4e%12.4e'%(i.BufSig[0], i.BufSig[1]))
                BondFlag = True
            elif i.Type in ["C2D"]:
                for j in i.BufSig: fp.write('%12.4e'%(j))
                mattype = True
            else:
                for j in i.BufSig: fp.write('%12.4e'%(j))
            fp.write('\n')
            
    if mattype:
        MatNamee = noTypesInc[0]+'_Mat'
        fp.write('        <DataArray type="Float64" Name="'+MatNamee+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in VTKNodeList:
            if i.Type in noTypesInc:
                if i.material == "IsoD" :    fp.write('%12.4e'%(0))
                else                    :    fp.write('%12.4e'%(1))
                fp.write('\n')
    if BondFlag:
        bondName = noTypesInc[0]+'_bond'
        fp.write('        <DataArray type="Float64" Name="'+bondName+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in VTKNodeList:
            if i.Type in noTypesInc:
#                fp.write('%12.4e'%(i.BufBond))
                fp.write('%12.4e'%(i.BufSig[5]))
                fp.write('\n')
        
        slipName= noTypesInc[0]+'_slip'
        fp.write('        <DataArray type="Float64" Name="'+slipName+'" NumberOfComponents="1" format="ascii"/>\n')
        for i in VTKNodeList:
            if i.Type in noTypesInc:
#                fp.write('%12.4e'%(i.Bufslip))
                fp.write('%12.4e'%(i.BufSig[2]))
#                    if i.Label == 641:  fslip.write('%12.4e\n'%(i.Bufslip))
                fp.write('\n')
#            fslip.flush()
        
    displName = noTypesInc[0]+'_displ'
    fp.write('        <DataArray type="Float64" Name="'+displName+'" NumberOfComponents=" 3" format="ascii"/>\n')
    for i in VTKNodeList:
        if i.Type in noTypesInc:
            if i.Type in ["BM2DE"]: fp.write("%12.4e%12.4e%12.4e"%(i.BufDis[0],i.BufDis[1],0))
            else:
                for j in i.BufDis: fp.write('%12.4e'%(j))
            fp.write('\n')
    fp.write('      </PointData>\n')
    # 
    fp.write('    </Piece>\n')
    fp.write('  </UnstructuredGrid>\n')
    fp.write('</VTKFile>\n')
    fp.close()
    return nn
#    return nn_
def PostElemVTK( FileDir, FileName, timeStr, ElList, NodeList,NoIndToCMInd, VecU, ElemResults,NodeResults, LinAlgFlag, Flag3D):
    def Distribute():
        for k,i in enumerate(el.Inzi):                                      # loop over element nodes
            no = NodeList[NoIndToCMInd[i]]
            for j in range(DatLen): no.BufSig[j] += X[j][k]                 # values of same type in "row" of X with "column" for node
            no.BufCount += 1
    # end def    
    for i in NodeList:                                                      # initialize and fill displacement buffer
#        dL = []
#        s = i.GlobDofStart
#        for j in range(len(i.DofT)): dL += [VecU[s+j]]
#        i.BufDis = dL                                                       # displacement buffer
        i.BufSig = zeros((6),dtype=float)                                   # buffer for nodal stress values used for VTK - be cautious with its length, must cover all elements
#        i.BufBond, i.Bufslip = 0., 0.
#        i.BufBondCount = 0
        i.BufCount = 0
        i.material = ""
        key = i.Label
        try:    i.BufDis = NodeResults[key][2]
        except: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 5")
    #
    for el in ElList:
        if el.Active:
            IntLen = el.nIntL                                               # number of integration points of element
            #
            if el.Type in ["C3D8"]:
                if IntLen != 8: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 3")
                DatLen = 6                                                  # corresponds to value types - actually 6 stress items
                Buf, X = zeros((DatLen,IntLen), dtype=float), []
                for i in range(IntLen):
                    key = str(el.Label)+str(el.Set)+str(i)
                    try:    ElRes = ElemResults[key][5][3:9]
                    except: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 1")
                    for j in range(DatLen): Buf[j,i] += [ElRes[j]]          # reorder rows of ElRes in columns of Buf so that rows of Buf hold values of same type
                for i in range(DatLen):                                     # loop over value types
                    X += [dot(C3D8.FormNI(el),Buf[i])]                      # multiply Buf-row with inverse of shape function to derive nodal values for each value type in row of
                Distribute()
            #
            elif el.Type in ["T3D3E"]:
                # special treatment to determine global displacement of center node
                no0 = NodeList[ NoIndToCMInd[el.Inzi[0]] ] 
                no1 = NodeList[ NoIndToCMInd[el.Inzi[1]] ] 
                no2 = NodeList[ NoIndToCMInd[el.Inzi[2]] ] 
                uEl = array([ no0.BufDis[0],no0.BufDis[1],no0.BufDis[2] , no1.BufDis[0],no1.BufDis[1],no1.BufDis[2], no2.BufDis[0]])
                NN = el.FormN( 0., 0., 0. )
                no2.BufDis = array([0.,0.,0.])                              # overwrite initial value
                for i, u in enumerate(dot(NN,uEl)): no2.BufDis[i] = u
                #
                if IntLen != 2: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 2")
                DatLen = 6                                                  # corresponds to value types - actually strain, stress for T3D3E, 6 bond items: bond item 1 for long slip, item 4 for long stress 
                ElRes, X = [], zeros((DatLen,3), dtype=float)  # 3 for three nodes of T3D3E
                for i in range(IntLen):
                    key = str(el.Label)+str(el.Set)+str(i)
                    try:    ElRes += [ElemResults[key][5][3:9]]
                    except: raise NameError("ConFemPostProcStandAalone::PostElemVTK: 4")
                for i in range(DatLen):                                     # loop over value types
                    f1 = ElRes[0][i]                                        # first integration point value 
                    f2 = ElRes[1][i]                                        # second integration point value
                    X[i,0] =  +1.366025404*f1 - .3660254039*f2              # extrapolation to left node
                    X[i,1] =   -.366025404*f1 + 1.366025404*f2              # extrapolation to right node
                    X[i,2] =   .4999999998*f1 + .5000000002*f2              # interpolation to center node
                Distribute()
    #
    for no in NodeList:
        if no.BufCount>0:
            no.BufSig = no.BufSig/no.BufCount
    # control
#    for el in ElList:
#        if el.Active and el.Type in ["T3D3E"]:
#            for i in el.Inzi:                       
#                no = NodeList[NoIndToCMInd[i]]
#                print(el.Label,no.Label,no.XCo,no.YCo,no.ZCo,no.BufCount,'__',no.BufSig[0],no.BufSig[1])
 
    off, numSDA = 0, 0
    ContinuumFirst = True
    fileName=FileDir+FileName+"_X"
    if LinAlgFlag: NodeList.sort(key=lambda t: t.Index)
    if ContinuumFirst:
        if Flag3D:
            #     VtkFile(fileName, NodeList, VTKNodeList,  marker,timeStr,   elTypesExc,                                        noTypesInc, offset, numSDA)
            off = VtkFile(fileName, ElList,NodeList, "C3D",  timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","T3D3E","T3D2E"], ["C3D"], off, numSDA)
            off = VtkFile(fileName, ElList,NodeList,    "T3D",  timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","C3D4","C3D8","C3D8_EFG","C3D8_IGA","C3D8_SDA"],  ["T3DE"], off, 0) # take care of collection of nodes belonging to a common type of elements in input file !!!
        else:
            off = VtkFile(fileName, ElList,NodeList, "CPS4", timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","T2D3E","T2D2E","T2D2", "BM2D2E", "BM2D2"], ["C2D"], off, numSDA)
            if False: #Beam:
                off = VtkFile(fileName, ElList,NodeList,    "BM2D",  timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","CPS4","CPS4_EFG","CPS4_IGA","CPS4_SDA","CPS3","CPS3_SDA"],  ["BM2DE"], off, 0) 
            else:
                off = VtkFile(fileName, ElList,NodeList,    "T2D",  timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","CPS4","CPS4_EFG","CPS4_IGA","CPS4_SDA","CPS3","CPS3_SDA"],  ["T2DE"], off, 0) # take care of collection of nodes belonging to a common type of elements in input file !!!
        Crack(ElList, fileName, timeStr)
        Principal( NodeList, ElList, fileName, timeStr)
    else:
        off = VtkFile(fileName, ElList,NodeList,    "T3D", timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","C3D8","C3D8_EFG","C3D8_IGA","C3D8_SDA"],  ["T3DE"], off, 0) # take care of collection of nodes belonging to a common type of elements in input file !!!
        off = VtkFile(fileName, ElList,NodeList, "C3D", timeStr, ["B2D2E","B2D3E","B3D2E","B3D3E","T3D3E","T3D2E"], ["C3D"], off, numSDA)
    if LinAlgFlag: NodeList.sort(key=lambda t: t.CMIndex)

def PostElem3D( ElList, NodeList,NoIndToCMInd, VecU, SN, mod_Z_limit):
    XList, YList, ZList, DisElemList, sig1POSMElemList, sig2POSMElemList, sig1NEGMElemList, sig2NEGMElemList, nxLi,nyLi,nxyLi,mxLi,myLi,mxyLi,qxLi,qyLi, aaLi = [], [], [], [], [], [], [], [],[],[],[],[],[],[],[],[],[]  # ts
    for i in NodeList:
        XList += [i.XCo]
        YList += [i.YCo]
        ZList += [i.ZCo]
    mX, mY, mZ = [np.min(XList),np.max(XList),np.mean(XList)], [np.min(YList),np.max(YList),np.mean(YList)], [np.min(ZList),np.max(ZList),np.mean(ZList)]
    for elem in ElList:
        DisList = []
        for o in elem.Inzi:
            ni = NoIndToCMInd[o]
            node = NodeList[ni]
            DofInd = node.GlobDofStart
            if   elem.Type in ['SH4','SH3']: DisList += [sqrt(VecU[DofInd]**2+VecU[DofInd+1]**2+VecU[DofInd+2]**2)]
            elif elem.Type=='SB3':           DisList += [VecU[DofInd]]       
        DisElemList += [np.mean(DisList)]
    plotDeformFigure('deform', DisElemList,VecU,ElList,NodeList,NoIndToCMInd,' deformed mesh', mX, mY, mZ)
    if mod_Z_limit==None: return 0
    ############################################################################################################################### # ts
    for elem in ElList: 
        if elem.Type in ['SH4','SH3']:
            offset, nRe = 0, 0
            if elem.ShellRCFlag: nRe = elem.Geom.shape[0]-2                     # number of reinforcement layers
            else:                nRe = 0
            Lis = elem.Lists1()                                                 # indices for first integration points in base area      
            LisLen = len(Lis)                                                   # number of integration points in base area
            nx, ny, nxy, qy, qx, mx, my, mxy, sig1POSSumme, sig2POSSumme, sig1NEGSumme, sig2NEGSumme, aa_ = 0., 0., 0., 0., 0., 0., 0., 0.,  0., 0., 0., 0., 0.
            for j in Lis:
                nx_, ny_, nxy_, qx_, qy_, mx_, my_, mxy_, aa = elem.StressIntegration( j, offset, nRe)
                nx += nx_  
                ny += ny_ 
                nxy+= nxy_  
                qx += qx_
                qy += qy_
                mx += mx_
                my += my_
                mxy+= mxy_  
                aa_ = aa_ + aa
            # end of base integration point loop - following for means of element 
    #        sig1NEGMElemList.append(sig1NEGSumme/LisLen)
    #        sig2NEGMElemList.append(sig2NEGSumme/LisLen)                                                                                               # ts sig1M = sig1Summe/(len(Lis))                                                                                                    # ts
    #        sig1POSMElemList.append(sig1POSSumme/LisLen)                                                                                               # ts 
    #        sig2POSMElemList.append(sig2POSSumme/LisLen)
            nxLi.append( nx/LisLen)
            nyLi.append( ny/LisLen)
            nxyLi.append(nxy/LisLen)
            mxLi.append( mx/LisLen)
            myLi.append( my/LisLen)
            mxyLi.append(mxy/LisLen)
            qxLi.append( qx/LisLen)
            qyLi.append( qy/LisLen)
            aaLi.append( aa_/LisLen)
    # end of element loop
    ############################################################################################################################### # ts                                 # ts
#        plotMainStressFigures2(sig1NEGMElemList,'S1 SNEG',XList,YList,ZList,ElList,NodeList)
#        plotMainStressFigures2(sig2NEGMElemList,'S2 SNEG',XList,YList,ZList,ElList,NodeList)
#        plotMainStressFigures2(sig1POSMElemList,'S1 SPOS',XList,YList,ZList,ElList,NodeList)
#        plotMainStressFigures2(sig2POSMElemList,'S2 SPOS',XList,YList,ZList,ElList,NodeList)
    plotMainStressFigures2(nxLi, 'n_x', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(nyLi, 'n_y', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(nxyLi,'n_xy',XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(mxLi, 'm_x', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(myLi, 'm_y', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(mxyLi,'m_xy',XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(qxLi, 'q_x', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(qyLi, 'q_y', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    plotMainStressFigures2(aaLi, 'thickness', XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit)
    #
    return 0

def plotMainStressFigures2(sig_MElemList,Sx_S,XList,YList,ZList,ElList,NodeList,NoIndToCMInd, mod_Z_limit):
    ColM = plt.get_cmap('brg')
    fig = pl.figure() #figsize=(15.,10.))
    fig.text( 0.05, 0.93, Sx_S, size='x-large')
    maxSig = max(sig_MElemList)
    minSig = min(sig_MElemList)
    if abs(maxSig-minSig)<1.e-6: return 0
    norm_ = colors.Normalize(vmin=minSig, vmax=maxSig)
    ax1, orie = fig.add_axes([0.90, 0.05, 0.02, 0.8]), 'vertical'
#    ax1, orie = fig.add_axes([0.05, 0.05, 0.9, 0.02]), 'horizontal' # left, bottom, width, height] where all quantities are in fractions of figure width and height
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=ColM,  norm=norm_, orientation=orie) # norm=norm_,
    ax = a3.Axes3D(fig)
#    GesamtDiff = abs(max(sig_MElemList) - min(sig_MElemList))
#    prozentualerNegAnteil = abs(min(sig_MElemList)) / GesamtDiff
#    prozentualerPosAnteil = abs(max(sig_MElemList)) / GesamtDiff
    for i in range(len(ElList)):
        elem, CoordList = ElList[i], []
        for j in range(len(elem.Inzi)):
            node = NodeList[elem.Inzi[j]]
            ni = NoIndToCMInd[elem.Inzi[j]]
            node = NodeList[ni]
            CoordList += [[node.XCo, node.YCo, node.ZCo]]
        tri = a3.art3d.Poly3DCollection([CoordList])
#        if sig_MElemList[i] <= 0.:
#            sig_MElemList[i] = minSig - sig_MElemList[i] # i-tes Element wird zu einem Differenzwert
#            sig_MElemList[i] = sig_MElemList[i] / minSig # i-tes Element wird zu einem Prozentwert, der den Spannungswert auf einer Skala von 0 bis 1 fuer den negativen Bereich abbildet
#            sig_MElemList[i] = sig_MElemList[i] * prozentualerNegAnteil # i-tes Element bildet den Prozentwert, der sich im vorigen Schritt noch nur auf den negativen Bereich bezogen hat, auf den gesamten Bereich ab
#        else:
#            sig_MElemList[i] = sig_MElemList[i] #sig_MElemList[i] = max(sig_MElemList) - sig_MElemList[i] # i-tes Element wird zu einem Differenzwert
#            sig_MElemList[i] = sig_MElemList[i] / maxSig # i-tes Element wird zu einem Prozentwert, welcher den Spannungswert auf einer Skala von 0 bis 1 fuer den positiven Bereich abbildet
#            sig_MElemList[i] = sig_MElemList[i] * prozentualerPosAnteil + prozentualerNegAnteil# i-tes Element bildet den Prozentwert, der sich im vorigen Schritt noch nur auf den positiven Bereich bezogen hat, auf den gesamten Bereich ab
#        tri.set_color(ColM(sig_MElemList[i]))
        tri.set_color(ColM((sig_MElemList[i] - minSig)/(maxSig-minSig)))
        tri.set_edgecolor('k')
        ax.add_collection3d(tri)
    ax.set_xlim3d(min(XList), max(XList))
    ax.set_ylim3d(min(YList), max(YList))
    ax.set_zlim3d(mod_Z_limit*min(ZList), mod_Z_limit*max(ZList))
    ax.set_xlabel('x',fontsize='x-large')
    ax.set_ylabel('y',fontsize='x-large')
    ax.set_zlabel('z',fontsize='x-large')
    ax.tick_params(axis='x', labelsize='x-large')
    ax.tick_params(axis='y', labelsize='x-large')
    ax.tick_params(axis='z', labelsize='x-large')
    ax.set_aspect('equal') #,'box', 'equal')
    ax.grid()
    return 0

def plotDeformFigure(flag, DisElemList,VecU,ElList,NodeList,NoIndToCMInd,text, mX, mY, mZ):
    xMin, xMax, XM = mX[0], mX[1], mX[2]
    yMin, yMax, YM = mY[0], mY[1], mY[2]
    zMin, zMax, ZM = mZ[0], mZ[1], mZ[2]
    dMax, dMin = max(DisElemList), min(DisElemList)
    dMax = max(abs(dMax),abs(dMin))
    scale_ = 10.                                                            # 20.  # displayed deformation should be roughly 1/scale_ of parts dimensions
    scale = max(XM,YM,ZM)/dMax/scale_
    ColM = plt.get_cmap('summer')
    fig = pl.figure()#(figsize=(15.,10.))
    fig.text( 0.05, 0.93, 'Deformed shape - scale: {:7.3f} '.format(scale)+text, size='x-large') # , ha='right')
    norm_ = colors.Normalize(vmin=0., vmax=dMax)
    ax1, orie = fig.add_axes([0.90, 0.05, 0.02, 0.8]), 'vertical'           # left, bottom, width, height] where all quantities are in fractions of figure width and height
#    ax1, orie = fig.add_axes([0.05, 0.05, 0.9, 0.02]), 'horizontal' # left, bottom, width, height] where all quantities are in fractions of figure width and height
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=ColM, norm=norm_, orientation=orie)
    ax = a3.Axes3D(fig)
    DisElemList = list(map(abs,DisElemList))/dMax
    for i in range(len(ElList)):
        elem, CoordList = ElList[i], []
        for j in elem.Inzi:
            ni = NoIndToCMInd[j]
            node = NodeList[ni]
            DofInd = node.GlobDofStart
            if   elem.Type in ['SH4','SH3']: CoordList += [[node.XCo-xMin+scale*VecU[DofInd], node.YCo-yMin+scale*VecU[DofInd+1], node.ZCo-zMin+scale*VecU[DofInd+2]]]
            elif elem.Type=='SB3':           CoordList += [[node.XCo-xMin, node.YCo-yMin, node.ZCo-zMin+scale*VecU[DofInd]]]
            else: 
                print(elem.Type)
                raise NameError ("XXX")
        tri = a3.art3d.Poly3DCollection([CoordList])
        tri.set_color(ColM(DisElemList[i]))
        tri.set_edgecolor('k')
        ax.add_collection3d(tri)
    limits = max(zMax - zMin, xMax -xMin, yMax - yMin)
    axesratio = 1.0                                                         # 0.3     # minimum of smallest dimension to largest dimension in 3d coordinate system displayed
    ax.set_zlim3d(0.,max( axesratio*limits, zMax - zMin))
    ax.set_xlim3d(0.,max( axesratio*limits, xMax - xMin))
    ax.set_ylim3d(0.,max( axesratio*limits, yMax - yMin))
    ax.set_xlabel('x',fontsize='x-large')
    ax.set_ylabel('y',fontsize='x-large')
    ax.set_zlabel('z',fontsize='x-large')
    ax.tick_params(axis='x', labelsize='x-large')
    ax.tick_params(axis='y', labelsize='x-large')
    ax.tick_params(axis='z', labelsize='x-large')
    ax.set_aspect('equal')
    ax.grid()
    return 0

def mainStressConverter(sig_x,sig_y,sig_xy): # Anwendung des Mohrschen Spannungskreises zur Ausgabe der Hauptspannungen
    point1 = [sig_x,sig_xy]
    #point2 = [sig_y,-sig_xy]
    xCenterPoint=min(sig_x,sig_y)+abs(sig_x - sig_y)/2
    centerPoint = [xCenterPoint,0] # entspricht dem Mittelpunkt des Morschenspannungskreises
    radiusMohr = sqrt((point1[0]-centerPoint[0])**2 + (point1[1]-centerPoint[1])**2)
    sigma1 = xCenterPoint + radiusMohr
    sigma2 = xCenterPoint - radiusMohr
    #kreisSehne = # fuer die Berechnung der Hauptspannungsorientierung
    #mittelPunktsWinkel = # fuer die Berechnung der Hauptspannungsorientierung
    return sigma1,sigma2

# store values of ndof per node 
def NodalDisplList( NodeList, NodeResults):
    DisplList = [] 
    for i in NodeList:
        key = i.Label
        nDof = NodeResults[key][1]
        Data = NodeResults[key][2]
        dL = []
        for j in range(nDof): dL += [Data[j]]
        DisplList += [dL]
        i.dL = dL
    return DisplList
# store values of ndof per element 
def ElementDisplList( i, DisplList, NoIndToCMInd, ScaleU):
    uEl, uxy = [], []                             
    for j in i.Inzi:                                                # loop over nodes affected by element -- equivalent to lines before
        j_ = NoIndToCMInd[j]
        for k in range(len(DisplList[j_])):
            uEl += [DisplList[j_][k]]                               # collects nodal displacements per element
            uxy += [ScaleU*DisplList[j_][k]]                        # collects scaled nodal displacements per element
    return uEl, uxy

def PrinC( xx, yy, xy ):                                 # calculation of principal values
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
    return( [s1, n11, n12, s2, n21, n22] )

def ReadLine(ff):
    z1 = ff.readline()
    z2 = z1.strip()
    return z2.split(',')

def ReadResultsFromFile( ff, ffN):
    ElemResults, NodeResults, EndFlag, EndFlagN = {}, {}, False, False
    # read from elemout_
    z3 = ReadLine(ff)
    while z3!="" and z3!=[""] and z3[0]!='Time':
        key = z3[0].strip()+z3[1].strip()+z3[4].strip()                     # key is label + elset + integration point
        Data = []
        for i in z3[5:-1]: Data += [float(i)]
        ElemResults[key] = [z3[0].strip(), z3[1].strip(), z3[2].strip(), z3[3].strip(), z3[4].strip(), Data]
        z3 = ReadLine(ff)
    if z3=="" or z3==[""]: EndFlag, Time = True,  None
    else:                  EndFlag, Time = False, float(z3[1])
    # read from nodeout_
    z3N = ReadLine(ffN)
    while z3N!="" and z3N!=[""] and z3N[0]!='Time':
        key  = int(z3N[0])                                                  # int node label, should be unique
        nDof = int(z3N[4])                                                  # number of degrees of freedom
        Data = []
        for i in z3N[5:5+nDof]: Data += [float(i)]
        NodeResults[key] = [ key, nDof, Data]
        z3N = ReadLine(ffN)
    if z3N=="" or z3N==[""]: EndFlagN, TimeN = True, None
    else:                    EndFlagN, TimeN = False, float(z3N[1])
        
    return EndFlag, EndFlagN, Time, TimeN, ElemResults, NodeResults

def PickleLoad( FilDir, FilName, StepCounter):
    import pickle
#    from scipy.spatial import KDTree
    import scipy.spatial as spatial
#    from os import path as pth
    if pth.isfile(FilDir+FilName+".pkl"):
        fd = open(FilDir+FilName+'.pkl', 'rb')
        NodeList=pickle.load(fd);ElList=pickle.load(fd);MatList=pickle.load(fd);StepList=pickle.load(fd);N=pickle.load(fd);WrNodes=pickle.load(fd);LineS=pickle.load(fd);FlElasticLT=pickle.load(fd);\
            VecU=pickle.load(fd);VecC=pickle.load(fd);VecI=pickle.load(fd);VecP=pickle.load(fd);VecP0=pickle.load(fd);VecP0old=pickle.load(fd);VecBold=pickle.load(fd);VecT=pickle.load(fd);VecS=pickle.load(fd);\
            VeaU=pickle.load(fd);VevU=pickle.load(fd);VeaC=pickle.load(fd);VevC=pickle.load(fd);VecY=pickle.load(fd);BCIn=pickle.load(fd);BCIi=pickle.load(fd);Time=pickle.load(fd);TimeOld=pickle.load(fd);\
            TimeEl=pickle.load(fd);TimeNo=pickle.load(fd);TimeS=pickle.load(fd);Step=pickle.load(fd);Mask=pickle.load(fd);Skyline=pickle.load(fd);SDiag=pickle.load(fd);SLen=pickle.load(fd);SymSys=pickle.load(fd);\
            NoLabToNoInd=pickle.load(fd);NoIndToCMInd=pickle.load(fd);ContinuumNodes=pickle.load(fd);CoNoToNoLi=pickle.load(fd);SecDic=pickle.load(fd);LinAlgFlag=pickle.load(fd);
        fd.close()
#        s = StepList[StepCounter]
        if len(ContinuumNodes)>0: CoorTree = spatial.cKDTree( ContinuumNodes ) # for search purposes, e.g. for EFG or aggregates or embedded truss elements
        else:                     CoorTree = None
    else: 
        raise NameError ("PickeLoad: cannot read data")
    return NodeList, ElList, MatList, CoorTree, CoNoToNoLi, VecU, NoIndToCMInd, NoLabToNoInd, SecDic, LinAlgFlag

def Main( FilDir,FilName,SubDir,StepCounter, PostFlagStress,PostFlagMaPlLi,Plot3D, ScU,ScaleS, TimeOut,VTK,VTK3D):
    # load system data
    NodeList,ElemList,MatList,CoorTree,CoNoToNoLi,VecU,NoIndToCMInd,NoLabToNoInd, SecDic, LinAlgFlag = PickleLoad( FilDir, FilName, StepCounter)
    # open read from elemout, nodeout
#    if SubDir != "": fileName = FilDir+FilName+"."+SubDir+"/"+FilName+".elemout"+"."+TimeEl+".txt"
#    else:            fileName = FilDir+"/"+FilName+".elemout_"+".txt"
    fileName = FilDir+"/"+FilName+".elemout_"+".txt"
    ff  = open(fileName,'r')
#    if SubDir != "": fileName = FilDir+FilName+"."+SubDir+"/"+FilName+".nodeout_"+"."+TimeEl+".txt"
#    else:            fileName = FilDir+"/"+FilName+".nodeout_"+".txt"
    fileName = FilDir+"/"+FilName+".nodeout_"+".txt"
    ffN = open(fileName,'r')
    # start reading
    z3  = ReadLine(ff)
    z3N = ReadLine(ffN)
    EndFlag, EndFlagN, Time, TimeN = False, False, float(z3[1]), float(z3N[1])
    if Time != TimeN: raise NameError("ConFemPostProcStandAlone: async time for elements and nodes",Time, TimeN)
    # loop over all times in result data sets
    while not EndFlag and not EndFlagN:
        EndFlag, EndFlagN, Time_, TimeN_, ElemResults, NodeResults = ReadResultsFromFile( ff, ffN) # reads until next time marker, will presumably work only for elemenout_, nodeout_
        if Time in TimeOut or len(TimeOut)==0:
            if not Plot3D:
                _ = Post2D( ElemList,NodeList,MatList, ScU,ScaleS, SecDic, PostFlagStress,PostFlagMaPlLi, NoIndToCMInd, ElemResults,NodeResults, Time, FilDir+"/")
            else:
                PostElem3D( ElemList, NodeList,NoIndToCMInd, VecU, 1.0, 1.0)
            if VTK:
                PostElemVTK(FilDir,FilName, str(int(1000*round(Time,5))), ElemList, NodeList,NoIndToCMInd, VecU, ElemResults,NodeResults, LinAlgFlag, VTK3D)
        Time, TimeN = Time_, TimeN_
        if Time != TimeN: raise NameError("ConFemPostProcStandAlone: async time for elements and nodes",Time, TimeN)
    # 
    ff.close()
    ffN.close()    
    WrNodes = None
    Name = FilDir+FilName
    if pth.isfile(Name+".opt.txt"):                                 # read options file if there is any
        f4=open( Name+".opt.txt", 'r')
        WrNodes, _, _, _ = ReadOptionsFile(f4, NodeList,NoLabToNoInd,NoIndToCMInd)
        f4.close()
        if WrNodes!=None and len(WrNodes)>0:
            f5=open( Name+".timeout.txt", 'r')
            print('PlotNodes')
            PlotNodes( f5, WrNodes )
            f5.close()
            return True
        return False
    return False
            
if __name__ == "__main__":
    StepCounter = 0
    PostFlagMaPlLi, PostFlagStress = True, True                             # flag for using motplotlib,  flag for plotting principal stresses 
    TimeOut = []
    Plot3D = False                                                          # for shells and slabs
    SubDir = ""
    VTK,VTK3D = False, False
    # ScP scaling factor principal stresses, ScU scaling factor displacements
    #                                                                  [output times], scaling factor displacements, [[SolidSection 1], [SolidSection 2], .... ]
    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataTmp/","phi-45coarse",          "", [],            5., [[100.],   [2.0e-3, 1.0, 0.1]] #
    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataTmp/","Deep_beam"   ,          "", [0.2,0.4,0.6], 5., [[0.0002], [0.0002] , [0.0005,0.0005,0.05,0.01]] #
    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataTmp/","Deep_beam_AacNL"   ,    "", [0.6], 5., [[0.0002], [0.0002], [0.002] , [0.0005,0.0005,0.05,0.01] ] #
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataBenchmarks/L_ShapedPanel/", "LSP-2079-7_6_DAMMK", "", [0.8], 5., [[0.001]]
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataTmp/","TwoElementCPS3_SDATest", "",[],          500., [[1.0]]
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataSpecimen/One2D/", "WillamsTest2D", "", [], 5.,  [[0.2]]
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataTmp/","Plane-el"   ,"", [0.5,1.0], 5.,[ [0.01] , [0.0005,0.05,1.0], [0.001]] #
#    FilDir,FilName, ScP,ScU, Plot3D ="C:/Users/uhc/Workspace/ConFemPy3/PreProc/dist/", "Hypar", 1., 10., True
#    FilDir,FilName, ScP,ScU, Plot3D ="C:/Users/uhc/Workspace/ConFemPy3/PreProc/", "test_007", 1., 10., True
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/","LSP_1","",[0.4],         200., [[1.0]] # 200.
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS,VTK, VTK3D = "../_DataDiv/","Cube8"      ,"",[],         200., [[1.0],[1.0]], True, True # 200.
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS = "../_DataRegularization/","TwoElementCPS3_SDATest_6n","",[], 200.,[[5.]]
#    FilDir,FilName,SubDir, TimeOut, ScU,ScaleS,VTK, VTK3D = "C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/","LSP_1_C3D8_","",[],200., [[1.0]], True, True # 200.
    Flag = Main( FilDir,FilName,SubDir,StepCounter, PostFlagStress,PostFlagMaPlLi,Plot3D, ScU,ScaleS, TimeOut,VTK,VTK3D)
    print('finish processing')
    if PostFlagMaPlLi or Flag: plt.show() 
    print('finish all')
