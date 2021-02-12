# ConFemAll -- 2014-01-13
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
'''
Created on 23.05.2011

@author: uhc
'''
import ConFem
import ConSimFem
import ConPlaD
import ConSimplex
import ConPostPlot

X = str(input('ConFem 0, SimFem 1, ConPlaD 2, ConSim 3, PostFem 4: '))
Name=str(input('Filename without extension: '))
PlotS=str(input('Post-Plots (y/n), default no: '))
if len(PlotS)>0 and PlotS.upper()[0]=='Y': PlotF = True
else:                                      PlotF = False 

if X=='0':
    LogName="../DataExamples/tmp"                             # to log temporary data
    ConFem_ = ConFem.ConFem()
    ConFem_.Run(Name, LogName, PlotF, True, False, "elemout", [None, None], [], [])
elif X=='1':
    SimFem_ = ConSimFem.ConSimFem()
    SimFem_.Run(Name, True, True, "elemout")
elif X=='2':
    ConPlaD_ = ConPlaD.ConPlaD()
    ElSet = str(input('Element Set ("EL1"): '))
    if ElSet.strip()=='': ElSet="EL1"
    Tim = float(input('Time (no default try 1.0): '))
    Tim = "%8.4f"%(Tim)
    Mat = str(input('Material Name (MAT1): '))
    if Mat.strip()=='': Mat="MAT1"
    Type = str(input('plate/slab: '))
    ConPlaD_.Run(Name, ElSet,Tim,Mat, True, Type)
elif X=='3':
    ConSim_ = ConSimplex.ConSimplex()
    ConSim_.Run(Name, True, True)
elif X=='4':
    PostFem_ = ConPostPlot.ConPostPlot()
    PostFem_.Run(Name)
else:
    print('Invalid index!')
