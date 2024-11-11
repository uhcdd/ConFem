# ConFemAll -- 2014-01-13
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
import ConFem
import ConSimFem
import ConPlaD
import ConSimplex
import ConFemPost
import ConStressStrain

# this script is for building a self-sufficient exe-version of the ConFem package by pyinstaller to run under windows
# tested for win 10, 11 and a limited selection of x64-processors only

Im      = True                                                              # hard coded ConFem only for batch computations - False for prompted selection out of package and plot
LogName = "./LogFiles"                                                      # to log temporary data
Name    = str(input('Filename without extension: '))

if Im:
    ConFem_ = ConFem.ConFem()
    ConFem_.Run(Name, LogName, False, True, False, "elemout", [None, None], [], [], False) # expects data set name in command line
else:
    X = str(input('ConFem 0, ConSimFem 1, ConPlaD 2, ConSimplex 3, ConFemPost 4, ConStressStrain 5: '))
    if X in ['0','1']:
        PlotS=str(input('ConFem / ConSimFem Post-Plots (y/n), default no: '))
        if len(PlotS)>0 and PlotS.upper()[0]=='Y': PlotF = True
        else:                                      PlotF = False
    if X=='0':
        ConFem_ = ConFem.ConFem()
        ConFem_.Run(Name, LogName, PlotF, True, False, "elemout", [None, None], [], [], False) # does not use DefData
    elif X=='1':
        SimFem_ = ConSimFem.ConSimFem()
        SimFem_.Run(Name, True, True, "elemout")
    elif X=='2':
        ConPlaD_ = ConPlaD.ConPlaD()
        ElSet = str(input('Element Set: '))
        Type = str(input('plate/slab: '))
        ConPlaD_.Run(Name, ElSet, True, Type)
    elif X=='3':
        ConSim_ = ConSimplex.ConSimplex()
        ConSim_.Run(Name, True, True)
    elif X=='4':
        ConFemPost_ = ConFemPost.ConFemPost()
        Flag = ConFemPost_.Run( "", Name,0, True,False,False, 1)
        if Flag:
            import matplotlib.pyplot as plt
            plt.show()
    elif X=='5':
        ElSet = str(input('ConMat requires extra data\nElement set name out of *.in: '))
        MatName = str(input('Material name corresponding to element set: '))
        NF = str(input('Normal forces separated by comma: '))
        NF_= NF.strip()
        NF = NF_.split(',')
        NormalForceList = []
        for i in NF: NormalForceList += [float(i)]
        ConMat = ConStressStrain.MatTest(Name + ".MatTester.txt")
        ConMat.Run(Name, MatName, ElSet, NormalForceList)
    else:
        print('Invalid index!')
