# ConFemBasics -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
def DefData():
    ElPlotTimes = []
    Plot = True
    StressStrainOut = []                                                    # stress strain output for selected elements labels with ip -- in a list [el,ip]
    VTK = False
    
    # DataExamples E03
#    Name="../DataExamples/E03/E3-01"                                        # 25-02-18
#    Name="../DataExamples/E03/E3-02C3"                                      # C1, C2, C3 25-02-18
    Name="../DataExamples/E03/E3-04"                                        # 25-02-18
    # E04 beams
#    Name="../DataExamples/E04/E3-02"                                        # 25-02-18
#    Name="../DataExamples/E04/E4-09"                                        # -03 , -04, -06, -08, -09 25-02-18
#    Name="../DataExamples/E04/E3-02_CircMises"                              # 25-02-18
    # E05 strut-and-tie
#    Name="../DataExamples/E05/E5-01plate"                                   # 25-02-18
#    Name="../DataExamples/E05/E5-02"                                        # 25-02-19
#    Name="../DataExamples/E05/E5-03"                                        # 25-02-19
    # E06 Multiaxial concrete
#    Name="../DataExamples/E06/E6-02"                                        # 25-02-19
    # E07 crack modeling and regularization
#    Name, Plot ="../DataExamples/E07/E7-01b", True                          # a, b, 25-02-20 c, premature end with a
#    Name, Plot ="../DataExamples/E07/E7-05", False                          # 25-02-20 1-D phase field
#    Name, Plot, StressStrainOut ="../DataExamples/E07-06/TensionEccTri2D", True, [] #[[1,0]]
#    Name, Plot ="../DataExamples/E07-06/E7-06c", False                      # 25-02-20 SDA
    # E08
#    Name, Plot ="../DataExamples/E08/E8-01", True                           # 25-02-20 deep beam elastic
#    Name, Plot = "../DataExamples/E08-02/E8-02", False                      # ok deep beam nonlinear
#    Name, Plot = "../DataExamples/E08/E8-03", True                          # 25-02-20 single fiber
#    Name, Plot = "../DataExamples/E08/E8-04", True                          # 25-02-20 plate with reinforcement
    # E09                             
#    Name, Plot ="../DataExamples/E09/E9-01", True                           # 25-03-03 elastic slab with opening and free edges
#    Name, Plot = "../DataExamples/E09/E9-04", True                          # 25-03-03 Elastoplastic slab (NLSLAB) with opening and free edges
#    Name, Plot = "../DataExamples/E09/E9-05", True                          # 25-03-03 simple Elastoplastic-extended (NLSLAB) slab under concentrated loading
    # E10
#    Name, Plot ="../DataExamples/E10/E10-01c", True                         # 25-03-04 flat shell convergence / shear locking
#    Name, Plot ="../DataExamples/E10/E10-02a", True                         # 25-03-12 one element RC shell
#    Name, Plot ="../DataExamples/E10-02/E10-02", False                      # simple slab nonlinear reinforcement
    # _DataAxisym
#    Name, Plot ="../_DataAxisym/BAX2E", True                                # 25-03-04
#    Name, Plot ="../_DataAxisym/CAX3", True                                 # 25-03-04
#    Name, Plot = "../_DataAxisym/CylAX", True                               # 25-03-04
#    Name, Plot = "../_DataAxisym/TAX23", True                               # 25-03-04
    # _DataBeams
#    Name, Plot = "../_DataBeams/E3-02_B23", True                            # 25-03-12
#    Name, Plot = "../_DataBeams/E3-02_B21E", True                           # 25-03-12
#    Name, Plot = "../_DataBeams/E3-02_B23E_PolyL", True                     # 23-03-12
    # _DataBenchmarks
#    Name, Plot ="../_DataBenchmarks/Arrea/arrea2D", True                    # 25-03-04
#    Name, Plot ="../_DataBenchmarks/Jenq/Jenq-L3-12122007-Final", True      # 25-03-04
#    Name, Plot ="../_DataBenchmarks/Nooru/Nooru-1550-5", True               # 25-03-04
#    Name, Plot, StressStrainOut ="../_DataBenchmarks/L_ShapedPanel/LSP_1", True, [[421,0]]     # 25-03-12k
#    Name, Plot, StressStrainOut ="../_DataBenchmarks/L_ShapedPanel/LSP-2079-7_6_DAMMK", True, []   # 25-03-12k
#    Name, Plot, StressStrainOut =" C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1_C3D8_", False, [[421,0]] # 25-03-12
#    Name, Plot, StressStrainOut = "../_DataBenchmarks/Aachentests/Deep_beam", True, []   # 25-03-13
#    Name, Plot, StressStrainOut = "../_DataBenchmarks/Aachentests/Deep_beam_AacNL", True, []  # 25-03-14
    # _DataBond
#    Name, Plot, StressStrainOut = "../_DataBond/bond_T2", True, []          # 25-03-15 TestSuite0
#    Name, Plot, StressStrainOut = "../_DataBond/Pullout", True, []          # 25-03-15 Pullout, Pullout2
#    Name, Plot, StressStrainOut = "../_DataBond/PulloutAxiSym", True, []    # 25-03-15 Isodamage
#    Name, Plot, StressStrainOut = "../_DataBond/PulloutLatPressure", True, []  # 25-03-15
    # _DataBuckling
#    Name, Plot ="../_DataBuckling/c_264_210000", True                       # 25-03-15 sphere cap
#    Name, Plot = "../_DataBuckling/E3-02V", True          # 25-03-15 2D Euler -- works only for 4 eigenmodes not for default 5
#    Name, Plot = "../_DataBuckling/fullCyl", True       # 25-03-15 full cylinder under compression - lower end horizontally free --works only for 1 eigenmodes not for default 5
#    Name, Plot = "../_DataBuckling/ShellBuckl", True                        # 25-03-15 simple flat shell with longitudinal compression
#    Name, Plot = "../_DataBuckling/a_640_210000", True                      # 25-03-15 full cylinder under compression -- lower end clamped
    # _DataC3D
#    Name, Plot, StressStrainOut, VTK ="../_DataC3D/Cube8", True, [], True   # 25-03-16
#    Name, Plot, StressStrainOut, VTK ="../_DataC3D/Deep_beam_SDA", True, [], True   # 25-03-16
    # _DataRegularization
#    Name, StressStrainOut = "../_DataRegularization/TwoElementCPS3_SDATest", [[1,0],[2,0]]   # 25-03-16
#    Name, Plot, StressStrainOut ="../_DataRegularization/TwoElementCPS3_SDATest_6n", True, [[1,0],[2,0]]   # 25-03-16
    # _DataShellsSlabs
#    Name, Plot ="../_DataShellsSlabs/HyparRand", True                       #  25-03-16 "Shell"  "bridge_el05m" "c_1461(0.08)_2.1e5_0.3_segment load", "HyparRand"
#    Name, Plot = "../_DataShellsSlabs/IuriisExamples/slab_sh4", True        # 4x4",
    # _DataSpecimen
#    Name, Plot, StressStrainOut = "../_DataSpecimen/UniaxTensionBar", True, []  # 25-03-18
#    Name, Plot, StressStrainOut ="../_DataSpecimen/One3D/Uniax3D", True, [[1,0]]# 25-03-18
#    Name, Plot ="../_DataSpecimen/One2D/WillamsTest2D", True                # 25-03-18
#    Name, Plot ="../_DataSpecimen/One3D/WillamsTest", False                 # 25-03-18 TestSuite0 microdamage
    # _DataTrusses all TestSuiteTrusses
#    Name, Plot = "../_DataTrusses/staebe_1_1d", True                        # 25-03-18 elastoplastic mises with un- and reloading
#    Name, Plot = "../_DataTrusses/staebe_4", True                           # 25-03-18 2D elastic large displacement
#    Name, Plot = "../_DataTrusses/staebe3D2", True                          # 25-03-18 3D elastic large displacement
    # _DataTmp
#    Name, Plot, StressStrainOut = "../_DataTmp/shellconfem", True, []
#    Name, Plot, StressStrainOut = "../_DataTmp/E7-06a", True, []
    #
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2023/Note_ShearPlateRandom/ConFem/ShearPanel/ShearPanelR_med", True, [] #[[1,0]]
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2023/Note_ShearPlateRandom/ConFem/ShearPanel/Explicit/ShearPanelR", False, []
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/OneElement/2D/Large/ConFem/Job-1", True, [[1,0]]
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/OneElement/3D/Large/ConFem/Set-2", True, [[1,0],[1,1],[1,2],[1,3],[1,4],[1,5],[1,6],[1,7]]
#    Name, Plot                  = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/NotchedBeam/2D/Fine/ConFem/notchedBeam_2d_fine", True
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/MixedMode/Specimen/Explicit/Willam/ConFem/Job-1", True, [[1,0]]
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Desktop/Note_FlatSlab_Comp/ExpDataSets/372-L1/372-L1", False, []
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Desktop/Note_FlatSlab_Comp/ExpDataLandler/1-M0-25-1.23/1-M0-25-1.23", True, []

    Restart = False
#    Restart = True
    Name = Name.strip()
    return Name, Plot, Restart, "elemout", ElPlotTimes, StressStrainOut, VTK
