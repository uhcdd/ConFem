# ConFemBasics -- 2022-09-27
# Copyright (C) [2014] [Ulrich Haeussler-Combe]
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (GNU GPLv3) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses
#
def DefData():
    nEig, niterEig, ElPlotTimes = None, None, []
    Plot = True #True 
    StressStrainOut = []                                                    # stress strain output for selected elements labels with ip -- in a list [el,ip]
    Restart = False
    VTK = False
    
    # DataExamples E03
#    Name="../DataExamples/E03/E3-01"                                        #  ok
#    Name="../DataExamples/E03/E3-02C3"                                      # C1, C2, C3 ok
    Name="../DataExamples/E03/E3-04"                                        # ok
    # E04 beams
#    Name="../DataExamples/E04/E3-02"                                        # ok
#    Name="../DataExamples/E04/E3-04"                                        # -03, -04, -06, -08, -09 ok
#    Name="../DataExamples/E04/E3-02_CircMises"                              # ok 2_CircMises"
    # E05 strut-and-tie
#    Name="../DataExamples/E05/E5-01plate"                                   # ok
#    Name="../DataExamples/E05/E5-02"                                        # ok
#    Name="../DataExamples/E05/E5-03"                                        # ok
    # E06 Multiaxial concrete
#    Name="../DataExamples/E06/E6-02"                                        # ok
    # E07 crack modeling and regularization
#    Name, Plot ="../DataExamples/E07/E7-01c", False                         # ok a, b, c, premature end with a
#    Name, Plot ="../DataExamples/E07/E7-05", False                          # ok 1D phase field
#    Name, Plot, StressStrainOut ="../DataExamples/E07-06/TensionEccTri2D", True, [] #[[1,0]]
#    Name, Plot ="../DataExamples/E07-06/E7-06a", False                      # ok SDA
    # E08
#    Name, Plot ="../DataExamples/E08/E8-01", True                           # ok deep beam elastic
#    Name, Plot = "../DataExamples/E08-02/E8-02", False                     # ok deep beam nonlinear
#    Name, Plot = "../DataExamples/E08/E8-03B23E", True                         # ok single fiber
#    Name, Plot = "../DataExamples/E08/E8-04", True                         # ok plate with reinforcement
    # E09                             
#    Name, Plot ="../DataExamples/E09/E9-01", True                           # ok elastic slab with opening and free edges
#    Name, Plot = "../DataExamples/E09/E9-04", True                          # ok Elastoplastic slab (NLSLAB) with opening and free edges
#    Name, Plot = "../DataExamples/E09/E9-05", True                           # ok simple Elastoplastic-extended (NLSLAB) slab under concentrated loading
    # E10
#    Name, Plot ="../DataExamples/E10/E10-01c", True                          # ok flat shell convergence / shear locking
#    Name, Plot ="../DataExamples/E10/E10-02a", True                           # one element RC shell
#    Name, Plot ="../DataExamples/E10-02/E10-02", False                     # simple slab nonlinear reinforcement
    # _DataAxisym
#    Name, Plot, StressStrainOut, VTK ="../_DataAxisym/BAX2E", True, [], False                   # ok
#    Name, Plot, StressStrainOut, VTK ="../_DataAxisym/CAX3", True, [], False                    # ok
#    Name, Plot, StressStrainOut, VTK = "../_DataAxisym/CylAX", True, [], False                  #ok
#    Name, Plot, StressStrainOut, VTK = "../_DataAxisym/TAX23", True, [], False                  # ok
    # _DataBeams
#    Name, Plot = "../_DataBeams/E3-02_B23", True                            # ok
#    Name, Plot = "../_DataBeams/E3-02_B21E", True                           # ok
#    Name, Plot = "../_DataBeams/E3-02_B23E_PolyL", True                     # ok                                                      # TestSuiteTrusses
    # _DataBenchmarks
#    Name, Plot ="../_DataBenchmarks/Arrea/arrea2D", True                           # ok
#    Name, Plot ="../_DataBenchmarks/Jenq/Jenq-L3-12122007-Final", True                 # ok
#    Name, Plot ="../_DataBenchmarks/Nooru/Nooru-1550-5", True                      # ok
#    Name, Plot, StressStrainOut ="../_DataBenchmarks/L_ShapedPanel/LSP_1", True, [[421,0]]     # ok
#    Name, Plot, StressStrainOut ="../_DataBenchmarks/L_ShapedPanel/LSP-2079-7_6_DAMMK", True, []   # ok
#    Name, Plot, StressStrainOut =" C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1_C3D8_", False, [[421,0]] # ok
#    Name, Plot, StressStrainOut = "../_DataBenchmarks/Aachentests/Deep_beam", True, []   # ok
#    Name, Plot, StressStrainOut = "../_DataBenchmarks/Aachentests/Deep_beam_AacNL", True, []   # ok
    # _DataBond
#    Name, Plot, StressStrainOut = "../_DataBond/bond_T2", True, []  # TestSuite0            # ok
#    Name, Plot, StressStrainOut = "../_DataBond/Pullout2", True, []   # Pullout, Pullout2           # ok
#    Name, Plot, StressStrainOut = "../_DataBond/PulloutAxiSym", True, []                                        # ok Isodamage
#    Name, Plot, StressStrainOut = "../_DataBond/PulloutLatPressure", True, []                    # ok
    # _DataBuckling
#    Name, nEig, niterEig, Plot ="../_DataBuckling/c_264_210000", 5, 200, True    # nEig, niterEig, sphere cap
#    Name, nEig, niterEig, Plot = "../_DataBuckling/E3-02V", 4, 200, True          # 2D Euler bat
#    Name, nEig, niterEig, Plot = "../_DataBuckling/fullCyl", 1, 200, True       # full cylinder under compression - lower end horizontally free
#    Name, nEig, niterEig, Plot = "../_DataBuckling/ShellBuckl", 5, 200, True    # simple flat shell with longitudinal compression
#    Name, nEig, niterEig, Plot = "../_DataBuckling/a_640_210000", 5, 200, True # full cylinder under compression -- lower end clamped
    # _DataC3D
#    Name, Plot, StressStrainOut, VTK ="../_DataC3D/Cube8", True, [], True                      # ok
#    Name, Plot, StressStrainOut, VTK ="../_DataC3D/Deep_beam_SDA", True, [], True              # ok
    # _DataRegularization
#    Name, StressStrainOut = "../_DataRegularization/TwoElementCPS3_SDATest", [[1,0],[2,0]]   # ok
#    Name, Plot, StressStrainOut ="../_DataRegularization/TwoElementCPS3_SDATest_6n", True, [[1,0],[2,0]]   # ok
    # _DataShellsSlabs
#    Name, Plot ="../_DataShellsSlabs/Shell", True #  ok "Shell"  "bridge_el05m" "c_1461(0.08)_2.1e5_0.3_segment load", "HyparRand.in"
    # _DataSpecimen
#    Name, Plot, StressStrainOut = "../_DataSpecimen/UniaxTensionBar", True, []  # ok
#    Name, Plot, StressStrainOut ="../_DataSpecimen/One3D/Uniax3D", True, [[1,0]]    # ok
#    Name, Plot ="../_DataSpecimen/One2D/WillamsTest2D", True                # ok
#    Name, Plot ="../_DataSpecimen/One3D/WillamsTest", False                 # ok microdamage, but no plotting avaliable for element type C3D8
    # _DataTrusses
#    Name, Plot = "../_DataTrusses/staebe_1_1d", True                        # ok looks like elastoplastic mises with un- and reloading
#    Name, Plot = "../_DataTrusses/staebe_4", True                                       # ok 2D elastic large displacement
#    Name, Plot = "../_DataTrusses/staebe3D2", True                             # ok 3D elastic large displacement
    # _DataTmp
#    Name, Plot, StressStrainOut = "../_DataTmp/shellconfem", True, []
#    Name, Plot, StressStrainOut = "../_DataTmp/BeamMises", True, []
#    Name, Plot = "../_DataTmp/staebe_4", True                                       # ok 2D elastic large displacement
    #
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2023/Note_ShearPlateRandom/ConFem/ShearPanel/ShearPanelR_med", True, [] #[[1,0]]
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2023/Note_ShearPlateRandom/ConFem/ShearPanel/Explicit/ShearPanelR", False, []
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2024/Note_FlatSlab/ConFem/PSlabN_1/Restart/PSlabN_1", True, []
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2024/Note_FlatSlab/ConFem/PSlabN_5/PSlabN_5__", True, []
    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/2024/Note_FlatSlab/ExpDataSets/_154-3F22/154-3F22", True, []
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/OneElement/2D/Large/ConFem/Job-1", True, [[1,0]]
#    Name, Plot, StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/OneElement/3D/Large/ConFem/Set-2", True, [[1,0],[1,1],[1,2],[1,3],[1,4],[1,5],[1,6],[1,7]]
#       C:\Users\uhc\Documents\Work\FoilPap\Abaqus\AbaqusModels\OneElement3D\Large\ConFem
#    Name, Plot = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/NotchedBeam/2D/Fine/ConFem/notchedBeam_2d_fine", True
#    Name, Plot,StressStrainOut = "C:/Users/uhc/Documents/Work/FoilPap/Abaqus/AbaqusModels/MixedMode/Specimen/Explicit/Willam/ConFem/Job-1", True, [[1,0]]

    Restart = False
#    Restart = True
    Name = Name.strip()
    return Name, Plot, Restart, "elemout", ElPlotTimes, StressStrainOut, VTK
