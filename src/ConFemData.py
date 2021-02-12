def DefData():
    nEig, niterEig, ElPlotTimes = None, None, []
    Plot = True #True 
    StressStrainOut = []
    Restart = False
    VTK = False
    
    # DataExamples
    Name="../DataExamples/E04/E2-01" #2C1"                     
#    Name="../DataExamples/E04/E2-04"                     
##    Name="../DataExamples/E05/E3-06" # 2_CircMises"
#    Name="../DataExamples/E06/E4-01" #02u03" #, -01plate, -01                 
##    Name="../DataExamples/E08/E6-02"
##    Name="../DataExamples/E08-03/E6-03"
##    Name="../DataExamples/E09/E7-01"                   
##    Name="../DataExamples/E09/E9-05c" 
##    Name="../DataExamples/E10/E8-01a"                   
##    Name, Plot ="../DataExamples/E10-02/E8-01du02", False
    # _DataTrusses
#    Name = "../_DataTrusses/staebe_4"
    # DataShellsSlabs                
#    Name="../_DataShellsSlabs/c_1461(0.08)_2.1e5_0.3_segment load"  #  "Shell"  "bridge_el05m"
#    _DataSpecimen
#    Name, Plot, StressStrainOut ="../_DataSpecimen/One2D/WillamsTest2D", False, [[1,0]]
#    Name, Plot ="../_DataSpecimen/One3D/WillamsTest", False
#    Name, Plot, StressStrainOut  ="../_DataSpecimen/One2D/Uniax2D", False, [[1,0]]
    # DataBuckling
#    Name, nEig, niterEig, Plot, StressStrainOut ="../_DataBuckling/c_264_210000", 5, 200, True, []
    # _DataBenchmarks
##    Name="../_DataBenchmarks/Arrea/arrea2D"
##    Name="../DataBenchmarks/Jenq/Jenq-L3-12122007-Final"
#    Name, Plot, StressStrainOut ="../_DataBenchmarks/L_ShapedPanel/LSP-2079-7_6_DAMMK", False, [[6,0]]
##    Name="../DataBenchmarks/Nooru/Nooru-1550-5"
#    Name, Plot, StressStrainOut ="C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1", False, [[421,0]]
#    Name, Plot, StressStrainOut ="C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1_C3D8_", False, [[421,0]]
    # _DataTemp
##    Name, Plot ="../_DataTmp/phi-45coarse", True
    Name, Plot, StressStrainOut ="../_DataTmp/Deep_beam_AacNL", False, []
    # _DataDiv
#    Name, Plot, StressStrainOut ="../_DataDiv/bond_T2", True, []
#    Name, Plot, StressStrainOut ="../_DataDiv/Cube8", False, []
    # Regularization
#    Name = "../_DataRegularization/TwoElementCPS3_SDATest"
#    Name, Plot, StressStrainOut ="../_DataRegularization/TwoElementCPS3_SDATest_6n", True, [[1,0],[2,0]]

    return Name, Plot, Restart, "elemout", [ nEig, niterEig ], ElPlotTimes, StressStrainOut, VTK
