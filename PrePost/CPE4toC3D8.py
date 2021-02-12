import sys
sys.path.insert(1, '../src')
import scipy.spatial as spatial
from ConFemInOut import ReadInputFile
from ConFemElem import Node
from numpy import array 

if __name__ == '__main__':
    OffsetZCoor = 1.0
    OffsetLabel = 2000
    Name, Plot ="C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_1", False
#    Name, Plot, StressStrainOut ="../_DataRegularization/TwoElementCPS3_SDATest", False, [[421,0]]
    f1=open( Name+".in.txt", 'r')
    f2=open( Name+"_C3D8.in.txt", 'w')
    NodeList, ElList, MatList, StepList, NoLabToNoInd, SecDic = ReadInputFile(f1, None, False) # read input file
    #
    print(f'*HEADING\n {Name:s} to C3D8', file=f2)
    print(f'*NODE', file=f2)
    for i, n in enumerate(NodeList):                                    # n is key of dictionary item
        print(f' {n.Label:10d}, {n.XCo:e}, {n.YCo:e}, 0.', file=f2)
        print(f' {n.Label+OffsetLabel:10d}, {n.XCo:e}, {n.YCo:e}, {OffsetZCoor:e}', file=f2)
    #
    print(f'***********Materials****', file=f2)
    #
    print(f'************************', file=f2)
    #
    print(f'*SOLID SECTION,ELSET=EL1,MATERIAL=???\n   1.0', file=f2)
    print(f'*ELEMENT, ELSET=EL1, TYPE=C3D8', file=f2)
    for el in ElList:
        if el.Type in ["CPS4","CPE4"]:
            n0 = NodeList[el.Inzi[0]]
            n1 = NodeList[el.Inzi[1]]
            n2 = NodeList[el.Inzi[2]]
            n3 = NodeList[el.Inzi[3]]
            n4 = n0.Label + OffsetLabel
            n5 = n1.Label + OffsetLabel
            n6 = n2.Label + OffsetLabel
            n7 = n3.Label + OffsetLabel
            print(f'{el.Label:d}, {n0.Label:d}, {n1.Label:d}, {n2.Label:d}, {n3.Label:d}, {n4:d}, {n5:d}, {n6:d}, {n7:d}', file=f2)
    #
    for st in StepList:
        AmpDic = {}
        print(f'****************************************************', file=f2)
        print(f'*STEP', file=f2)
        print(f'*CONTROLS, PARAMETERS=FIELD\n   ???\n*CONTROLS, PARAMETERS=TIME INCREMENTATION\n ???', file=f2)
        print(f'*STATIC\n   1.0, 1.0', file=f2)
        if len(st.BoundList) > 0:
            for bo in st.BoundList:
                li = [bo.NodeLabel,bo.Dof,bo.Val]
                Ad = bo.AddVal
                if AmpDic.get(bo.Amplitude,0) == 0:
                    AmpDic[bo.Amplitude] = [li]
                else:
                    AmpDic[bo.Amplitude] += [li]
        for am in AmpDic:
            print(f'*BOUNDARY, AMPLITUDE={am:s}, OP=???', file=f2)
            li = AmpDic[am]
            for bo in li:
                print(f' {bo[0]:d},  {bo[1]:d}, {bo[1]:d}, {bo[2]:f}', file=f2)
                print(f' {bo[0]+OffsetLabel:d},  {bo[1]:d}, {bo[1]:d}, {bo[2]:f}', file=f2)
        #
        print(f'*EL FILE, FREQUENCY=1.0\nS\nE\n*NODE FILE, FREQUENCY=1.0\nU', file=f2)
        print(f'*END STEP', file=f2)
    
    print(f'****************************************************', file=f2)
    print('finished')

