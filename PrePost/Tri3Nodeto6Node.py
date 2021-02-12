import sys
sys.path.insert(1, '../src')
import scipy.spatial as spatial
from ConFemInOut import ReadInputFile
from ConFemElem import Node
from numpy import array 

if __name__ == '__main__':
    Name, Plot, StressStrainOut ="C:/Users/uhc/Documents/Work/FoilPap/2020/Note_ConFemBenchmarks/LShapedP/ConFem/LSP_4", False, [[421,0]]
#    Name, Plot, StressStrainOut ="../_DataRegularization/TwoElementCPS3_SDATest", False, [[421,0]]
    f1=open( Name+".in.txt", 'r')
    f2=open( Name+"_6n.in.txt", 'w')
    NodeList, ElList, MatList, StepList, NoLabToNoInd, SecDic = ReadInputFile(f1, None, False) # read input file
    # build tree of new nodes
    newNodes = []
    for el in ElList:
        if el.Type in ["CPS3","CPE3"]:
            for i in range(3):
                i0 = el.Inzi[i]
                try:    i1 = el.Inzi[i+1]
                except: i1 = el.Inzi[0]
                no0 = NodeList[i0] 
                no1 = NodeList[i1]
                newNodes += [[0.5*(no0.XCo+no1.XCo), 0.5*(no0.YCo+no1.YCo), 0.]]
    if len(newNodes)>0: CoorTree = spatial.cKDTree( newNodes )
    else:               CoorTree = None
    newUsedNodes, extendInzi = {}, {}                                       # key of newUsedNodes is given by index of item in CoorTree
    # consolidate list of new nodes and extend Inzlist
    for el in ElList:
        if el.Type in ["CPS3","CPE3"]:
            Inzi_ = []
            for i in range(3):
                i0 = el.Inzi[i]
                try:    i1 = el.Inzi[i+1]
                except: i1 = el.Inzi[0]
                no0 = NodeList[i0] 
                no1 = NodeList[i1]
                x = 0.5*(no0.XCo+no1.XCo)
                y = 0.5*(no0.YCo+no1.YCo)
                res = CoorTree.query([array([x, y, 0.])],k=1)
                j = res[1][0]                                               # index of node in CoorTree
                if j not in newUsedNodes: newUsedNodes[j] = [x,y]
                Inzi_ += [j]
            extendInzi[el.Label] = Inzi_ 
    # determine largest node label
    maxLab = 0
    for no in NodeList:
        if no.Label>maxLab: maxLab = no.Label
    maxLab += 1
    # add labels to new nodes and write nodal data
    newUseNodesL = {}
    for i, n in enumerate(newUsedNodes):                                    # n is key of dictionary item
        j = newUsedNodes[n]                                                 # data of dictionary item 
        newUseNodesL[n] = i+maxLab
        print(f' {i+maxLab:d}, {j[0]:e}, {j[1]:e}, 0.', file=f2)
    # collect and write element data
    for el in ElList:
        if el.Type in ["CPS3","CPE3"]:
            n1 = NodeList[el.Inzi[0]]
            n2 = NodeList[el.Inzi[1]]
            n3 = NodeList[el.Inzi[2]]
            In = extendInzi[el.Label]                                       # list of indices from CoorTree
            n4 = newUseNodesL[In[0]]
            n5 = newUseNodesL[In[1]]
            n6 = newUseNodesL[In[2]]
            print(f'{el.Label:d}, {n1.Label:d}, {n2.Label:d}, {n3.Label:d}, {n4:d}, {n5:d}, {n6:d}', file=f2)
    #
    print('finished')

