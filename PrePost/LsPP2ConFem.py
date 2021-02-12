import numpy as np
#import math


class Node(object):
    def __init__(self, Label, XCo, YCo, ZCo):
        self.Label = Label
        self.XCo = XCo
        self.YCo = YCo
        self.ZCo = ZCo
        self.DirL= []
        self.Director = []
        self.CLoadX = 0.
        self.CLoadY = 0.
        self.CLoadZ = 0.
        self.Used = False

class Element(object):
    def __init__(self, Label, ElType, Inzi, partId):
        self.Label  = Label
        self.ElType = ElType
        self.Inzi   = Inzi
        self.partId = partId
        
class NodeSet(object):
    def __init__(self, Label ):
        self.Label = Label
        self.Nodes = []

class BoundaryCondition(object):
    def __init__(self, Label ):
        self.Label = Label
        self.NSetInd = 0                                        # integer label of node set
        self.DofList = []                                       # markers for 6 dofs

class SLoad(object):                                            # load for node set
    def __init__(self, NSetInd, NSetDof, xyz, NSetVal ):
        self.NSetInd = NSetInd
        self.NSetDof = NSetDof                                  # 
        self.NSetxyz = xyz                                      # 
        self.NSetVal = NSetVal                                  # 
        self.Nodes  = []                                        # 

class PointLoad(object):
    def __init__(self, NoInd, Dof, Val ):
        self.NodeIndex = NoInd
        self.Dof = Dof
        self.Val = Val
        
class Material(object):
    def __init__(self, idMat, Emod, nu ):
        self.Type = "elastic"
        self.idMat= idMat
        self.Emod = Emod                                        # Young's modulus
        self.nu   = nu                                          # Poisson's ratio

def FindIndexByLabel(List, Key):                                # find index from label
    if len(List)==0: return None
    for i in range( len(List)):                                 # loop over all entries
        if List[i].Label == Key: return i
    return None                                                 # no entry found

class ShellSection(object):                         
    def __init__(self, Val, IdS ):
        self.Val = Val                                          # shell thickness
        self.Label = IdS

class BeamSection(object):                         
    def __init__(self, SecId, Val ):
        self.Label = SecId
        self.Val   = Val                                          # cross section

class LoadSegment(object):
    def __init__(self, Label, Val ):
        self.Label = Label
        self.Val   = Val

class Segments(object):
    def __init__(self, Ref ):
        self.Ref = Ref
        self.Segment = []

class Part(object):
    def __init__(self, Id, idSec, idMat ):
        self.Id       = Id
        self.idMat    = idMat
        self.idSec    = idSec
        self.ElementsSH4 = []
        self.ElementsSH3 = []
        self.ElementsB   = []                               # beam / truss elements
        self.ElementsT3D2= []
        self.ElementsC3D4= []
        self.ElementsC3D8= []
        self.Material    = []
        if idSec==0: self.SectionType = "solid"
        else:        self.SectionType = "structure"

def ReadkFile( fileIn ):
    NodeList, ElListAll, NoSet, BCond, PLoad, Mat, LoadSeg, SegmentList, ShellSec, BeamSecList = [],[],[],[],[],[],[],[],[],[]
    Parts, ElListB, ShellSecAll, SetLoadAll = [], [], [], []
    f1 = open(fileIn,'r')
    z1 = f1.readline()
    IType, IType2 = "", ""
    while z1!="":
        z2 = z1.strip()
        if z2:
            z3 = z2.split()
            print(z3)
            if z3[0]  =="*ELEMENT_SHELL":       IType = "element"
            elif z3[0]=="*ELEMENT_BEAM":        IType = "element_beam"
            elif z3[0]=="*ELEMENT_SOLID":       IType = "element_solid"
            elif z3[0]=="*NODE":                IType = "node"
            elif z3[0]=="*LOAD_NODE_POINT":     IType = "loadpoint"
            elif z3[0]=="*LOAD_NODE_SET":       IType, Flag = "setload", True
            elif z3[0]=="*BOUNDARY_SPC_SET":    IType = "pointc"                        # single point constraint without ID
            elif z3[0]=="*BOUNDARY_SPC_SET_ID": IType, Flag = "pointc1", True           # single point constraint with ID
            elif z3[0]=="*BOUNDARY_SPC_NODE":   IType = "pointcN"                       # single point constraint
            elif z3[0]=="*SET_NODE_LIST_TITLE": IType, Flag, LCounter = "nodeset", True, 0
            elif z3[0]=="*SET_NODE_LIST":       IType, Flag, LCounter = "nodeset", True, 1
            elif z3[0]=="*MAT_ELASTIC":         IType = "matElastic"                    # material without ID
            elif z3[0]=="*MAT_ELASTIC_TITLE":   IType, LCounter = "matElastic1", 0      # material without ID
            elif z3[0]=="*SECTION_SHELL":       IType = "shellsection"
            elif z3[0]=="*SECTION_SHELL_TITLE": IType, LCounter = "shellsection1", 0
            elif z3[0]=="*SECTION_BEAM_TITLE":  IType, LCounter = "beamsection1", 0
            elif z3[0]=="*LOAD_SEGMENT_SET":    IType = "loadSegment"
            elif z3[0]=="*LOAD_SEGMENT_SET_ID": IType, LCounter = "loadSegment1", 0
            elif z3[0]=="*SET_SEGMENT_TITLE":   IType, LCounter = "setSegment", 0
            elif z3[0]=="*SET_SEGMENT":         IType, LCounter = "setSegment", 1
            elif z3[0]=="*PART":                IType, LCounter = "part", 0
            elif z3[0] in ["*KEYWORD","*SECTION_SHELL_TITLE","*DEFINE_CURVE_TITLE"]: IType ="" 
            elif z3[0].find("$#")>-1:   pass                                    # comment line

            elif z3[0].find("*END")>-1: break
            
            elif IType!="":
                if IType  =="node":      NodeList += [Node( int(z3[0]), float(z3[1]), float(z3[2]), float(z3[3]))]
                elif IType=="element":
                    Inzi, InziL = [ int(z3[2]), int(z3[3]), int(z3[4]), int(z3[5])], []
                    for i in range(len(Inzi)): 
                        if Inzi[i]!=Inzi[i-1]: InziL += [Inzi[i]]               # to avoid multiple entries
                    if len(InziL)==4: ElType = 'SH4'
                    if len(InziL)==3: ElType = 'SH3'
                    ElListAll += [Element( int(z3[0]), ElType, InziL, int(z3[1])) ]
                elif IType=="element_solid":
                    Inzi, InziL = [ int(z3[2]), int(z3[3]), int(z3[4]), int(z3[5]), int(z3[6]), int(z3[7]), int(z3[8]), int(z3[9])], []
                    for i in range(len(Inzi)): 
                        if Inzi[i]!=Inzi[i-1]: InziL += [Inzi[i]]               # to avoid multiple entries
                    if    len(InziL)==8: ElType = 'C3D8'
                    elif  len(InziL)==4: ElType = 'C3D4'
                    else: raise NameError("ReadkFile:: something wrong with solid elements")
                    ElListAll += [Element( int(z3[0]), ElType, InziL, int(z3[1])) ]
                elif IType=="element_beam": ElListAll += [Element( int(z3[0]), 'T2D2', [int(z3[2]), int(z3[3])], int(z3[1]) )]
                elif IType=="loadpoint":    PLoad += [PointLoad( int(z3[0]), int(z3[1]), float(z3[3])) ]
                elif IType=="setload":      
                    SetLoad = SLoad(int(z3[0]), int(z3[1]), int(z3[2]), float(z3[3]))
                    SetLoadAll += [SetLoad]
                elif IType=="pointc":    
                        BCond += [BoundaryCondition(None) ]
                        BCond[-1].NSetInd = int(z3[0])
                        BCond[-1].DofList = [ int(z3[2]), int(z3[3]), int(z3[4]), int(z3[5]), int(z3[6]), int(z3[7]) ]
                elif IType=="pointc1":    
                    if Flag:
                        BCond += [BoundaryCondition(z3[0]) ]
                        Flag = False
                    else:
                        BCond[-1].NSetInd = int(z3[0])
                        BCond[-1].DofList = [ int(z3[2]), int(z3[3]), int(z3[4]), int(z3[5]), int(z3[6]), int(z3[7]) ]
                elif IType=="pointcN":
                    BCond += [BoundaryCondition(int(z3[0])) ]
                    BCond[-1].DofList = [ int(z3[2]), int(z3[3]), int(z3[4]), int(z3[5]), int(z3[6]), int(z3[7]) ]
                elif IType=="nodeset":
                    if LCounter==0: LCounter += 1                               # 1st "useless" line  
                    elif LCounter==1:           
                        NoSet += [ NodeSet(int(z3[0])) ]                        # 2nd line: index of node set
                        LCounter += 1
                    else: 
                        for i in z3:
                            if(int(i)>0): NoSet[-1].Nodes += [int(i)]           # nodes of node set
                elif IType=="matElastic": 
                    Mat += [Material( int(z3[0]), float(z3[2]), float(z3[3]) )]
                elif IType=="matElastic1":
                    if LCounter == 0:
                        LCounter += 1
                    else:
                        Mat += [Material( int(z3[0]), float(z3[2]), float(z3[3]) )]
                elif IType=="shellsection": 
                    if LCounter<1: 
                        IdS = int(z3[0])                                          # 1st line: section id
                        LCounter += 1  
                    else:
                        ShellSec = ShellSection( float(z3[0]), Ids) #  2nd line: 1st of 4 values used! 
                        ShellSecAll += [ShellSec]  
                elif IType=="shellsection1":
                    if LCounter<1: 
                        LCounter += 1                                           # 1st line: title line, not used here
                    elif LCounter<2:
                        IdS = int(z3[0])                                        # 2nd line section id
                        LCounter += 1
                    else: 
                        ShellSec = ShellSection( float(z3[0]), IdS)                # 3rd line: 1st of 4 values used!
                        ShellSecAll += [ShellSec]
                elif IType=="beamsection1":
                    if   LCounter == 0: LCounter += 1                                           # 1st line: title line, not used
                    elif LCounter == 1:                            
                        SecId = z3[0]
                        LCounter += 1
                    else: BeamSecList += [BeamSection( int(SecId), float(z3[0]))]        # 2nd line 1st value cross section
                elif IType=="loadSegment": LoadSeg += [ LoadSegment(int(z3[0]), float(z3[2])) ]
                elif IType=="setSegment":
                    if    LCounter<1: LCounter += 1                             # useless line
                    elif  LCounter<2: 
                        SegmentList += [Segments( int(z3[0]) )]                 # to which load segment id does this refer to
                        LCounter += 1
                    else: SegmentList[-1].Segment += [[ int(z3[0]),int(z3[1]),int(z3[2]),int(z3[3]) ]]
                elif IType=="part":
                    if LCounter==0: LCounter += 1                               # 1st "useless" line  
                    elif LCounter==1:           
                        Parts += [ Part( int(z3[0]), int(z3[1]), int(z3[2]) ) ] # 2nd line: part data
        z1 = f1.readline()
    f1.close()
    return NodeList, ElListAll, NoSet, BCond, PLoad, Mat, LoadSeg, SegmentList, Parts, BeamSecList, ShellSecAll, SetLoadAll

def NodalDirectors( NoList, ElList, fileLog ):
    f1 = open(fileLog,'w')
    for i in ElList:
        if i.ElType in ["SH4","SH3"]:
            NoCo, NoNo, NoIn = [], [], []
            for j in i.Inzi:
                k = FindIndexByLabel(NoList, j)
                if k in NoIn: continue                                          # for identical nodes
                no = NoList[k]
                NoCo += [[no.XCo,no.YCo,no.ZCo]]
                NoNo += [no]
                NoIn += [k]
            lNC = len(NoCo)
            for j in range(lNC):
                XX, YY, ZZ = NoCo[j][0],       NoCo[j][1],       NoCo[j][2]
                X1, Y1, Z1 = NoCo[j-1][0],     NoCo[j-1][1],     NoCo[j-1][2]
                X2, Y2, Z2 = NoCo[j-lNC+1][0], NoCo[j-lNC+1][1], NoCo[j-lNC+1][2]
                d1 = [X1-XX, Y1-YY, Z1-ZZ]
                d2 = [X2-XX, Y2-YY, Z2-ZZ]
                NN = np.cross( d2, d1 )
                lN = np.sqrt(NN[0]**2+NN[1]**2+NN[2]**2)
                NoNo[j].DirL += [NN/lN]                                         # NoNo holds handles for NoList, NN/lN is a single unit vector
    for i in NoList:
        print(i.Label, i.DirL, file=f1)
        l = len(i.DirL)
        if l==0: 
#            print("node ", i.Label, " seems not to belong to an element", file=f1)
#            print("node ", i.Label, " seems not to belong to an element")
            continue
        for k in range(1,l):                                                # for control purposes
            for j in range(l-k):
                cc = np.dot(i.DirL[k-1],i.DirL[j+k])
                if cc<0.9: 
                    print('YYY', i.Label, i.DirL, '__', k-1, j+k, cc)
                    raise NameError("Seems to be a folded shell, not yet covered!")
        xD, yD, zD = 0., 0., 0.
        for j in i.DirL:
            xD += j[0]
            yD += j[1]
            zD += j[2]
        i.Director = [ xD/l, yD/l, zD/l]                                    # mean of all element-directors on this node
        print(i.Label, i.Director, file=f1)
    f1.close()
    return 0

def AssignElementsMaterialsToParts( ElListAll, Materials, Parts):
    for i in range(len(Parts)):
        part  = Parts[i]
        Id    = part.Id
        for j in range(len(ElListAll)):
            elem = ElListAll[j]
            if elem.partId==Id: 
                if    elem.ElType=='SH4':  part.ElementsSH4 += [j]
                elif  elem.ElType=='SH3':  part.ElementsSH3 += [j]
                elif  elem.ElType=='T2D2': part.ElementsB   += [j]
                elif  elem.ElType=='C3D4': part.ElementsC3D4 += [j]
                elif  elem.ElType=='C3D8': part.ElementsC3D8 += [j]
                elif  elem.ElType=='BEAM': part.ElementsT3D2 += [j]
                else: raise NameError("unknown element type")
        for m in Materials:
            if m.idMat == part.idMat: part.Material += [m]
        if len(part.Material) != 1: raise NameError("wrong material assignment for part")     
    return 0 

def UsedNodes( NoList, ElemList):
    for el in ElemList:
        for i in el.Inzi:
            j = FindIndexByLabel( NoList, i)
            NoList[j].Used = True
    for no in NoList:
        if not no.Used: print('node ',no.Label,' seems not to be used')

def SegmentLoad( NoList, ElListAll, LoadSeg, SegmentList ):
    for i in SegmentList:
        j = FindIndexByLabel(LoadSeg, i.Ref)
        val = LoadSeg[j].Val
        totArea = 0.
        counter = 0
        for k in i.Segment:                                                 # k -> list of 3/4 coordinates --> single segments                      
            if k[-1]==k[-2]: del k[-1]                                      # to avoid duplicate entries
            CoorLi, NodeLi = [], []
            for l in k:                                                     # loop over nodes of segment - first determine total load of segment
                n = FindIndexByLabel(NoList, l)
                try:    node = NoList[n]
                except: raise NameError("node of SET_SEGMENT not found")
                CoorLi += [[node.XCo,node.YCo,node.ZCo]]
                NodeLi += [node]
            if len(CoorLi)==3:                                              # segment surface area - triangle
                d1 = [ CoorLi[1][0]-CoorLi[0][0], CoorLi[1][1]-CoorLi[0][1], CoorLi[1][2]-CoorLi[0][2] ] 
                d2 = [ CoorLi[2][0]-CoorLi[0][0], CoorLi[2][1]-CoorLi[0][1], CoorLi[2][2]-CoorLi[0][2] ]
                dd = np.cross(d1,d2) 
                AA = 0.5*np.sqrt(dd[0]*dd[0]+dd[1]*dd[1]+dd[2]*dd[2])       # area of triangle in space
                totArea += AA
                counter += 1
#                print 'XXX', counter, AA, totArea
                loadval = AA*val/3.
            elif len(CoorLi)==4:                                            # segment surface area - quadrilateral
                d1 = [ CoorLi[1][0]-CoorLi[0][0], CoorLi[1][1]-CoorLi[0][1], CoorLi[1][2]-CoorLi[0][2] ] 
                d2 = [ CoorLi[2][0]-CoorLi[0][0], CoorLi[2][1]-CoorLi[0][1], CoorLi[2][2]-CoorLi[0][2] ]
                d3 = [ CoorLi[3][0]-CoorLi[0][0], CoorLi[3][1]-CoorLi[0][1], CoorLi[3][2]-CoorLi[0][2] ]
                dd1= np.cross(d1,d2) 
                A1 = 0.5*np.sqrt(dd1[0]*dd1[0]+dd1[1]*dd1[1]+dd1[2]*dd1[2]) # area of triangle in space
                dd2= np.cross(d2,d3) 
                A2 = 0.5*np.sqrt(dd2[0]*dd2[0]+dd2[1]*dd2[1]+dd2[2]*dd2[2]) # area of triangle in space
                AA = A1+A2
                totArea += AA
                counter += 1
            else: raise NameError("to much/less nodes for segments")
            n = len(NodeLi)
            for no in NodeLi:                                               # loop over nodes of segment - second determine concentrated loads for each node of segment 
                CL = AA*val/n                                               # concentrated load per node
                no.CLoadX = no.CLoadX + CL*no.Director[0]
                no.CLoadY = no.CLoadY + CL*no.Director[1]
                no.CLoadZ = no.CLoadZ + CL*no.Director[2]
        totXYZ = 0.
        for no in NoList:                                                   # sum up concentrated load for control
            totXYZ += np.sqrt( no.CLoadX*no.CLoadX + no.CLoadY*no.CLoadY + no.CLoadZ*no.CLoadZ )
        print("segment "+str(i.Ref)+" surface area "+ str(totArea) + ", load from pressure: "+ str(val*totArea)+ ", load from C-loads: "+ str(totXYZ))

def WriteConFile( fileOut, NoList, NoSet, BCond, PLoad, Materials, Parts, ElListAll, BeamSecList, ShellSecAll, SetLoadAll, SegmentList, Flag3D, CaeFem ):
    el = ["EL1", "EL2", "EL3", "EL4", "EL5" ]                               # should be no more then 5 element sets in this data set
    mat= ["mat1","mat2","mat3","mat4","mat5"]                               # should be no more then 5 element sets in this data set
    f1 = open(fileOut,'w')
    f1.write('*HEADING\n shell\n')
    #
    f1.write('*NODE\n')
    for i in NoList:
        if len(i.Director)>0: f1.write('%5i,%12.4e,%12.4e,%12.4e,    %12.4e,%12.4e,%12.4e\n'%(i.Label,i.XCo,i.YCo,i.ZCo,i.Director[0],i.Director[1],i.Director[2]))
        else:                 f1.write('%5i,%12.4e,%12.4e,%12.4e\n'                         %(i.Label,i.XCo,i.YCo,i.ZCo))
    #
    print('Nodes written')
    for i in range(len(Materials)):
        mat = Materials[i]
        if mat.Type == "elastic":
            f1.write('*MATERIAL, NAME=%s\n*ELASTIC\n'%("MAT"+str(mat.idMat)))
            f1.write(' %12.2e,  %8.4f\n'%( mat.Emod, mat.nu))
    #
    print('Materials written')
    i_ = 1
    if Flag3D: 
        SecType, ElType4, ElType3 = '*SHELL SECTION', 'SH4', 'SH3'
    else:      
        SecType, ElType4, ElType3 = '*SOLID SECTION','CPS4', 'CPS3'
    for i in range(len(Parts)):
        part = Parts[i]
        mat  = part.Material[0]
        if len(part.ElementsSH4)>0:
            j = FindIndexByLabel(ShellSecAll, part.idSec)
            f1.write('%s, ELSET=%s, MATERIAL=%s\n%10.3e\n'%(SecType,"EL"+str(i_),"MAT"+str(mat.idMat),ShellSecAll[j].Val))
            f1.write('*ELEMENT, ELSET=%s, TYPE=%s\n'%("EL"+str(i_), ElType4))
            for j in part.ElementsSH4:
                el =  ElListAll[j]
                f1.write('%6i,  %10i,%10i,%10i,%10i\n'%(el.Label,el.Inzi[0],el.Inzi[1],el.Inzi[2],el.Inzi[3]))
            i_ += 1
        if len(part.ElementsSH3)>0:
            j = FindIndexByLabel(ShellSecAll, part.idSec)
            f1.write('%s, ELSET=%s, MATERIAL=%s\n%10.3e\n'%(SecType,"EL"+str(i_),"MAT"+str(mat.idMat),ShellSecAll[j].Val))
            f1.write('*ELEMENT, ELSET=%s, TYPE=%s\n'%("EL"+str(i_), ElType3))
            for j in part.ElementsSH3: 
                el =  ElListAll[j]
                f1.write('%6i,  %10i,%10i,%10i\n'     %(el.Label,el.Inzi[0],el.Inzi[1],el.Inzi[2]))
            i_ += 1
        if len(part.ElementsB)>0:
            j = FindIndexByLabel(BeamSecList, part.idSec)
            f1.write('*SOLID SECTION, ELSET=%s, MATERIAL=%s\n%10.3e\n'%("EL"+str(i_),"MAT"+str(mat.idMat),BeamSecList[j].Val))
            f1.write('*ELEMENT, ELSET=%s, TYPE=%s\n'%("EL"+str(i_), 'T2D2'))
            for j in part.ElementsB: 
                el =  ElListAll[j]
                f1.write('%6i,  %10i,%10i\n'     %(el.Label,el.Inzi[0],el.Inzi[1]))
            i_ += 1
        if len(part.ElementsC3D4)>0:
            f1.write('*SOLID SECTION, ELSET=%s, MATERIAL=%s\n%10.3e\n'%("EL"+str(i_),"MAT"+str(mat.idMat),1.0))
            f1.write('*ELEMENT, ELSET=%s, TYPE=%s\n'%("EL"+str(i_), 'C3D4'))
            for j in part.ElementsC3D4: 
                el =  ElListAll[j]
                f1.write('%6i,  %10i,%10i,%10i,%10i\n'     %(el.Label,el.Inzi[0],el.Inzi[1],el.Inzi[2], el.Inzi[3]))
            i_ += 1
        if len(part.ElementsC3D8)>0:
            raise NameError("not yet implemented")
    print('Sections and elements written')
        
    #
    f1.write('*STEP\n')
    if CaeFem:  f1.write('*CONTROLS, ITOL= 1.e-5, NITER= 20\n')
    else:       f1.write('*CONTROLS, PARAMETERS=FIELD\n    1.e-5\n*CONTROLS, PARAMETERS=TIME INCREMENTATION\n   100\n')
    f1.write('*STATIC\n   1.0,   1.0\n')
    f1.write('*AMPLITUDE, NAME=AMPLB\n 0.0, 0.0,  1.0, 1.0\n')
    f1.write('*BOUNDARY, OP=NEW\n')
    if Flag3D: DofMax = 5                                   # ConFem shell has 5 dofs per node
    else:      DofMax = 2                                   # ConDem plate has 2 dofs per node

###
    for i in BCond:
        DofC = []
        for j in range(len(i.DofList)):
            if i.DofList[j] == 1 and j<DofMax:              # 1 --> dof constrained  
                DofC += [j+1]                               # build list of constrained dof markers
        k = FindIndexByLabel(NoSet, i.NSetInd)              # find index of assigned node set
        if k!=None:
            for j in NoSet[k].Nodes:
                for jj in DofC:
                    if CaeFem:  f1.write('%6i,  %10i,  0.\n'%(j,jj))
                    else:       f1.write('%6i,  %10i,%10i,  0.\n'%(j,jj,jj))
        else:
            for jj in DofC:
                if CaeFem:  f1.write('%6i,  %10i,  0.\n'%(i.Label,jj))
                else:       f1.write('%6i,  %10i,%10i,  0.\n'%(i.Label,jj,jj))
 
    # for point load 
    if len(PLoad)>0:                                                        # point load
        f1.write('*AMPLITUDE, NAME=AMPLP\n 0.0, 0.0,  1.0, 1.0\n')
        f1.write('*CLOAD, AMPLITUDE=AMPLP\n')
        xyz = 0.
        for i in PLoad: 
            f1.write('%6i,  %10i,%12.4e\n'%(i.NodeIndex,i.Dof,i.Val))
            xyz += i.Val
        print('total point load for ',len(PLoad),'nodes:',xyz)
    # for point load segments
    for SLoad in SetLoadAll:
        k = FindIndexByLabel(NoSet, SLoad.NSetInd)              # find index of assigned node set
        f1.write('*AMPLITUDE, NAME=AMPL'+str(SLoad.NSetInd)+'\n 0.0, 0.0,  1.0, 1.0\n')
        f1.write('*CLOAD, AMPLITUDE=AMPL'+str(SLoad.NSetInd)+'\n')
        for j in NoSet[k].Nodes:
             f1.write('%6i,  %10i,%12.4e\n'%(j,SLoad.NSetDof,SLoad.NSetVal))
        print('total load in node-set',SLoad.NSetInd,':  ',len(NoSet[k].Nodes)*SLoad.NSetVal)
    # for segment load - pressure
    if len(SegmentList) > 0:
        f1.write('*AMPLITUDE, NAME=AMPLS\n 0.0, 0.0,  1.0, 1.0\n')
        f1.write('*CLOAD, AMPLITUDE=AMPLS\n')
        for i in NoList:
            if abs(i.CLoadX)>0.: f1.write('%6i,  %10i,%12.4e\n'%(i.Label,1,i.CLoadX))
            if abs(i.CLoadY)>0.: f1.write('%6i,  %10i,%12.4e\n'%(i.Label,2,i.CLoadY))
            if abs(i.CLoadZ)>0.: f1.write('%6i,  %10i,%12.4e\n'%(i.Label,3,i.CLoadZ))
        print('generated segment pressure')

    if CaeFem:  f1.write('*EL FILE, FREQUENCY=1.0\n*NODE FILE, FREQUENCY=1.0\n*END STEP')
    else:       f1.write('*EL FILE, FREQUENCY=1\nS\nE\n*NODE FILE, FREQUENCY=1\nU\n*EL PRINT, FREQUENCY=1\nS\n*NODE PRINT, FREQUENCY=1\nU\nRF\n*END STEP')
    print('Step written')
    f1.close()
    return 0

if __name__ == '__main__':
 #   Name, Flag3D, CaeFem = "a_6.9e4_625", True, False
#    Name, Flag3D = "a_SECTION_SHELL_TITLE", True
#    Name, Flag3D = "a_SECTION_SHELL", True
#    Name, Flag3D, CaeFem = "c_1461(0.08)_2.1e5_0.3_segment load", True, False
#    Name, Flag3D, CaeFem = "temp_model_", True, True
#    Name, Flag3D = "c_5503_210000", True
#    Name, Flag3D = "HighBeam2", False
    Name, Flag3D, CaeFem = "SHB", False, True
    Name, Flag3D, CaeFem = "Hypar", True, False
    Name, Flag3D, CaeFem = "temp_model2", False, False
    
    Name=str(input('Filename without extension: '))
    
    fileIn = Name+".k"
    if CaeFem:  fileOut= Name+".txt"
    else:       fileOut= Name+".in.txt"
    fileLog= Name+".log.txt"
    NodeList, ElemListAll, NodeSet, BoundC, CLoad, Materials, LoadSeg, SegmentList, Parts, BeamSecList, ShellSecAll, SetLoadAll = ReadkFile( fileIn )
    RC = AssignElementsMaterialsToParts( ElemListAll, Materials, Parts)
    RC = NodalDirectors( NodeList, ElemListAll, fileLog )
    RC = UsedNodes( NodeList, ElemListAll )
    RC = SegmentLoad( NodeList, ElemListAll, LoadSeg, SegmentList )
    RC = WriteConFile( fileOut, NodeList, NodeSet, BoundC, CLoad, Materials, Parts, ElemListAll, BeamSecList, ShellSecAll, SetLoadAll, SegmentList, Flag3D, CaeFem )
    print('finish')
