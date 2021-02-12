
def ReadLine( ff):
    z1 = ff.readline()
    if z1 != "":
        z2 = z1.strip()
        z3 = z2.split()
        return z3
    else: 
        return None
def ReadDataElemOut(fName, Time, column):
    rList, counter = [], 0
    ff = open(fName,'r')
    while True:
        zz = ReadLine(ff)
        if zz == None: 
            break
        else:
            if zz[0] == "Time": 
                Time_ = float(zz[1])
                continue
            if zz[0] == "El":
                Elem = int(zz[1])
                ElSet= zz[3]
                counter += 1
                continue
            if Time_ == Time:
                rList += [[Elem,float(zz[0]),float(zz[1]),float(zz[2]),float(zz[column-1]),ElSet]]
    ff.close()
    return rList, counter
def WriteData( fName, rList ):
    ff = open(fName,'w')
    for i in rList:
        if   len(i)==5: print(f'{i[0]:6d}, {i[1]:8.3f}, {i[2]:8.3f}, {i[3]:8.3f}, {i[4]:12.4e}', file=ff )
        elif len(i)==6: print(f'{i[0]:6d}, {i[1]:8.3f}, {i[2]:8.3f}, {i[3]:8.3f}, {i[4]:12.4e}, {i[5]:s}', file=ff )
    ff.close()   
def ReadDataNodeOut(fName, Time, column):
    rList, counter = [], 0
    ff = open(fName,'r')
    while True:
        zz = ReadLine(ff)
        if zz == None: 
            break
        else:
            if len(zz) == 1: 
                Time_ = float(zz[0])
                continue
            if Time_ == Time:
                counter += 1
                rList += [[int(zz[0]),float(zz[1]),float(zz[2]),float(zz[3]),float(zz[column-1])]]
    ff.close()
    return rList, counter

if __name__ == "__main__":
    Dir  = "C:/Users/uhc/ConFemExe/ConFem/"
    Name = "HyparRand"
    ElColumn = 5
    NoColumn = 7
    Time = 1.0
    # exe versiond
    Dir = ""
    Name=str(input('Filename without extension: '))
    ElColumn=int(input('Column for elemeout: '))
    NoColumn=int(input('Column for nodeout: '))
    #
    fName = Dir+Name+".elemout.txt"
    rList, counter = ReadDataElemOut(fName, Time, ElColumn)
    rList.sort(key=lambda t: t[4])
    print("checked ",counter," elements in column ",ElColumn)
    print("min ",rList[0:5])
    print("max ",rList[-5:-1])
    fName = Dir+Name+".elemout.col" + str(ElColumn)+".txt"
    WriteData( fName, rList )
    #
    fName = Dir+Name+".nodeout.txt"
    rList, counter = ReadDataNodeOut(fName, Time, NoColumn)
    rList.sort(key=lambda t: t[4])
    print("checked ",counter," nodes in column ",NoColumn)
    print("min ",rList[0:5])
    print("max ",rList[-5:-1])
    fName = Dir+Name+".nodeout.col" + str(NoColumn)+".txt"
    WriteData( fName, rList )
    #
    print("finished")
