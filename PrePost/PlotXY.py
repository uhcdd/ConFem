import matplotlib.pyplot as plt

def DefPlot(Label, ResList):
    P0 = plt.figure()
    p0 = P0.add_subplot(111)
    p0.set_title(Label,fontsize='x-large')
    p0.tick_params(axis='x', labelsize='large')
    p0.tick_params(axis='y', labelsize='large')
    p0.grid()
    p1 = None
    for i in ResList:
        if i[0]==1: 
            p1 = p0.twinx()
            break
    return P0, p0, p1 

def PlotXY_ser(fileNames, pp, pt, ResList, ColList, x_Factor, y_Factor, x_Range):
    PlotList, XList = {}, []
    for i in ResList:
        PlotList[i[1]] = []
    for i in fileNames:
        ff = open(i,'r')
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split()
        while z1!="":
            x_Val = float(z3[0])
            if len(x_Range)>0 and x_Val<x_Range[0]: continue
            XList += [x_Factor*x_Val]
            for j in ResList:
                ind = j[1]
                PlotList[ind] += [float(z3[ind])]
            if len(x_Range)>0 and x_Val>x_Range[1]: break
            z1 = ff.readline()
            z2 = z1.strip()
            z3 = z2.split()

    for i, j in enumerate(ResList):
        if j[0] == 0: pp.plot(XList,PlotList[j[1]], color=ColList[i])
        else:         pt.plot(XList,PlotList[j[1]], color=ColList[i])

def PlotXY_par(fileNames, pp, pt, ResList, ColList, x_Factor, y_Factor):
    PlotList, XList = {}, {}
    for i in FileNames:
        PlotList[i] = []
        XList[i]    = []
    for i in fileNames:
        ff = open(i,'r')
        z1 = ff.readline()
        z2 = z1.strip()
        z3 = z2.split()
        while z1!="":
            XList[i] += [x_Factor*float(z3[0])]
            for j in ResList:
                ind = j[1]
                PlotList[i] += [y_Factor*float(z3[ind])]
            z1 = ff.readline()
            z2 = z1.strip()
            z3 = z2.split()

    for i, s in enumerate(FileNames):
        for j in ResList:
            if j[0] == 0: pp.plot(XList[s],PlotList[s], color=ColList[i])
            else:         pt.plot(XList[s],PlotList[s], color=ColList[i])

if __name__ == "__main__":
    
#    FileNames, ResList, PType, ColorList, x_Factor, y_Factor, x_Range = ['../Data/_tmp/dumbbell-test-random-fiber.ResS0/0_OutCra.txt',
#                 '../Data/_tmp/dumbbell-test-random-fiber.restart.ResS-1/0_OutCra.txt',
#                 '../Data/_tmp/dumbbell-test-random-fiber.restart2.S-1/0_OutCra.txt',
#                 '../Data/_tmp/dumbbell-test-random-fiber.restart3.S-1/0_OutCra.txt'],\
#                 [[0,5],[0,6],[1,7]], "ser",  ['darkorange','forestgreen','steelblue'], 0.1, 1.0, [0.0, 25]
##                 [[0,9],[0,10],[0,11],[0,12]], "ser", ['darkorange','darkorange','forestgreen','forestgreen'], 0.1, 1.0, [0, 25]
    FileNames, ResList, x_Factor, y_Factor, PType = ['../DataExamples/E09/E9-05c.timeout_50.txt',
                             '../DataExamples/E09/E9-05c.timeout_25.txt',
                             '../DataExamples/E09/E9-05c.timeout_00.txt'],\
                 [[0,1]], -1., -4., "par"
    
    P0, p0, p1 = DefPlot('Label', ResList)
    if   PType == 'par': PlotXY_par(FileNames, p0, p1, ResList, ['steelblue','darkorange','forestgreen','forestgreen'], x_Factor, y_Factor)
    elif PType == 'ser': PlotXY_ser(FileNames, p0, p1, ResList, ColorList, x_Factor, y_Factor, x_Range )
    P0.autofmt_xdate()
    plt.show()
    print('finish')