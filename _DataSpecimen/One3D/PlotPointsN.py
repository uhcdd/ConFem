import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("P_u")

#Plot0(p3,'UniaxTension.timeout.txt','red', 0, 1, 2.) # large element
Plot0(p3,'../One2D/References/Uniax2D_MIPL.timeout.txt','magenta', 0, 1, 100.,'dashdot') # small
Plot0(p3,'../One2D/References/Uniax2D_LUBL.timeout.txt','blue', 0, 1, 100.,'dashdot') # small
Plot0(p3,'../One2D/References/Uniax2D_ISOD.timeout.txt','green', 0, 1, 100.,'dashdot') # small
Plot0(p3,'Uniax3D.timeout.txt','red', 0, 1, 200.,'solid') # small
Plot0(p3,'./References/Uniax3D_LUBL.timeout.txt','blue', 0, 1, 200.,'solid') # small
Plot0(p3,'./References/Uniax3D_ISOD.timeout.txt','green', 0, 1, 200.,'solid') # small
Plot0(p3,'./References/Uniax3D_MIPL.timeout.txt','magenta', 0, 1, 200.,'solid') # small

P3.autofmt_xdate()
plt.show()
