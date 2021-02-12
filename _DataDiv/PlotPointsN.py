import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem/PrePost')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("P_u")
Plot0(p3,'C:/Users\uhc/eclipse-workspace/CaeFem/Data/3D/Cube8.timeout.txt','green', 0, 1, 100.,'solid')
Plot0(p3,'Cube8.timeout.txt','red', 0, 1, 100.,'solid')

P3.autofmt_xdate()
plt.show()
