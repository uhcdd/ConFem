import sys
sys.path.insert(1, '../')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("P_u")
#Plot0(p3,'Plane-el.timeout - Kopie.txt','blue', 0, 1)
#Plot0(p3,'Plane-el.timeout - Kopie (2).txt','magenta', 0, 1)
#Plot0S('C:\Users\uhc\Workspace\CaeFem3\Data\_tmp\Ahmad\Deep_beam_AacNL.timeout.txt', p3, 1., 1.,'green')
#Plot0(p3,'Plane-el.timeout.txt','red', 0, 1)
#Plot0(p3,'Deep_beam_AacNL.timeout.txt','red', 0, 1)
#Plot0(p3,'Deep_beam_AacNL.timeout - Kopie.txt','green', 0, 1)
#Plot0(p3,'Uniax2D.timeout.txt','green', 0, 1)
#Plot0(p3,'Uniax2D_Lubl.timeout.txt','red', 0, 1)
Plot0(p3,'TwoElementCPS3_SDATest_6n.timeout.txt','red', 0, 1)
P3.autofmt_xdate()
plt.show()
