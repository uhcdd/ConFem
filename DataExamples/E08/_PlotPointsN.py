import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem/PrePost')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("L-D")
# 'solid', 'dashdot', 

sc_x, sc_y = -1000.0, -1000.0
Plot0(p3,'..\..\DataExamplesRef\E08\E8-04.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','-')
#Plot0(p3,'E8-03B23E.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','-')
Plot0(p3,'E8-04.timeout.txt','red',    0,1, sc_x,sc_y,'solid','-')

#p3.legend(loc='upper right', prop={'size': 10})
p3.set_xlabel('deflection [mm]')
p3.set_ylabel('loading [kN]')
P3.autofmt_xdate()
plt.show()
