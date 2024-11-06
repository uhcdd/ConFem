import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem/PrePost')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("Uniax Tension")
# 'solid', 'dashdot', 

sc_x, sc_y = 1, 1
Plot0(p3,'BAX2.timeout_BAX23.txt','green',    0,1, sc_x,sc_y,'solid','-')
Plot0(p3,'BAX2E.timeout_BAX21E.txt','blue',    0,1, sc_x,sc_y,'solid','-')
Plot0(p3,'BAX2.timeout_BAX21.txt','red',    0,1, sc_x,sc_y,'solid','-')
Plot0(p3,'BAX2.timeout.txt','magenta',    0,1, sc_x,sc_y,'solid','-')

p3.legend(loc='upper right', prop={'size': 10})
p3.set_xlabel('long strain %')
p3.set_ylabel('long stress MN/m^2')
P3.autofmt_xdate()
plt.show()
