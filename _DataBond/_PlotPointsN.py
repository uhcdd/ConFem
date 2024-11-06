import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem/PrePost')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("Uniax Tension")
# 'solid', 'dashdot', 

#ll = 0.004 #0.01
#sc_x, sc_y = 100./ll, 1./ll
#Plot0(p3,'Uniax2D.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','ISOD')
#Plot0(p3,'Uniax2D.timeout - Kopie.txt','red',    0,1, sc_x,sc_y,'solid','ISOD')

sc_x, sc_y = 1, 1
Plot0(p3,'PulloutAxisym.timeout - Kopie.txt','red',    0,1, sc_x,sc_y,'solid','-')
Plot0(p3,'PulloutAxiSym.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','-')

#sc_x, sc_y = 1, 1
#Plot0(p3,'phi-45coarse.timeout - Kopie.txt','red',    0,1, sc_x,sc_y,'solid','-')
#Plot0(p3,'phi-45coarse.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','-')

p3.legend(loc='upper right', prop={'size': 10})
p3.set_xlabel('long strain %')
p3.set_ylabel('long stress MN/m^2')
P3.autofmt_xdate()
plt.show()
