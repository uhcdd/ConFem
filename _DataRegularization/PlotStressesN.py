import sys
sys.path.insert(1, '../')
import imp
import PlotStressesT
imp.reload(PlotStressesT)
from PlotStressesT import *

P3,p3, P4,p4, P5,p5, P6,p6 = Assemble(' ')
Plot0(p3,p4,p5,p6,'TwoElementCPS3_SDATest.elemout.1_0.txt','','red',200.)
P3.autofmt_xdate()
plt.show()
