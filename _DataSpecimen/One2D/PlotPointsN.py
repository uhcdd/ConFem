import sys
sys.path.insert(1, 'C:/Users/uhc/eclipse-workspace/ConFem/PrePost')
import imp
import PlotPointsT
imp.reload(PlotPointsT)
from PlotPointsT import *

P3, p3 = DefPlot("Willams Test - size 0.1 m")
# 'solid', 'dashdot', 

case = 'WT'
###### Uniax2D
#Plot0(p3,'.\References\Uniax2D_LUBL.timeout.txt','green', 0, 1, 100.,'solid','LUBL')
#Plot0(p3,'.\References\Uniax2D_ISOD.timeout.txt','blue', 0, 1, 100.,'solid','ISOD')
#Plot0(p3,'.\References\Uniax2D_MIPL.timeout.txt','red', 0, 1, 100.,'solid','MIPL')
#Plot0(p3,'Uniax2D.timeout - Kopie.txt','magenta', 0, 1, 100.,'solid')
#Plot0(p3,'Uniax2D.timeout.txt','red', 0, 1, 100.,'solid')
#sc_x, sc_y = 10*100., 100.      # scaling factors strain %/displ large element
#Plot0(p3,'Uniax2D_ISOD.timeout.txt','blue',    0,1, sc_x,sc_y,'solid','ISOD')
#Plot0(p3,'Uniax2D_MIPL.timeout.txt','magenta', 0,1, sc_x,sc_y,'solid','MIPL')
#Plot0(p3,'Uniax2D_LUBL.timeout.txt','green',     0,1, sc_x,sc_y,'solid','LUBL')
#Plot0(p3,'Uniax2D_ELLT.timeout.txt','red',       0,1, sc_x,sc_y,'solid','ELLT')
#sc_x, sc_y = 100*100., 1000.      # scaling factors strain small element
#sc_x, sc_y = 1000., 1000.      # scaling factors displ small element
#Plot0(p3,'Uniax2D_ISOD_s.timeout.txt','blue',    0,1, sc_x,sc_y,'dashed','ISOD_s')
#Plot0(p3,'Uniax2D_MIPL_s.timeout.txt','magenta', 0,1, sc_x,sc_y,'dashed','MIPL_s')
#Plot0(p3,'Uniax2D_LUBL_s.timeout.txt','green',     0,1, sc_x,sc_y,'dashed','LUBL_s')
#Plot0(p3,'Uniax2D_ELLT_d.timeout.txt','red',       0,1, sc_x,sc_y,'dashed','ELLT_s')

##### TensionCompression
if case == 'TC':
    sc_x, sc_y = 1000., 200.      # scaling factors large element
    Plot0(p3,'TensionCompression2D_ISOD_2_c.timeout.txt','blue', 0, 1, sc_x,sc_y,'solid','ISOD_2')
    Plot0(p3,'TensionCompression2D_ISOD_2_t.timeout.txt','blue', 0, 1, sc_x,sc_y,'solid','ISOD_2')
    Plot0(p3,'TensionCompression2D_MIPL_2_t.timeout.txt','magenta', 0, 1, sc_x,sc_y,'solid','MIPL_2')
    Plot0(p3,'TensionCompression2D_MIPL_2_c.timeout.txt','magenta', 0, 1, sc_x,sc_y,'solid','MIPL_2')
    Plot0(p3,'TensionCompression2D_LUBL_t.timeout.txt','green', 0, 1, sc_x,sc_y,'solid','LUBL')
    Plot0(p3,'TensionCompression2D_LUBL_c.timeout.txt','green', 0, 1, sc_x,sc_y,'solid','LUBL')
    sc_x, sc_y = 100*100., 10*200.      # scaling factors small element
    Plot0(p3,'TensionCompression2D_ISOD_2_t_s.timeout.txt','blue', 0, 1, sc_x,sc_y,'dashed','ISOD_2_s')
    Plot0(p3,'TensionCompression2D_ISOD_2_c_s.timeout.txt','blue', 0, 1, sc_x,sc_y,'dashed','ISOD_2_s')
    Plot0(p3,'TensionCompression2D_MIPL_2_t_s.timeout.txt','magenta', 0, 1, sc_x,sc_y,'dashed','MIPL_2_s')
    Plot0(p3,'TensionCompression2D_MIPL_2_c_s.timeout.txt','magenta', 0, 1, sc_x,sc_y,'dashed','MIPL_2_s')
    Plot0(p3,'TensionCompression2D_LUBL_t_s.timeout.txt','green',     0,1, sc_x,sc_y,'dashed','LUBL_s')
    Plot0(p3,'TensionCompression2D_LUBL_c_s.timeout.txt','green',     0,1, sc_x,sc_y,'dashed','LUBL_s')

#### Biax2D a is for eps_x/eps_y=1, b for eps_x/eps_y=3/2
elif case == 'BI':
    sc_x, sc_y = 10*100., 100.      # scaling factors large element
    Plot0(p3,'Biax2D_ISOD_2_a.timeout.txt','blue', 0, 1, sc_x,sc_y,'solid','ISOD_2_a')
    Plot0(p3,'Biax2D_ISOD_2_b.timeout.txt','blue', 0, 1, sc_x,sc_y,'dashed','ISOD_2_b')
    Plot0(p3,'Biax2D_MIPL_2_a.timeout.txt','magenta', 0, 1, sc_x,sc_y,'solid','MIPL_2_a')
    Plot0(p3,'Biax2D_MIPL_2_b.timeout.txt','magenta', 0, 1, sc_x,sc_y,'dashed','MIPL_2_b')
    Plot0(p3,'Biax2D_LUBL_a.timeout.txt','green', 0, 1, sc_x,sc_y,'solid','LUBL_a')
    Plot0(p3,'Biax2D_LUBL_b.timeout.txt','green', 0, 1, sc_x,sc_y,'dashed','LUBL_b')
#    Plot0(p3,'Biax2D.timeout.txt','red', 0, 1, sc_x,sc_y,'solid','ISOD_2_a')

# y for tension, s for shear
elif case == 'TS':
    sc_x, sc_y = 10*100., 100.      # scaling factors large element
    Plot0(p3,'TensionShear2D_ISOD_2_y.timeout.txt','blue', 0, 1, sc_x,sc_y,'dashed','ISOD_2_y')
    Plot0(p3,'TensionShear2D_ISOD_2_s.timeout.txt','blue', 0, 1, sc_x,sc_y,'solid','ISOD_2_s')
    Plot0(p3,'TensionShear2D_MIPL_2_y.timeout.txt','magenta', 0, 1, sc_x,sc_y,'dashed','MIPL_2_y')
    Plot0(p3,'TensionShear2D_MIPL_2_s.timeout.txt','magenta', 0, 1, sc_x,sc_y,'solid','MIPL_2_s')
    Plot0(p3,'TensionShear2D_LUBL_y.timeout.txt','green', 0, 1, sc_x,sc_y,'dashed','LUBL_y')
    Plot0(p3,'TensionShear2D_LUBL_s.timeout.txt','green', 0, 1, sc_x,sc_y,'solid','LUBL_s')
#    Plot0(p3,'TensionShear2D.timeout.txt','red', 0, 1, sc_x,sc_y,'solid','ISOD')

##### Willams Test
elif case == 'WT':
    sc_x, sc_y = 10*100., 100.
#    Plot0(p3,'.\References\WillamsTest2D_ISOD_2.timeout.txt','blue', 0, 1, 1.)
#    Plot0(p3,'.\References\WillamsTest2D_LUBL.timeout.txt','green', 0, 1, 1.)
##    Plot0(p3,'WillamsTest2D.timeout - Kopie.txt','magenta', 0, 1, 1.)
    Plot0(p3,'WillamsTest2D_ISOD_2.timeout.txt','blue', 0, 1,sc_x,sc_y,'solid','ISOD_2')
    Plot0(p3,'WillamsTest2D_MIPL_2.timeout.txt','magenta', 0, 1,sc_x,sc_y,'solid','MIPL_2')
    Plot0(p3,'WillamsTest2D_LUBL.timeout.txt','green', 0, 1,sc_x,sc_y,'solid','LUBL')
    Plot0(p3,'WillamsTest2D_ELLT.timeout.txt','red', 0, 1,sc_x,sc_y,'solid','ELLT')
#    Plot0(p3,'WillamsTest2D.timeout.txt','black', 0, 1,sc_x,sc_y,'solid','-')
    sc_x, sc_y = 10*100., 1000.
    Plot0(p3,'WillamsTest2D_ISOD_2_s.timeout.txt','blue', 0, 1,sc_x,sc_y,'dashed','ISOD_2_s')
    Plot0(p3,'WillamsTest2D_MIPL_2_s.timeout.txt','magenta', 0, 1,sc_x,sc_y,'dashed','MIPL_2_s')
    Plot0(p3,'WillamsTest2D_LUBL_s.timeout.txt','green', 0, 1,sc_x,sc_y,'dashed','LUBL_s')
    Plot0(p3,'WillamsTest2D_ELLT_s.timeout.txt','red', 0, 1,sc_x,sc_y,'dashed','ELLT_s')
#    Plot0(p3,'WillamsTest2D.timeout.txt','black', 0, 1,sc_x,sc_y,'solid','-')

p3.legend(loc='upper right', prop={'size': 10})
p3.set_xlabel('shear strain % / displ mm')
p3.set_ylabel('stress tau_xy MN/m^2')
P3.autofmt_xdate()
plt.show()
