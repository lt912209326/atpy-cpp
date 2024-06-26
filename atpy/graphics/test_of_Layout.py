from layout import Layout
from data import data
from data2 import data2
import matplotlib.pyplot as plt
lat=Layout("STCF")
lat.add_beamline(data)
lat.add_beamline(data2,drift_color="red")
# lat.ax.axis([-110,110,-110,110])
lat.fig.set_facecolor("white")
# lat.fig.savefig("stcf_layout.svg",dpi=80)
plt.show()
