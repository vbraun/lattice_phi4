
from sage.all_cmdline import *

from phi4 import Simulation
from phase_plot import PhasePlot



# values = phase_plot_mass_squared(1.0)
# print values




plt = PhasePlot(32, 32, equilibrate=1512, iterations=1512)
#plt.mass_squared = srange(-2, 0, 0.2)
#plt.coupling_constant = 1.0
plt.mass_squared = srange(-1.5, 0, 0.1)
plt.coupling_constant = srange(0.01, 1, 0.1)
# print plt.plot()


# plt.plot(vmax=5).show()
#  plt.plot3d().show()

gfx = \
      plt.plot3d() + \
      text3d(r'lambda', [0, 0.75, 1]) + \
      text3d(r'mu squared', [-1, 0, 1]) + \
      text3d(r'susceptibility', [-1.5, 0, 25])


gfx.show()

