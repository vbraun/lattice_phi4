
from sage.all_cmdline import *

from phi4 import Simulation
from phase_plot import PhasePlot

plt = PhasePlot(32, 32, equilibrate=1512, iterations=1512, observable=Simulation.abs_phi)
plt.mass_squared = srange(-1.5, 0, 0.1)
plt.coupling_constant = srange(0.1, 1, 0.1)

gfx = \
      plt.plot3d() + \
      text3d(r'lambda', [0, 0.75, 1]) + \
      text3d(r'mu squared', [-0.5, 0, 0.2]) + \
      text3d(r'abs(phi)', [-1.5, 0, 2])

gfx.show()

