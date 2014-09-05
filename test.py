from sage.all_cmdline import *

from phi4 import Simulation


sim = Simulation(-1.26, 1.00, equilibrate=1000)
sim.run(1000)

phi = sim.phi()
print len(phi), phi.average(), phi.abs_average(), phi.pow_average(0.3)

print phi.minmax()
print phi.histogram()

abs_phi = sim.abs_phi()
print len(abs_phi), abs_phi.average(), abs_phi.abs_average(), abs_phi.pow_average(0.3)

energy = sim.energy()
print len(energy), energy.average(), energy.abs_average(), energy.pow_average(0.3)

auto_corr = abs_phi.autocorrelation()
# list_plot(auto_corr, scale='semilogy')


t_corr = sim.autocorrelation_time()

print '---'
print sim.specific_heat(), 1.016279
print sim.Binder_cumulant()
print sim.susceptibility(), 25.724055

print sim.matrix()
