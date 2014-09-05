

from sage.all_cmdline import *

from phi4 import Simulation



class PhasePlot(object):
    
    def __init__(self, Nx, Ny, equilibrate=4096, iterations=4096, observable=None):
        self._Nx = Nx
        self._Ny = Ny
        self._equilibrate = equilibrate
        self._iterations = iterations
        self._mu2 = -1.26
        self._mu2_is_range = False
        self._lambda = 1.0
        self._lambda_is_range = False
        if observable is None:
            self._obs = Simulation.susceptibility
        else:
            self._obs = observable
        self._title = r'Lattice simulation of $\phi^4$ theory in 2d'
        self._legend = None
        self._observable = None
        self._color = 'blue'

    @property
    def mass_squared(self):
        return self._mu2

    @mass_squared.setter
    def mass_squared(self, value):
        try:
            self._mu2 = tuple(sorted(value))
            self._mu2_is_range = True
        except TypeError:
            self._mu2 = RDF(value)
            self._mu2_is_range = False

    @property
    def coupling_constant(self):
        return self._lambda

    @coupling_constant.setter
    def coupling_constant(self, value):
        try:
            self._lambda = tuple(sorted(value))
            self._lambda_is_range = True
        except TypeError:
            self._lambda = RDF(value)
            self._lambda_is_range = False

    @parallel
    def _eval(self, mass_squared, coupling_constant, **kwds):
        sim = Simulation(mass_squared, coupling_constant, 
                         self._Nx, self._Ny, equilibrate=self._equilibrate)
        sim.run(self._iterations)
        return self._obs(sim)

    def _parallel_eval(self, points):
        args = [((mass_squared, coupling_constant), dict(axis=axis))
                for mass_squared, coupling_constant, axis in points]
        result = self._eval(args)
        result = list(result)
        result = [tuple(list(r[0][1]['axis']) + [r[1]]) for r in result]
        return sorted(result)

    def plot(self, *args, **kwds):
        if not(self._mu2_is_range or self._lambda_is_range):
            s =  'Plotting requires a range for mass_squared and/or coupling_constant.'
            value = self._eval(self._mu2, self._lambda)
            s += 'Value of observable at single point is ' + str(value)
        elif not self._mu2_is_range:
            return self._plot_lambda(*args, **kwds)
        elif not self._lambda_is_range:
            return self._plot_mu2(*args, **kwds)
        else:
            return self._plot_2d(*args, **kwds)
            
    def _plot_1d(self, points, x_axis_label, y_axis_label, legend, title, color):
        result = self._parallel_eval(points)
        return list_plot(result, plotjoined=True, color=color) + \
            list_plot(result, title=title, 
                      axes_labels=[x_axis_label, y_axis_label], 
                      legend_label=legend, color=color)
        
    @property
    def legend(self):
        if self._legend is not None:
            return self._legend
        return self._default_legend

    @legend.setter
    def legend(self, name):
        self._legend = name

    def set_default_legend(self, legend, add_lattice=True, add_samples=True):
        result = [legend]
        if add_lattice:
            result.append(r'${0}\times{1}$'.format(self._Nx, self._Ny))
        if add_samples:
            result.append(r'{0} samples'.format(self._iterations))
        self._default_legend = ', '.join(result)

    @property
    def observable(self):
        if self._observable is not None:
            return self._observable
        return self._obs.__name__

    @observable.setter
    def observable(self, name):
        self._observable = name

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, title):
        self._title = title

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, name):
        self._color = name
    
    def _plot_lambda(self):
        points = [(self._mu2, l, [l]) for l in self._lambda]
        self.set_default_legend(r'$\mu^2 = {0}$'.format(self._mu2))
        return self._plot_1d(points, r'$\lambda$', self.observable,
                             self.legend, self.title, self.color)

    def _plot_mu2(self):
        points = [(m, self._lambda, [m]) for m in self._mu2]
        self.set_default_legend(r'$\lambda = {0}$'.format(self._lambda))
        return self._plot_1d(points, r'$\mu^2$', self.observable,
                             self.legend, self.title, self.color)
        
    def _plot_2d(self, cmap='winter', vmin=None, vmax=None):
        points = [(m, l, [i,j]) 
                  for i, m in enumerate(self._mu2)
                  for j, l in enumerate(self._lambda)]
        result = self._parallel_eval(points)
        value_min = vmin if vmin is not None else min(r[2] for r in result)
        value_max = vmax if vmax is not None else max(r[2] for r in result)
        
        delta_mu2 = [(self._mu2[i] - self._mu2[i-1])/2.0
                     for i in range(1, len(self._mu2))]
        delta_lambda = [(self._lambda[i] - self._lambda[i-1])/2.0
                        for i in range(1, len(self._lambda))]

        from sage.plot.colors import get_cmap
        cmap = get_cmap(cmap)
        def tile(i, j, value):
            m = self._mu2[i]
            if i == 0:
                dm = (delta_mu2[0], delta_mu2[0])
            elif i == len(delta_mu2):
                dm = (delta_mu2[-1], delta_mu2[-1])
            else:
                dm = (delta_mu2[i-1], delta_mu2[i])
            l = self._lambda[j]
            if j == 0:
                dl = (delta_lambda[0], delta_lambda[0])
            elif j == len(delta_lambda):
                dl = (delta_lambda[-1], delta_lambda[-1])
            else:
                dl = (delta_lambda[j-1], delta_lambda[j])
            norm_value = (value - value_min) / (value_max - value_min)
            print value, norm_value
            c = cmap(norm_value)[:3]
            v0 = (l-dl[0], m-dm[0])
            v1 = (l+dl[1], m-dm[0])
            v2 = (l+dl[1], m+dm[1])
            v3 = (l-dl[0], m+dm[1])
            return polygon2d([v0, v1, v2, v3], rgbcolor=c, fill=True)
        g = Graphics()
        g.axes_labels([r'$\lambda$', r'$\mu^2$'])
        g._set_extra_kwds({
            'title': self.title,
            'axes_labels': [r'$\lambda$', r'$\mu^2$'],
            'title_pos': (0.5, 1.1),
        })
        for i, j, value in result:
            g = tile(i, j, value) + g
        return g

    def plot3d(self):
        points = [(m, l, [m, l]) 
                  for i, m in enumerate(self._mu2)
                  for j, l in enumerate(self._lambda)]
        result = self._parallel_eval(points)
        return list_plot3d(result, point_list=True)
        m = matrix(RDF, len(self._mu2), len(self._lambda))
        for i_m, i_l, obs in result:
            m[i_m, i_l] = obs
        return list_plot3d(m)


# values = phase_plot_mass_squared(1.0)
# print values


def susceptibility(mass_squared):
    N = 32
    coupling_constant = 1.0
    sim = Simulation(mass_squared, coupling_constant, N, N, equilibrate=500)
    sim.run(500)
    return sim.susceptibility()



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
