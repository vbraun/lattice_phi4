import math

include "phi4.pxi"


# TODO
def cached_method(x):
    return x


cdef class Observable(object):

    cdef CppObservable* thisptr

    def __cinit__(self):
        self.thisptr = NULL

    def __init__(self):
        raise RuntimeError('You cannot instantiate Observable by hand')

    def __dealloc__(self):
        del self.thisptr

    def __len__(self):
        return self.thisptr.get().size()
        
    def __getitem__(self, unsigned int i):
        """
        Return the ``i``-th element.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()
            sage: v[1]
            456
        """
        if not (i >= 0 and i < self.thisptr.get().size()):
            raise ValueError('index is out of bounds')
        return self.thisptr.get().at(i)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def samples(self):
        return self.thisptr.get()

    def _check_nonempty(self):
        if len(self) == 0:
            raise RuntimeError('no samples present, run the simulation first')

    @cached_method
    def average(self):
        self._check_nonempty()
        return self.thisptr.average()

    @cached_method
    def abs_average(self):
        self._check_nonempty()
        return self.thisptr.abs_average()

    @cached_method
    def square_average(self):
        self._check_nonempty()
        return self.thisptr.square_average()

    @cached_method
    def pow_average(self, exponent):
        self._check_nonempty()
        return self.thisptr.pow_average(exponent)

    @cached_method
    def autocorrelation(self):
        return self.thisptr.autocorrelation()
        
    @cached_method
    def minmax(self):
        self._check_nonempty()
        return self.thisptr.minmax()

    def standard_deviation(self, autocorrelation_time):
        """
        Return the standard deviation

        INPUT:

        - ``autocorrelation_time`` -- float. Autocorrelation time
          measured in units of iterations of the simulation. See
          :meth:`Simulation.autocorrelation_time`.

        OUTPUT:

        Float. The standard deviation, where we only treat each
        ``autocorrelation_time`` interval as statistically
        independent.
        """
        s1 = self.average()
        s2 = self.square_average()
        t = float(autocorrelation_time) / len(self)
        return math.sqrt(2 * t * (s2 - s1**2))
    
    def histogram(self, bins=21, lower_bound=None, upper_bound=None):
        minmax = self.minmax()
        if lower_bound is None:
            lower_bound = minmax[0]
        if upper_bound is None:
            upper_bound = minmax[1]
        bin_width = (upper_bound - lower_bound) / float(bins)
        result = [0 for i in range(bins)]
        for x in self:
            pos = int((x - lower_bound) / bin_width)
            if pos < 0 or pos >= bins:
                continue
            result[pos] += 1
        return result


cdef wrap_Observable(CppObservable observable):
    cdef Observable result = Observable.__new__(Observable)
    result.thisptr = new CppObservable(observable)
    return result




cdef class Simulation(object):

    cdef CppLattice* lattice
    cdef CppSimulation* sim

    def __cinit__(self):
        self.lattice = NULL
        self.sim = NULL

    def __dealloc__(self):
        del self.sim
        del self.lattice

    def __init__(self, mass_squared, coupling_constant, 
                 n_x=32, n_y=32, random_seed=0, equilibrate=16384):
        r"""
        Construct a Lattice QFT simulation

        This is the main object holding the data for a lattice
        simulation of the 2-dimensional `\phi^4` theory.

        INPUT:

        - ``mass_squared`` -- float. The mass-squared bare parameter
          in the Lagrangian (can be negative, its just a name).

        - ``coupling_constant`` -- float. The coupling constant
          multiplying the `\phi^4` interaction.

        - ``n_x`` -- positive integer (default: 32). The number of
          lattice points in the x-direction.

        - ``n_y`` -- positive integer (default: 32). The number of
          lattice points in the y-direction.

        - ``random_seed`` -- positive integer (default: 0). The random
          seed to use for the simulation. A value of 0 means to not
          randomize, use this for exactly reproducable computations.
        
        - ``equilibrate`` -- positive integer (default: 16384). The
          number of lattice iterations to perform from the random
          starting point to get close to an equilibrium.

        EXAMPLES::

            sage: from phi4 import Simulation
            sage: sim = Simulation(-0.16, 1.0)
            sage: sim.lattice()
        """
        self.lattice = new CppLattice(mass_squared, coupling_constant, n_x, n_y)
        self.lattice.randomize(random_seed)
        self.lattice.equilibrate(1000)
        self.sim = new CppSimulation(self.lattice[0])
        
    def run(self, iterations=16384):
        """
        Run the simulation.

        INPUT:
        
        - ``iterations`` -- positive integer (default: 16384). The
          number of iterations to run for.

        This method may be called more than once. Every call to
        :meth:`run` will add ``iterations`` samples to the
        observables.
        """
        self.sim.run(iterations)
    
    def phi(self):
        return wrap_Observable(self.sim.get_phi())

    def abs_phi(self):
        return wrap_Observable(self.sim.get_abs_phi())

    def energy(self):
        return wrap_Observable(self.sim.get_energy())
    
    def autocorrelation_time(self):
        correlation = self.abs_phi().autocorrelation()
        if len(correlation) == 1:
            return 0
        result = 0
        for i in range(1, len(correlation)):
            result += i / math.log(correlation[0] / correlation[i])
        return result / (len(correlation) - 1)

    def matrix(self):
        """
        Return the final field configuration as matrix.
        """
        Nx, Ny = self.shape()
        from sage.rings.real_double import RDF
        from sage.matrix.constructor import matrix
        return matrix(RDF, Nx, Ny, self.lattice.get_phi())
        
    _matrix_ = matrix

    def shape(self):
        """
        Return the shape of the lattice
        """
        return (self.lattice.Nx, self.lattice.Ny)

    def specific_heat(self):
        Nx, Ny = self.shape()
        energy = self.energy()
        return (energy.square_average() - energy.average() ** 2) * (Nx * Ny)
        
    def Binder_cumulant(self):
        r"""
        Return the fourth-order cumulant.

        .. math::

            U = 1 - \frac{\langle \phi^4 \rangle}{3 \langle \phi^2 \rangle^2}
        """
        phi = self.phi()
        phi_squared = phi.square_average()
        phi_fourth = phi.pow_average(4)
        return 1 - phi_fourth / (3 * phi_squared * phi_squared);

    def susceptibility(self):
        Nx, Ny = self.shape()
        abs_phi = self.abs_phi()        
        return (abs_phi.square_average() - abs_phi.average() ** 2) * (Nx * Ny)



