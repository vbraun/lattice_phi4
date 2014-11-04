
#include 'libcpp/vector.pxd'
#include 'libcpp/pair.pxd'

from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef extern from 'lattice.hpp':

    cdef cppclass CppLattice "Lattice":
        CppLattice(double mass_squared, double coupling_constant, int Nx, int  Ny)
        void randomize(int seed)
        void equilibrate(int iterations)
        int Nx, Ny
        vector[double] get_phi()
        double total_energy()
        double average_phi()
        double average_abs_phi()

cdef extern from 'simulation.hpp':

    cdef cppclass CppObservable "Observable<Lattice::RR>":
          CppObservable(CppObservable copy)
          double average()
          double abs_average()
          double square_average()
          double pow_average(double exponent)
          pair[double, double] minmax()
          vector[double] get()
          vector[double] autocorrelation()

    cdef cppclass CppSimulation "Simulation":
        CppSimulation(CppLattice)
        void run(int)
        CppObservable get_phi()
        CppObservable get_abs_phi()
        CppObservable get_energy()




