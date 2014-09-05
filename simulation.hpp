#ifndef SIMULATION__HPP
#define SIMULATION__HPP

#include <vector>

#include "lattice.hpp"


class Simulation;


template<typename Value>
class Observable
{
private:
  friend class Simulation;

protected:
  std::vector<Value> data;

public:
  const std::vector<Value> get() const { return data; };
  Value average() const;
  Value abs_average() const;
  Value square_average() const;
  Value pow_average(const Value exponent) const;
  std::pair<Value, Value> minmax() const;
  std::vector<Value> autocorrelation() const;
};



class Simulation
{
public:
  typedef Lattice::RR RR;

protected:
  Lattice& lattice;
  Observable<RR> phi;
  Observable<RR> abs_phi;
  Observable<RR> energy;  
  
public:
  Simulation(Lattice& lattice);
 
  // Run the simulation for a given number of iterations
  /* There is no warmup / equilibration performed. It is your
   * responsibility to construct with a suitable lattice.
   */
  void run(unsigned int iterations);

  const Observable<RR> get_phi() const { return phi; };
  const Observable<RR> get_abs_phi() const { return abs_phi; };
  const Observable<RR> get_energy() const { return energy; };
};



#endif // SIMULATION__HPP
