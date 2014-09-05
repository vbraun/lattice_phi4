#include <numeric>
#include <cmath>
#include <iterator>
#include <limits>
#include <algorithm>
#include <iostream>

#include "simulation.hpp"




template<typename Value>
Value Observable<Value>::average() const
{
  const Value sum = std::accumulate(data.begin(), data.end(), Value(0));
  return sum / data.size();
}

template<typename Value>
Value Observable<Value>::abs_average() const
{
  const Value sum = std::accumulate
    (data.begin(), data.end(), Value(0),
     [] (double partial_sum, double x){
      return partial_sum + std::abs(x); 
    });
  return sum / data.size();
}

template<typename Value>
Value Observable<Value>::square_average() const
{
  const Value sum = std::accumulate
    (data.begin(), data.end(), Value(0),
     [] (double partial_sum, double x){
      return partial_sum + x * x; 
    });
  return sum / data.size();
}

template<typename Value>
Value Observable<Value>::pow_average(const Value exponent) const
{
  const Value sum = std::accumulate
    (data.begin(), data.end(), Value(0),
     [&exponent] (double partial_sum, double x) {
      return partial_sum + std::pow(x, exponent); 
    });
  return sum / data.size();
}

template<typename Value>
std::pair<Value, Value> Observable<Value>::minmax() const
{
  auto pos = std::minmax_element(data.begin(), data.end());
  return std::make_pair(*(pos.first), *(pos.second));
}

template<typename Value>
std::vector<Value> Observable<Value>::autocorrelation() const
{
  std::vector<Value> result;
  const Value avg = average();
  Value previous_correlation = std::numeric_limits<Value>::max();
  for (auto i = data.begin(); i != data.end(); i++) {
    const auto separation = std::distance(data.begin(), i);
    const Value overlap = std::distance(i, data.end());
    Value sum = 0;
    for (auto j = data.begin(); j != data.begin() + overlap; j++)
      sum += (*j) * *(j + separation);
    const Value correlation = 
      (sum / overlap) - avg * avg;

    if (0 > correlation || correlation > previous_correlation)
      return result;

    result.push_back(correlation);
    previous_correlation = correlation;
  }
  return result;
}





Simulation::Simulation(Lattice& lattice_)
  : lattice(lattice_)
{}


void Simulation::run(unsigned int iterations)
{
  for (unsigned int i = 0; i < iterations; i++) {
    lattice.iteration();
    energy.data.push_back(lattice.total_energy());
    const RR average_phi = lattice.average_phi();
    phi.data.push_back(average_phi);
    abs_phi.data.push_back(std::abs(average_phi));
  }
}




template class Observable<double>;
