#ifndef LATTICE__HH
#define LATTICE__HH

#include "lattice_iterator.hpp"

#include <vector>
#include <unordered_set>
#include <cstdint>

#include <gsl/gsl_rng.h>     // Random number generators



template<typename Site> 
class LatticeNeighborhoodXY {
public:
  Site site;
  Site prev_x;
  Site next_x;
  Site prev_y;
  Site next_y;
  LatticeNeighborhoodXY(Site site_, 
                        Site prev_x_, Site next_x_, 
                        Site prev_y_, Site next_y_)
    : site(site_), prev_x(prev_x_), next_x(next_x_), prev_y(prev_y_), next_y(next_y_) {};
};


template<typename RR>
class LatticeFieldXY
{
public:
  RR& value;
  RR& prev_x;
  RR& next_x;
  RR& prev_y;
  RR& next_y;
  LatticeFieldXY(RR& value_, 
                 RR& prev_x_, RR& next_x_, 
                 RR& prev_y_, RR& next_y_)
    : value(value_), prev_x(prev_x_), next_x(next_x_), prev_y(prev_y_), next_y(next_y_) {};
};




class Lattice {
public:

  // Index of a lattice point
  typedef std::size_t ZZ;

  // Field values
  typedef double RR;

private:
  // The range of "random" fields: [-PHI_RANGE, PHI_RANGE]
  /* This is used in the initial randomization and in field updates
   * during the iteration.
   */
  const RR PHI_RANGE = 1.5;

  // The number of Metropolis sweeps per iteration
  const ZZ METROPOLIS_SWEEPS_PER_ITERATION = 5;

  // The GNU Scientific Library random number generator 
  gsl_rng* rng;  

protected:
  // Return a uniformly-distributed random number in [0, 1]
  inline RR random() const {  return gsl_rng_uniform(rng); };
  inline RR random_phi() const {  return PHI_RANGE * (2*random() - 1); };

  // Index of a single site without considering neigbors
  typedef const ZZ Site;
  // Index of a lattice site together with indices of its 4 nearest neigbors
  typedef LatticeNeighborhoodXY<const Site> SiteXY;

  // Field value at a single Site
  typedef RR& Phi;
  typedef const RR& ConstPhi;
  // Field value at site together with its 4 nearest neigbors
  typedef LatticeFieldXY<RR> PhiXY;
  typedef LatticeFieldXY<const RR> ConstPhiXY;

  // Construct a SiteXY instance from a Site
  inline SiteXY get_site_xy(Site site) const;

  // Construct a Phi instance from a Site
  inline Phi      get_phi(Site site);
  inline ConstPhi get_phi(Site site) const;
  // Construct a PhiXY instance  from a Site
  inline PhiXY      get_phi_xy(Site site);
  inline ConstPhiXY get_phi_xy(Site site) const;

  // Return a random site;
  inline Site random_site() const;
  inline SiteXY random_site_xy() const;

public:
  friend class ConstLatticeSiteIterator<Lattice const>;
  typedef ConstLatticeSiteIterator<Lattice const> const_iterator;
  const_iterator begin() const { return const_iterator(*this, 0); };
  const_iterator end() const { return const_iterator(*this, size); };

  friend class ConstLatticePhiIterator<Lattice const>;
  typedef ConstLatticePhiIterator<Lattice const> const_phi_iterator;
  const_phi_iterator begin_phi() const { return const_phi_iterator(*this, 0); };
  const_phi_iterator end_phi() const { return const_phi_iterator(*this, size); };

  friend class LatticePhiIterator<Lattice>;
  typedef LatticePhiIterator<Lattice> phi_iterator;
  phi_iterator begin_phi() { return phi_iterator(*this, 0); };
  phi_iterator end_phi() { return phi_iterator(*this, size); };

  friend class ConstLatticePhiXYIterator<Lattice const>;
  typedef ConstLatticePhiXYIterator<Lattice const> const_phi_xy_iterator;
  const_phi_xy_iterator begin_phi_xy() const { return const_phi_xy_iterator(*this, 0); };
  const_phi_xy_iterator end_phi_xy() const { return const_phi_xy_iterator(*this, size); };

private:
  // The actual field data
  RR* phi_data;

  // Implementation of the Metropolis update
  inline void step_Metropolis(const PhiXY& phi);

  // Cluster membership for Wolff algorithm
  std::unordered_set<ZZ> cluster;
  std::vector<ZZ> cluster_queue;

  // Implementation of the Wolff update step
  inline bool cluster_try_add(bool sign, Phi phi, Site new_site);
  void flip_cluster();

public:
  // The two bare parameters of the Lagrangian
  const RR mu2;
  const RR lambda;

  // Rescaled parameters to simplify formulae
  const RR mu2_tilde;
  const RR lambda_tilde;

  // Number of lattice points in x and y direction
  const ZZ Nx;
  const ZZ Ny;

  // Number of sites in lattice
  const ZZ size;
  // short for size - Nx
  const ZZ size_minus_Nx;

  Lattice(const RR mu2_, const RR lambda_, ZZ Nx_, ZZ Ny_);
  ~Lattice();

  // Access to the field data
  std::vector<RR> get_phi() const 
  { return std::vector<RR>(phi_data, phi_data + size); };

  // Randomize the lattice variables
  void randomize(unsigned long int seed = 0);

  // Perform a single lattice iteration
  /* An iteration consists of a number of Metropolis sweeps plus one
   * Wolf cluster flip.
   */
  void iteration();
  
  // Equilibrate for N iterations
  void equilibrate(const ZZ N);

  // Metropolis algorithm step
  /* Randomly modifies a single lattice site
   */
  void sweep_Metropolis();

  // Wolff algorithm step
  /* Randomly flips the sign of multiple adjacent lattice sites (a
   * cluster). Returns cluster size.
   */
  ZZ step_Wolff(const Site& site);


  ////////////////////////////////////////////////////////
  //
  // Calculation of observables
  // 
  ////////////////////////////////////////////////////////
  
  // Energy contribution from a lattice site
  inline RR energy(const ConstPhiXY& phi) const;

  RR total_energy() const;
  RR average_phi() const;   
  RR average_abs_phi() const;   
};






#endif // LATTICE__HH
