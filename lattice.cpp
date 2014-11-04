#include <iostream>
#include <cmath>
#include <cassert>
#include <ctime>
#include <gsl/gsl_sf_exp.h>  // Exponential functions

#include "lattice.hpp"
#include "util.hpp"


Lattice::Lattice(RR mu2_, RR lambda_, ZZ Nx_, ZZ Ny_) 
  : mu2(mu2_), lambda(lambda_), 
    mu2_tilde(2 + (mu2_ / 2)), lambda_tilde(lambda_ / 4), 
    Nx(Nx_), Ny(Ny_), 
    size(Nx * Ny),
    size_minus_Nx(Nx * (Ny-1))
{
  phi_data = new RR[size];
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  randomize();
}


Lattice::~Lattice() 
{
  gsl_rng_free(rng);
  delete [] phi_data;
}


void Lattice::randomize(unsigned long int seed)
{
  if (seed != 0) 
    gsl_rng_set(rng, seed);
  for (phi_iterator pi = begin_phi(); pi != end_phi(); pi++)
    *pi = random_phi();
}


void Lattice::iteration()
{
  sweep_Metropolis();
  step_Wolff(random_site());
}


void Lattice::equilibrate(const ZZ N)
{
  const std::clock_t t0 = std::clock();
  for (ZZ i = 0; i < N; i++) 
    iteration();
  const std::clock_t t1 = std::clock();
  std::cout << "Equilibrate: " 
            << 1000 * (t1 - t0) / CLOCKS_PER_SEC 
            << "ms" << std::endl;
}



////////////////////////////////////////////////////////
//
// Accessing Lattice Sites
// 
////////////////////////////////////////////////////////

inline Lattice::SiteXY Lattice::get_site_xy(Site site) const
{
  ZZ prev_x, next_x, prev_y, next_y;
  if (likely(site < size_minus_Nx)) {
    next_x = site + 1;
    next_y = site + Nx;
  }
  else if (likely(site < size - 1)) {
    next_x = site + 1;
    next_y = site + Nx - size;
  }
  else {
    assert(site == size - 1);
    next_x = 0;
    next_y = Nx - 1;
  }

  if (likely(site >= Nx)) {
    prev_x = site - 1;
    prev_y = site - Nx;
  }
  else if (likely(site > 0)) {
    prev_x = site - 1;
    prev_y = site + size - Nx;
  }
  else {
    assert(site == 0);
    prev_x = size - 1;
    prev_y = size - Nx;
  }
  return SiteXY(site, prev_x, next_x, prev_y, next_y);
}


inline Lattice::Phi Lattice::get_phi(Site site)
{
  return phi_data[site];
}


inline Lattice::ConstPhi Lattice::get_phi(Site site) const
{
  return phi_data[site];
}


inline Lattice::PhiXY Lattice::get_phi_xy(Site site)
{
  SiteXY s = get_site_xy(site);
  return PhiXY(phi_data[s.site], 
               phi_data[s.prev_x], phi_data[s.next_x], 
               phi_data[s.prev_y], phi_data[s.next_y]);
}


inline Lattice::ConstPhiXY Lattice::get_phi_xy(Site site) const
{
  SiteXY s = get_site_xy(site);
  return ConstPhiXY(phi_data[s.site], 
                    phi_data[s.prev_x], phi_data[s.next_x], 
                    phi_data[s.prev_y], phi_data[s.next_y]);
}


inline Lattice::Site Lattice::random_site() const
{
  return Site(floor(size * random()));
}


////////////////////////////////////////////////////////
//
// Metropolis Algorithm
// 
////////////////////////////////////////////////////////

void noinline Lattice::sweep_Metropolis() {
  for (ZZ i = 0; i < METROPOLIS_SWEEPS_PER_ITERATION * size; i++) {
    step_Metropolis(get_phi_xy(random_site()));
  }
}


inline void Lattice::step_Metropolis(const PhiXY& phi) {
  const double phi_2 = phi.value * phi.value;
  const double newphi = phi.value + random_phi();
  const double newphi_2 = newphi * newphi;
  // Calculate energy difference
  const double delta = 
    (phi.value - newphi) * (phi.next_x + phi.next_y + phi.prev_x + phi.prev_y)
    +
    mu2_tilde * (newphi_2 - phi_2)
    + 
    lambda_tilde * (newphi_2 * newphi_2 - phi_2 * phi_2);
  // Flip if difference <= 0 or probabilistic acceptance
  if (delta <= 0 || random() < gsl_sf_exp(-delta))
    phi.value = newphi;
}


////////////////////////////////////////////////////////
//
// Wolff Algorithm
// 
////////////////////////////////////////////////////////

// Conditionally add new_site to the cluster.
/* Acceptance of new_site depends on the sign at the lattice site and
 * a probabilistic factor 
 */
inline bool Lattice::cluster_try_add(bool positive, Phi phi, Site new_site) {
  if ((get_phi(new_site) > 0) != positive)
    return false;
  if (cluster.find(new_site) != cluster.end())
    return false;
  const RR probability = 
    1 - gsl_sf_exp(-2 * phi * get_phi(new_site));
  if (random() < probability) {
    cluster.insert(new_site);
    cluster_queue.push_back(new_site);
    return true;
  }
  return false;
}


void Lattice::flip_cluster() {
  assert(cluster_queue.empty());
  for (auto ci = cluster.begin(); ci != cluster.end(); ci++) {
    get_phi(*ci) *= -1;
  }
  cluster.clear();
}


Lattice::ZZ noinline Lattice::step_Wolff(const Site& site) {
  const bool positive = (get_phi(site) > 0);
  cluster.clear();
  cluster.insert(site);
  cluster_queue.clear();
  cluster_queue.push_back(site);
  while (!cluster_queue.empty()) {
    Site site = cluster_queue.back();
    cluster_queue.pop_back();
    Phi phi = get_phi(site);
    SiteXY nbhd = get_site_xy(site);
    cluster_try_add(positive, phi, nbhd.prev_x);
    cluster_try_add(positive, phi, nbhd.next_x);
    cluster_try_add(positive, phi, nbhd.prev_y);
    cluster_try_add(positive, phi, nbhd.next_y);
  }
  const ZZ cluster_size = cluster.size();
  flip_cluster();
  return cluster_size;
}


////////////////////////////////////////////////////////
//
// Calculation of observables
// 
////////////////////////////////////////////////////////

inline Lattice::RR Lattice::energy(const ConstPhiXY& phi) const
{
  const RR phi_2 = phi.value * phi.value;
  return 
    - phi.value * (phi.next_x + phi.next_y)
    + mu2_tilde * phi_2
    + lambda_tilde * phi_2 * phi_2;
}


Lattice::RR Lattice::total_energy() const 
{
  RR result = 0;
  for (const_phi_xy_iterator phi = begin_phi_xy(); phi != end_phi_xy(); phi++)
    result += energy(*phi);
  return result / size;
}


Lattice::RR Lattice::average_phi() const
{
  RR sum = 0;
  for (const_phi_iterator phi_iter = begin_phi(); phi_iter != end_phi(); phi_iter++)
    sum += *phi_iter;
  return sum / size;
}

Lattice::RR Lattice::average_abs_phi() const
{
  RR sum = 0;
  for (const_phi_iterator phi_iter = begin_phi(); phi_iter != end_phi(); phi_iter++)
    sum += abs(*phi_iter);
  return sum / size;
}
