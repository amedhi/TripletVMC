/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-20 09:43:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 15:03:31
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <stdexcept>
#include "./groundstate.h"

namespace var {

void GroundState::update(const input::Parameters& inputs)
{
  throw std::runtime_error("GroundState::update_parameters: function must be overriden");
}

void GroundState::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  throw std::runtime_error("GroundState::update_parameters: function must be overriden");
}

void GroundState::get_wf_amplitudes(Matrix& psi)
{
  throw std::runtime_error("GroundState::get_wf_amplitudes: function must be overriden");
}

void GroundState::get_wf_gradient(std::vector<Matrix>& psi_gradient)
{
  throw std::runtime_error("GroundState::get_wf_gradients: function must be overriden");
}

void GroundState::set_particle_num(const input::Parameters& inputs)
{
  hole_doping_ = inputs.set_value("hole_doping", 0.0);
  if (std::abs(last_hole_doping_-hole_doping_)<1.0E-15) {
    // no change in hole doping, particle number remails same
    return;
  }
  last_hole_doping_ = hole_doping_;
  band_filling_ = 1.0-hole_doping_;
  int num_sites = static_cast<int>(num_sites_);
  if (pairing_type_) {
    int n = static_cast<int>(std::round(0.5*band_filling_*num_sites));
    if (n<0 || n>num_sites) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_upspins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_upspins_;
    num_spins_ = num_upspins_ + num_dnspins_;
    band_filling_ = static_cast<double>(2*n)/num_sites;
  }
  else{
    int n = static_cast<int>(std::round(band_filling_*num_sites));
    if (n<0 || n>2*num_sites) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_spins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_ - num_dnspins_;
    band_filling_ = static_cast<double>(n)/num_sites;
  }
  hole_doping_ = 1.0 - band_filling_;
}

double GroundState::get_noninteracting_mu(void)
{
  std::vector<double> ek;
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    ek.insert(ek.end(),es_k_up.eigenvalues().data(),
      es_k_up.eigenvalues().data()+kblock_dim_);
  }
  std::sort(ek.begin(),ek.end());
  //for (const auto& e : ek) std::cout << e << "\n";
  //double e = 0.0;
  //for (unsigned i=0; i<num_upspins_; ++i) e += ek[i];
  //std::cout << "energy = " << 2*e/num_sites_ << "\n";
  //std::cout << "upspins = " << num_upspins_ << "\n";
  if (num_upspins_ < num_sites_) {
    //std::cout << 0.5*(ek[num_upspins_-1]+ek[num_upspins_]) << "\n";
    return 0.5*(ek[num_upspins_-1]+ek[num_upspins_]);
  }
  else
    return ek[num_upspins_-1];
}

void GroundState::set_ft_matrix(const lattice::LatticeGraph& graph)
{
  // matrix for transformation from site-basis to k-basis
  FTU_.resize(num_kpoints_,num_kpoints_);
  double one_by_sqrt_nk = 1.0/std::sqrt(static_cast<double>(num_kpoints_));
  unsigned i = 0;
  for (unsigned n=0; n<num_kpoints_; ++n) {
    auto Ri = graph.site_cellcord(i);
    //std::cout << Ri << "\n"; getchar();
    for (unsigned k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      FTU_(n,k) = std::exp(ii()*kvec.dot(Ri)) * one_by_sqrt_nk;
    }
    i += kblock_dim_;
    // i is first basis site in next unitcell
  }
  //std::cout << FTU_ << "\n"; 
  //getchar();
}


} // end namespace var
