/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-09 22:48:45
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-22 00:00:22
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <algorithm>
#include "wavefunction.h"

namespace var {

void Wavefunction::bcs_init(void)
{
  // resizing
  mat_work.resize(block_dim_,block_dim_);
  mat_delta_k.resize(block_dim_,block_dim_);
  mat_dphi_k.resize(block_dim_,block_dim_);
  cphi_k.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) cphi_k[k].resize(block_dim_,block_dim_);
  bcs_large_number_ = 1.0E+2;
}

void Wavefunction::bcs_oneband(void)
{
  // BCS wave function for one band system 
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    double ek = std::real(mf_model_.quadratic_spinup_block()(0,0)); 
    auto delta_k = mf_model_.pairing_part()(0,0);
    mf_model_.construct_kspace_block(-kvec);
    ek += std::real(mf_model_.quadratic_spinup_block()(0,0));;
    delta_k += mf_model_.pairing_part()(0,0);
    delta_k *= 0.5;
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (deltak_sq<1.0E-12 && ek<0.0) {
      cphi_k[k](0,0) = bcs_large_number_ * std::exp(ii()*std::arg(delta_k));
    }
    else {
      cphi_k[k](0,0) = 2.0*delta_k/ek_plus_Ek;
    }
    //std::cout << ek << "\n";
    //std::cout << delta_k << "\n";
    //std::cout << k << " " << cphi_k[k](0,0) << "\n";
    //getchar();
  }
}

// Assumption: TR symmetry of the Hamiltonian
void Wavefunction::bcs_multiband(void)
{
  // BCS wave function for multiband system (INTRABAND pairing only)
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    hk.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    mat_delta_k = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    hminusk.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    mat_work = 0.5*(mat_delta_k + mf_model_.pairing_part().transpose());
    // transform pairing part
    mat_delta_k = hk.eigenvectors().adjoint() * mat_work * 
      hminusk.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    mat_dphi_k.setZero();
    for (unsigned i=0; i<block_dim_; ++i) {
      double ek = hk.eigenvalues()[i] + hminusk.eigenvalues()[i];
      double deltak_sq = std::norm(mat_delta_k(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      if (deltak_sq<1.0E-12 && ek<0.0) {
        mat_dphi_k(i,i) = bcs_large_number_ * std::exp(ii()*std::arg(mat_delta_k(i,i)));
      }
      else {
        mat_dphi_k(i,i) = 2.0*mat_delta_k(i,i)/ek_plus_Ek;
      }
    }
    // bcs ampitudes in original basis 
    for (unsigned i=0; i<block_dim_; ++i) 
      mat_work.col(i) = hk.eigenvectors().col(i) * mat_dphi_k(i,i);
    cphi_k[k] = mat_work * hminusk.eigenvectors().transpose();
    //std::cout << delta_k << "\n";
  } 
}

void Wavefunction::bcs_disordered(const lattice::LatticeGraph& graph)
{
  assert (cphi_k.size() == 1);  // only k=0 point
  // diagonalize quadratic part
  hk.compute(mf_model_.quadratic_up_matrix(graph));
  mf_model_.get_pairing_varparms(mat_delta_k); // assuming only diagonal elems
  mat_dphi_k.setZero();
  for (unsigned i=0; i<block_dim_; ++i) {
    double ek = hk.eigenvalues()[i] + hminusk.eigenvalues()[i];
    double deltak_sq = std::norm(mat_delta_k(i,i));
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (deltak_sq<1.0E-12 && ek<0.0) {
      mat_dphi_k(i,i) = bcs_large_number_ * std::exp(ii()*std::arg(mat_delta_k(i,i)));
    }
    else {
      mat_dphi_k(i,i) = 2.0*mat_delta_k(i,i)/ek_plus_Ek;
    }
  }
  // bcs ampitudes in original basis 
  for (unsigned i=0; i<block_dim_; ++i) 
    mat_work.col(i) = hk.eigenvectors().col(i) * mat_dphi_k(i,i);
  cphi_k[0] = mat_work * hk.eigenvectors().transpose();
  //std::cout << delta_k << "\n";
}

double Wavefunction::get_noninteracting_mu(void)
{
  std::vector<double> ek;
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    hk.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    ek.insert(ek.end(),hk.eigenvalues().data(),hk.eigenvalues().data()+hk.eigenvalues().size());
  }
  std::sort(ek.begin(),ek.end());
  //for (const auto& e : ek) std::cout << e << "\n";
  //std::cout << "spins = " << num_spins_ << "\n";
  if (num_upspins_ < num_sites_)
    return 0.5*(ek[num_upspins_-1]+ek[num_upspins_]);
  else
    return ek[num_upspins_-1];
}

} // end namespace var











