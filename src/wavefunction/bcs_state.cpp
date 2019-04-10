/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 23:06:41
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 11:32:25
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./bcs_state.h"

namespace var {

BCS_State::BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph) 
  : GroundState(true)
{
  init(order_type, inputs, graph);
}

int BCS_State::init(const bcs& order_type, const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph)
{
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 5.0E+4;

  // build MF Hamiltonian
  order_type_ = order_type;
  varparms_.clear();
  mf_model_.init(graph.lattice());
  std::string name;
  double defval, lb, ub;
  using namespace model;
  model::CouplingConstant cc;
  if (order_type_==bcs::swave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else if (order_type_==bcs::dwave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
    mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else {
    throw std::range_error("BCS_State::BCS_State: unidefined bcs order");
  }
  // chemical potential
  int info;
  mf_model_.add_parameter(name="mu", defval=0.0, inputs, info);
  if (info == 0) noninteracting_mu_ = false;
  else noninteracting_mu_ = true;
  if (inputs.set_value("mu_variational", false, info)) 
    varparms_.add("mu", defval=0.0, -2.0, +1.0);
  // finalize MF Hamiltonian
  mf_model_.finalize(graph);
  num_varparms_ = varparms_.size();

  // bloch basis
  blochbasis_.construct(graph);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(graph);
  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  delta_k_.resize(kblock_dim_,kblock_dim_);
  dphi_k_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  work_k_.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
    work_k_[k].resize(kblock_dim_,kblock_dim_);
  } 
  return 0;
}

void BCS_State::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
  // update MF model
  mf_model_.update(inputs);
  // update variational parameters
  for (auto& p : varparms_) 
    p.change_value(mf_model_.get_parameter_value(p.name()));
  // chemical potential
  if (noninteracting_mu_) {
    // the next line is 'really' needed 
    mf_model_.update_site_parameter("mu", 0.0);
    double mu_0 = get_noninteracting_mu();
    //std::cout << "mu = " << mu_0 << "\n";
    mf_model_.update_site_parameter("mu", mu_0);
  }
  // check MF energy
  //std::cout << " MF energy = " << get_mf_energy() << "\n"; 
}

void BCS_State::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // update from variational parameters
  unsigned i = 0;
  for (auto& p : varparms_) {
    auto x = pvector[start_pos+i];
    p.change_value(x);
    mf_model_.update_parameter(p.name(), x);
    ++i;
  }
  mf_model_.update_terms();
}

void BCS_State::get_wf_amplitudes(Matrix& psi) 
{
  // k-space pair amplitudes
  if (kblock_dim_==1) {
    get_pair_amplitudes_oneband(phi_k_);
  }
  else {
    get_pair_amplitudes_multiband(phi_k_);
  }
  // 'lattice space' pair amplitudes
  get_pair_amplitudes_sitebasis(phi_k_, psi);
}

void BCS_State::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  unsigned i=0; 
  for (const auto& p : varparms_) {
    double h = p.diff_h();
    double inv_2h = 0.5/h;
    double x = p.value();
    mf_model_.update_parameter(p.name(), x+h);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(phi_k_);
    else get_pair_amplitudes_multiband(phi_k_);
    mf_model_.update_parameter(p.name(), x-h);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(work_k_);
    else get_pair_amplitudes_multiband(work_k_);
    // model to original state
    mf_model_.update_parameter(p.name(), x);
    mf_model_.update_terms();
    for (unsigned k=0; k<num_kpoints_; ++k) {
      phi_k_[k] -= work_k_[k];
      phi_k_[k] *= inv_2h;
    }
    //std::cout << phi_k_[0] << "\n"; getchar();
    // phi_k_ is now the derivative wrt i-th parameter
    // wave function gradients
    get_pair_amplitudes_sitebasis(phi_k_, psi_gradient[i]);
    ++i;
  }
}

void BCS_State::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  unsigned p = 0;
  for (unsigned i=0; i<num_kpoints_; ++i) {
    unsigned q = 0;
    for (unsigned j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (unsigned k=0; k<num_kpoints_; ++k) {
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // copy transformed block
      //psi.block(p,q,kblock_dim_,kblock_dim_) = 
      for (unsigned m=0; m<kblock_dim_; ++m) {
        for (unsigned n=0; n<kblock_dim_; ++n) 
          psi(p+m,q+n) = ampl_part(work_(m,n));
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  /*
  for (unsigned i=0; i<num_sites_; ++i) {
    for (unsigned j=0; j<num_sites_; ++j) {
      std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      getchar();
    }
  }*/

}

void BCS_State::get_pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for one band system 
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    double ek = std::real(mf_model_.quadratic_spinup_block()(0,0)); 
    auto delta_k = mf_model_.pairing_part()(0,0);
    mf_model_.construct_kspace_block(-kvec);
    ek += std::real(mf_model_.quadratic_spinup_block()(0,0));;
    delta_k += mf_model_.pairing_part()(0,0);
    delta_k *= 0.5;
    //----------------------------------
    //std::cout << "** Hack at BCS_State\n";
    //ek = -4.0*(std::cos(kvec[0])+std::cos(kvec[1]));
    //delta_k = 0.5 * (std::cos(kvec[0])-std::cos(kvec[1])); 
    //----------------------------------
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
      phi_k[k](0,0) = large_number_ * std::exp(ii()*std::arg(delta_k));
    }
    else {
      phi_k[k](0,0) = 2.0*delta_k/ek_plus_Ek;
    }
    //std::cout << "k = " << k << "\n";
    //std::cout << "ek = " << 0.5*ek << "\n";
    //std::cout << "delta_k = " << delta_k << "\n";
    //std::cout << "phi_k = " << phi_k[k](0,0) << "\n";
    //getchar();
  }
}

void BCS_State::get_pair_amplitudes_multiband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for multi-band system (INTRABAND pairing only)
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    dphi_k_.setZero();
    for (unsigned i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i] + es_minusk_up.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
        dphi_k_(i,i) = large_number_ * std::exp(ii()*std::arg(delta_k_(i,i)));
      }
      else {
        dphi_k_(i,i) = 2.0*delta_k_(i,i)/ek_plus_Ek;
      }
    }
    // bcs ampitudes in original basis 
    for (unsigned i=0; i<kblock_dim_; ++i) 
      work_.col(i) = es_k_up.eigenvectors().col(i) * dphi_k_(i,i);
    phi_k[k] = work_ * es_minusk_up.eigenvectors().transpose();
    //std::cout << delta_k << "\n";
  } 
}

double BCS_State::get_mf_energy(void)
{
  double mf_energy = 0.0;
  double delta = mf_model_.get_parameter_value("delta_sc");
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    /*
    double ek = -2.0*(std::cos(kvec[0])+std::cos(kvec[1]));
    double deltak = delta * (std::cos(kvec[0])-std::cos(kvec[1]));
    double Ek = std::sqrt(ek*ek + deltak*deltak);
    mf_energy += (ek - Ek);
    */
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    //std::cout << work_ << "\n"; getchar();
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    for (unsigned i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double Ek = std::sqrt(ek*ek + deltak_sq);
      //std::cout << ek << " " << Ek << "\n"; getchar();
      mf_energy += (ek - Ek);
    }
  }
  // constant term
  double J = 0.35;  // tJ model
  mf_energy += num_bonds_*delta*delta/(2.0*J);  
  //std::cout << delta << " ";
  return mf_energy/num_sites_;
}



} // end namespace var
