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
  num_states_ = 2*num_sites_;
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
    mf_model_.add_bondterm(name="pairing", cc, op::create_singlet());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else if (order_type_==bcs::pwave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_00", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_uu", defval=0.0, inputs);
    mf_model_.add_parameter(name="delta_ud", defval=0.0, inputs);
    mf_model_.add_parameter(name="delta_dd", defval=0.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    cc = CouplingConstant({0,"delta_00"}, {1,"-delta_00"});
    mf_model_.add_bondterm(name="pairing", cc, op::create_singlet());
    mf_model_.add_bondterm(name="pairing", cc="delta_uu", op::create_triplet_uu());
    mf_model_.add_bondterm(name="pairing", cc="delta_ud", op::create_triplet_ud());
    mf_model_.add_bondterm(name="pairing", cc="delta_dd", op::create_triplet_dd());
    // variational parameters
    varparms_.add("delta_00", defval=1.0, lb=0.0, ub=2.0);
    varparms_.add("delta_uu", defval=1.0, lb=0.0, ub=2.0);
    varparms_.add("delta_ud", defval=1.0, lb=0.0, ub=2.0);
    varparms_.add("delta_dd", defval=1.0, lb=0.0, ub=2.0);
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
  dim_ = kblock_dim_;
  dim2_ = 2*dim_;
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(graph);
  // work arrays
  bdg_mat_.resize(4,4);
  uk_mat_.resize(2,2);
  vk_mat_.resize(2,2);
  work_.resize(dim2_,dim2_);
  delta_k_.resize(dim2_,dim2_);
  dphi_k_.resize(dim2_,dim2_);
  phi_k_.resize(num_kpoints_);
  work_k_.resize(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(dim2_,dim2_);
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
  int i = 0;
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
  if (dim_==1) {
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
  int p = 0;
  for (int i=0; i<num_kpoints_; ++i) {
    int q = 0;
    for (int j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (int k=0; k<num_kpoints_; ++k) {
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // copy transformed block
      // UP-UP block
      for (int m=0; m<dim_; ++m) {
        for (int n=0; n<dim_; ++n) 
          psi(p+m,q+n) = ampl_part(work_(2*m,2*n));
      }
      // UP-DN block
      for (int m=0; m<dim_; ++m) {
        for (int n=0; n<dim_; ++n) 
          psi(num_sites_+p+m,q+n) = ampl_part(work_(2*m,2*n+1));
      }
      // UP-DN block
      for (int m=0; m<dim_; ++m) {
        for (int n=0; n<dim_; ++n) 
          psi(p+m,num_sites_+q+n) = ampl_part(work_(2*m+1,2*n));
      }
      // DN-DN block
      for (int m=0; m<dim_; ++m) {
        for (int n=0; n<dim_; ++n) 
          psi(num_sites_+p+m,num_sites_+q+n) = ampl_part(work_(2*m+1,2*n+1));
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  /* 
  for (int i=0; i<num_states_; ++i) {
    for (int j=0; j<num_states_; ++j) {
      std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      getchar();
    }
  }
  getchar();*/
}

void BCS_State::get_pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for one-band system 
  dphi_k_.setZero();
  RealVector eps_k(2);
  RealVector Ek(2);
  RealVector Eqp(2);
  RealVector Dk(2);
  ComplexMatrix Delta_sq(2,2);
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //std::cout << "kvec = " << kvec.transpose() << "\n";
    //-------------'+k block'-------------
    mf_model_.construct_kspace_block(kvec);
    delta_k_ = mf_model_.pairing_part();
    es_k_up.compute(mf_model_.quadratic_block());

    /* 
      Analytical expression for uk & vk (See Vollhardt, Woelfle, 
      "The superfluid phases of Helium 2", p69)
    */
    Delta_sq = delta_k_*delta_k_.adjoint();
    for (int i=0; i<2; ++i) {
      eps_k(i) = es_k_up.eigenvalues()(i);
      Ek(i) = std::sqrt(eps_k(i)*eps_k(i)+std::real(Delta_sq(i,i)));
      Eqp(i) = eps_k(i) + Ek(i);
      Dk(i) = std::sqrt(2.0*Ek(i)*(eps_k(i)+Ek(i)));
    }
    uk_mat_(0,0) = Eqp(0)/Dk(0);
    uk_mat_(1,1) = Eqp(1)/Dk(1);
    uk_mat_(0,1) = 0.0;
    uk_mat_(1,0) = 0.0;
    /*
      In the 'vk' expression, '-' sign as in above reference not needed
      as here we define the paring term with a '-' sign in front.
    */
    vk_mat_(0,0) = delta_k_(0,0)/Dk(0);  
    vk_mat_(0,1) = delta_k_(0,1)/Dk(0);
    vk_mat_(1,0) = delta_k_(1,0)/Dk(1);
    vk_mat_(1,1) = delta_k_(1,1)/Dk(1);

    // Pair amplitudes 'dphi(k)' matrix.
    for (int i=0; i<2; ++i) {
      dphi_k_(i,0) = vk_mat_(i,0)/uk_mat_(i,i);
      dphi_k_(i,1) = vk_mat_(i,1)/uk_mat_(i,i);
      /*
      if (std::abs(uk_mat_(i,i))<1.0E-12) {
        dphi_k_(i,0) = large_number_*std::arg(vk_mat_(i,0));
        dphi_k_(i,1) = large_number_*std::arg(vk_mat_(i,1));
      }
      else {
        dphi_k_(i,0) = vk_mat_(i,0)/uk_mat_(i,i);
        dphi_k_(i,1) = vk_mat_(i,1)/uk_mat_(i,i);
      }*/
    }
    // Pair amplitudes 'phi(k)' matrix.
    phi_k[k] = dphi_k_;
  } 
/*
    // bdg matrix
    int n = 0;
    for (int i=0; i<dim_; ++i) {
      bdg_mat_.setZero();
      bdg_mat_(0,0) = 1.0*es_k_up.eigenvalues()[i];
      bdg_mat_(1,1) = 1.0*es_k_up.eigenvalues()[i+dim_];
      bdg_mat_(2,2) = -1.0*es_minusk_up.eigenvalues()[i];
      bdg_mat_(3,3) = -1.0*es_minusk_up.eigenvalues()[i+dim_];
      // upper diagonal non-zero elements
      //bdg_mat_(0,2) = -delta_k_(i,i);
      bdg_mat_(2,0) = -std::conj(delta_k_(i,i));
      //bdg_mat_(0,3) = -delta_k_(i,i+dim_);
      bdg_mat_(3,0) = -std::conj(delta_k_(i,i+dim_));
      //bdg_mat_(1,2) = -delta_k_(i+dim_,i);
      bdg_mat_(2,1) = -std::conj(delta_k_(i+dim_,i));
      //bdg_mat_(1,3) = -delta_k_(i+dim_,i+dim_);
      bdg_mat_(3,1) = -std::conj(delta_k_(i+dim_,i+dim_));
      //std::cout << "bdg_mat = \n" << bdg_mat_ << "\n"; getchar();
      // diagonalize
      es_bdg_.compute(bdg_mat_);
      // uk & vk matrix. By definition
      //           |uk    -vk|
      //  U^\dag = |         |
      //           |vk^* uk^*|
      //
      //std::cout << "U = \n" << es_bdg_.eigenvectors() << "\n"; getchar();
      auto eigvec1 = es_bdg_.eigenvectors().col(3).conjugate();
      auto eigvec2 = es_bdg_.eigenvectors().col(2).conjugate();
      if (std::abs(eigvec1(0))>= std::abs(eigvec2(0))) {
        uk_mat_(0,0) = eigvec1(0);
        uk_mat_(0,1) = eigvec2(0);
        uk_mat_(1,0) = eigvec1(1);
        uk_mat_(1,1) = eigvec2(1);

        vk_mat_(0,0) = -eigvec1(2);
        vk_mat_(0,1) = -eigvec1(3);
        vk_mat_(1,0) = -eigvec2(2);
        vk_mat_(1,1) = -eigvec2(3);
      }
      else {
        uk_mat_(0,0) = eigvec2(0);
        uk_mat_(0,1) = eigvec1(0);
        uk_mat_(1,0) = eigvec2(1);
        uk_mat_(1,1) = eigvec1(1);

        vk_mat_(0,0) = -eigvec2(2);
        vk_mat_(0,1) = -eigvec2(3);
        vk_mat_(1,0) = -eigvec1(2);
        vk_mat_(1,1) = -eigvec1(3);
      }
      //bdg_mat_ = es_bdg_.eigenvectors(); //.adjoint();
      //std::cout << es_bdg_.eigenvalues().transpose() << "\n"; 
      //std::cout << "U = \n" << bdg_mat_ << "\n"; getchar();
      //uk_mat_ = bdg_mat_.block(0,0,2,2); 
      //vk_mat_ = -bdg_mat_.block(0,2,2,2); 
      //uk_mat_(0,0) = std::conj(bdg_mat_(0,3));
      //uk_mat_(0,1) = std::conj(bdg_mat_(0,2));
      //uk_mat_(1,0) = std::conj(bdg_mat_(2,3));
      //uk_mat_(1,1) = std::conj(bdg_mat_(1,2));
      //vk_mat_(0,0) = std::conj(bdg_mat_(2,3));
      //vk_mat_(0,1) = std::conj(bdg_mat_(2,2));
      //vk_mat_(1,0) = std::conj(bdg_mat_(3,3));
      //vk_mat_(1,1) = std::conj(bdg_mat_(3,2));
      //
      // Pair amplitudes 'dphi(k)' matrix.
      for (int j=0; j<2; ++j) {
        if (std::abs(uk_mat_(j,j))<1.0E-12) {
          dphi_k_(n+j,n) = large_number_ * std::arg(vk_mat_(j,0));
          dphi_k_(n+j,n+1) = large_number_ * std::arg(vk_mat_(j,1));
        }
        else {
          dphi_k_(n+j,n) = vk_mat_(j,0)/uk_mat_(j,j);
          dphi_k_(n+j,n+1) = vk_mat_(j,1)/uk_mat_(j,j);
        }
      }
      // offset index for the next band
      n += 2;
    }
    std::cout << "dphi_k = \n" << dphi_k_ << "\n\n\n"; getchar();
    //
    // Pair ampitudes in original 'ck' basis 
    // (yes one-band system, but the following still might be needed)
    //
    phi_k[k] = es_k_up.eigenvectors() * dphi_k_ * 
               es_minusk_up.eigenvectors().transpose();
    //std::cout << "phi_k = \n" << phi_k[k] << "\n"; getchar();
  } 
*/
}

void BCS_State::get_pair_amplitudes_multiband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for multi-band system (INTRABAND pairing only)
  bdg_mat_.setZero();
  dphi_k_.setZero();
  RealVector eps_k(2);
  RealVector Ek(2);
  RealVector Eqp(2);
  RealVector Dk(2);
  ComplexMatrix Delta_k(2,2);
  ComplexMatrix Delta_sq(2,2);
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_block());
    work_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_block());
    // pairing part in band-operator basis
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();

    //------Assuming INTRABAND pairing only----------
    //  Analytical expression for uk & vk (See comments in 'oneband' case).
    int band = 0;
    for (int n=0; n<dim_; ++n) {
      eps_k(0) = es_k_up.eigenvalues()[n];
      eps_k(1) = es_k_up.eigenvalues()[n+dim_];
      Delta_k(0,0) = delta_k_(n,n);
      Delta_k(0,1) = delta_k_(n,n+dim_);
      Delta_k(1,0) = delta_k_(n+dim_,n);
      Delta_k(1,1) = delta_k_(n+dim_,n+dim_);
      Delta_sq = Delta_k*Delta_k.adjoint();
      for (int i=0; i<2; ++i) {
        Ek(i) = std::sqrt(eps_k(i)*eps_k(i)+std::real(Delta_sq(i,i)));
        Eqp(i) = eps_k(i) + Ek(i);
        Dk(i) = std::sqrt(2.0*Ek(i)*(eps_k(i)+Ek(i)));
      }
      uk_mat_(0,0) = Eqp(0)/Dk(0);
      uk_mat_(1,1) = Eqp(1)/Dk(1);
      uk_mat_(0,1) = 0.0;
      uk_mat_(1,0) = 0.0;
      vk_mat_(0,0) = delta_k_(0,0)/Dk(0);  
      vk_mat_(0,1) = delta_k_(0,1)/Dk(0);
      vk_mat_(1,0) = delta_k_(1,0)/Dk(1);
      vk_mat_(1,1) = delta_k_(1,1)/Dk(1);
      // Pair amplitudes 'dphi(k)' matrix.
      for (int i=0; i<2; ++i) {
        dphi_k_(band+i,band) = vk_mat_(i,0)/uk_mat_(i,i);
        dphi_k_(band+i,band+1) = vk_mat_(i,1)/uk_mat_(i,i);
        /*
        if (std::abs(uk_mat_(i,i))<1.0E-12) {
          dphi_k_(band+i,band) = large_number_*std::arg(vk_mat_(i,0));
          dphi_k_(band+i,band+1) = large_number_*std::arg(vk_mat_(i,1));
        }
        else {
          dphi_k_(band+i,band) = vk_mat_(i,0)/uk_mat_(i,i);
          dphi_k_(band+i,band+1) = vk_mat_(i,1)/uk_mat_(i,i);
        }*/
      }
      // offset index for the next band
      band += 2;
    }
    //std::cout << delta_k << "\n";
    // Pair ampitudes in original 'ck' basis 
    phi_k[k] = es_k_up.eigenvectors() * dphi_k_ * 
               es_minusk_up.eigenvectors().transpose();
  }
  /*
    // bdg matrix
    int n = 0;
    for (int i=0; i<dim_; ++i) {
      bdg_mat_(0,0) = 0.5*es_k_up.eigenvalues()[i];
      bdg_mat_(1,1) = 0.5*es_k_up.eigenvalues()[i+dim_];
      bdg_mat_(2,2) = -0.5*es_minusk_up.eigenvalues()[i];
      bdg_mat_(3,3) = -0.5*es_minusk_up.eigenvalues()[i+dim_];
      // upper diagonal non-zero elements
      bdg_mat_(0,2) = -delta_k_(i,i);
      bdg_mat_(0,3) = -delta_k_(i,i+dim_);
      bdg_mat_(1,2) = -delta_k_(i+dim_,i);
      bdg_mat_(1,3) = -delta_k_(i+dim_,i+dim_);
      // diagonalize
      es_bdg_.compute(bdg_mat_);
      // uk & vk matrix. By definition
      //           |uk    -vk|
      //  U^\dag = |         |
      //           |vk^* uk^*|
      //
      bdg_mat_ = es_bdg_.eigenvectors().adjoint();
      uk_mat_ = bdg_mat_.block(0,0,2,2); 
      vk_mat_ = -bdg_mat_.block(0,2,2,2); 

      // Pair amplitudes 'dphi(k)' matrix.
      //  For INTRABAND pairing only, it is block diagonal.
      //
      for (int j=0; j<2; ++j) {
        if (std::abs(uk_mat_(j,j))<1.0E-12) {
          auto large_ampl = large_number_ * std::arg(vk_mat_(j,j));
          dphi_k_(n+j,n) = large_ampl;
          dphi_k_(n+j,n+1) = large_ampl;
        }
        else {
          dphi_k_(n+j,n) = vk_mat_(j,0)/uk_mat_(j,j);
          dphi_k_(n+j,n+1) = vk_mat_(j,1)/uk_mat_(j,j);
        }
      }
      // offset index for the next band
      n += 2;
    }
    //std::cout << delta_k << "\n";
    // Pair ampitudes in original 'ck' basis 
    phi_k[k] = es_k_up.eigenvectors() * dphi_k_ * 
               es_minusk_up.eigenvectors().transpose();
  } 
  */
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
