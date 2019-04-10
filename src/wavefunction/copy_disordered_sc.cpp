/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:46:55
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-29 21:12:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./disordered_sc.h"
#include <fstream>

namespace var {

DisorderedSC::DisorderedSC(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
  : GroundState(true)
{
  init(inputs, graph);
}

int DisorderedSC::init(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 1.0E+2;

  // the hamiltonian to be stored in Sparse matrix formats
  quadratic_ham_.resize(num_sites_,num_sites_);
  pairing_ham_.resize(num_sites_,num_sites_);
  quadratic_coeffs_.resize(num_sites_+num_bonds_);
  pairing_coeffs_.resize(num_bonds_);
  delta_.resize(num_sites_,num_sites_);
  phi_k_.resize(num_sites_);
  work_.resize(num_sites_,num_sites_);
  psi_work_.resize(num_sites_,num_sites_);
  //psi_work2_.resize(num_sites_,num_sites_);

  // variational parameters
  varparms_.clear();
  std::string name;
  double defval, lb, ub;
  // site potentials
  mu_start_ = varparms_.size();
  for (unsigned i=0; i<num_sites_; ++i) {
    name = "mu_" + std::to_string(i);
    //varparms_.add(name, defval=0.0, lb=-4.0, ub=4.0);
    varparms_.add(name, defval=-0.50, lb=-2.0, ub=2.0);
  }
  t_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "t_" + std::to_string(i);
    //varparms_.add(name, defval=-1.0, lb=-2.0, ub=2.0);
    varparms_.add(name, defval=1.0, lb=0.5, ub=1.5);
  }
  delta_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "delta_" + std::to_string(i);
    //varparms_.add(name, defval=0.5, lb=-4.0, ub=4.0);
    varparms_.add(name, defval=0.4, lb=-1.0, ub=1.0);
  }
  num_varparms_ = varparms_.size();

  // Hamiltonian coefficients
  // mu-values
  for (unsigned i=0; i<num_sites_; ++i) {
    // mu_term:  -mu\sum_n_{i\sigma}
    double mu = varparms_[mu_start_+i].value();
    quadratic_coeffs_[i] = EigenTriplet(i,i,-mu);
  }
  // t & delta values (only lower triangular elements)
  unsigned i = num_sites_;
  unsigned j = 0;
  unsigned it = t_start_;
  unsigned id = delta_start_;
  std::cout << "** DisorderedSC::init: ---TEMPORARY HACK-----\n";
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    // t-values
    double t = varparms_[it++].value();
    double delta = varparms_[id++].value();
    if (graph.bond_type(b)==1) delta = -delta;
    if (site_i >= site_j) {
      quadratic_coeffs_[i++] = EigenTriplet(site_i,site_j,-t);
      pairing_coeffs_[j++] = EigenTriplet(site_i,site_j,-delta);
    }
    else {
      quadratic_coeffs_[i++] = EigenTriplet(site_j,site_i,-t);
      pairing_coeffs_[j++] = EigenTriplet(site_j,site_i,-delta);
    }
  }
  /*for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    double t = varparms_[it++].value();
    quadratic_coeffs_[i++] = EigenTriplet(site_i,site_j,t);
    quadratic_coeffs_[i++] = EigenTriplet(site_j,site_i,t);
    double delta = varparms_[id++].value();
    pairing_coeffs_[j++] = EigenTriplet(site_i,site_j,delta);
    pairing_coeffs_[j++] = EigenTriplet(site_j,site_i,delta);
  }*/
  // set the sparse matrices
  //quadratic_ham_.setFromTriplets(quadratic_coeffs_.begin(), quadratic_coeffs_.end());
  //pairing_ham_.setFromTriplets(pairing_coeffs_.begin(), pairing_coeffs_.end());
  quadratic_ham_.setZero();
  pairing_ham_.setZero();
  for (auto& coeff : quadratic_coeffs_) {
    quadratic_ham_(coeff.row(), coeff.col()) = coeff.value();
    //quadratic_ham_(coeff.col(), coeff.row()) = coeff.value();
  }
  for (auto& coeff : pairing_coeffs_) {
    pairing_ham_(coeff.row(), coeff.col()) = coeff.value();
    pairing_ham_(coeff.col(), coeff.row()) = coeff.value();
  }

  // check whether sparse matrix diagonalization is okay
  /*delta_.setZero();
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int i = graph.source(b);
    int j = graph.target(b);
    delta_(i,j) = -1.0;
    delta_(j,i) = -1.0;
  }
  es_quad_.compute(quadratic_ham_);
  es_pair_.compute(delta_);
  std::cout << es_quad_.eigenvectors() - es_pair_.eigenvectors() << "\n\n";
  getchar();*/

  // HACK
  num_varparms_ = 1;
  varparms_.resize(1);

  return 0;
}

void DisorderedSC::update(const input::Parameters& inputs)
{
  // HACK
  double x = inputs.set_value("varp", 1.0);
  for (unsigned i=0; i<1; i++) {
    auto coeff = quadratic_coeffs_[i];
    quadratic_ham_(coeff.row(),coeff.col()) = -x;
    //quadratic_ham_(coeff.col(),coeff.row()) = -x;
    //std::cout << coeff.row() << " " << coeff.col() << -x << "\n";
  }
  //std::cout << "x="<<quadratic_ham_.coeffRef(coeff.row(),coeff.col()) << "\n";
  //std::ofstream of("res_varp.txt", std::fstream::app);
  //of << x << "\n";
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
}

void DisorderedSC::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // update matrix elements
  // '-' sign because H = -\sum t_{ij} cicj - mu\sum_ni -sum_{ij}Delta_{ij}cicj etc
  unsigned i = 0;
  for (auto& coeff : quadratic_coeffs_) {
    coeff = EigenTriplet(coeff.row(), coeff.col(), -pvector[start_pos+i]);
    ++i;
  }
  for (auto& coeff : pairing_coeffs_) {
    coeff = EigenTriplet(coeff.row(), coeff.col(), -pvector[start_pos+i]);
    ++i;
  }
  // new Hamiltonian matrices
  //quadratic_ham_.setFromTriplets(quadratic_coeffs_.begin(), quadratic_coeffs_.end());
  //pairing_ham_.setFromTriplets(pairing_coeffs_.begin(), pairing_coeffs_.end());
  for (auto& coeff : quadratic_coeffs_) {
    quadratic_ham_(coeff.row(), coeff.col()) = coeff.value();
    //quadratic_ham_(coeff.col(), coeff.row()) = coeff.value();
  }
  for (auto& coeff : pairing_coeffs_) {
    pairing_ham_(coeff.row(), coeff.col()) = coeff.value();
    pairing_ham_(coeff.col(), coeff.row()) = coeff.value();
  }
}

void DisorderedSC::get_wf_amplitudes(Matrix& psi) 
{
  // diagonalize the quadratic part
  es_quad_.compute(quadratic_ham_);
  // work = \Delta + Transpose(Delta)
  //work_ = pairing_ham_ + Eigen::SparseMatrix<double>(pairing_ham_.transpose());
  //std::cout << work_ << "\n"; getchar();
  // transform the pairing part
  delta_ = es_quad_.eigenvectors().adjoint() * pairing_ham_ * 
      es_quad_.eigenvectors().conjugate();
  //delta_ = es_quad_.eigenvectors().transpose() * work_ * 
  //    es_quad_.eigenvectors();
  // wf amplitudes
  get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
    delta_, psi);

  //std::cout << psi_work_-psi  << "\n"; getchar();
  //psi_work_ = psi;
  //std::cout << es_quad_.eigenvalues() << "\n"; getchar();
  //std::cout << delta_.diagonal() << "\n"; getchar();
}

void DisorderedSC::hack_gradient(std::vector<Matrix>& psi_gradient) 
{
  //ComplexMatrix tmp(num_sites_,num_sites_);
  unsigned i = 0;
  double h = 1.0E-1;// varparms_[0].diff_h();
  double inv_2h = 0.5/h;
  auto coeff = quadratic_coeffs_[i];
  double x = -quadratic_ham_(coeff.row(),coeff.col());
  std::cout << "x = " << x << "\n";
  std::cout << "h = " << h << "\n";
  std::cout << "1/2h = " << inv_2h << "\n";
  quadratic_ham_(coeff.row(),coeff.col()) = -(x+h);
  //quadratic_ham_(coeff.col(),coeff.row()) = -(x+h);
  //std::cout << quadratic_ham_ << "\n"; getchar();
  get_wf_amplitudes(psi_gradient[i]);
  quadratic_ham_(coeff.row(),coeff.col()) = -(x-h);
  //quadratic_ham_(coeff.col(),coeff.row()) = -(x-h);
  //std::cout << quadratic_ham_ << "\n"; getchar();
  get_wf_amplitudes(psi_work_);
  // restore the hamtonian
  quadratic_ham_(coeff.row(),coeff.col()) = -x;
  //quadratic_ham_(coeff.col(),coeff.row()) = -x;
  // gradient
  psi_gradient[i] -= psi_work_;
  //std::cout << psi_gradient[i] << "\n"; getchar();
  //std::cout << "-----hi-------" << "\n"; getchar();
  psi_gradient[i] *= inv_2h;
}

void DisorderedSC::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  return hack_gradient(psi_gradient);
  std::cout << "---NO GO HERE-------" << "\n";

/*
  // derivative wrt quadratic parameters
  unsigned i = 0;
  for (auto& coeff : quadratic_coeffs_) {
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    quadratic_ham_.coeffRef(coeff.row(),coeff.col()) = coeff.value()+h;
    get_wf_amplitudes(psi_gradient[i]);
    quadratic_ham_.coeffRef(coeff.row(),coeff.col()) = coeff.value()-h;
    get_wf_amplitudes(psi_work_);
    // restore the hamtonian
    quadratic_ham_.coeffRef(coeff.row(),coeff.col()) = coeff.value();
    // gradient
    psi_gradient[i] -= psi_work_;
    psi_gradient[i] *= inv_2h;
    ++i;
  }
  // derivative wrt pairing parameters
  // quadratic part remains constant
  es_quad_.compute(quadratic_ham_);
  // paring matrix elements
  work_ = pairing_ham_ + Eigen::SparseMatrix<double>(pairing_ham_.transpose());
  for (auto& coeff : pairing_coeffs_) {
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    // wavefunction 'at (x+h)'
    work_.coeffRef(coeff.row(),coeff.col()) = coeff.value()+h;
    work_.coeffRef(coeff.col(),coeff.row()) = coeff.value()+h;
    //delta_ = es_quad_.eigenvectors().adjoint() * work_ * 
    //  es_quad_.eigenvectors().conjugate();
    delta_ = es_quad_.eigenvectors().transpose() * work_ * 
      es_quad_.eigenvectors();
    get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
      delta_, psi_gradient[i]);

    // wavefunction 'at (x-h)'
    work_.coeffRef(coeff.row(),coeff.col()) = coeff.value()-h;
    work_.coeffRef(coeff.col(),coeff.row()) = coeff.value()-h;
    //delta_ = es_quad_.eigenvectors().adjoint() * work_ * 
    //  es_quad_.eigenvectors().conjugate();
    delta_ = es_quad_.eigenvectors().transpose() * work_ * 
      es_quad_.eigenvectors();
    get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
      delta_, psi_work_);
    // restore the elements
    work_.coeffRef(coeff.row(),coeff.col()) = coeff.value();
    work_.coeffRef(coeff.col(),coeff.row()) = coeff.value();
    // gradient
    psi_gradient[i] -= psi_work_;
    psi_gradient[i] *= inv_2h;
    ++i;
  }
  */
}

void DisorderedSC::get_pair_amplitudes_sitebasis(const Eigen::VectorXd& es_eigenvalues, 
  const RealMatrix& es_eigenvectors, const RealMatrix& delta, Matrix& psi)
{
  // for INTRABAND pairing, phi_a = v_a/u_a for each band a
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_eigenvalues[k];
    auto delta_k = delta(k,k);
    double deltak_sq = delta_k * delta_k;
    if (deltak_sq<1.0E-12 && ek<0.0) {
      phi_k_[k] = large_number_;
    }
    else {
      phi_k_[k] = delta_k/(ek+std::sqrt(ek*ek+deltak_sq));
    }
    //std::cout << ek << "\n";
    //std::cout << delta_k << "\n";
    //std::cout << phi_k_[k] << "\n";
    //getchar();
  }
  //getchar();
  // pair ampitudes in site basis 
  for (unsigned i=0; i<num_sites_; ++i) {
    work_.col(i) = es_eigenvectors.col(i) * phi_k_[i];
  }
  work_ *= es_eigenvectors.transpose();
  double one_by_nk = 1.0/static_cast<double>(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) {
    for (unsigned j=0; j<num_sites_; ++j) {
      psi(i,j) = ampl_part(work_(i,j)) * one_by_nk;
      //std::cout << psi_work2_(i,j) << "\n"; getchar();
    }
  }
}



} // end namespace var