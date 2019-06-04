/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:47:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 00:13:28
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./energy.h"

namespace vmc {

//-------------------ENERGY--------------------------------
void Energy::setup(const lattice::LatticeGraph& graph, 
  const model::Hamiltonian& model)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = graph.num_sites();
  std::vector<std::string> elem_names;
  model.get_term_names(elem_names);
  this->resize(elem_names.size(), elem_names);
  this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}
  
void Energy::measure(const lattice::LatticeGraph& graph, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  using op_id = model::op_id;
  //for (auto& elem : config_value_) elem = 0.0;
  config_value_.setZero();
  // bond energies
  if (model.have_bondterm()) {
    Matrix matrix_elem(model.num_bondterms(),graph.num_bond_types());
    matrix_elem.setZero();
    for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
      unsigned type = graph.bond_type(b);
      unsigned site_i = graph.source(b);
      unsigned site_j = graph.target(b);
      int bc_phase = graph.bond_sign(b);
      // matrix elements each term & bond type
      unsigned term = 0;
      for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
        matrix_elem(term,type) += config.apply(it->qn_operator(),site_i,site_j,bc_phase);
        term++;
      }
    }
    unsigned i = 0;
    for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
      for (unsigned btype=0; btype<graph.num_bond_types(); ++btype) {
        config_value_(i) += std::real(it->coupling(btype)*matrix_elem(i,btype));
      }
      i++;
    }
  }

  // site energies
  if (model.have_siteterm()) {
    Matrix matrix_elem(model.num_siteterms(),graph.num_site_types());
    matrix_elem.setZero();
    Eigen::VectorXi hubbard_nd(graph.num_site_types());
    hubbard_nd.setZero();
    // do it little differently
    unsigned term = model.num_bondterms();
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
          unsigned site = graph.site(s);
          unsigned type = graph.site_type(s);
          hubbard_nd(type) += config.apply_ni_updn(site);
        }
      }
      else {
        for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
          unsigned site = graph.site(s);
          unsigned type = graph.site_type(s);
          matrix_elem(term,type) += config.apply(it->qn_operator(), site);
        }
      }
      term++;
    }
    unsigned i = 0;
    unsigned n = model.num_bondterms();
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          config_value_(n+i) += std::real(it->coupling(stype)*hubbard_nd(stype));
        }
      }
      else {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          config_value_(n+i) += std::real(it->coupling(stype)*matrix_elem(i,stype));
        }
      }
      i++;
    }
  }

  // energy per site
  config_value_ /= num_sites_;
  // add to databin
  *this << config_value_;
}

//-------------------ENERGY GRADIENT--------------------------------
void EnergyGradient::setup(const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_varp_ = config.num_varparms();
  this->resize(config.num_varparms(),config.varp_names());
  grad_logpsi_.resize(config.num_varparms());
  config_value_.resize(2*config.num_varparms());
  grad_terms_.resize(2*config.num_varparms());
  setup_done_ = true;
}
  
void EnergyGradient::measure(const SysConfig& config, const double& config_energy)
{
  config.get_grad_logpsi(grad_logpsi_);
  unsigned n = 0;
  for (unsigned i=0; i<num_varp_; ++i) {
    config_value_[n] = config_energy*grad_logpsi_[i];
    config_value_[n+1] = grad_logpsi_[i];
    n += 2;
  }
  grad_terms_ << config_value_;
}

void EnergyGradient::finalize(const double& mean_energy)
{
  config_value_ = grad_terms_.mean_data(); 
  unsigned n = 0;
  mcdata::data_t energy_grad(num_varp_);
  for (unsigned i=0; i<num_varp_; ++i) {
    energy_grad(i) = 2.0*(config_value_[n] - mean_energy*config_value_[n+1]);
    n += 2;
  }
  //-------------------------------------------------
  //std::cout << mean_energy << "\n";
  //std::cout << energy_grad2_ << "\n";
  //for (int i=0; i<energy_grad_.size(); ++i) std::cout << energy_grad_[i] << "\n";
  //-------------------------------------------------
  *this << energy_grad;
}

//-------------------SR_Matrix (Stochastic Reconfiguration)---------------
void SR_Matrix::setup(const lattice::LatticeGraph& graph, const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = graph.num_sites();
  num_varp_ = config.num_varparms();
  // '\del(ln(psi))' plus upper triangular part of the sr_matrix 
  unsigned n = num_varp_ + num_varp_*(num_varp_+1)/2;
  this->resize(n);
  config_value_.resize(n);
  setup_done_ = true;
}

void SR_Matrix::measure(const RealVector& grad_logpsi)
{
  assert(grad_logpsi.size()==num_varp_);
  // operator 'del(ln(psi))' terms
  for (unsigned i=0; i<num_varp_; ++i) config_value_[i] = grad_logpsi[i];
  // flatten the upper triangular part to a vector
  unsigned k = num_varp_;
  for (unsigned i=0; i<num_varp_; ++i) {
    double x = grad_logpsi[i];
    for (unsigned j=i; j<num_varp_; ++j) {
      double y = grad_logpsi[j];
      config_value_[k] = x * y;
      ++k;
    }
  }
  *this << config_value_;
}

void SR_Matrix::get_matrix(Eigen::MatrixXd& sr_matrix) const
{
  // 'config_value' is used as temporary storage
  config_value_ = MC_Observable::mean_data();
  unsigned k = num_varp_;
  for (unsigned i=0; i<num_varp_; ++i) {
    double x = config_value_[i];
    for (unsigned j=i; j<num_varp_; ++j) {
      double y = config_value_[j];
      sr_matrix(i,j) = (config_value_[k] - x*y)/num_sites_;
      sr_matrix(j,i) = sr_matrix(i,j);
      ++k;
    }
  }
}


} // end namespace vmc
