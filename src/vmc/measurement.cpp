/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-17 23:30:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:27:40
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include "vmc.h"

namespace vmc {

void VMC::init_measurements(void)
{
  //if (observables.sccf()) {
  //}
}

void VMC::do_measurements(const observable_set& obs_set)
{
  if (observables.sc_corr()) {
    observables.sc_corr().measure(graph, model, config);
  } 
  if (observables.eenergy()) {
    observables.eenergy().measure(graph, model, config);
  } 

  // observable_set::normal (as specied in input)
  if (obs_set==observable_set::normal) {
    if (observables.need_energy()) {
      // energy
      term_energy_ = get_energy();
      if (observables.energy()) observables.energy() << term_energy_;
      // energy gradients
      if (observables.energy_grad()) {
        double total_energy = term_energy_.sum();
        observables.total_energy() << total_energy;
        measure_energy_grad(total_energy);
      } 
      // total energy  
      else if (observables.total_energy()) {
        double total_energy = term_energy_.sum();
        observables.total_energy() << total_energy;
      }
    }
    return;
  }

  // observable_set::sr_matrix
  else if (obs_set==observable_set::sr_coeffs) {
    // total energy
    double total_energy = get_energy().sum();
    observables.total_energy() << total_energy;
    // energy gradient
    measure_energy_grad(total_energy);
    // sr matrix
    // operator 'del(ln(psi))' terms
    for (unsigned i=0; i<num_varparms_; ++i) sr_coeffs_[i] = grad_logpsi_[i];
    // flatten the upper triangular part to a vector
    unsigned k = num_varparms_;
    for (unsigned i=0; i<num_varparms_; ++i) {
      double x = grad_logpsi_[i];
      for (unsigned j=i; j<num_varparms_; ++j) {
        double y = grad_logpsi_[j];
        sr_coeffs_[k] = x * y;
        ++k;
      }
    }
    observables.sr_coeffs() << sr_coeffs_;
    return;
  }

  // observable_set::energy_grad
  else if (obs_set==observable_set::energy_grad) {
    double total_energy = get_energy().sum();
    observables.total_energy() << total_energy;
    if (observables.energy_grad()) {
      measure_energy_grad(total_energy);
    }
    return;
  }
  else {
    throw std::range_error("VMC::do_measurements: unknown observables set");
  }
}

void VMC::measure_energy_grad(const double& total_en)
{
  config.get_grad_logpsi(grad_logpsi_);
  unsigned n = 0;
  for (unsigned i=0; i<num_varparms_; ++i) {
    energy_grad2_[n] = total_en * grad_logpsi_[i];
    energy_grad2_[n+1] = grad_logpsi_[i];
    n += 2;
  }
  observables.energy_grad2() << energy_grad2_;
}

int VMC::finalize_energy_grad(void)
{
  double mean_energy = observables.total_energy().mean();
  energy_grad2_ = observables.energy_grad2().mean_data(); 
  unsigned n = 0;
  for (unsigned i=0; i<num_varparms_; ++i) {
    energy_grad_(i) = 2.0*(energy_grad2_[n] - mean_energy * energy_grad2_[n+1]);
    n += 2;
  }
  //-------------------------------------------------
  //std::cout << mean_energy << "\n";
  //std::cout << energy_grad2_ << "\n";
  //for (int i=0; i<energy_grad_.size(); ++i) std::cout << energy_grad_[i] << "\n";
  //-------------------------------------------------
  observables.energy_grad() << energy_grad_;
  return 0;
}


Observable::data_t VMC::get_energy(void) const
{
  using op_id = model::op_id;
  //for (auto& elem : term_energy_) elem = 0.0;
  term_energy_.resize(model.num_terms());
  term_energy_.setZero();
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
        term_energy_(i) += std::real(it->coupling(btype)*matrix_elem(i,btype));
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
          hubbard_nd(type) += config.apply_niup_nidn(site);
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
    /*for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
      unsigned type = graph.site_type(s);
      // matrix elements each term & site type
      unsigned term = model.num_bondterms();
      for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
        //matrix_elem(term,type) += config.apply(it->qn_operator(), s);
        term++;
      }
    }*/
    unsigned i = 0;
    unsigned n = model.num_bondterms();
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          term_energy_(n+i) += std::real(it->coupling(stype)*hubbard_nd(stype));
        }
      }
      else {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          term_energy_(n+i) += std::real(it->coupling(stype)*matrix_elem(i,stype));
        }
      }
      i++;
    }
  }

  // disorder energies
  if (site_disorder_) {
    unsigned n = model.num_bondterms()+model.num_siteterms();
    double disorder_en = 0.0;
    for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
      unsigned site = graph.site(s);
      int n_i = config.apply(model::op::ni_sigma(), site);
      disorder_en += std::real(n_i * site_disorder_.potential(site));
    }
    term_energy_(n) = disorder_en;
  }

  // energy per site
  return term_energy_/num_sites_;
}



} // end namespace vmc
