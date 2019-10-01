/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-09-26 22:16:30
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./particle.h"

namespace vmc {

//-------------------Site Occupancy--------------------------------
void SiteOccupancy::setup(const lattice::LatticeGraph& graph, 
	const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = graph.num_sites();
  num_basis_sites_ = graph.lattice().num_basis_sites();
  num_particles_ = config.num_particles();
  std::vector<std::string> elem_names(num_basis_sites_);
  for (int i=0; i<num_basis_sites_; ++i) {
  	std::ostringstream ss;
  	ss << "site-"<<i;
  	elem_names[i] = ss.str();
  }
  this->resize(elem_names.size(), elem_names);
  this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void SiteOccupancy::measure(const lattice::LatticeGraph& graph, 
	const SysConfig& config) 
{
  IntVector matrix_elem(num_basis_sites_);
  matrix_elem.setZero();
  for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
    int site = graph.site(s);
    int basis = graph.site_uid(s);
    matrix_elem(basis) += config.apply(model::op::ni_sigma(), site);
  }
  for (int i=0; i<num_basis_sites_; ++i) {
  	config_value_[i] = static_cast<double>(matrix_elem[i])/num_particles_;
  }
  // add to databin
  *this << config_value_;
}




} // end namespave vmc








