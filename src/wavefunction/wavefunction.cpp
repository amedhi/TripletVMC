/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-11 10:28:20
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "wavefunction.h"
#include <boost/algorithm/string.hpp>

namespace var {

Wavefunction::Wavefunction(const lattice::LatticeGraph& graph,
  const input::Parameters& inputs)
  : num_sites_(graph.num_sites())
{
  num_states_ = 2*num_sites_;
  name_ = inputs.set_value("wavefunction", "NONE");
  boost::to_upper(name_);
  if (name_ == "FERMISEA") {
    groundstate_.reset(new Fermisea(inputs,graph));
  }
  else if (name_ == "SWAVE_SC") {
    groundstate_.reset(new BCS_State(bcs::swave,inputs,graph));
  }
  else if (name_ == "DWAVE_SC") {
    groundstate_.reset(new BCS_State(bcs::dwave,inputs,graph));
  }
  else if (name_ == "SC_P+IP") {
    groundstate_.reset(new BCS_State(bcs::p_plus_ip,inputs,graph));
  }
  else {
    throw std::range_error("Wavefunction::Wavefunction: unidefined wavefunction");
  }
  // resize
  psi_up_.resize(num_states_,num_states_);
  psi_gradient_.resize(varparms().size());
  for (int i=0; i<varparms().size(); ++i)
    psi_gradient_[i].resize(num_states_,num_states_);
}

std::string Wavefunction::signature_str(void) const
{
  // signature string
  std::ostringstream signature;
  signature << "wf_N"; 
  signature << std::setfill('0'); 
  signature << std::setw(3) << groundstate_->num_upspins(); 
  signature << std::setw(3) << groundstate_->num_dnspins(); 
  return signature.str();
}

int Wavefunction::compute(const lattice::LatticeGraph& graph, 
  const input::Parameters& inputs, const bool& psi_gradient)
{
  groundstate_->update(inputs);
  groundstate_->get_wf_amplitudes(psi_up_);
  if (psi_gradient) {
    groundstate_->get_wf_gradient(psi_gradient_);
    have_gradient_ = true;
  }
  else have_gradient_ = false;
  return 0;
}

int Wavefunction::compute(const lattice::LatticeGraph& graph, const var::parm_vector& pvector,
  const unsigned& start_pos, const bool& psi_gradient)
{
  groundstate_->update(pvector,start_pos);
  groundstate_->get_wf_amplitudes(psi_up_);
  if (psi_gradient) {
    groundstate_->get_wf_gradient(psi_gradient_);
    have_gradient_ = true;
  }
  else have_gradient_ = false;
  return 0;
}

void Wavefunction::get_amplitudes(Matrix& psi, const std::vector<int>& spin_states) const 
{
  // Pfafian matrix for the state
  for (int i=0; i<spin_states.size(); ++i) {
    for (int j=i; j<spin_states.size(); ++j) {
      psi(i,j) = psi_up_(spin_states[i], spin_states[j]);
      psi(j,i) = -psi(i,j);
    }
  }
}

void Wavefunction::get_amplitudes(ColVector& psi_row, RowVector& psi_col,
  const int& spin, const int& new_state, const std::vector<int>& spin_states) const
{
  for (int i=0; i<spin; ++i) {
    psi_col(i) = psi_up_(spin_states[i],new_state);
    psi_row(i) = -psi_col(i);
  }
  for (int i=spin+1; i<spin_states.size(); ++i) {
    psi_row(i) = psi_up_(new_state,spin_states[i]);
    psi_col(i) = -psi_row(i);
  }
  psi_row(spin) = 0.0;
  psi_col(spin) = 0.0;
}

void Wavefunction::get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const
{
  for (int j=0; j<col.size(); ++j)
    psi_vec[j] = psi_up_(irow,col[j]);
}

void Wavefunction::get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const
{
  for (int j=0; j<row.size(); ++j)
    psi_vec[j] = psi_up_(row[j],icol);
}

void Wavefunction::get_amplitudes(amplitude_t& elem, const int& irow, 
  const int& jcol) const
{
  elem = psi_up_(irow,jcol);
}

void Wavefunction::get_gradients(Matrix& psi_grad, const int& n, 
  const std::vector<int>& spin_states) const
{
  if (!have_gradient_) 
    throw std::logic_error("Wavefunction::get_gradients: gradients were not computed");

  for (int i=0; i<spin_states.size(); ++i) {
    for (int j=i; j<spin_states.size(); ++j) {
      psi_grad(i,j) = psi_gradient_[n](spin_states[i], spin_states[j]);
      psi_grad(j,i) = -psi_grad(i,j);
    }
  }
}

void Wavefunction::get_vparm_names(std::vector<std::string>& vparm_names, 
  int start_pos) const
{
  int i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_names[start_pos+i] = p.name(); ++i;
  }
}

void Wavefunction::get_vparm_values(var::parm_vector& vparm_values, 
  int start_pos)
{
  int i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_vector(std::vector<double>& vparm_values, 
  int start_pos)
{
  int i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_lbound(var::parm_vector& vparm_lb, 
  int start_pos) const
{
  int i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_lb[start_pos+i] = p.lbound(); ++i;
  }
}

void Wavefunction::get_vparm_ubound(var::parm_vector& vparm_ub, 
  int start_pos) const
{
  int i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_ub[start_pos+i] = p.ubound(); ++i;
  }
}

} // end namespace var











