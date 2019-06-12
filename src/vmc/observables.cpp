/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace vmc {

ObservableSet::ObservableSet() 
  : energy_("Energy")
  , energy_grad_("EnergyGradient")
  , sc_corr_("SC_Correlation")
  , sr_matrix_("SR_Matrix")
{
}

//void ObservableSet::init(const input::Parameters& inputs, 
//    void (&print_copyright)(std::ostream& os), const lattice::LatticeGraph& graph, 
//    const model::Hamiltonian& model, const SysConfig& config)
void ObservableSet::init(const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph, const model::Hamiltonian& model, 
  const SysConfig& config)
{
  // file open mode
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  // check which observables to calculate
  //for (auto& obs : *this) obs.get().check_on(inputs, replace_mode_);
  // heading message
  //print_copyright(headstream_);
  //model.print_info(headstream_);
  num_xvars_ = 0; 

  // switch on required observables
  energy_.check_on(inputs, replace_mode_);
  energy_grad_.check_on(inputs, replace_mode_);
  if (energy_grad_) energy_.switch_on();
  sc_corr_.check_on(inputs, replace_mode_);
  sr_matrix_.check_on(inputs, replace_mode_);


  // set up observables
  if (energy_) energy_.setup(graph, model);
  if (energy_grad_) energy_grad_.setup(config);
  if (sc_corr_) sc_corr_.setup(graph);
  if (sr_matrix_) sr_matrix_.setup(graph, config);
}

void ObservableSet::reset(void)
{
  if (energy_) energy_.reset();
  if (energy_grad_) energy_grad_.reset();
  if (sc_corr_) sc_corr_.reset();
  if (sr_matrix_) sr_matrix_.reset();
}

int ObservableSet::do_measurement(const lattice::LatticeGraph& graph, 
    const model::Hamiltonian& model, const SysConfig& config)
{
  if (energy_) energy_.measure(graph,model,config);
  if (energy_grad_) {
    if (!energy_) 
      throw std::logic_error("ObservableSet::measure: dependency not met for 'energy'");
    energy_grad_.measure(config, energy_.config_value().sum());
  }
  if (sc_corr_) sc_corr_.measure(graph,model,config);
  if (sr_matrix_) {
    if (!energy_grad_) 
      throw std::logic_error("ObservableSet::measure: dependency not met for 'sr_matrix_'");
    sr_matrix_.measure(energy_grad_.grad_logpsi());
  }
  return 0;
}

void ObservableSet::finalize(void)
{
  if (energy_grad_) {
    energy_grad_.finalize(energy_.mean_data().sum());
  }
}


void ObservableSet::as_functions_of(const std::vector<std::string>& xvars)
{
  xvars_ = xvars;
  num_xvars_ = xvars_.size();
}

void ObservableSet::as_functions_of(const std::string& xvar)
{
  xvars_ = {xvar};
  num_xvars_ = 1;
}

void ObservableSet::switch_off(void) {
  energy_.switch_off();
  energy_grad_.switch_off();
  sc_corr_.switch_off();
  sr_matrix_.switch_off();
}

void ObservableSet::print_heading(void)
{
  energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
  energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
  sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
  sr_matrix_.print_heading(headstream_.rdbuf()->str(),xvars_);
}

void ObservableSet::print_results(const std::vector<double>& xvals) 
{
  if (num_xvars_ == xvals.size()) 
    throw std::invalid_argument("Observables::print_result: 'x-vars' size mismatch");
  if (energy_) {
    energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_.print_result(xvals);
  }
  if (energy_grad_) {
    energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_grad_.print_result(xvals);
  }
  if (sc_corr_) {
    sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    sc_corr_.print_result(xvals);
  }
}

void ObservableSet::print_results(const double& xval) 
{
  if (num_xvars_ != 1) 
    throw std::invalid_argument("ObservableSet::print_result: 'x-vars' size mismatch");
  std::vector<double> xvals{xval};
  if (energy_) {
    energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_.print_result(xvals);
  }
  if (energy_grad_) {
    energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_grad_.print_result(xvals);
  }
  if (sc_corr_) {
    sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    sc_corr_.print_result(xvals);
  }
}


} // end namespace vmc

