/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MC_OBSERVABLES_H
#define MC_OBSERVABLES_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./sysconfig.h"
#include "./energy.h"
#include "./sccorr.h"

namespace vmc {

class ObservableSet 
{
public:
  ObservableSet();
  ~ObservableSet() {}
  std::stringstream& headstream(void) { return headstream_; }
  //void init(const input::Parameters& inputs, 
  //  void (&print_copyright)(std::ostream& os), const lattice::LatticeGraph& graph, 
  //  const model::Hamiltonian& model, const SysConfig& config);
  void init(const input::Parameters& inputs, const lattice::LatticeGraph& graph, 
    const model::Hamiltonian& model, const SysConfig& config);
  void as_functions_of(const std::vector<std::string>& xvars=std::vector<std::string>());
  void as_functions_of(const std::string& xvar);
  void switch_off(void);
  void reset(void); 
  int do_measurement(const lattice::LatticeGraph& graph, 
    const model::Hamiltonian& model, const SysConfig& config);
  inline Energy& energy(void) { return energy_; }
  inline EnergyGradient& energy_grad(void) { return energy_grad_; }
  inline SC_Correlation& sc_corr(void) { return sc_corr_; }
  inline SR_Matrix& sr_matrix(void) { return sr_matrix_; }
  void finalize(void);
  void print_heading(void);
  void print_results(const std::vector<double>& xvals=std::vector<double>()); 
  void print_results(const double& xval); 
private:
  bool replace_mode_{true};
  std::stringstream headstream_;
  std::vector<std::string> xvars_;
  unsigned num_xvars_{0};
  Energy energy_;
  EnergyGradient energy_grad_;
  SC_Correlation sc_corr_;
  SR_Matrix sr_matrix_;
};





} // end namespace vmc

#endif
