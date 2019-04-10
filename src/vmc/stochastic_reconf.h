/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:07:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-18 12:33:52
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef STOCHASTIC_RECONF_H
#define STOCHASTIC_RECONF_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "../utils/utils.h"
#include "./vmc.h"

namespace vmc {

class StochasticReconf 
{
public:
  StochasticReconf() {} 
  StochasticReconf(const input::Parameters& parms); 
  ~StochasticReconf() {}
  int init(const input::Parameters& parms, const VMC& vmc);
  int optimize(VMC& vmc);
  const var::parm_vector& optimal_parms(void) const { return vparms_; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  mcdata::MC_Observable optimal_parms_;
  unsigned num_parms_;
  var::parm_vector vparms_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  Eigen::MatrixXd sr_matrix_;
  Eigen::VectorXd grad_;
  // Mann-Kendall trend test for converegence
  util::MK_Statistic mk_statistic_;
  // optimization parameters
  unsigned num_sim_samples_{1000};
  unsigned num_opt_samples_{30};
  unsigned max_iter_{500};
  unsigned refinement_cycle_{100};
  unsigned mk_series_len_{40};
  double start_tstep_{0.05};
  double stabilizer_{1.0E-4};
  double grad_tol_{0.01};
  double mk_thresold_{0.30};
  bool print_progress_{false};
  bool print_log_{true};

  // progress file
  std::ofstream logfile_;
};


} // end namespace vmc

#endif