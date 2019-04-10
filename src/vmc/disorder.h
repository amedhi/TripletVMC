/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-13 11:22:16
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 15:45:32
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef DISORDER_H
#define DISORDER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <Eigen/Core>
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./random.h"
#include "./sysconfig.h"

namespace vmc {

class SiteDisorder
{
public:
  SiteDisorder() {}
  /*SiteDisorder(const input::Parameters& inputs);
  SiteDisorder(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model, RandomGenerator& rng)
  { init(inputs, graph, model, rng); }
  */
  ~SiteDisorder() {}
  const bool& check(const input::Parameters& inputs);
  int init(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model, const SysConfig& config,
    RandomGenerator& rng); 
  void set_current_config(const unsigned& current_config); 
  void save_optimal_parms(const var::parm_vector& optimal_parms); 
  const double& potential(const unsigned& site) const
    { return (*current_config_it_)(site); }
  operator bool(void) const { return exists_; }
  const bool& exists(void) const { return exists_; }
  const unsigned& num_configs(void) const { return num_configs_; }
  const var::parm_vector& get_optimal_parms(void); 
  bool optimal_parms_exists(const unsigned& config); 
private:
  bool exists_{false};
  unsigned num_sites_;
  unsigned num_bonds_;
  unsigned num_configs_;
  int current_config_{-1};
  unsigned num_opt_parms_;
  double bandwidth_{1.0};
  //bool potential_set_{false};
  std::vector<Eigen::VectorXd> disorder_pot_;
  std::vector<Eigen::VectorXd>::const_iterator current_config_it_;
  std::string datafile_prefix_;
  std::string optparm_path_;
  var::parm_vector optimal_parms_;
};





} // end namespace vmc

#endif