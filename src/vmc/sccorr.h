/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:32:31
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:07:57
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_SCCORR_H
#define OBS_SCCORR_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../mcdata/mc_observable.h"
#include "./sysconfig.h"

namespace vmc {

class SC_Correlation : public mcdata::MC_Observable
{
public:
  using site_t = lattice::LatticeGraph::site_descriptor;
  using MC_Observable::MC_Observable;
  void setup(const lattice::LatticeGraph& graph);
  void measure(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
    const SysConfig& config);
  const unsigned& num_site_pairs(void) const { return src_pairs_size_; }
  const std::pair<site_t,site_t>& site_pair(const unsigned& i) const 
    { return src_pairs_[i]; }
  void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  int max_dist_{0};
  unsigned num_bond_types_{0};
  unsigned src_pairs_size_{0};
  std::vector<std::pair<site_t,site_t> > src_pairs_;
  std::vector<int> pair_distance_;
  std::vector<Eigen::MatrixXd> bond_pair_corr_;
  std::vector<Eigen::MatrixXi> num_symm_pairs_;
};


} // end namespace vmc

#endif
