/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-09-26 13:53:41
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_PARTICLE_H
#define OBS_PARTICLE_H

#include "../mcdata/mc_observable.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./sysconfig.h"

namespace vmc {

class SiteOccupancy : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::LatticeGraph& graph, const SysConfig& config);
  void measure(const lattice::LatticeGraph& graph, const SysConfig& config);
  const mcdata::data_t& config_value(void) const { return config_value_; }
private:
  bool setup_done_{false};
  int num_sites_{0};
  int num_basis_sites_{0};
  int num_particles_{0};
  mcdata::data_t config_value_;
};

} // end namespave vmc



#endif