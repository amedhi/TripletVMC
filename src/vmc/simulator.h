/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include "./vmc.h"
#include "./stochastic_reconf.h"
//#include "../optimizer/optimizer.h"

namespace vmc {

class Simulator : public scheduler::Worker
{
public:
  Simulator(const input::Parameters& parms); 
  ~Simulator() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  int run(const input::Parameters& parms, 
    const scheduler::mpi_communicator& mpi_comm) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  const VMC& sim(void) { return vmc; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  VMC vmc;
  StochasticReconf sreconf;
  //optimizer::Optimizer nlopt_;
  bool optimization_mode_{false};

  //static double enfunc(const std::vector<double>& x, std::vector<double>& grad, 
  //  void *my_func_data);
};

//double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data);

} // end namespace vmc

#endif
