/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-19 07:07:40
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:12:26
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "./cppoptlib/boundedproblem.h"
#include "./cppoptlib/solver/lbfgsbsolver.h"
#include "../scheduler/task.h"

namespace optimizer {

class Problem : public cppoptlib::BoundedProblem<double> 
{
  public:
  using super_type = BoundedProblem<double>;
  using VectorXd = super_type::TVector;
  using super_type::super_type;
  Problem(int RunDim=Eigen::Dynamic) : super_type(RunDim) {} 
  virtual ~Problem() {} 
  virtual double operator()(const TVector &x, TVector &grad) = 0;
  double operator()(const TVector &x) { return 0; }
  double value(const TVector &x) { return 0; }
};

class Optimizer
{
public:
  Optimizer() {}
  ~Optimizer() {}
  int init(const input::Parameters& parms, Problem& problem);
  bool optimize(Problem& problem);
  //{ solver.minimize(problem, x); }
private:
  cppoptlib::LbfgsbSolver<Problem> solver_;
};


} // end namespace vmc

#endif