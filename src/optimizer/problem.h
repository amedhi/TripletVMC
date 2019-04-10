/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-19 16:22:11
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-19 16:26:30
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H



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


} // end namespace vmc

#endif