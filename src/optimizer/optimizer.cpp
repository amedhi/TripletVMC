/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-19 16:03:53
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-19 23:14:39
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./optimizer.h"

namespace optimizer {

int Optimizer::init(const input::Parameters& parms, Problem& problem)
{
  //problem.set_box_constraints();
  return 0;
}


bool Optimizer::optimize(Problem& problem)
{
  Problem::TVector x = 0.5*(problem.lowerBound() + problem.upperBound()); 
  solver_.minimize(problem, x);

  std::cout << "optimal parms = " << x.transpose() << "\n";

  return true;
}



} // end name space optimizer
