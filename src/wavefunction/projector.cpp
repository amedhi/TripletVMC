/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:17:49
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-05 11:47:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./projector.h"

namespace var {

void WavefunProjector::init(const input::Parameters& parms) 
{
  int info;
  gutzwiller_proj_ = parms.set_value("gutzwiller_proj", false, info);
  double g = 1.0;
  if (gutzwiller_proj_) {
    g = parms.set_value("gfactor", 1.0);
    if (g<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'g'-value.");
    varparms_.add("gfactor", g, 1.0E-2, 1.0);
  }
  //pfactors_.insert({"g", g});
  // gw ratio
  gw_ratio_ = {1.0,1.0,1.0};
  //gw_ratio_[0] = 1.0/g;  // nd_increament = -1
  //gw_ratio_[1] = 1.0;    // nd_increament = 0
  //gw_ratio_[2] = g;      // nd_increament = 1
}

void WavefunProjector::update(const input::Parameters& inputs) 
{ 
  for (auto& p : varparms_) {
    double x = inputs.set_value(p.name(), p.value());
    p.change_value(x);
  }
  //varparms_.update(inputs); 
  if (gutzwiller_proj_) set_gw_ratio();
} 

void WavefunProjector::update(const var::parm_vector& pvector, const unsigned& start_pos)
{ 
  //varparms_.update(vparms,begin,end); 
  unsigned i = 0;
  for (auto& p : varparms_) {
    p.change_value(pvector[start_pos+i]);
    i++;
  }
  if (gutzwiller_proj_) set_gw_ratio();
}

void WavefunProjector::get_vparm_names(std::vector<std::string>& vparm_names, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_names[start_pos+i] = p.name(); ++i;
  }
}

void WavefunProjector::get_vparm_values(var::parm_vector& vparm_values, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void WavefunProjector::get_vparm_vector(std::vector<double>& vparm_values, unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void WavefunProjector::get_vparm_lbound(var::parm_vector& vparm_lb, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_lb[start_pos+i] = p.lbound(); ++i;
  }
}

void WavefunProjector::get_vparm_ubound(var::parm_vector& vparm_ub, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_ub[start_pos+i] = p.ubound(); ++i;
  }
}


double WavefunProjector::gw_factor(void) const 
{
  //return varparms_.values()[varparms_.at("gfactor")];
  return varparms_["gfactor"].value();
}

void WavefunProjector::set_gw_ratio(void) 
{ 
  double g = gw_factor();
  //std::cout << "### g = " << g << "\n";
  if (g<0.0) throw std::range_error("WavefunProjector::set_gw_ratio: out-of-range 'g'-value");
  gw_ratio_[0] = 1.0/g;  // nd_increament = -1
  gw_ratio_[1] = 1.0;    // nd_increament = 0
  gw_ratio_[2] = g;      // nd_increament = 1
} 




} // end namespace var
