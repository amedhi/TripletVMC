/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-20 10:47:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-12 22:52:23
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./varparm.h"

namespace var {


VariationalParms::VariationalParms(void)
{
  clear();
  num_parms_ = 0;
}

int VariationalParms::add(const std::string& name, const double& val, 
    const double& lb, const double& ub, const double& diff_h)
{
  if (val<lb || val>ub) {
    throw std::invalid_argument("VariationalParm::add_parameter: invalid parameter bound");
  }
  name_id_map::const_iterator it;
  bool status;
  std::tie(it, status) = name_id_map_.insert({name,num_parms_});
  if (!status) throw std::logic_error("VariationalParm::add_parameter: parameter already exists");
  push_back({val, lb, ub, diff_h});
  back().set_name(it);
  ++num_parms_;
  return num_parms_;
}

/*
unsigned VariationalParms::num_sets_ = 0;

VariationalParms::VariationalParms(void)
{
  this->clear();
  my_id_ = num_sets_++;
  set_idx_.insert({my_id_,std::make_pair(0,0)});
}

void VariationalParms::clear(void)
{
  std::map<std::string,unsigned>::clear();
  val_.clear();
  lower_.clear();
  upper_.clear();
  num_parms_ = 0;
} 

int VariationalParms::add(const std::string& name, const double& value, 
    const double& lower, const double& upper)
{
  if (set_idx_.size() > 1) {
    throw std::logic_error("VariationalParm::add_parameter: parameter can't be added, new set exist");
  }
  if (!insert({name,num_parms_}).second) {
    throw std::logic_error("VariationalParm::add_parameter: parameter with same name exists");
  }
  val_.push_back(value);
  lower_.push_back(lower);
  upper_.push_back(upper);
  ++num_parms_;
  set_idx_.at(my_id_).second = num_parms_;
  return num_parms_;
}

int VariationalParms::append_set(const VariationalParms& other_set)
{
  if (other_set.size()==0) return num_parms_;
  // new set to be appended
  set_idx_.insert({other_set.id(),std::make_pair(num_parms_,num_parms_)});
  for (unsigned i=0; i<other_set.size(); ++i) {
    append(other_set[i].name,other_set[i].val,other_set[i].lower,other_set[i].upper);
  }
  set_idx_.at(other_set.id()).second = num_parms_;
  return num_parms_;
}

parm_t VariationalParms::operator[](const unsigned& i) const
{
  parm_t p;
  auto it = this->cbegin();
  for (unsigned n=0; n<i; ++n) {
    ++it;
    if (it == this->cend()) 
      throw std::range_error("VariationalParms::operator[]: out-of-range element");
  }
  p.name = it->first;
  p.val = val_[i];
  p.lower = lower_[i];
  p.upper = upper_[i];
  return p;
}

void VariationalParms::update(const input::Parameters& inputs)
{
  unsigned begin = set_idx_.at(my_id_).first;
  unsigned end = set_idx_.at(my_id_).second;
  auto it = this->begin();
  for (unsigned i=begin; i<end; ++i) {
    val_[i] = inputs.set_value(it->first,val_[i]);
    ++it;
  }
}

void VariationalParms::update(const std::vector<double>& pvec, const unsigned& begin, 
  const unsigned& end)
{
  if ((end-begin) != size()) 
    throw std::invalid_argument("VariationalParm::update: vector size mismatch");
  unsigned i = 0;
  for (unsigned j=begin; j<end; ++j) {
    val_[i] = pvec[j];
    ++i;
  }
}

void VariationalParms::append(const std::string& name, const double& value, 
    const double& lower, const double& upper)
{
  if (!insert({name,num_parms_}).second) {
    throw std::logic_error("VariationalParm::append: parameter with same name exists");
  }
  val_.push_back(value);
  lower_.push_back(lower);
  upper_.push_back(upper);
  ++num_parms_;
}
*/


} // end namespace var


