/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-20 10:18:07
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-30 22:12:35
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef VARPARM_H
#define VARPARM_H

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <Eigen/Core>
#include "../scheduler/task.h"

namespace var {
//using parm_vector = std::vector<double>;
using parm_vector = Eigen::VectorXd;
using name_id_map = std::map<std::string,unsigned>; 
const double default_h = 0.005; // step size for numerical differentiation

class varparm_t
{
public:
  varparm_t() : val_{0.0}, lb_{0.0}, ub_{0.0}, diff_h_{0.0} {}
  varparm_t(const double& val, const double& lb, const double& ub, 
    const double& diff_h=default_h)
  : val_{val}, lb_{lb}, ub_{ub}, diff_h_{diff_h} 
    {
      // actual bounds to allow for differentiation 
      lb_ += 2*diff_h_;  
      ub_ -= 2*diff_h_;  
    }
  ~varparm_t() {}
  bool change_value(const double& newval) 
  {
    if (newval<lb_ || newval>ub_) return false;
    val_ = newval; return true;
  }
  void set_name(const name_id_map::const_iterator& it ) { map_it_=it; }
  const double& value(void) const { return val_; }
  const double& lbound(void) const { return lb_; }
  const double& ubound(void) const { return ub_; }
  const double& diff_h(void) const { return diff_h_; }
  const std::string& name(void) const { return map_it_->first; }
private:
  double val_;
  double lb_;
  double ub_;
  double diff_h_; // step size for finite differential
  name_id_map::const_iterator map_it_;
};

class VariationalParms : public std::vector<varparm_t>
{
public:
  VariationalParms(); 
  ~VariationalParms() {} 
  int add(const std::string& name, const double& val, const double& lb, 
    const double& ub, const double& diff_h=default_h);
  using std::vector<varparm_t>::operator[];
  const varparm_t& operator[](const std::string& pname) const 
    { return std::vector<varparm_t>::operator[](name_id_map_.at(pname)); } 
private:
  unsigned num_parms_{0};
  name_id_map name_id_map_;
};

/*
struct parm_t
{
  std::string name;
  double val;
  double lower;
  double upper;
};

class VariationalParms : public std::map<std::string,unsigned>
{
public:
  VariationalParms(); 
  ~VariationalParms() { --num_sets_; }
  void clear(void);
  int add(const std::string& name, const double& value, 
    const double& lower, const double& upper);
  int append_set(const VariationalParms& other_set);
  void update(const input::Parameters& inputs);
  void update(const std::vector<double>& pvec, const unsigned& begin, const unsigned& end);
  const unsigned& id(void) const { return my_id_; }
  parm_t operator[](const unsigned& i) const;
  const std::vector<double>& values(void) const { return val_; }
  const std::vector<double>& lbounds(void) const { return lower_; }
  const std::vector<double>& ubounds(void) const { return upper_; }
  //int update_set(VariationalParms& other_set);
private:

  using idx_pair = std::pair<unsigned,unsigned>;
  std::vector<double> val_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::map<unsigned,idx_pair> set_idx_;
  //std::vector<std::pair<> > set_offset_;
  unsigned num_parms_{0};
  unsigned my_id_{0};
  //std::vector<std::pair<it_begin, it_end> > v;
  static unsigned num_sets_; 

  void append(const std::string& name, const double& value, const double& lower, 
    const double& upper);
};
*/

} // end namespace var

#endif