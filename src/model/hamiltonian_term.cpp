/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-13 22:10:41
*----------------------------------------------------------------------------*/
#include "hamiltonian_term.h"

namespace model {

/*--------------------------CouplingConstant--------------------*/
const int CouplingConstant::global_type = -1;

CouplingConstant::CouplingConstant(const std::string& expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1;
  valid_ = true;
}

std::pair<CouplingConstant::iterator, bool> CouplingConstant::insert(const value_type& val) 
{
  std::pair<iterator, bool> res = super_type::insert(val);
  num_types_ = super_type::size();
  valid_ = true;
  return res;
}

CouplingConstant& CouplingConstant::operator=(const std::string expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

CouplingConstant::CouplingConstant(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  create(type0, type1, type2, type3, type4, type5); 
}

void CouplingConstant::clear(void)
{
  super_type::clear();
  num_types_ = 0;
  valid_ = false;
}

void CouplingConstant::create(const unsigned& num_types) 
{
  super_type::clear();
  num_types_ = num_types;
  valid_ = false;
}

void CouplingConstant::create(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  super_type::insert(type0);
  num_types_=1;
  if (type1.second != "_null_") {
    super_type::insert(type1); num_types_++;
  }
  if (type2.second != "_null_") {
    super_type::insert(type2); num_types_++;
  }
  if (type3.second != "_null_") {
    super_type::insert(type3); num_types_++;
  }
  if (type4.second != "_null_") {
    super_type::insert(type4); num_types_++;
  }
  if (type5.second != "_null_") {
    super_type::insert(type5); num_types_++;
  }
  valid_=true;
}

void CouplingConstant::add_type(const unsigned& type, const std::string& expr)
{
  super_type::insert({type, expr});
  valid_ = (num_types_==static_cast<int>(size()));
}

void CouplingConstant::add_type(const value_type& val)
{
  super_type::insert(val);
  valid_ = (num_types_==static_cast<int>(size()));
}



//-----------------------HamiltonianTerm-------------------------
void HamiltonianTerm::construct(const std::string& name, const op::quantum_op& op, 
  const CouplingConstant& cc, const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("HamiltonianTerm:: Invalid CouplingConstant");
  name_ = name;
  op_ = op;
  cc_ = cc;
  max_operand_types_ = size;

  cc_values_.resize(max_operand_types_);
  is_defined_.resize(max_operand_types_);
  for (unsigned i=0; i<max_operand_types_; ++i) {
    cc_values_[i] = 0.0;
    is_defined_[i] = false;
  }

  // if the 'cc' is implicitly defined for all types 
  if (cc_.size()==1 && cc_.begin()->first==CouplingConstant::global_type) {
    for (unsigned i=0; i<max_operand_types_; ++i) is_defined_[i] = true;
  } 
  else {
    // operator is defined only for those types for which 'cc' is set explicitly
    for (const auto& p : cc) is_defined_[p.first] = true;
  }
}

void HamiltonianTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& c : cvals) vars[c.first] = c.second;
  for (const auto& p : pvals) vars[p.first] = p.second;
  try { 
    // if the 'cc' is implicitly defined for all types 
    if (cc_.size()==1 && cc_.begin()->first==CouplingConstant::global_type) {
      double val = expr.evaluate(cc_.begin()->second, vars);
      for (auto& v : cc_values_) v = val;
    }
    else {
      for (const auto& p : cc_) 
        cc_values_[p.first] = expr.evaluate(p.second, vars); 
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "BondOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}



} // end namespace model
