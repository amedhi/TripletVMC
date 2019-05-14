/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:20:28
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-20 04:54:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "basisstate.h"
#include <stdexcept>
#include <algorithm>

namespace vmc {

void FockBasis::init(const int& num_sites, const bool& allow_dbl) 
{
  num_sites_ = num_sites;
  num_states_ = 2*num_sites_;
  state_.resize(num_states_);
  state_.setZero();
  spin_id_.resize(num_states_);
  spin_id_.setConstant(-1); 
  double_occupancy_ = allow_dbl;
  spin_states_.clear();
  up_states_.clear();
  dn_states_.clear();
  uphole_states_.clear();
  dnhole_states_.clear();
  proposed_move_ = move_t::null;
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
}

void FockBasis::init_spins(const int& num_upspins, const int& num_dnspins)
{
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  num_spins_ = num_upspins_+num_dnspins_;
  if (num_upspins_>num_sites_ || num_dnspins_>num_sites_)
    throw std::range_error("* FockBasis::init_spins: spin number exceeds capacity");
  if (!double_occupancy_ && num_spins_>num_sites_)
    throw std::range_error("* FockBasis::init_spins: spin number exceeds capacity");
  num_upholes_ = num_sites_ - num_upspins;
  num_dnholes_ = num_sites_ - num_dnspins;
  // resizing
  spin_states_.resize(num_spins_);
  up_states_.resize(num_upspins_);
  dn_states_.resize(num_dnspins_);
  dnspin_sites_.resize(num_dnspins_);
  uphole_states_.resize(num_upholes_);
  dnhole_states_.resize(num_dnholes_);
  // random generator
  if (num_upspins_>0) rng_.set_upspin_generator(0,num_upspins_-1);
  else rng_.set_upspin_generator(0,0);
  if (num_dnspins_>0) rng_.set_dnspin_generator(0,num_dnspins_-1);
  else rng_.set_dnspin_generator(0,0);
  if (num_upholes_>0) rng_.set_uphole_generator(0,num_upholes_-1);
  else rng_.set_uphole_generator(0,0);
  if (num_dnholes_>0) rng_.set_dnhole_generator(0,num_dnholes_-1);
  else rng_.set_dnhole_generator(0,0);
  // random initial configuration
  set_random();
}

void FockBasis::set_random(void)
{
  proposed_move_ = move_t::null;
  state_.setZero();
  std::vector<int> all_up_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_up_states[i] = i;
  std::shuffle(all_up_states.begin(),all_up_states.end(),rng_);
  for (int i=0; i<num_upspins_; ++i) {
    int state = all_up_states[i];
    state_[state] = 1;
    spin_id_[state] = i;
    up_states_[i] = state;
    spin_states_[i] = state;
  }
  int j=0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    uphole_states_[j++] = all_up_states[i];
  }

  // DN spins & holes
  if (double_occupancy_) {
    std::vector<int> all_dn_states(num_sites_);
    for (int i=0; i<num_sites_; ++i) all_dn_states[i] = num_sites_+i;
    std::shuffle(all_dn_states.begin(),all_dn_states.end(),rng_);
    int j = num_upspins_;
    for (int i=0; i<num_dnspins_; ++i) {
      int state = all_dn_states[i];
      state_[state] = 1;
      spin_id_[state] = i;
      dn_states_[i] = state;
      spin_states_[j++] = state;
    }
    j = 0;
    for (int i=num_dnspins_; i<num_sites_; ++i) {
      dnhole_states_[j++] = all_dn_states[i];
    }
  }
  else {
    // DN spins
    int j = 0;
    for (int i=num_upspins_; i<num_spins_; ++i) {
      int state = num_sites_+all_up_states[i];
      state_[state] = 1;
      spin_id_[state] = j;
      dn_states_[j] = state;
      spin_states_[i] = state;
      ++j;
    }
    // DN holes
    j = 0;
    for (int i=0; i<num_upspins_; ++i) {
      int state = num_sites_+all_up_states[i];
      dnhole_states_[j++] = state;
    }
    for (int i=num_spins_; i<num_sites_; ++i) {
      int state = num_sites_+all_up_states[i];
      dnhole_states_[j++] = state;
    }
  }
  // number of doubly occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (int i=0; i<num_sites_; ++i) {
      if (state_[i]==1 && state_[i+num_sites_]==1)
        num_dblocc_sites_++;
    }
  }
}

void FockBasis::set_custom(void)
{
  proposed_move_ = move_t::null;
  state_.setZero();
  std::vector<int> all_up_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_up_states[i] = i;
  //std::shuffle(all_up_states.begin(),all_up_states.end(),rng_);
  for (int i=0; i<num_upspins_; ++i) {
    int state = all_up_states[i];
    state_[state] = 1;
    spin_id_[state] = i;
    up_states_[i] = state;
    spin_states_[i] = state;
  }
  int j=0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    uphole_states_[j++] = all_up_states[i];
  }

  // DN spins & holes
  std::vector<int> all_dn_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_dn_states[i] = num_sites_+i;
  //std::shuffle(all_dn_states.begin(),all_dn_states.end(),rng_);
  int last_site = num_sites_-1;
  j = num_upspins_;
  for (int i=0; i<num_dnspins_; ++i) {
    int state = all_dn_states[last_site-i];
    state_[state] = 1;
    spin_id_[state] = i;
    dn_states_[i] = state;
    spin_states_[j++] = state;
  }
  j = 0;
  for (int i=num_dnspins_; i<num_sites_; ++i) {
    dnhole_states_[j++] = all_dn_states[last_site-i];
  }

  // number of doublely occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (int i=0; i<num_sites_; ++i) {
      if (state_[i]==1 && state_[i+num_sites_]==1)
        num_dblocc_sites_++;
    }
  }
}

bool FockBasis::gen_upspin_hop(void)
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (num_upholes_==0 || num_upspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }
  mv_upspin_ = rng_.random_upspin();
  mv_uphole_ = rng_.random_uphole();
  //std::cout << " rng test = " << spin_site_pair.first << "\n";
  up_fr_state_ = up_states_[mv_upspin_]; 
  up_to_state_ = uphole_states_[mv_uphole_]; 
  if (!double_occupancy_ && state_[num_sites_+up_to_state_]) {
    proposed_move_ = move_t::null;
    return false;
  }
  else {
    proposed_move_=move_t::upspin_hop;
    state_[up_fr_state_] = 0;
    state_[up_to_state_] = 1;
    dblocc_increament_ = state_[num_sites_+up_to_state_]; // must be 0 or 1
    dblocc_increament_ -= state_[num_sites_+up_fr_state_];
    //fr_state = up_fr_state_;
    //to_state = up_to_state_;
    return true;
  }
}

bool FockBasis::gen_dnspin_hop(void)
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (num_dnholes_==0 || num_dnspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }
  mv_dnspin_ = rng_.random_dnspin();
  mv_dnhole_ = rng_.random_dnhole();
  //std::cout << " rng test = " << spin_site_pair.first << "\n";
  dn_fr_state_ = dn_states_[mv_dnspin_]; 
  dn_to_state_ = dnhole_states_[mv_dnhole_]; 
  if (!double_occupancy_ && state_[dn_to_state_-num_sites_]) {
    proposed_move_ = move_t::null;
    return false;
  }
  else {
    proposed_move_=move_t::dnspin_hop;
    state_[dn_fr_state_] = 0;
    state_[dn_to_state_] = 1;
    dblocc_increament_ = state_[dn_to_state_-num_sites_]; // must be 0 or 1
    dblocc_increament_ -= state_[dn_fr_state_-num_sites_];
    return true;
  }
}

bool FockBasis::gen_exchange_move(void)
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (num_upholes_==0 || num_upspins_==0) return false;
  if (num_dnholes_==0 || num_dnspins_==0) return false;
  mv_upspin_ = rng_.random_upspin();
  mv_dnspin_ = rng_.random_dnspin();
  up_fr_state_ = up_states_[mv_upspin_]; 
  dn_fr_state_ = dn_states_[mv_dnspin_]; 
  up_to_state_ = dn_fr_state_-num_sites_; 
  dn_to_state_ = num_sites_+up_fr_state_; 
  mv_uphole_ = -1;
  for (int i=0; i<num_upholes_; ++i) {
    if (uphole_states_[i]==up_to_state_) {
      mv_uphole_ = i;
      break;
    }
  }
  if (mv_uphole_<0) return false;
  mv_dnhole_ = -1;
  for (int i=0; i<num_dnholes_; ++i) {
    if (dnhole_states_[i]==dn_to_state_) {
      mv_dnhole_ = i;
      break;
    }
  }
  if (mv_dnhole_<0) return false;
  // valid move
  proposed_move_ = move_t::exchange;
  state_[up_fr_state_] = 0;
  state_[up_to_state_] = 1;
  state_[dn_fr_state_] = 0;
  state_[dn_to_state_] = 1;
  return true;
}

const int& FockBasis::which_upspin(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return mv_upspin_;
  }
  else if (proposed_move_==move_t::exchange) {
    return mv_upspin_;
  }
  else {
    throw std::logic_error("FockBasis::which_upspin: no upspin move exists");
  }
}

const int& FockBasis::which_dnspin(void) const
{
  if (proposed_move_==move_t::dnspin_hop) {
    return mv_dnspin_;
  }
  else if (proposed_move_==move_t::exchange) {
    return mv_dnspin_;
  }
  else {
    throw std::logic_error("FockBasis::which_dnspin: no dnspin move exists");
  }
}

int FockBasis::which_spin(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return mv_upspin_;
  }
  else if (proposed_move_==move_t::dnspin_hop) {
    return num_upspins_+mv_dnspin_;
  }
  else {
    throw std::logic_error("FockBasis::which_spin: no upspin move exists");
  }
}

int FockBasis::which_site(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return up_to_state_;
  }
  else if (proposed_move_==move_t::dnspin_hop) {
    return dn_to_state_-num_sites_;
  }
  else {
    throw std::logic_error("FockBasis::which_site: no existing move");
  }
}

const int& FockBasis::which_state(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return up_to_state_;
  }
  else if (proposed_move_==move_t::dnspin_hop) {
    return dn_to_state_;
  }
  else {
    throw std::logic_error("FockBasis::which_state: no existing move");
  }
}

int FockBasis::op_ni_up(const int& site) const
{
  return state_[site];
}

int FockBasis::op_ni_dn(const int& site) const
{
  return state_[num_sites_+site];
}

int FockBasis::op_ni_updn(const int& site) const
{
  if (state_[site] && state_[num_sites_+site]) return 1;
  else return 0;
}

bool FockBasis::op_cdagc_up(const int& site_i, const int& site_j) const
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (state_[site_i]==0 && state_[site_j]==1) {
    up_fr_state_ = site_j;
    up_to_state_ = site_i;
  }
  else if (state_[site_i]==1 && state_[site_j]==0) {
    up_fr_state_ = site_i;
    up_to_state_ = site_j;
  }
  else return false;
  mv_upspin_ = spin_id_[up_fr_state_];
  op_sign_ = 1;
  dblocc_increament_ = 0;
  if (up_fr_state_==up_to_state_ && state_[up_fr_state_]) return true;
  // actual move now
  proposed_move_ = move_t::upspin_hop;
  state_[up_fr_state_] = 0;
  state_[up_to_state_] = 1;
  // change in no of doubly occupied sites
  dblocc_increament_ = state_[num_sites_+up_to_state_]; // must be 0 or 1
  dblocc_increament_ -= state_[num_sites_+up_fr_state_];
  // sign (considered that the state is aready changed above)
  for (int i=up_to_state_+1; i<up_fr_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  for (int i=up_fr_state_+1; i<up_to_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  return true;
}

bool FockBasis::op_cdagc_dn(const int& site_i, const int& site_j) const
{
  if (proposed_move_!=move_t::null) undo_last_move();
  int idx_i = num_sites_+site_i;
  int idx_j = num_sites_+site_j;
  if (state_[idx_i]==0 && state_[idx_j]==1) {
    dn_fr_state_ = idx_j; 
    dn_to_state_ = idx_i; 
  }
  else if (state_[idx_i]==1 && state_[idx_j]==0) {
    dn_fr_state_ = idx_i;
    dn_to_state_ = idx_j;
  }
  else return false;
  mv_dnspin_ = spin_id_[dn_fr_state_];
  op_sign_ = 1;
  dblocc_increament_ = 0;
  if (dn_fr_state_==dn_to_state_ && state_[dn_fr_state_]) return true;
  // actual move now
  proposed_move_ = move_t::dnspin_hop;
  state_[dn_fr_state_] = 0;
  state_[dn_to_state_] = 1;
  // change in no of doubly occupied sites
  dblocc_increament_ = state_[dn_to_state_-num_sites_]; // must be 0 or 1
  dblocc_increament_ -= state_[dn_fr_state_-num_sites_];
  // sign (considered that the state is aready changed above)
  for (int i=dn_to_state_+1; i<dn_fr_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  for (int i=dn_fr_state_+1; i<dn_to_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  return true;
}

int FockBasis::op_exchange_ud(const int& site_i, const int& site_j) const
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (site_i == site_j) return 1;
  auto* ni_up = &state_[site_i];
  auto* nj_up = &state_[site_j];
  auto* ni_dn = &state_[num_sites_+site_i];
  auto* nj_dn = &state_[num_sites_+site_j];
  if (*ni_up==1 && *nj_up==0 && *ni_dn==0 && *nj_dn==1) {
    *ni_up = 0;
    *nj_up = 1;
    *ni_dn = 1;
    *nj_dn = 0;
    up_fr_state_ = site_i;
    up_to_state_ = site_j;
    dn_fr_state_ = num_sites_+site_j;
    dn_to_state_ = num_sites_+site_i;
  }
  else if (*ni_up==0 && *nj_up==1 && *ni_dn==1 && *nj_dn==0) {
    *ni_up = 1;
    *nj_up = 0;
    *ni_dn = 0;
    *nj_dn = 1;
    up_fr_state_ = site_j;
    up_to_state_ = site_i;
    dn_fr_state_ = num_sites_+site_i;
    dn_to_state_ = num_sites_+site_j;
  }
  else {
    return 0;
  }
  // sign (considered that the state is aready changed above)
  op_sign_ = 1;
  for (int i=up_to_state_; i<up_fr_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  for (int i=up_fr_state_; i<up_to_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  for (int i=dn_to_state_; i<dn_fr_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  for (int i=dn_fr_state_; i<dn_to_state_; ++i) {
    if (state_[i]) op_sign_ = -op_sign_;
  }
  proposed_move_ = move_t::exchange;
  return op_sign_;
}

void FockBasis::commit_last_move(void)
{
  // double occupancy count
  switch (proposed_move_) {
    case move_t::upspin_hop:
      num_dblocc_sites_ += dblocc_increament_;
      //operator[](up_fr_state_) = 0;
      //operator[](up_to_state_) = 1;
      spin_id_[up_fr_state_] = null_id_;
      spin_id_[up_to_state_] = mv_upspin_;
      up_states_[mv_upspin_] = up_to_state_;
      spin_states_[mv_upspin_] = up_to_state_;
      uphole_states_[mv_uphole_] = up_fr_state_;
      proposed_move_ = move_t::null;
      break;
    case move_t::dnspin_hop:
      num_dblocc_sites_ += dblocc_increament_;
      //operator[](dn_fr_state_) = 0;
      //operator[](dn_to_state_) = 1;
      spin_id_[dn_fr_state_] = null_id_;
      spin_id_[dn_to_state_] = mv_dnspin_;
      dn_states_[mv_dnspin_] = dn_to_state_;
      spin_states_[num_upspins_+mv_dnspin_] = dn_to_state_;
      dnhole_states_[mv_dnhole_] = dn_fr_state_;
      proposed_move_ = move_t::null;
      break;
    case move_t::exchange:
      spin_id_[up_fr_state_] = null_id_;
      spin_id_[up_to_state_] = mv_upspin_;
      spin_id_[dn_fr_state_] = null_id_;
      spin_id_[dn_to_state_] = mv_dnspin_;
      up_states_[mv_upspin_] = up_to_state_;
      spin_states_[mv_upspin_] = up_to_state_;
      uphole_states_[mv_uphole_] = up_fr_state_;
      dn_states_[mv_dnspin_] = dn_to_state_;
      spin_states_[num_upspins_+mv_dnspin_] = dn_to_state_;
      dnhole_states_[mv_dnhole_] = dn_fr_state_;
      proposed_move_ = move_t::null;
      break;
    case move_t::null:
      break;
  }
  // check
  /*
  int m = 0;
  int n = 0;
  for (int i=0; i<num_sites_; ++i) m += operator[](i);
  for (int i=num_sites_; i<num_states_; ++i) n += operator[](i);
  if (m!= num_upspins_ || n!= num_dnspins_) {
    throw std::logic_error("FockBasis::commit_last_move");
  }*/
}


void FockBasis::undo_last_move(void) const
{
  // double occupancy count
  switch (proposed_move_) {
    case move_t::upspin_hop:
      state_[up_fr_state_] = 1;
      state_[up_to_state_] = 0;
      break;
    case move_t::dnspin_hop:
      state_[dn_fr_state_] = 1;
      state_[dn_to_state_] = 0;
      break;
    case move_t::exchange:
      state_[up_fr_state_] = 1;
      state_[up_to_state_] = 0;
      state_[dn_fr_state_] = 1;
      state_[dn_to_state_] = 0;
      break;
    case move_t::null:
      break;
  }
  dblocc_increament_ = 0;
  proposed_move_ = move_t::null;
}

std::ostream& operator<<(std::ostream& os, const FockBasis& bs)
{
  os << "state: |" << bs.state_.transpose() << ">\n";
  return os;
}


/*
unsigned Hole::total_num_ = 0;
unsigned SpinUp::total_num_ = 0;
unsigned SpinDn::total_num_ = 0;
unsigned Doublon::total_num_ = 0;
*/

BasisState::BasisState(void) 
{
  set_vaccuum(0);
} 

BasisState::BasisState(const unsigned& num_sites) 
{
  set_vaccuum(num_sites);
} 

BasisState::BasisState(const unsigned& num_sites, const bool& allow_dbl) 
{
  set_vaccuum(num_sites, allow_dbl);
} 

void BasisState::set_vaccuum(const unsigned& num_sites, const bool& allow_dbl) 
{
  num_sites_ = num_sites;
  resize(num_sites_);
  double_occupancy_ = allow_dbl;
  if (double_occupancy_) site_capacity_ = 2;
  else site_capacity_ = 1;
  for (unsigned i=0; i<num_sites_; ++i) {
    operator[](i).put_uphole(i);
    operator[](i).put_dnhole(i);
  }
  upspin_sites_.clear();
  dnspin_sites_.clear();
  uphole_sites_.clear();
  dnhole_sites_.clear();
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
}

void BasisState::clear(void)
{
  upspin_sites_.clear();
  dnspin_sites_.clear();
} 

void BasisState::init_spins(const unsigned& num_upspins, const unsigned& num_dnspins)
{
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  if (num_upspins_>num_sites_ || num_dnspins_>num_sites_)
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  if (!double_occupancy_ && (num_upspins_+num_dnspins_)>num_sites_)
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  num_upholes_ = num_sites_ - num_upspins;
  num_dnholes_ = num_sites_ - num_dnspins;
  // resizing
  upspin_sites_.resize(num_upspins_);
  dnspin_sites_.resize(num_dnspins_);
  uphole_sites_.resize(num_upholes_);
  dnhole_sites_.resize(num_dnholes_);
  // random generator
  // assuming 'num_upspins>0', 'num_dnspins>0'
  rng_.set_upspin_generator(0,num_upspins_-1);
  rng_.set_dnspin_generator(0,num_dnspins_-1);
  // hole numbers may be zero
  int m = std::max(static_cast<int>(num_upholes_),1);
  rng_.set_uphole_generator(0,m-1);
  int n = std::max(static_cast<int>(num_dnholes_),1);
  if (num_dnholes_>0)
    rng_.set_dnhole_generator(0,n-1);
} 

void BasisState::allow_double_occupancy(const bool& allow)
{
  double_occupancy_ = allow;
  if (double_occupancy_) site_capacity_ = 2;
  else site_capacity_ = 1;
} 

void BasisState::set_random(void)
{
  std::vector<unsigned> all_sites(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) all_sites[i] = i;
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);

  // UP spins & holes
  for (unsigned i=0; i<num_upspins_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_upspin(i);
    upspin_sites_[i] = site;
  }
  unsigned uh = 0;
  for (unsigned i=num_upspins_; i<num_sites_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_uphole(uh);
    uphole_sites_[uh] = site;
    uh++;
  }
  // DN spins & holes
  if (double_occupancy_) {
    std::shuffle(all_sites.begin(),all_sites.end(),rng_);
    for (unsigned i=0; i<num_dnspins_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(i);
      dnspin_sites_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=num_dnspins_; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(dh);
      dnhole_sites_[dh] = site;
      dh++;
    }
  }
  else {
    unsigned total_spins = num_upspins_+num_dnspins_;
    // DN spins
    unsigned ds = 0;
    for (unsigned i=num_upspins_; i<total_spins; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(ds);
      dnspin_sites_[ds] = site;
      ds++;
    }
    // DN holes
    for (unsigned i=0; i<num_upspins_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(i);
      dnhole_sites_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=total_spins; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(num_upspins_+dh);
      dnhole_sites_[num_upspins_+dh] = site;
      dh++;
    }
  }
  // in addition, in case of 'no double occupancy' 
  // (put block on singly occupied sites):
  if (!double_occupancy_) {
    for (unsigned i=0; i<num_upspins_; ++i) {
      if (operator[](uphole_sites_[i]).have_dnspin())
        uphole_sites_[i] = -1;
    }
    for (unsigned i=0; i<num_dnspins_; ++i) {
      if (operator[](dnhole_sites_[i]).have_upspin())
        dnhole_sites_[i] = -1;
    }
  }
  // in case there are no holes
  if (num_upholes_==0) {
    uphole_sites_.resize(1);
    upspin_sites_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_sites_.resize(1);
    dnspin_sites_[0] = -1;
  }
  // number of doublely occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (const auto& s : *this) 
      if (s.count()==2) num_dblocc_sites_++;
  }
}

void BasisState::set_custom(void)
{
  std::vector<unsigned> all_sites(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) all_sites[i] = i;
  //std::shuffle(all_sites.begin(),all_sites.end(),rng_);

  // UP spins & holes
  for (unsigned i=0; i<num_upspins_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_upspin(i);
    upspin_sites_[i] = site;
  }
  unsigned uh = 0;
  for (unsigned i=num_upspins_; i<num_sites_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_uphole(uh);
    uphole_sites_[uh] = site;
    uh++;
  }

  // DN spins & holes
  unsigned last_site = num_sites_-1;
  for (unsigned i=0; i<num_dnspins_; ++i) {
    unsigned site = all_sites[last_site-i];
    operator[](site).put_dnspin(i);
    dnspin_sites_[i] = site;
  }
  unsigned dh = 0;
  for (unsigned i=num_dnspins_; i<num_sites_; ++i) {
    unsigned site = all_sites[last_site-i];
    operator[](site).put_dnhole(dh);
    dnhole_sites_[dh] = site;
    dh++;
  }

  // in addition, in case of 'no double occupancy' 
  // (put block on singly occupied sites):
  if (!double_occupancy_) {
    for (unsigned i=0; i<num_upspins_; ++i) {
      if (operator[](uphole_sites_[i]).have_dnspin())
        uphole_sites_[i] = -1;
    }
    for (unsigned i=0; i<num_dnspins_; ++i) {
      if (operator[](dnhole_sites_[i]).have_upspin())
        dnhole_sites_[i] = -1;
    }
  }
  // in case there are no holes
  if (num_upholes_==0) {
    uphole_sites_.resize(1);
    upspin_sites_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_sites_.resize(1);
    dnspin_sites_[0] = -1;
  }
  // number of doublely occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (const auto& s : *this) 
      if (s.count()==2) num_dblocc_sites_++;
  }
}


std::pair<int,int> BasisState::gen_upspin_hop(void)
{
  mv_upspin_ = rng_.random_upspin();
  mv_uphole_ = rng_.random_uphole();
  //std::cout << " rng test = " << spin_site_pair.first << "\n";
  up_tosite_ = uphole_sites_[mv_uphole_]; 
  if (up_tosite_ < 0) { // in case 'no dbl occupancy'
    proposed_move_ = move_t::null;
    dblocc_increament_ = 0;
  }
  else {
    proposed_move_=move_t::upspin_hop;
    // double occupancy count
    dblocc_increament_ = operator[](up_tosite_).count(); // must be 0 or 1
    if (operator[](upspin_sites_[mv_upspin_]).count() == 2)
      dblocc_increament_--;
  }
  return std::make_pair(mv_upspin_,up_tosite_);
}

std::pair<int,int> BasisState::gen_dnspin_hop(void)
{
  mv_dnspin_ = rng_.random_dnspin();
  mv_dnhole_ = rng_.random_dnhole();
  dn_tosite_ = dnhole_sites_[mv_dnhole_]; 
  if (dn_tosite_ < 0) { // in case 'no dbl occupancy'
    proposed_move_ = move_t::null;
    dblocc_increament_ = 0;
  }
  else {
    proposed_move_ = move_t::dnspin_hop;
    // double occupancy count
    dblocc_increament_ = operator[](dn_tosite_).count(); // must be 0 or 1
    if (operator[](dnspin_sites_[mv_dnspin_]).count() == 2)
      dblocc_increament_--;
  }
  return std::make_pair(mv_dnspin_,dn_tosite_);
}

std::pair<int,int> BasisState::gen_exchange_move(void)
{
  mv_upspin_ = rng_.random_upspin();
  mv_dnspin_ = rng_.random_dnspin();
  up_tosite_ = dnspin_sites_[mv_dnspin_];
  dn_tosite_ = upspin_sites_[mv_upspin_];
  if (operator[](up_tosite_).have_upspin()) {
    proposed_move_ = move_t::null;
    up_tosite_ = -1;
    return std::make_pair(mv_upspin_,up_tosite_);
  }
  if (operator[](dn_tosite_).have_dnspin()) {
    proposed_move_ = move_t::null;
    up_tosite_ = -1;
    return std::make_pair(mv_upspin_,up_tosite_);
  }
  proposed_move_=move_t::exchange;
  dblocc_increament_ = 0;
  return std::make_pair(mv_upspin_,up_tosite_);
}

std::pair<int,int> BasisState::exchange_move_part(void)
{
  return std::make_pair(mv_dnspin_,dn_tosite_);
}

void BasisState::accept_last_move(void)
{
  // double occupancy count
  int site;
  SiteState* src_state;
  SiteState* tgt_state;
  num_dblocc_sites_ += dblocc_increament_;
  switch (proposed_move_) {
    case move_t::upspin_hop:
      // source state
      site = upspin_sites_[mv_upspin_];
      src_state = &operator[](site);
      src_state->put_uphole(mv_uphole_);
      // put dblocc block for 'uphole'?
      if (src_state->count()!=site_capacity_)
        uphole_sites_[mv_uphole_] = site;
      else
        uphole_sites_[mv_uphole_] = -1;
      // remove dblocc block for 'dnhole'?
      if (src_state->have_dnhole())
        dnhole_sites_[src_state->dnhole_id()] = site;
      // target state
      tgt_state = &operator[](up_tosite_);
      tgt_state->put_upspin(mv_upspin_);
      upspin_sites_[mv_upspin_] = up_tosite_;
      // put dblocc block for 'dnhole'?
      if (!double_occupancy_) {
        if (tgt_state->have_dnhole())
          dnhole_sites_[tgt_state->dnhole_id()] = -1;
      }
      proposed_move_ = move_t::null;
      break;
    case move_t::dnspin_hop:
      // source state
      site = dnspin_sites_[mv_dnspin_];
      src_state = &operator[](site);
      src_state->put_dnhole(mv_dnhole_);
      // put dblocc block for 'dnhole'?
      if (src_state->count()!=site_capacity_)
        dnhole_sites_[mv_dnhole_] = site;
      else
        dnhole_sites_[mv_dnhole_] = -1;
      // remove dblocc block for 'uphole'?
      if (src_state->have_uphole())
        uphole_sites_[src_state->uphole_id()] = site;
      // target state
      tgt_state = &operator[](dn_tosite_);
      tgt_state->put_dnspin(mv_dnspin_);
      dnspin_sites_[mv_dnspin_] = dn_tosite_;
      // put dblocc block for 'uphole'?
      if (!double_occupancy_) {
        if (tgt_state->have_uphole())
          uphole_sites_[tgt_state->uphole_id()] = -1;
      }
      proposed_move_ = move_t::null;
      break;
    case move_t::exchange:
      // source & target states
      src_state = &operator[](dn_tosite_);
      tgt_state = &operator[](up_tosite_);
      mv_uphole_ = tgt_state->uphole_id();
      mv_dnhole_ = src_state->dnhole_id();
      // spins move
      tgt_state->put_upspin(mv_upspin_);
      upspin_sites_[mv_upspin_] = up_tosite_;
      src_state->put_dnspin(mv_dnspin_);
      dnspin_sites_[mv_dnspin_] = dn_tosite_;
      // holes move
      if (double_occupancy_) {
        src_state->put_uphole(mv_uphole_);
        uphole_sites_[mv_uphole_] = dn_tosite_;
        tgt_state->put_dnhole(mv_dnhole_);
        dnhole_sites_[mv_dnhole_] = up_tosite_;
      }
      else {
        src_state->put_uphole(mv_uphole_);
        uphole_sites_[mv_uphole_] = -1;
        tgt_state->put_dnhole(mv_dnhole_);
        dnhole_sites_[mv_dnhole_] = -1;
      }
      proposed_move_ = move_t::null;
      break;
    case move_t::null:
      break;
  }
}

std::ostream& operator<<(std::ostream& os, const BasisState& bs)
{
  unsigned len = 12;
  unsigned i = 0;
  os << std::string(60,'-') << std::endl;
  os << "Basis state:\n";
  for (const auto& s : bs) {
    os << "(" << s.to_string() << ") ";
    i++;
    if (i==len) {
      os << std::endl; i = 0;
    }
  }
  os << std::endl; 
  os << std::string(60,'-') << std::endl;
  return os;
}



} // end namespace vmc
