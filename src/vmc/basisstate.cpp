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
  occu_n_.resize(num_states_);
  occu_n_.setZero();
  spin_id_.resize(num_states_);
  spin_id_.setConstant(-1); 
  double_occupancy_ = allow_dbl;
  spin_states_.clear();
  hole_states_.clear();

#ifdef CONSERVE_SPIN_NUM
  up_states_.clear();
  dn_states_.clear();
  uphole_states_.clear();
  dnhole_states_.clear();
#endif

  proposed_move_ = move_t::null;
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
  else rng_.set_site_generator(0,0);
}

void FockBasis::init_spins(const int& num_upspins, const int& num_dnspins)
{
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  num_spins_ = num_upspins_+num_dnspins_;
  num_holes_ = num_states_-num_spins_;
  if (num_spins_>num_states_)
    throw std::range_error("* FockBasis::init_spins: spin number exceeds capacity");
  if (!double_occupancy_ && num_spins_>num_sites_)
    throw std::range_error("* FockBasis::init_spins: spin number exceeds capacity");

#ifdef CONSERVE_SPIN_NUM
  num_upholes_ = num_sites_-num_upspins;
  num_dnholes_ = num_sites_-num_dnspins;
#endif

  spin_states_.resize(num_spins_);
  hole_states_.resize(num_holes_);
  if (num_spins_>0) rng_.set_spin_generator(0,num_spins_-1);
  else rng_.set_spin_generator(0,0);
  if (num_holes_>0) rng_.set_hole_generator(0,num_holes_-1);
  else rng_.set_hole_generator(0,0);

#ifdef CONSERVE_SPIN_NUM
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
#endif
  // random initial configuration
  set_random();
}

void FockBasis::set_random(void)
{
  proposed_move_ = move_t::null;
  occu_n_.setZero();

#ifndef CONSERVE_SPIN_NUM
  std::vector<int> all_states(num_states_);
  for (int i=0; i<num_states_; ++i) all_states[i] = i;
  std::shuffle(all_states.begin(),all_states.end(),rng_);
  if (double_occupancy_) {
    for (int i=0; i<num_spins_; ++i) {
      int state = all_states[i];
      spin_states_[i] = state;
      spin_id_[state] = i;
      occu_n_[state] = 1;
    }
    int j = 0;
    for (int i=num_spins_; i<num_states_; ++i) {
      hole_states_[j++] = all_states[i];
    }
  }
  else {
    // no double occupancy constraint
    int i = 0;
    for (auto& state : all_states) {
      if (state<num_sites_) {
        if (occu_n_[state+num_sites_]) continue;
      } 
      else {
        if (occu_n_[state-num_sites_]) continue;
      }
      spin_states_[i] = state;
      spin_id_[state] = i;
      occu_n_[state] = 1;
      i++;
      if (i == num_spins_) break;
    }
    int j = 0;
    for (auto& state : all_states) {
      if (occu_n_[state]==0) hole_states_[j++] = state;
    }
  }
#else
  // no of UP & DN spins conserved separately
  std::vector<int> all_up_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_up_states[i] = i;
  std::shuffle(all_up_states.begin(),all_up_states.end(),rng_);
  for (int i=0; i<num_upspins_; ++i) {
    int state = all_up_states[i];
    spin_states_[i] = state;
    spin_id_[state] = i;
    occu_n_[state] = 1;
  }
  int j=0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    hole_states_[j++] = all_up_states[i];
  }

  // DN spins & holes
  if (double_occupancy_) {
    std::vector<int> all_dn_states(num_sites_);
    for (int i=0; i<num_sites_; ++i) all_dn_states[i] = num_sites_+i;
    std::shuffle(all_dn_states.begin(),all_dn_states.end(),rng_);
    int j = num_upspins_;
    for (int i=0; i<num_dnspins_; ++i) {
      int state = all_dn_states[i];
      spin_states_[j] = state;
      spin_id_[state] = j;
      occu_n_[state] = 1;
      j++;
    }
    j = num_sites_-num_upspins;
    for (int i=num_dnspins_; i<num_sites_; ++i) {
      hole_states_[j++] = all_dn_states[i];
    }
  }
  else {
    // DN spins
    int j = num_upspins_;
    for (int i=num_upspins_; i<num_spins_; ++i) {
      int state = num_sites_+all_up_states[i];
      spin_states_[j] = state;
      spin_id_[state] = j;
      occu_n_[state] = 1;
      j++;
    }
    // holes
    j = num_sites_-num_upspins;
    for (int i=0; i<num_upspins_; ++i) {
      int state = num_sites_+all_up_states[i];
      hole_states_[j++] = state;
    }
    for (int i=num_spins_; i<num_sites_; ++i) {
      int state = num_sites_+all_up_states[i];
      hole_states_[j++] = state;
    }
  }
#endif
  // number of doubly occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (int i=0; i<num_sites_; ++i) {
      if (occu_n_[i]==1 && occu_n_[i+num_sites_]==1)
        num_dblocc_sites_++;
    }
  }
}

void FockBasis::set_custom(void)
{
  proposed_move_ = move_t::null;
  occu_n_.setZero();
  std::vector<int> all_up_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_up_states[i] = i;
  for (int i=0; i<num_upspins_; ++i) {
    int state = all_up_states[i];
    spin_states_[i] = state;
    spin_id_[state] = i;
    occu_n_[state] = 1;
  }
  int j=0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    hole_states_[j++] = all_up_states[i];
  }
  // DN spins & holes
  std::vector<int> all_dn_states(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_dn_states[i] = num_sites_+i;
  int last_site = num_sites_-1;
  j = num_upspins_;
  for (int i=0; i<num_dnspins_; ++i) {
    int state = all_dn_states[last_site-i];
    spin_states_[j] = state;
    spin_id_[state] = j;
    occu_n_[state] = 1;
    j++;
  }
  j = num_sites_-num_upspins_;
  for (int i=num_dnspins_; i<num_sites_; ++i) {
    hole_states_[j++] = all_dn_states[last_site-i];
  }
  // number of doubly occupied sites
  num_dblocc_sites_ = 0;
  if (double_occupancy_) {
    for (int i=0; i<num_sites_; ++i) {
      if (occu_n_[i]==1 && occu_n_[i+num_sites_]==1) 
        num_dblocc_sites_++;
    }
  }
}

bool FockBasis::gen_spin_hop(void)
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (num_holes_==0 || num_spins_==0) return false;
  mv_spin_ = rng_.random_spin();
  mv_hole_ = rng_.random_hole();
  fr_state_ = spin_states_[mv_spin_]; 
  to_state_ = hole_states_[mv_hole_]; 
  if (!double_occupancy_) {
    if (to_state_<num_sites_) {
      if (occu_n_[to_state_+num_sites_]) return false;
    }
    else {
      if (occu_n_[to_state_-num_sites_]) return false;
    }
  }
  proposed_move_=move_t::spin_hop;
  occu_n_[fr_state_] = 0;
  occu_n_[to_state_] = 1;
  if (to_state_<num_sites_) {
    delta_nd_ = occu_n_[to_state_+num_sites_]; // must be 0 or 1
  }
  else {
    delta_nd_ = occu_n_[to_state_-num_sites_]; // must be 0 or 1
  }
  if (fr_state_<num_sites_) {
    delta_nd_ -= occu_n_[fr_state_+num_sites_]; 
  }
  else {
    delta_nd_ -= occu_n_[fr_state_-num_sites_];
  }
  return true;
}

bool FockBasis::gen_exchange_move(void)
{
  return false;
#ifdef CONSERVE_SPIN_NUM
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
#endif
}

int FockBasis::which_spin(void) const
{
  if (proposed_move_==move_t::spin_hop) {
    return mv_spin_;
  }
  else {
    throw std::logic_error("FockBasis::which_spin: no spin move exists");
  }
}

const int& FockBasis::which_state(void) const
{
  if (proposed_move_==move_t::spin_hop) {
    return to_state_;
  }
  else {
    throw std::logic_error("FockBasis::which_state: no existing move");
  }
}

int FockBasis::op_ni_up(const int& site) const
{
  return occu_n_[site];
}

int FockBasis::op_ni_dn(const int& site) const
{
  return occu_n_[num_sites_+site];
}

int FockBasis::op_ni_updn(const int& site) const
{
  if (occu_n_[site] && occu_n_[num_sites_+site]) return 1;
  else return 0;
}

bool FockBasis::op_cdagc_up(const int& site_i, const int& site_j) const
{
  if (proposed_move_!=move_t::null) undo_last_move();
  if (occu_n_[site_i]==0 && occu_n_[site_j]==1) {
    fr_state_ = site_j;
    to_state_ = site_i;
  }
  else if (occu_n_[site_i]==1 && occu_n_[site_j]==0) {
    fr_state_ = site_i;
    to_state_ = site_j;
  }
  else return false;
  mv_spin_ = spin_id_[fr_state_];
  op_sign_ = 1;
  delta_nd_ = 0;
  if (fr_state_==to_state_ && occu_n_[fr_state_]) return true;
  // actual move now
  proposed_move_ = move_t::spin_hop;
  occu_n_[fr_state_] = 0;
  occu_n_[to_state_] = 1;
  // change in no of doubly occupied sites
  delta_nd_ = occu_n_[num_sites_+to_state_]; // must be 0 or 1
  delta_nd_ -= occu_n_[num_sites_+fr_state_];
  // sign (considered that the state is aready changed above)
  for (int i=to_state_+1; i<fr_state_; ++i) {
    if (occu_n_[i]) op_sign_ = -op_sign_;
  }
  for (int i=fr_state_+1; i<to_state_; ++i) {
    if (occu_n_[i]) op_sign_ = -op_sign_;
  }
  return true;
}

bool FockBasis::op_cdagc_dn(const int& site_i, const int& site_j) const
{
  if (proposed_move_!=move_t::null) undo_last_move();
  int idx_i = num_sites_+site_i;
  int idx_j = num_sites_+site_j;
  if (occu_n_[idx_i]==0 && occu_n_[idx_j]==1) {
    fr_state_ = idx_j; 
    to_state_ = idx_i; 
  }
  else if (occu_n_[idx_i]==1 && occu_n_[idx_j]==0) {
    fr_state_ = idx_i;
    to_state_ = idx_j;
  }
  else return false;
  mv_spin_ = spin_id_[fr_state_];
  op_sign_ = 1;
  delta_nd_ = 0;
  if (fr_state_==to_state_ && occu_n_[fr_state_]) return true;
  // actual move now
  proposed_move_ = move_t::spin_hop;
  occu_n_[fr_state_] = 0;
  occu_n_[to_state_] = 1;
  // change in no of doubly occupied sites
  delta_nd_ = occu_n_[to_state_-num_sites_]; // must be 0 or 1
  delta_nd_ -= occu_n_[fr_state_-num_sites_];
  // sign (considered that the state is aready changed above)
  for (int i=to_state_+1; i<fr_state_; ++i) {
    if (occu_n_[i]) op_sign_ = -op_sign_;
  }
  for (int i=fr_state_+1; i<to_state_; ++i) {
    if (occu_n_[i]) op_sign_ = -op_sign_;
  }
  return true;
}

int FockBasis::op_exchange_ud(const int& site_i, const int& site_j) const
{
  return 0;
#ifdef CONSERVE_SPIN_NUM
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
#endif
}

void FockBasis::commit_last_move(void)
{
  switch (proposed_move_) {
    case move_t::spin_hop:
      num_dblocc_sites_ += delta_nd_;
      spin_id_[fr_state_] = null_id_;
      spin_id_[to_state_] = mv_spin_;
      spin_states_[mv_spin_] = to_state_;
      hole_states_[mv_hole_] = fr_state_;
      if (fr_state_<num_sites_) num_upspins_--;
      else num_dnspins_--;
      if (to_state_<num_sites_) num_upspins_++;
      else num_dnspins_++;
      proposed_move_ = move_t::null;
      break;
    case move_t::exchange:
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
    case move_t::spin_hop:
      occu_n_[fr_state_] = 1;
      occu_n_[to_state_] = 0;
      break;
    case move_t::exchange:
      break;
    case move_t::null:
      break;
  }
  delta_nd_ = 0;
  proposed_move_ = move_t::null;
}

std::ostream& operator<<(std::ostream& os, const FockBasis& bs)
{
  os << "state: |" << bs.occu_n_.transpose() << ">\n";
  return os;
}


} // end namespace vmc
