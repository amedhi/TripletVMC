/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:16:02
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-12 21:26:24
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BASISSTATE_H
#define BASISSTATE_H

#include <iostream>
#include <vector>
#include <utility>
#include <bitset>
#include <array>
#include "../wavefunction/matrix.h"
#include "./random.h"

namespace vmc {

//#define CONSERVE_SPIN_NUM

enum class spin {UP, DN};

enum class move_t {spin_hop, exchange, null};

class FockBasis 
{
public:
  FockBasis() { init(0); }
  FockBasis(const int& num_sites, const bool& allow_dbl=true) 
    { init(num_sites, allow_dbl); }
  ~FockBasis() {}
  RandomGenerator& rng(void) const { return rng_; }
  void init(const int& num_sites, const bool& allow_dbl=true);
  void init_spins(const int& num_upspins, const int& num_dnspins);
  void set_double_occupancy(const bool& yesno) { double_occupancy_=yesno; }
  const ivector& state(void) const { return occu_n_; }
  const std::vector<int>& spin_states(void) const { return spin_states_; }
  const bool& double_occupancy(void) const { return double_occupancy_; }
  void set_random(void);
  void set_custom(void);
  bool gen_spin_hop(void);
  bool gen_exchange_move(void);
  void do_exchange_first_move(void);
  void undo_exchange_first_move(void);
  const int& which_spin(void) const;
  const int& which_state(void) const;
  const int& which_second_spin(void) const;
  const int& which_second_state(void) const;
  void commit_last_move(void);
  void undo_last_move(void) const;
  int op_ni_up(const int& site) const;
  int op_ni_dn(const int& site) const;
  int op_ni_updn(const int& site) const;
  bool op_cdagc_up(const int& fr_site, const int& to_site) const;
  bool op_cdagc_dn(const int& fr_site, const int& to_site) const;
  bool op_sisj(const int& site_i, const int& site_j) const;
  int op_exchange_ud(const int& fr_site, const int& to_site) const;
  const int op_sign(void) const { return op_sign_; }
  const int delta_nd(void) const { return delta_nd_; }
  const int dbocc_count(void) const { return num_dblocc_sites_; }
  friend std::ostream& operator<<(std::ostream& os, const FockBasis& bs);
private:
  mutable RandomGenerator rng_;
  mutable ivector occu_n_;
  ivector spin_id_;  // store which UP-spin for a given state index
  int num_sites_{0};
  int num_states_{0};
  int num_upspins_{0};
  int num_dnspins_{0};
  int num_spins_{0};
  int num_holes_{0};
  int num_dblocc_sites_{0};
  bool double_occupancy_{true};
  std::vector<int> spin_states_;
  std::vector<int> hole_states_;
  mutable std::vector<int> dnspin_sites_;
  // update moves
  mutable move_t proposed_move_;
  mutable int delta_nd_{0};
  //move_type accepted_move;
  mutable int mv_spin_;
  mutable int mv_hole_;
  mutable int fr_state_;
  mutable int to_state_;
  mutable int mv_spin2_;
  mutable int mv_hole2_;
  mutable int fr_state2_;
  mutable int to_state2_;
  mutable int op_sign_;
  int null_id_{-1};
  void clear(void); 
};
 




} // end namespace vmc

#endif
