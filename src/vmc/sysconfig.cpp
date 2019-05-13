/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-18 14:01:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:16:31
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./sysconfig.h"
#include <Eigen/SVD>

namespace vmc {

SysConfig::SysConfig(const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph, const model::Hamiltonian& model)
  : BasisState(graph.num_sites(), model.double_occupancy())
  , basis_state_(graph.num_sites(), model.double_occupancy())
  , wf_(graph, inputs)
  , pj_(inputs)
  , num_sites_(graph.num_sites())
{
  // variational parameters
  num_pj_parms_ = pj_.varparms().size();
  num_wf_parms_ = wf_.varparms().size();
  num_varparms_ = (num_pj_parms_ + num_wf_parms_);
  vparm_names_.resize(num_varparms_);
  vparm_lbound_.resize(num_varparms_);
  vparm_ubound_.resize(num_varparms_);
  // names
  pj_.get_vparm_names(vparm_names_,0);
  wf_.get_vparm_names(vparm_names_,num_pj_parms_);
  // values are not static and may change
  // bounds
  pj_.get_vparm_lbound(vparm_lbound_,0);
  wf_.get_vparm_lbound(vparm_lbound_,num_pj_parms_);
  pj_.get_vparm_ubound(vparm_ubound_,0);
  wf_.get_vparm_ubound(vparm_ubound_,num_pj_parms_);
}

const var::parm_vector& SysConfig::vparm_values(void) 
{
  // values as 'var::parm_vector'
  vparm_values_.resize(num_varparms_);
  pj_.get_vparm_values(vparm_values_,0);
  wf_.get_vparm_values(vparm_values_,num_pj_parms_);
  return vparm_values_;
}

const std::vector<double>& SysConfig::vparm_vector(void) 
{
  // values as 'std::double'
  vparm_vector_.resize(num_varparms_);
  pj_.get_vparm_vector(vparm_vector_,0);
  wf_.get_vparm_vector(vparm_vector_,num_pj_parms_);
  return vparm_vector_;
}

int SysConfig::build(const lattice::LatticeGraph& graph, const input::Parameters& inputs,
    const bool& with_gradient)
{
  if (num_sites_==0) return -1;
  pj_.update(inputs);
  wf_.compute(graph, inputs, with_gradient);
  return init_config();
}

int SysConfig::build(const lattice::LatticeGraph& graph, const var::parm_vector& pvector,
  const bool& need_psi_grad_)
{
  if (num_sites_==0) return -1;
  pj_.update(pvector, 0);
  unsigned start_pos = pj_.varparms().size();
  wf_.compute(graph, pvector, start_pos, need_psi_grad_);
  return init_config();
}

int SysConfig::init_config(void)
{
  num_upspins_ = wf_.num_upspins();
  num_dnspins_ = wf_.num_dnspins();
  num_spins_ = num_upspins_ + num_dnspins_;
  if (num_upspins_==0 && num_dnspins_==0) return -1;
  if (num_upspins_ != num_dnspins_) 
    throw std::range_error("*SysConfig::init_config: unequal UP & DN spin case not implemented");
  // small 'gfactor' caution
  bool tmp_restriction = false;
  bool original_state = BasisState::double_occupancy();
  if (pj_.have_gutzwiller()) {
    if (pj_.gw_factor()<gfactor_cutoff()) {
      BasisState::allow_double_occupancy(false);
      tmp_restriction = true;
    }
  }
  BasisState::init_spins(num_upspins_, num_dnspins_);
  basis_state_.init_spins(num_upspins_, num_dnspins_);
  psi_mat_.resize(num_spins_,num_spins_);
  psi_inv_.resize(num_spins_,num_spins_);

  // try for a well condictioned amplitude matrix
  //BasisState::set_random();
  basis_state_.set_random();
  int num_attempt = 0;
  while (true) {
    wf_.get_amplitudes(psi_mat_,basis_state_.up_states(), basis_state_.dn_states());
    //std::cout << basis_state_ << "\n"; getchar();
    //std::cout << psi_mat_ << "\n"; getchar();
    // reciprocal conditioning number
    Eigen::JacobiSVD<ComplexMatrix> svd(psi_mat_);
    // reciprocal cond. num = smallest eigenval/largest eigen val
    double rcond = svd.singularValues()(svd.singularValues().size()-1)/svd.singularValues()(0);
    if (std::isnan(rcond)) rcond = 0.0; 
    //std::cout << "rcondition number = "<< rcond << "\n";
    if (rcond>1.0E-15) break;
    // try new basis state
    basis_state_.set_random();
    if (++num_attempt > 1000) {
      throw std::underflow_error("*SysConfig::init: configuration wave function ill conditioned.");
    }
  }
  if (tmp_restriction) allow_double_occupancy(original_state);

  // amplitude matrix invers
  psi_inv_ = psi_mat_.inverse();
  // run parameters
  set_run_parameters();
  return 0;
}

int SysConfig::set_run_parameters(void)
{
  num_updates_ = 0;
  refresh_cycle_ = 100;
  // number of moves per mcstep
  int n_up = static_cast<int>(num_upspins_);
  int n_dn = static_cast<int>(num_dnspins_);
  if (double_occupancy()) {
    num_uphop_moves_ = num_upspins_;
    num_dnhop_moves_ = num_dnspins_;
    num_exchange_moves_ = std::min(n_up, n_dn);
    //num_exchange_moves_ = 2*std::min(n_up, n_dn);
  }
  else {
    int num_holes = num_sites_-(num_upspins_+num_dnspins_);
    num_uphop_moves_ = std::min(n_up,num_holes);
    num_dnhop_moves_ = std::min(n_dn,num_holes);
    num_exchange_moves_ = std::min(n_up, n_dn);
    //num_exchange_moves_ = 4*std::min(n_up, n_dn);
  }
  for (int i=0; i<move_t::end; ++i) {
    num_proposed_moves_[i] = 0;
    num_accepted_moves_[i] = 0;
  }
  last_proposed_moves_ = 1;
  last_accepted_moves_ = 1;

  // work arrays 
  psi_row_.resize(num_spins_);
  psi_col_.resize(num_spins_);
  inv_row_.resize(num_spins_);
  psi_grad_.resize(num_spins_,num_spins_);
  return 0;
}

int SysConfig::update_state(void)
{
  for (int n=0; n<num_uphop_moves_; ++n) do_upspin_hop();
  for (int n=0; n<num_dnhop_moves_; ++n) do_dnspin_hop();
  //for (int n=0; n<num_exchange_moves_; ++n) do_spin_exchange();
  num_updates_++;
  if (num_updates_ % refresh_cycle_ == 0) {
    psi_inv_ = psi_mat_.inverse();
  }
  return 0;
}

int SysConfig::do_upspin_hop(void)
{
  if (basis_state_.gen_upspin_hop()) {
    int upspin = basis_state_.which_upspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_row_,to_site,basis_state_.dnspin_sites());
    amplitude_t det_ratio = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
    if (std::abs(det_ratio) < 1.0E-12) {
      // for safety
      basis_state_.undo_last_move();
      return 0; 
    } 
    auto weight_ratio = det_ratio;
    double transition_proby = std::norm(weight_ratio);
    num_proposed_moves_[move_t::uphop]++;
    last_proposed_moves_++;
    if (basis_state_.rng().random_real()<transition_proby) {
      num_accepted_moves_[move_t::uphop]++;
      last_accepted_moves_++;
      // upddate state
      basis_state_.commit_last_move();
      // update amplitudes
      inv_update_upspin(upspin,psi_row_,det_ratio);
    }
    else {
      basis_state_.undo_last_move();
    }
  } 
  return 0;
}

int SysConfig::do_dnspin_hop(void)
{
  if (basis_state_.gen_dnspin_hop()) {
    int dnspin = basis_state_.which_dnspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_col_, basis_state_.upspin_sites(), to_site);
    amplitude_t det_ratio = psi_col_.cwiseProduct(psi_inv_.row(dnspin)).sum();
    if (std::abs(det_ratio) < 1.0E-12) { // for safety
      basis_state_.undo_last_move();
      return 0; 
    } 
    amplitude_t weight_ratio = det_ratio;
    double transition_proby = std::norm(weight_ratio);
    num_proposed_moves_[move_t::dnhop]++;
    last_proposed_moves_++;
    if (basis_state_.rng().random_real()<transition_proby) {
      num_accepted_moves_[move_t::dnhop]++;
      last_accepted_moves_++;
      // upddate state
      basis_state_.commit_last_move();
      // update amplitudes
      inv_update_dnspin(dnspin,psi_col_,det_ratio);
    }
    else {
      basis_state_.undo_last_move();
    }
  } 
  return 0;
}

/*
int SysConfig::do_spin_exchange(void)
{
  int upspin, up_tosite;
  std::tie(upspin, up_tosite) = gen_exchange_move();
  if (up_tosite < 0) return 0; // valid move not found
  num_proposed_moves_[move_t::exch]++;
  last_proposed_moves_++;
  // for upspin hop forward
  wf_.get_amplitudes(psi_row_, up_tosite, dnspin_sites());
  amplitude_t det_ratio1 = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
  if (std::abs(det_ratio1) < dratio_cutoff()) return 0; // for safety

  // now for dnspin hop backward
  int dnspin, dn_tosite;
  std::tie(dnspin, dn_tosite) = exchange_move_part();
  // new col for this move
  wf_.get_amplitudes(psi_col_, upspin_sites(), dn_tosite);
  // since the upspin should have moved
  wf_.get_amplitudes(psi_col_(upspin), up_tosite, dn_tosite);
  // updated 'dnspin'-th row of psi_inv_
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;
  // elements other than 'upspin'-th
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
    inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
  }
  for (int i=upspin+1; i<static_cast<int>(num_upspins_); ++i) {
    amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
    inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
  }
  inv_row_(upspin) = ratio_inv * psi_inv_(dnspin,upspin);
  // ratio for the dnspin hop
  amplitude_t det_ratio2 = psi_col_.cwiseProduct(inv_row_).sum();
  if (std::abs(det_ratio2) < dratio_cutoff()) return 0; // for safety

  // weight ratio (gw pj_ does not play here)
  amplitude_t weight_ratio = det_ratio1 * det_ratio2;
  double transition_proby = std::norm(weight_ratio);
  //std::cout << "W = " << transition_proby << "\n";
  if (rng().random_real()<transition_proby) {
    num_accepted_moves_[move_t::exch]++;
    last_accepted_moves_++;
    // upddate state
    accept_last_move();
    // update amplitudes
    inv_update_upspin(upspin,psi_row_,det_ratio1);
    inv_update_dnspin(dnspin,psi_col_,det_ratio2);
  }
  return 0;
}
*/

int SysConfig::inv_update_upspin(const int& upspin, const ColVector& psi_row, 
  const amplitude_t& det_ratio)
{
  psi_mat_.row(upspin) = psi_row;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv_.col(i)).sum();
    psi_inv_.col(i) -= beta * psi_inv_.col(upspin);
  }
  for (int i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv_.col(i)).sum();
    psi_inv_.col(i) -= beta * psi_inv_.col(upspin);
  }
  psi_inv_.col(upspin) *= ratio_inv;
  return 0;
}

int SysConfig::inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
  const amplitude_t& det_ratio)
{
  psi_mat_.col(dnspin) = psi_col;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<dnspin; ++i) {
    amplitude_t beta = ratio_inv*psi_col_.cwiseProduct(psi_inv_.row(i)).sum();
    psi_inv_.row(i) -= beta * psi_inv_.row(dnspin);
  }
  for (int i=dnspin+1; i<num_dnspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_col_.cwiseProduct(psi_inv_.row(i)).sum();
    psi_inv_.row(i) -= beta * psi_inv_.row(dnspin);
  }
  psi_inv_.row(dnspin) *= ratio_inv;
  return 0;
}

amplitude_t SysConfig::apply(const model::op::quantum_op& qn_op, const unsigned& site_i, 
    const unsigned& site_j, const int& bc_phase) const
{
  amplitude_t term(0); 
  switch (qn_op.id()) {
    case model::op_id::cdagc_sigma:
      term = apply_upspin_hop(site_i,site_j,bc_phase);
      term+= apply_dnspin_hop(site_i,site_j,bc_phase);
      break;
    case model::op_id::sisj_plus:
      term = apply_sisj_plus(site_i,site_j); 
      break;
    default: 
      throw std::range_error("SysConfig::apply: undefined bond operator.");
  }
  return term;
}

int SysConfig::apply(const model::op::quantum_op& qn_op, const int& i) const
{
  switch (qn_op.id()) {
    case model::op_id::ni_sigma:
      return basis_state_.op_ni_up(i)+basis_state_.op_ni_up(i);
    /*case model::op_id::ni_up:
      return basis_state_.op_ni_up(i);
    case model::op_id::ni_dn:
      return basis_state_.op_ni_dn(i);
    */
    case model::op_id::niup_nidn:
      return basis_state_.op_ni_updn(i);
    default: 
      throw std::range_error("SysConfig::apply: undefined site operator");
  }
}

int SysConfig::apply_niup_nidn(const unsigned& site_i) const
{
  if (operator[](site_i).count()==2) return 1;
  else return 0;
}

amplitude_t SysConfig::apply_upspin_hop(const int& src, const int& tgt,
  const int& bc_phase) const
{
  if (src == tgt) return ampl_part(basis_state_.op_ni_up(src));
  if (basis_state_.op_cdagc_up(src,tgt)) {
    int upspin = basis_state_.which_upspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_row_,to_site,basis_state_.dnspin_sites());
    amplitude_t det_ratio = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
    int delta_nd = basis_state_.delta_nd();
    return amplitude_t(bc_phase) * det_ratio * pj_.gw_ratio(delta_nd);
  }
  else {
    return amplitude_t(0.0);
  }
}

amplitude_t SysConfig::apply_dnspin_hop(const int& src, const int& tgt,
  const int& bc_phase) const
{
  if (src == tgt) return ampl_part(basis_state_.op_ni_dn(src));
  if (basis_state_.op_cdagc_dn(src,tgt)) {
    int dnspin = basis_state_.which_dnspin();
    int to_site = basis_state_.which_site();
    wf_.get_amplitudes(psi_col_,basis_state_.upspin_sites(),to_site);
    amplitude_t det_ratio = psi_col_.cwiseProduct(psi_inv_.row(dnspin)).sum();
    int delta_nd = basis_state_.delta_nd();
    return amplitude_t(bc_phase) * det_ratio * pj_.gw_ratio(delta_nd);
  }
  else {
    return amplitude_t(0.0);
  }
}

amplitude_t SysConfig::apply_sisj_plus(const unsigned& i, const unsigned& j) const
{
/* It evaluates the following operator:
 !   O = (S_i.S_j - (n_i n_j)/4)
 ! The operator can be cast in the form,
 !   O = O_{ud} + O_{du}
 ! where,
 !   O_{ud} = 1/2*(- c^{\dag}_{j\up}c_{i\up} c^{\dag}_{i\dn}c_{j\dn}
 !                 - n_{i\up} n_{j_dn})
 ! O_{du} is obtained from O_{ud} by interchanging spin indices. */

  const SiteState* state_i = &operator[](i);
  const SiteState* state_j = &operator[](j);
  // ni_nj term
  double ninj_term;
  if (state_i->have_upspin() && state_j->have_dnspin()) 
    ninj_term = -0.5;
  else if (state_i->have_dnspin() && state_j->have_upspin()) 
    ninj_term = -0.5;
  else ninj_term = 0.0;

  // spin exchange term
  // if any of the two sites doubly occupied, no exchange possible
  if (state_i->count()==2 || state_j->count()==2) return amplitude_t(ninj_term);

  unsigned upspin, up_tosite;
  unsigned dnspin, dn_tosite;
  if (state_i->have_upspin() && state_j->have_dnspin()) {
    upspin = state_i->upspin_id();
    up_tosite = j;
    dnspin = state_j->dnspin_id();
    dn_tosite = i;
  }
  else if (state_i->have_dnspin() && state_j->have_upspin()) {
    upspin = state_j->upspin_id();
    up_tosite = i;
    dnspin = state_i->dnspin_id();
    dn_tosite = j;
  }
  else return amplitude_t(ninj_term);

  // det_ratio for the term
  wf_.get_amplitudes(psi_row_, up_tosite, dnspin_sites());
  amplitude_t det_ratio1 = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
  // now for dnspin hop 
  wf_.get_amplitudes(psi_col_, upspin_sites(), dn_tosite);
  // since the upspin should have moved
  wf_.get_amplitudes(psi_col_(upspin), up_tosite, dn_tosite);
  // updated 'dnspin'-th row of psi_inv_
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;

  // for safety: if 'det_ratio1 == 0', result is zero
  if (std::isinf(std::abs(ratio_inv))) {
    return amplitude_t(ninj_term);
  }

  // elements other than 'upspin'-th
  for (unsigned i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
    inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
  }
  for (unsigned i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
    inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
  }
  inv_row_(upspin) = ratio_inv * psi_inv_(dnspin,upspin);
  // ratio for the dnspin hop
  amplitude_t det_ratio2 = psi_col_.cwiseProduct(inv_row_).sum();
  amplitude_t det_ratio = ampl_part(std::conj(det_ratio1*det_ratio2));
  /*
  if (std::isnan(det_ratio)) {
    std::cout << std::scientific<< det_ratio1 << "\n\n";
    std::cout << std::scientific<< ratio_inv << "\n\n";
    std::cout << std::scientific<< det_ratio2 << "\n\n";
    std::cout << "NaN detected\n"; getchar();
  }*/
  return -0.5 * det_ratio + amplitude_t(ninj_term);
}

amplitude_t SysConfig::apply_bondsinglet_hop(const unsigned& i_dag, 
  const unsigned& ia_dag, const int& bphase_i, const unsigned& j, 
  const unsigned& jb, const int& bphase_j) const
{
  // Evaluates the following operator:
  //   F_{ab}(i,j) = (c^{\dag}_{i\up}c^{\dag}_{i+a\dn} -  c^{\dag}_{i\dn}c^{\dag}_{i+a,\up})/sqrt(2)
  //          x (c_{j+b\dn}c_{j\up} - c_{j+b\up}c_{j\dn})/sqrt(2)
  //          = 0.5 * [c^{\dag}_{i\up}c_{j\up} x c^{\dag}_{i+a\dn}c_{j+b\dn}
  //                 + c^{\dag}_{i+a\up}c_{j\up} x c^{\dag}_{i\dn}c_{j+b\dn}
  //                 + c^{\dag}_{i+a\up}c_{j+b\up} x c^{\dag}_{i\dn}c_{j\dn}
  //                 + c^{\dag}_{i\up}c_{j+b\up} x c^{\dag}_{i+a\dn}c_{j\dn}]

  int num_terms = 4;
  unsigned up_fromsite, up_tosite;
  unsigned dn_fromsite, dn_tosite;
  amplitude_t net_ratio(0.0);
  for (int iterm=0; iterm<num_terms; ++iterm) {
    switch (iterm) {
      case 0:
        up_fromsite = j; up_tosite = i_dag;
        dn_fromsite = jb; dn_tosite = ia_dag;
        break;
      case 1:
        up_fromsite = j; up_tosite = ia_dag;
        dn_fromsite = jb; dn_tosite = i_dag;
        break;
      case 2:
        up_fromsite = jb; up_tosite = ia_dag;
        dn_fromsite = j; dn_tosite = i_dag;
        break;
      case 3:
        up_fromsite = jb; up_tosite = i_dag;
        dn_fromsite = j; dn_tosite = ia_dag;
        break;
    }
    const SiteState* state_j = &operator[](up_fromsite);
    const SiteState* state_i = &operator[](up_tosite);
    const SiteState* state_jb = &operator[](dn_fromsite);
    const SiteState* state_ia = &operator[](dn_tosite);

    //std::cout << " up_from = " << up_fromsite << "\n";
    //std::cout << " up_to   = " << up_tosite << "\n";
    //std::cout << " dn_from = " << dn_fromsite << "\n";
    //std::cout << " dn_to   = " << dn_tosite << "\n\n";
    // non-zero contribution from term only if followings
    if (!state_j->have_upspin()) continue;
    if (state_i->have_upspin() && up_fromsite!=up_tosite) continue;
    if (!state_jb->have_dnspin()) continue;
    if (state_ia->have_dnspin() && dn_fromsite!=dn_tosite) continue;

    unsigned upspin, dnspin;
    int delta_nd = 0;
    amplitude_t det_ratio1, det_ratio2;

    // first hop the up-spin
    upspin = state_j->upspin_id();
    if (up_fromsite == up_tosite) {
      det_ratio1 = amplitude_t(1.0);
    }
    else {
      wf_.get_amplitudes(psi_row_, up_tosite, dnspin_sites());
      det_ratio1 = psi_row_.cwiseProduct(psi_inv_.col(upspin)).sum();
    }

    // next hop the dn-spin
    dnspin = state_jb->dnspin_id();
    if (dn_fromsite == dn_tosite) {
      det_ratio2 = amplitude_t(1.0);
    }
    else {
      wf_.get_amplitudes(psi_col_, upspin_sites(), dn_tosite);
      // since one upspin have moved
      wf_.get_amplitudes(psi_col_(upspin), up_tosite, dn_tosite);
      // updated 'dnspin'-th row of psi_inv_
      amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;
      // elements other than 'upspin'-th
      for (unsigned i=0; i<upspin; ++i) {
        amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
        inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
      }
      for (unsigned i=upspin+1; i<num_upspins_; ++i) {
        amplitude_t beta = ratio_inv*psi_row_.cwiseProduct(psi_inv_.col(i)).sum();
        inv_row_(i) = psi_inv_(dnspin,i) - beta * psi_inv_(dnspin,upspin);
      }
      inv_row_(upspin) = ratio_inv * psi_inv_(dnspin,upspin);
      // ratio for the dnspin hop
      det_ratio2 = psi_col_.cwiseProduct(inv_row_).sum();
    }
    // net ratio for up & dn spin hop
    amplitude_t det_ratio = ampl_part(std::conj(det_ratio1*det_ratio2));

    if (pj_.have_gutzwiller()) {
      // change in double occupancy for up-spin hop
      delta_nd = 0;
      if (up_fromsite != up_tosite) {
        delta_nd = state_i->count(); // must be 0 or one
        if (state_j->count()==2) delta_nd--;
      }
      // for dn-spin hop
      if (dn_fromsite != dn_tosite) {
        delta_nd += state_ia->count(); 
        if (state_jb->count()==2) delta_nd--;
      }
      det_ratio *= pj_.gw_ratio(delta_nd);
    }
    // contribution from this term
    net_ratio += det_ratio;
    //std::cout << " det_ratio1 = " << det_ratio1 << "\n";
    //std::cout << " det_ratio2 = " << det_ratio2 << "\n";
    //std::cout << " delta_nd = " << delta_nd << "\n";
    //getchar();
  }

  int bc_phase = bphase_i * bphase_j;
  return 0.5 * bc_phase * net_ratio;
}

void SysConfig::get_grad_logpsi(RealVector& grad_logpsi) const
{
  // grad_logpsi wrt pj_ parameters
  int p = pj_.varparms().size();
  for (int n=0; n<p; ++n) {
    if (pj_.varparms()[n].name()=="gfactor") {
      double g = pj_.varparms()[n].value();
      grad_logpsi(n) = static_cast<double>(basis_state_.dbocc_count())/g;
    }
    else {
      throw std::range_error("SysConfig::get_grad_logpsi: this pj_ parameter not implemented\n");
    }
  }
  // grad_logpsi wrt wf_ parameters
  for (int n=0; n<wf_.varparms().size(); ++n) {
    wf_.get_gradients(psi_grad_,n,basis_state_.upspin_sites(), basis_state_.dnspin_sites());
    grad_logpsi(p+n) = std::real(psi_grad_.cwiseProduct(psi_inv_.transpose()).sum());
  }
}

double SysConfig::accept_ratio(void)
{
  // acceptance ratio wrt particle number
  return static_cast<double>(last_accepted_moves_)/
         static_cast<double>(num_upspins_+num_dnspins_); 
  //return static_cast<double>(last_accepted_moves_)/
  //       static_cast<double>(last_proposed_moves_); 
}

void SysConfig::reset_accept_ratio(void)
{
  last_proposed_moves_ = 0;
  last_accepted_moves_ = 0;
}

void SysConfig::print_stats(std::ostream& os) const
{
  long proposed_hops = num_proposed_moves_[move_t::uphop]
                         + num_proposed_moves_[move_t::dnhop];
  long proposed_exch = num_proposed_moves_[move_t::exch];
  long accepted_hops = num_accepted_moves_[move_t::uphop] 
                     + num_accepted_moves_[move_t::dnhop];
  long accepted_exch = num_accepted_moves_[move_t::exch];
  double accept_ratio = 100.0*double(accepted_hops+accepted_exch)/(proposed_hops+proposed_exch);
  double hop_ratio = double(100.0*accepted_hops)/(proposed_hops);
  double exch_ratio = double(100.0*accepted_exch)/(proposed_exch);
  os << "--------------------------------------\n";
  os << " total mcsteps = " << num_updates_ <<"\n";
  os << " total accepted moves = " << (accepted_hops+accepted_exch)<<"\n";
  os << " acceptance ratio = " << accept_ratio << " %\n";
  os << " hopping = " << hop_ratio << " %\n";
  os << " exchange = " << exch_ratio << " %\n";
  os << "--------------------------------------\n";
}


} // end namespace vmc


