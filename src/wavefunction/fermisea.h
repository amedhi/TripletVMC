/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-02-20 12:21:42
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-02-20 12:21:42
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef FERMISEA_H
#define FERMISEA_H

#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {


class Fermisea : public GroundState
{
public:
  Fermisea() : GroundState(true) {}
  Fermisea(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph); 
  ~Fermisea() {} 
  int init(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph);
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  bool noninteracting_mu_{true};
  // ground state
  bool have_TP_symmetry_{true};
  double fermi_energy_;
  double total_energy_;
  bool degeneracy_warning_{false};
  struct kshell_t {unsigned k; unsigned nmin; unsigned nmax;};
  std::vector<kshell_t> kshells_up_;
  std::vector<kshell_t> kshells_dn_;

  // matrices
  ComplexMatrix work_;
  std::vector<ComplexMatrix> phi_k_;

  void construct_groundstate(void);
  void get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, Matrix& psi);
  double get_mf_energy(void);
};


} // end namespace var


#endif