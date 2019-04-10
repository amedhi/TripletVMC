/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:32:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 10:47:01
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef GROUNDSTATE_H
#define GROUNDSTATE_H

#include <vector>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class GroundState
{
public:
  GroundState(const bool& pairing_type)
    : pairing_type_{pairing_type} {}
  virtual ~GroundState() {} 
  virtual void update(const input::Parameters& inputs);
  virtual void update(const var::parm_vector& pvector, const unsigned& start_pos=0);
  virtual void get_wf_amplitudes(Matrix& psi);
  virtual void get_wf_gradient(std::vector<Matrix>& psi_gradient); 
  const VariationalParms& varparms(void) { return varparms_; }
  const unsigned& num_varparms(void) const { return num_varparms_; }
  const bool& pairing_type(void) const { return pairing_type_; }
  const unsigned& num_upspins(void) const { return num_upspins_; }
  const unsigned& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
protected:
  unsigned num_sites_{0};
  unsigned num_bonds_{0};
  unsigned num_kpoints_{0};
  unsigned kblock_dim_{0};
  unsigned num_varparms_{0};
  basis::BlochBasis blochbasis_;
  MF_Model mf_model_;
  VariationalParms varparms_;
  // fourier transform matrix
  ComplexMatrix FTU_;
  // solvers
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_minusk_up;

  void set_pairing_type(const bool& pairing_type) { pairing_type_=pairing_type; }
  void set_particle_num(const input::Parameters& inputs);
  void set_ft_matrix(const lattice::LatticeGraph& graph);
  double get_noninteracting_mu(void);
private:
  bool pairing_type_{false};
  unsigned num_spins_{0};
  unsigned num_upspins_{0};
  unsigned num_dnspins_{0};
  double last_hole_doping_{2.5}; // unlikely input
  double hole_doping_{0.0};
  double band_filling_{1.0};
  //MF_Model mf_model_;
  //basis::BlochBasis blochbasis_;
};


} // end namespace var

#endif