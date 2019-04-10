/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-01 21:13:21
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-22 22:41:27
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BLOCHBASIS_H
#define BLOCHBASIS_H
#include <iostream>
#include <vector>
#include <set>
#include <stdexcept>
#include "../lattice/graph.h"
#include "matrix.h"

namespace basis {

using basis_state = lattice::LatticeGraph::site_descriptor;
using kpoint = Vector3d;

class BlochBasis : public std::vector<kpoint>
{
public:
  // ctors
  BlochBasis() {}
  BlochBasis(const lattice::LatticeGraph& graph);
  ~BlochBasis() {}
  int construct(const lattice::LatticeGraph& graph);
  const unsigned& num_kpoints(void) const { return num_kpoint_; }
  const unsigned& subspace_dimension(void) const { return subspace_dimension_; }
  kpoint kvector(const unsigned& k) const { return operator[](k); }
  const basis_state& site_state(const unsigned& idx) const 
    { return subspace_basis_[idx]; }
  const unsigned& representative_state_idx(const basis_state& s) const 
    { return representative_state_idx_[s]; }
  //const Vector3d& translation_vector(const basis_state& s) const 
  //  { return translation_vectors_[s]; }
  //const basis_state& representative_state(const unsigned& site);
  /*
  unsigned dimension(void) const { return basis_states.size(); }
  basis_state representative_state(const basis_state& s, const lattice::LatticeGraph& graph,
    Vector3d& R) const;
  basis_state site_basis(const unsigned& idx) const 
    { return basis_states[idx]; }
  unsigned state_idx(const basis_state& s) const 
    {  auto pos = state_indices.find(s); return pos != state_indices.end() ? pos->second : null_index; }
  unsigned null_idx(void) const { return null_index; }
  */
  // friends
private:
  Vector3d b1;
  Vector3d b2;
  Vector3d b3;
  unsigned num_kpoint_;
  unsigned subspace_dimension_;
  //std::vector<Vector3d> translation_vectors_;
  std::vector<basis_state> subspace_basis_;
  std::vector<unsigned> representative_state_idx_;
  unsigned null_idx_;

  // helper functions
  void make_kpoints(const lattice::Lattice& lattice);
  void make_subspace_basis(const lattice::LatticeGraph& graph);
  //void make_site_basis(const lattice::Lattice& lattice, const lattice::LatticeGraph& graph);
};


} // end namespace basis

#endif
