/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-01 21:13:27
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-22 22:41:24
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "blochbasis.h"

namespace basis {

BlochBasis::BlochBasis(const lattice::LatticeGraph& graph) 
{
  construct(graph);
}

int BlochBasis::construct(const lattice::LatticeGraph& graph)
{
  /*
  if (disorded_system) {
    // no translational symmetry
    // only k=0 point
    num_kpoint_ = 1;
    this->clear();
    this->push_back(Vector3d(0.0,0.0,0.0));
    // subspace basis
    subspace_dimension_ = graph.num_sites();
    subspace_basis_.resize(subspace_dimension_);
    for (unsigned i=0; i<subspace_dimension_; ++i) {
      basis_state s = graph.site(i);
      subspace_basis_[i] = s; 
    }
    null_idx_ = subspace_basis_.size();
    // index of the 'representative state' of a site
    representative_state_idx_.resize(graph.num_sites());
    for (unsigned i=0; i<graph.num_sites(); ++i) {
      representative_state_idx_[i] = i;
    }
    return 0;
  }
  */
  make_kpoints(graph.lattice());
  make_subspace_basis(graph);
  return 0;
}

void BlochBasis::make_kpoints(const lattice::Lattice& lattice)
{
  Vector3d a1, a2, a3;
  double v;
  using bc = lattice::boundary_type;

  // reciprocal lattice vectors 
  a1 = lattice.vector_a1();
  a2 = lattice.vector_a2();
  a3 = lattice.vector_a3();
  b1 = Vector3d(0.0,0.0,0.0);
  b2 = Vector3d(0.0,0.0,0.0);
  b3 = Vector3d(0.0,0.0,0.0);

  unsigned symmetry_type = 0;
  if (lattice.bc1() == bc::periodic) {
    b1 = two_pi() * a1 / a1.dot(a1); 
    symmetry_type = symmetry_type + 1;
  } 

  if (lattice.bc2() == bc::periodic) {
    switch (symmetry_type) {
      case 0:
        b2 = two_pi() * a2 / a2.dot(a2); 
        break;
      case 1:
        a3 = a1.cross(a2);
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b2 = two_pi() * a3.cross(a1) / v;
        break;
      default: break;
    }
    symmetry_type = symmetry_type + 2;
  } 

  if (lattice.bc3() == bc::periodic) {
    switch (symmetry_type) {
      case 0:
        b3 = two_pi() * a3 / a3.dot(a3); 
        break;
      case 1:
        a2 = a3.cross(a1);
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      case 2:
        a1 = a2.cross(a3);
        v = a1.dot(a2.cross(a3));
        b2 = two_pi() * a3.cross(a1) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      case 3:
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b2 = two_pi() * a3.cross(a1) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      default: break;
    }
  }

  // antiperiodic boundary condition
  Vector3d antipb_shift(0.0, 0.0, 0.0);
  if (lattice.bc1_periodicity()==bc::antiperiodic) antipb_shift(0) = 0.5/lattice.size1();
  if (lattice.bc2_periodicity()==bc::antiperiodic) antipb_shift(1) = 0.5/lattice.size2();
  if (lattice.bc3_periodicity()==bc::antiperiodic) antipb_shift(2) = 0.5/lattice.size3();

  // k-points & translation vectors
  double x1, x2, x3;
  Vector3i n = {0,0,0};
  Vector3i m = {-lattice.size1()/2, -lattice.size2()/2, -lattice.size3()/2};
  this->clear();
  num_kpoint_ = lattice.num_unitcells();
  for (unsigned i=0; i<num_kpoint_; i++) {
    x1 = static_cast<double>(m(0)+n(0))/lattice.size1() + antipb_shift(0);
    x2 = static_cast<double>(m(1)+n(1))/lattice.size2() + antipb_shift(1);
    x3 = static_cast<double>(m(2)+n(2))/lattice.size3() + antipb_shift(2);
    this->push_back(x1*b1 + x2*b2 + x3*b3);
    //auto kvec = x1 * b1 + x2 * b2 + x3 * b3;
    //std::cout << i << ": " << kvec.transpose() << "\n";
    //translation_vectors.push_back(n);
    n = lattice.get_next_bravindex(n);
  }
}

void BlochBasis::make_subspace_basis(const lattice::LatticeGraph& graph)
{
  subspace_dimension_ = graph.lattice().num_basis_sites();
  subspace_basis_.resize(subspace_dimension_);
  for (unsigned i=0; i<subspace_dimension_; ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    if (s != graph.site(uid))
      throw std::logic_error("BlochBasis::make_site_basis: unexpected graph property.");
    subspace_basis_[i] = s; 
  }
  null_idx_ = subspace_basis_.size();
  // index of the 'representative state' of a site
  representative_state_idx_.resize(graph.num_sites());
  //translation_vectors_.resize(graph.num_sites());
  for (auto& idx : representative_state_idx_) idx = null_idx_;
  for (unsigned i=0; i<graph.num_sites(); ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    representative_state_idx_[s] = uid;
    //translation_vectors_[s] = graph.site_cellcord(i);
    //std::cout << translation_vectors_[s] << "\n"; getchar();
  }
}


#ifdef ON
int BlochBasis::construct(const lattice::Lattice& lattice, const lattice::LatticeGraph& graph)
{
  // reset
  kpoints.resize(lattice.num_unitcells());
  //translation_vectors.clear();
  basis_states.resize(lattice.num_basis_sites());
  state_indices.clear();
  //null_state = lattice.num_sites();
  make_bloch_vectors(lattice);
  make_site_basis(lattice, graph);

  return 0;
}


void BlochBasis::make_site_basis(const lattice::Lattice& lattice, const lattice::LatticeGraph& graph)
{
  /*
  * The basis states are ! combination of site states connected by lattice 
  * translational symmetry. Among these states in a 'bloch cycle', only the state
  * (site) with smallest 'id' (representative state=rs) is stored. 
  * The basis states are to be regarded as: 
  *   |s> = 1/sqrt(N) \sum_R e^{ik.R} |rs> 
  */

  // basis states & indexing
  lattice::vertex_iterator vi, vi_end;
  basis_state s;
  unsigned i;
  for (i=0; i<basis_states.size(); ++i) {
    s = graph.vertex(i);
    //std::cout << graph.vertex_uid(i) << "\n";
    if (graph.vertex_id(s) != i) {
      throw std::logic_error("*error: blochbasis:: unexpected graph property");
    }
    basis_states[i] = s;
    state_indices.insert({s,i});
    //std::cout << state_index[s] << "\n";
  }
  null_index = basis_states.size();
}

basis_state BlochBasis::representative_state(const basis_state& s, const lattice::LatticeGraph& graph,
    Vector3d& R) const
{
  unsigned rs = graph.vertex_uid(s);
  R = graph.vertex_cellcord(s);
  return graph.vertex(rs);
}

#endif

} // end namespace basis
