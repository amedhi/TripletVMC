/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-02-20 12:21:42
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-05-24 22:45:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <numeric>
#include "./fermisea.h"

namespace var {

Fermisea::Fermisea(const input::Parameters& inputs, const lattice::LatticeGraph& graph) 
  : GroundState(true)
{
  init(inputs, graph);
}

int Fermisea::init(const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph)
{
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);

  // build MF Hamiltonian
  varparms_.clear();
  mf_model_.init(graph.lattice());
  std::string name;
  double defval;
  using namespace model;
  model::CouplingConstant cc;
  mf_model_.add_parameter(name="t", defval=1.0, inputs);
  mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
  // chemical potential
  noninteracting_mu_ = true;
  // finalize MF Hamiltonian
  mf_model_.finalize(graph);
  num_varparms_ = varparms_.size();

  // bloch basis
  blochbasis_.construct(graph);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(graph);
  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
  } 
  return 0;
}

void Fermisea::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
  // update MF model
  mf_model_.update(inputs);
}

void Fermisea::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
}

void Fermisea::get_wf_amplitudes(Matrix& psi) 
{
  construct_groundstate();
  get_pair_amplitudes(phi_k_);
  get_pair_amplitudes_sitebasis(phi_k_, psi);
}

void Fermisea::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  for (auto& mat : psi_gradient) mat.setZero();
}

void Fermisea::get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k)
{
  for (int i=0; i<kshells_up_.size(); ++i) {
    int k = kshells_up_[i].k;
    int m = kshells_up_[i].nmax+1;
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_block_up());
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_block_up());
    phi_k[k] = es_k_up.eigenvectors().block(0,0,kblock_dim_,m)
      		 * es_minusk_up.eigenvectors().conjugate().block(0,0,m,kblock_dim_);
    //std::cout << kvec.transpose() << "\n"; 
    //std::cout << phi_k[k] << "\n"; getchar();
  }
}

void Fermisea::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  int p = 0;
  for (int i=0; i<num_kpoints_; ++i) {
    int q = 0;
    for (int j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (int ks=0; ks<kshells_up_.size(); ++ks) {
    	int k = kshells_up_[ks].k;
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // copy transformed block
      //psi.block(p,q,kblock_dim_,kblock_dim_) = 
      for (int m=0; m<kblock_dim_; ++m) {
        for (int n=0; n<kblock_dim_; ++n) {
          psi(p+m,q+n) = ampl_part(work_(m,n));
      	}
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  //std::cout << psi << "\n";
  //getchar();
}

double Fermisea::get_mf_energy(void)
{
  double mf_energy = 0.0;
  return mf_energy/num_sites_;
}

void Fermisea::construct_groundstate(void)
{
  int num_upspins_ = num_upspins();
  int num_dnspins_ = num_dnspins();
  if (have_TP_symmetry_ || num_upspins_!=num_dnspins_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
    std::vector<std::pair<int,int>> qn_list; // list of (k,n)
    std::vector<double> ek;
    for (int k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      mf_model_.construct_kspace_block(kvec);
      es_k_up.compute(mf_model_.quadratic_block_up(), Eigen::EigenvaluesOnly);
      ek.insert(ek.end(),es_k_up.eigenvalues().data(),
        es_k_up.eigenvalues().data()+kblock_dim_);
      //std::cout << kvec.transpose() << " " << es_k_up_.eigenvalues() << "\n"; getchar();
      for (int n=0; n<kblock_dim_; ++n) {
        qn_list.push_back({k, n});
      }
    }
    //std::sort(ek.begin(),ek.end());
    //for (const auto& e : ek) std::cout << e << "\n";

    // Indices in the original ek-array is sorted according to increasing ek
    std::vector<int> idx(ek.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),[&ek](const int& i1, const int& i2) 
      {return ek[i1]<ek[i2];});
    /*for (int i=0; i<ek.size(); ++i) {
      //std::cout << i << "  " << idx[i] << "  " << ek[idx[i]] << "\n";
      std::cout << i+1 << " " << blochbasis_.kvector(idx[i]).transpose()  << "\n";
    }
    getchar();
    */
    
    // mean energy
    /*double e0 = 0.0;
    for (int i=0; i<num_upspins_; ++i) {
      e0 += ek[idx[i]];
    }
    e0 = 2 * e0 / num_sites_;
    std::cout << "e0 = " << e0 << "\n"; getchar();
    */

    // check for degeneracy 
    double degeneracy_tol = 1.0E-12;
    int top_filled_level = num_upspins()-1;
    int num_degen_states = 1;
    int num_valence_particle = 1;
    fermi_energy_ = ek[idx[top_filled_level]];
    total_energy_ = 0.0;
    for (int i=0; i<=top_filled_level; ++i) {
      total_energy_ += ek[idx[i]];
    }
    total_energy_ = 2.0*total_energy_/num_sites_;
    std::cout << "Total KE =" << total_energy_ << "\n";

    // look upward in energy
    for (int i=top_filled_level+1; i<ek.size(); ++i) {
      if (std::abs(fermi_energy_-ek[idx[i]])>degeneracy_tol) break;
      num_degen_states++;
    }
    // look downward in energy
    if (num_degen_states>1) {
      for (int i=top_filled_level-1; i>=0; --i) {
        if (std::abs(fermi_energy_ - ek[idx[i]])>degeneracy_tol) break;
        num_degen_states++;
        num_valence_particle++;
      }
      // warn
      if (degeneracy_warning_) {
        std::cout << ">> warning: groundstate degeneracy: " << 2*num_valence_particle <<
          " particles in " << 2*num_degen_states << " states." << "\n";
      }
    }
    /* 
      Filled k-shells. A k-shell is a group of energy levels having same 
      value of quantum number k.
    */
    // find 'nmax' values of filled k-shells
    std::vector<int> shell_nmax(num_kpoints_);
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    int k, n;
    for (int i=0; i<num_upspins_; ++i) {
      int state = idx[i]; 
      std::tie(k,n) = qn_list[state];
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    // store the filled k-shells
    kshells_up_.clear();
    for (unsigned k=0; k<num_kpoints_; ++k) {
      int nmax = shell_nmax[k];
      if (nmax != -1) kshells_up_.push_back({k,0,static_cast<unsigned>(nmax)});
    }
    
    /*
    for (unsigned k=0; k<kshells_up_.size(); ++k) {
      std::cout << kshells_up_[k].k << " " << kshells_up_[k].nmin << "  "
          << kshells_up_[k].nmax << "\n";
    }
    getchar();
    */
    
  }
  else {
    /* 
      No T.P (Time Reversal * Inversion) symmetry. 
      Need to consider both UP and DN states.
    */
    throw std::range_error("Fermisea::construct_groundstate: case not implemented\n");
  }
}



} // end namespace var
