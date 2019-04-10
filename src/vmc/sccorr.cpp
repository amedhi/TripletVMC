/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:35:10
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:08:00
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./sccorr.h"

namespace vmc {

void SC_Correlation::setup(const lattice::LatticeGraph& graph)
{
  max_dist_ = graph.lattice().size1()/2+1;
  num_bond_types_ = graph.num_bond_types();
  bond_pair_corr_.resize(max_dist_);
  num_symm_pairs_.resize(max_dist_);
  for (int d=0; d<max_dist_; ++d) {
    bond_pair_corr_[d].resize(num_bond_types_,num_bond_types_);
    num_symm_pairs_[d].resize(num_bond_types_,num_bond_types_);
    bond_pair_corr_[d].setZero();
    num_symm_pairs_[d].setZero();
  }
  // all the 'source site' pairs (along 'x'-direction) & their distances
  src_pairs_.clear();
  pair_distance_.clear();
  for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
    auto si = *s;
    for (int d=0; d<max_dist_; ++d) {
      auto sj = graph.translated_site(si, Eigen::Vector3i(d,0,0)); 
      //std::cout << si << "   " << sj << "\n";
      src_pairs_.push_back({si, sj});
      pair_distance_.push_back(d);
    }
  }
  src_pairs_size_ = src_pairs_.size();
  // number of symmetrical bond pairs at a given distance
  lattice::LatticeGraph::out_bond_iterator b1, b1_end, b2, b2_end;
  for (int i=0; i<src_pairs_.size(); ++i) {
    for (std::tie(b1,b1_end)=graph.out_bonds(src_pairs_[i].first); b1!=b1_end; ++b1) {
      for (std::tie(b2,b2_end)=graph.out_bonds(src_pairs_[i].second); b2!=b2_end; ++b2) {
        int d = pair_distance_[i];  
        num_symm_pairs_[d](graph.bond_type(b1),graph.bond_type(b2)) += 1;
      }
    }
  }
  /*for (int d=0; d<max_dist_; ++d) {
    for (int i=0; i<num_bond_types_; ++i) {
      for (int j=0; j<num_bond_types_; ++j) {
        std::cout <<"d= "<<d<<" t1= "<<i<<" t2= "<<j<<" n= "<<
        num_symm_pairs_[d](i,j)<<"\n";
      }
    }
  }*/
  // data size & element names
  std::vector<std::string> elem_names;
  for (int i=0; i<num_bond_types_; ++i)
    for (int j=0; j<num_bond_types_; ++j)
      elem_names.push_back(std::to_string(i)+"-"+std::to_string(j));
  unsigned n = max_dist_*num_bond_types_*num_bond_types_;
  this->resize(n, elem_names);
}

void SC_Correlation::measure(const lattice::LatticeGraph& graph, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  lattice::LatticeGraph::out_bond_iterator b1, b1_end, b2, b2_end;
  for (int d=0; d<max_dist_; ++d) {
    bond_pair_corr_[d].setZero();
  }
  for (int i=0; i<src_pairs_.size(); ++i) {
    int d = pair_distance_[i];  
    for (std::tie(b1,b1_end)=graph.out_bonds(src_pairs_[i].first); b1!=b1_end; ++b1) {
      unsigned i_cdag = graph.source(b1);
      unsigned ia_cdag = graph.target(b1);
      unsigned type_i = graph.bond_type(b1);
      int i_phase = graph.bond_sign(b1);
      for (std::tie(b2,b2_end)=graph.out_bonds(src_pairs_[i].second); b2!=b2_end; ++b2) {
        unsigned j_c = graph.source(b2);
        unsigned ja_c = graph.target(b2);
        unsigned type_j = graph.bond_type(b2);
        int j_phase = graph.bond_sign(b2);
        double term = std::real(config.apply_bondsinglet_hop(i_cdag,ia_cdag,
          i_phase,j_c,ja_c,j_phase));
        // pair corr
        bond_pair_corr_[d](type_i,type_j) += term;
        //std::cout << "d = "<<d << " term = "<<term << "\n"; 
      }
    }
  }
  int n = max_dist_ * num_bond_types_ * num_bond_types_;
  mcdata::data_t corr_data_(n);
  n = 0;
  for (int d=0; d<max_dist_; ++d) {
    for (int i=0; i<num_bond_types_; ++i) {
      for (int j=0; j<num_bond_types_; ++j) {
        int multi = std::max(1,num_symm_pairs_[d](i,j));
        corr_data_[n] = bond_pair_corr_[d](i,j)/multi;
        ++n;
      }
    }
  }
  // add to databin
  *this << corr_data_;
  //std::cout << "SC_Correlation::measure\n"; getchar();
}

void SC_Correlation::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!is_open()) open_file();
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  // total value
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  fs_ << std::setw(6)<<"d";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  /*if (xvars.size()>0) {
    fs_ << "# ";
    for (const auto& p : xvars) fs_<<std::setw(14)<<p.substr(0,14);
    fs_ << std::endl;
  }*/
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void SC_Correlation::print_result(const std::vector<double>& xvals) 
{
  if (!is_on()) return;
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  /*if (xvals.size()>0) {
    fs_ << "#";
    for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::endl;
  }*/
  std::vector<double> odlro;
  int n = 0;
  for (int d=0; d<max_dist_; ++d) {
    if (d>=2) {
      for (const auto& p : xvals) fs_ << std::setw(14) << p;
      fs_ << std::setw(6) << d; 
    }
    for (int i=0; i<num_bond_types_; ++i) {
      for (int j=0; j<num_bond_types_; ++j) {
        if (d>=2) fs_ << MC_Data::result_str(n);
        // value at max distance
        if (d==(max_dist_-1)) {
          odlro.push_back(MC_Data::mean(n));
          odlro.push_back(MC_Data::stddev(n));
        }
        ++n;
      }
    }
    if (d>=2) {
      fs_ << MC_Data::conv_str(0); 
      fs_ << std::endl; 
    }
  }
  // print ODLRO
  fs_ << "#odlro "; 
  for (const auto& p : xvals) fs_ << std::setw(14) << p;
  //fs_ << std::setw(6) << max_dist_-1; 
  for (int i=0; i<odlro.size(); i+=2) {
    fs_ << std::scientific << std::setw(15) << std::sqrt(std::abs(odlro[i]));
    fs_ << std::fixed << std::setw(10) << std::sqrt(std::abs(odlro[i+1]));
  }
  fs_ << MC_Data::conv_str(0); 
  fs_ << std::endl; 

  //fstream() << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::endl;
  fs_ << std::endl;
  fs_ << std::flush;
  close_file();
}

} // end namespace vmc

