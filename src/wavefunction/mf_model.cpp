/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-30 09:14:19
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_model.h"
#include <boost/algorithm/string.hpp>
#include "../expression/expression.h"
//#include "../xml/pugixml.hpp"

namespace var {

int MF_Model::init(const lattice::Lattice& lattice)
{
  return Model::init(lattice);
}

int MF_Model::finalize(const lattice::LatticeGraph& graph)
{
  Model::finalize(graph.lattice());
  dim_ = graph.lattice().num_basis_sites();
  dim2_ = 2*dim_;
  work_.resize(dim_,dim_);
  work2_.resize(dim_,dim_);
  quadratic_block_up_.resize(dim_,dim_);
  quadratic_block_dn_.resize(dim_,dim_);
  quadratic_block_.resize(dim2_,dim2_);
  pairing_block_.resize(dim2_,dim2_);
  build_unitcell_terms(graph);

  quadratic_block_.setZero();
  pairing_block_.setZero();
  return 0;
}

void MF_Model::update(const input::Parameters& inputs)
{
  Model::update_parameters(inputs);
  update_terms();
}

void MF_Model::update_terms(void)
{
  update_unitcell_terms();
}

void MF_Model::update_site_parameter(const std::string& pname, const double& pvalue)
{
  Model::update_parameter(pname, pvalue);
  for (int i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}

void MF_Model::build_unitcell_terms(const lattice::LatticeGraph& graph)
{
  /*uc_siteterms_.resize(Model::num_siteterms());
  uc_bondterms_.resize(Model::num_bondterms());
  unsigned i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    uc_siteterms_[i].build_siteterm(*sterm, graph);
    i++;
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    uc_bondterms_[i].build_bondterm(*bterm, graph);
    i++;
  }*/
  usite_terms_.resize(Model::num_siteterms());
  ubond_terms_.resize(Model::num_bondterms());
  int i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    usite_terms_[i].build_siteterm(*sterm, graph);
    i++;
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    ubond_terms_[i].build_bondterm(*bterm, graph);
    i++;
  }
}

void MF_Model::construct_kspace_block(const Vector3d& kvec)
{
  //work.setZero(); 
  quadratic_block_up_.setZero();
  quadratic_block_dn_.setZero();
  //quadratic_block_.setZero();
  pairing_block_.setZero();
  //work2 = Matrix::Zero(dim_,dim_);
  // bond terms
  //for (const auto& term : uc_bondterms_) {
  for (const auto& term : ubond_terms_) {
    for (int i=0; i<term.num_out_bonds(); ++i) {
      //----------------HOPPING terms--------------
      Vector3d delta = term.bond_vector(i);
      work_ = term.coeff_matrix(i)*std::exp(ii()*kvec.dot(delta));
      if (term.qn_operator().is_quadratic()) {
        if (term.qn_operator().spin_up()) {
          quadratic_block_up_ += work_;
        }
        if (term.qn_operator().spin_dn()) {
          quadratic_block_dn_ += work_;
        }
      }

      if (term.qn_operator().is_pairing()) {
        //----------------PAIRING terms--------------
        /* 
          In the paring term, sum over 'delta' should include all 
          the nearest neighbours, which is effectively obtained
          by adding the 'adjoint' term. The factor of 0.5 is taken
          as a convention.
        */
        work2_ = 0.5*(work_+work_.adjoint());
        /*
          The quantum number ordering scheme is: band index varies
          first and them spin index.  
        */
        if (term.qn_operator().id()==model::op_id::create_singlet) {
          for (int m=0; m<dim_; ++m) {
            for (int n=0; n<dim_; ++n) {
              pairing_block_(m,dim_+n) += work2_(m,n);
              pairing_block_(dim_+m,n) += -work2_(m,n);
            }
          }
          //pairing_block_.block(0,dim_,dim_,dim_) += term.coeff_matrix(i)*std::cos(kvec.dot(delta)); 
          //pairing_block_.block(dim_,0,dim_,dim_) += -term.coeff_matrix(i)*std::cos(kvec.dot(delta));
          //std::cout << "delta = " << delta.transpose() << "\n";
          //std::cout << "work = \n" << work << "\n";
          //std::cout << "delta_k = \n" << pairing_block_.block(0,dim_,dim_,dim_) << "\n\n";
          //getchar();
        }
        if (term.qn_operator().id()==model::op_id::create_triplet_uu) {
          //pairing_block_.block(0,0,dim_,dim_) += work_;
          pairing_block_.block(0,0,dim_,dim_) += work2_;
        }
        if (term.qn_operator().id()==model::op_id::create_triplet_ud) {
          //pairing_block_.block(0,dim_,dim_,dim_) += work_;
          //pairing_block_.block(dim_,0,dim_,dim_) += work_;
          for (int i=0; i<dim_; ++i) {
            for (int j=0; j<dim_; ++j) {
              pairing_block_(i,dim_+j) += work2_(i,j);
              pairing_block_(i+dim_,j) += work2_(i,j);
            }
          }
        }
        if (term.qn_operator().id()==model::op_id::create_triplet_dd) {
          //pairing_block_.block(dim_,dim_,dim_,dim_) += work_;
          pairing_block_.block(dim_,dim_,dim_,dim_) += work2_;
        }
      }
    }
  }

  // Final quadratic spin-UP block + SITE terms
  work_ = quadratic_block_up_ + quadratic_block_up_.adjoint();
  // site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_up()) {
      work_ += term.coeff_matrix();
    }
  }
  quadratic_block_.block(0,0,dim_,dim_) = work_;
  // quadratic spin-DN block
  work_ = quadratic_block_dn_ + quadratic_block_dn_.adjoint();
  // site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_dn()) {
      work_ += term.coeff_matrix();
    }
  }
  quadratic_block_.block(dim_,dim_,dim_,dim_) = work_;
  //for (const auto& term : uc_siteterms_) {
  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "ek = " << quadratic_block_(0,0) << "\n";
}

void MF_Model::update_unitcell_terms(void)
{
  for (int i=0; i<ubond_terms_.size(); ++i) 
    ubond_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  for (int i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}

/* Write a bond term like,
 H = \sum_{Ia,Jb}c^{\dag}_{Ia} t_{Ia,Jb} c_{Jb}
 for lattices with multiple sites per unit cell as
 H = \sum_{I,delta} Psi^{\dag}_{I} M^{\delta} Psi_{I+delta}
 Assumption: 'site's in the Graph are numbered contigously. For
 sites in the unitcell, the 'sl number' is same as 'uid'.
*/
void UnitcellTerm::build_bondterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  dim_ = graph.lattice().num_basis_sites();
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  // get number of unique 'cell bond vectors'
  num_out_bonds_ = 0;
  for (int i=0; i<dim_; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      int id = graph.vector_id(ei);
      if (id > num_out_bonds_) num_out_bonds_ = id;
    }
  }
  num_out_bonds_++;
  //std::cout << "num_out_bonds_ = " << num_out_bonds_ << "\n";
  bond_vectors_.resize(num_out_bonds_);
  coeff_matrices_.resize(num_out_bonds_);
  for (auto& M : coeff_matrices_) {
    M.resize(dim_, dim_);
    M.setZero();
  }
  expr_matrices_.clear();
  expr_matrices_.resize(num_out_bonds_);
  for (auto& M : expr_matrices_) {
    M.resize(dim_);
    for (int i=0; i<dim_; ++i) M[i].resize(dim_);
  }

  // operator
  op_ = hamterm.qn_operator();
  // build the matrices (for each 'bond vector')
  for (int i=0; i<dim_; ++i) {
    for (std::tie(ei,ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      int id = graph.vector_id(ei);
      auto t = graph.target(ei);
      auto j = graph.site_uid(t);
      int btype = graph.bond_type(ei);
      // expression
      std::string cc_expr(hamterm.coupling_expr(btype));
      boost::trim(cc_expr);
      if (cc_expr.size()>0) {
        if (cc_expr[0]!='-') expr_matrices_[id][i][j] += "+";
        expr_matrices_[id][i][j] += cc_expr;
      }
      // values
      coeff_matrices_[id](i,j) += hamterm.coupling(btype);
      //std::cout << id << " " << coeff_matrices_[id](i,j) << "\n";
      bond_vectors_[id] = graph.vector(ei);
    }
  }
}

void UnitcellTerm::build_siteterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  dim_ = graph.lattice().num_basis_sites();
  num_out_bonds_ = 1; // dummy, no real contrib
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim_,dim_);
  coeff_matrices_[0].setZero();
  expr_matrices_.resize(1);
  expr_matrices_[0].resize(dim_);
  for (int i=0; i<dim_; ++i) expr_matrices_[0][i].resize(dim_);
  // operator
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (int i=0; i<dim_; ++i) {
    int stype = graph.site_type(i);
    coeff_matrices_[0](i,i) = hamterm.coupling(stype);
    // expression
    std::string cc_expr(hamterm.coupling_expr(stype));
    boost::trim(cc_expr);
    if (cc_expr.size()>0) {
      if (cc_expr[0]!='-') expr_matrices_[0][i][i] += "+";
      expr_matrices_[0][i][i] += cc_expr;
    }
  }
  bond_vectors_[0] = Vector3d(0.0,0.0,0.0);
}

void UnitcellTerm::eval_coupling_constant(const model::ModelParams& pvals, const model::ModelParams& cvals)
{
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& p : pvals) {
    vars[p.first] = p.second;
    //std::cout << p.first << " = " << p.second << "\n"; getchar();
  }
  for (const auto& c : cvals) vars[c.first] = c.second;
  try { 
    for (int n=0; n<num_out_bonds_; ++n) {
      for (int i=0; i<dim_; ++i) {
        for (int j=0; j<dim_; ++j) {
          std::string cc_expr(expr_matrices_[n][i][j]);
          if (cc_expr.size()>0) {
            coeff_matrices_[n](i,j) = expr.evaluate(cc_expr, vars); 
            //std::cout << "cc = " << coeff_matrices_[n](i,j) << "\n"; getchar();
          }
          else
            coeff_matrices_[n](i,j) = 0.0;
        }
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "UnitcellTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}

/*void MF_Model::check_xml(void)
{
  std::cout << "Checking XML parser\n";
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("model.xml", pugi::parse_trim_pcdata);
  //std::cout << "Load result: " << result.description() << ", mesh name: " << "\n"; 
  pugi::xml_node model = doc.child("model");
  //std::cout << model.child("parameter").attribute("default").value() << std::endl;
  for (pugi::xml_node p = model.child("parameter"); p; p = p.next_sibling())
  {
    std::cout << "Parameter: ";
    for (pugi::xml_attribute attr = p.first_attribute(); attr; attr = attr.next_attribute())
    {
      std::cout << attr.name() << " = " << attr.as_double() << "\n";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------------------------\n\n";
}*/



} // end namespace var











