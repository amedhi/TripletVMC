/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 17:40:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <Eigen/SVD>
#include <boost/algorithm/string.hpp>
#include "./stochastic_reconf.h"

namespace vmc {

int StochasticReconf::init(const input::Parameters& inputs, const VMC& vmc) 
{
  // problem size
  num_parms_ = vmc.num_varp();
  vparms_.resize(num_parms_);
  lbound_ = vmc.varp_lbound();
  ubound_ = vmc.varp_ubound();
  grad_.resize(num_parms_);
  sr_matrix_.resize(num_parms_,num_parms_);

  // optimization parameters
  int nowarn;
  num_sim_samples_ = inputs.set_value("sr_measure_steps", 1000, nowarn);
  num_opt_samples_ = inputs.set_value("sr_opt_samples", 30, nowarn);
  max_iter_ = inputs.set_value("sr_max_iter", 500, nowarn);
  start_tstep_ = inputs.set_value("sr_start_tstep", 0.05, nowarn);
  mk_series_len_ = inputs.set_value("sr_series_len", 40, nowarn);
  stabilizer_ = inputs.set_value("sr_stabilizer", 1.0E-4, nowarn);
  mk_thresold_ = inputs.set_value("sr_fluctuation_tol", 0.30, nowarn);
  grad_tol_ = inputs.set_value("sr_grad_tol", 0.05, nowarn);
  print_progress_ = inputs.set_value("sr_progress_stdout", false, nowarn);
  print_log_ = inputs.set_value("sr_progress_log", true, nowarn);
  refinement_cycle_ = max_iter_/5;
  if (mk_series_len_ > refinement_cycle_) {
    throw std::domain_error("StochasticReconf::init: sr_series_len > sr_max_iter/5");
  }
  // Mann-Kendall statistic
  mk_statistic_.resize(num_parms_);
  mk_statistic_.set_maxlen(mk_series_len_);
  // optimal parameter values
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  bool replace_mode = true;
  if (mode=="APPEND") replace_mode = false;
  optimal_parms_.init("OptParams", replace_mode);
  optimal_parms_.resize(vmc.num_varp(), vmc.varp_names());
  // observable file header
  std::stringstream heading;
  vmc.copyright_msg(heading);
  vmc.print_info(heading);
  std::vector<std::string> as_funct_of{"x"};
  optimal_parms_.print_heading(heading.rdbuf()->str(), as_funct_of);
  // progress file
  if (print_log_) {
    logfile_.open("log_optimization.txt");
    if (!logfile_.is_open())
      throw std::runtime_error("StochasticReconf::init: file open failed");
    vmc.copyright_msg(logfile_);
    vmc.print_info(logfile_);
    logfile_ << "#" << std::string(72, '-') << std::endl;
    logfile_ << "Stochastic Reconfiguration" << std::endl;
    logfile_ << "max_iter = " << max_iter_ << std::endl;
    logfile_ << "start_tstep = " << start_tstep_ << std::endl;
    logfile_ << "mk_series_len = " << mk_thresold_ << std::endl;
    logfile_ << "stabilizer = " << stabilizer_ << std::endl;
    logfile_ << "grad_tol = " << grad_tol_ << std::endl;
    logfile_ << "fluctuation_tol = " << mk_thresold_ << std::endl;
    logfile_ << "optimization samples = " << num_opt_samples_ << std::endl;
    logfile_ << "#" << std::string(72, '-') << std::endl;
  }
  if (print_progress_) {
    std::cout << "#" << std::string(72, '-') << std::endl;
    std::cout << "Stochastic Reconfiguration" << std::endl;
    std::cout << "max_iter = " << max_iter_ << std::endl;
    std::cout << "start_tstep = " << start_tstep_ << std::endl;
    std::cout << "mk_series_len = " << mk_thresold_ << std::endl;
    std::cout << "stabilizer = " << stabilizer_ << std::endl;
    std::cout << "grad_tol = " << grad_tol_ << std::endl;
    std::cout << "fluctuation_tol = " << mk_thresold_ << std::endl;
    std::cout << "optimization samples = " << num_opt_samples_ << std::endl;
    std::cout << "#" << std::string(72, '-') << std::endl;
  }
  return 0;
}

int StochasticReconf::optimize(VMC& vmc)
{
  // start optimization
  optimal_parms_.reset();
  for (unsigned n=0; n<num_opt_samples_; ++n) {
    //std::cout << " optimal sample = " << n << "\n";
    if (print_log_) {
      logfile_ << "Starting sample " << n << " of " 
        << num_opt_samples_ << " ... " << std::flush;
    }
    if (print_progress_) {
      std::cout << "Starting sample " << n << " of " 
        << num_opt_samples_ << " ... " << std::flush;
    }
    // starting value of variational parameters
    vparms_ = vmc.varp_values();
    //vparms_ = lbound_+(ubound_-lbound_)*vmc.rng().random_real();
    // Stochastic reconfiguration iterations
    mk_statistic_.reset();
    double search_tstep = start_tstep_;
    int mc_samples = num_sim_samples_;
    unsigned iter;
    for (iter=1; iter<=max_iter_; ++iter) {
      double en = vmc.sr_function(vparms_, grad_, sr_matrix_, mc_samples);
      // apply to stabilizer to sr matrix 
      for (unsigned i=0; i<num_parms_; ++i) sr_matrix_(i,i) += stabilizer_;
      //for (unsigned i=0; i<num_parms_; ++i) 
       // sr_matrix_(i,i) += sr_matrix_(i,i) * stabilizer_;

      /*
      // stabilization by truncation of redundant direaction
      // reciprocal conditioning number
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(sr_matrix_,
        Eigen::ComputeFullU);
      Eigen::VectorXd lambda_inv(num_parms_);
      double lambda0_inv = 1.0/svd.singularValues()[0];
      lambda_inv[0] = lambda0_inv;
      unsigned num_kept = num_parms_;
      for (int i=1; i<num_parms_; ++i) {
        double lambdai = svd.singularValues()[i];
        if (lambdai * lambda0_inv < 1.0E-4) {
          num_kept = i; break;
        }
        lambda_inv[i] = 1.0/lambdai;
      }
      Eigen::VectorXd del_x(num_parms_);
      for (unsigned k=0; k<num_parms_; ++k) {
        double isum = 0.0;
        for (unsigned i=0; i<num_parms_; ++i) {
          double jsum = 0.0;
          for (unsigned j=0; j<num_kept; ++j) {
            jsum += lambda_inv[j] * svd.matrixU()(k,j) * svd.matrixU()(i,j);  
          }
          isum += jsum * grad_[i];
        }
        del_x[k] = -search_tstep * isum;
      }
      // update variables
      vparms_ += del_x;
      //std::cout << "singular values\n" << svd.singularValues() << "\n";
      //std::cout << "num_kept \n" << num_kept << "\n";
      //getchar();
      */

      
      // search direction
      Eigen::VectorXd search_dir = sr_matrix_.fullPivLu().solve(-grad_);
      //Eigen::VectorXd search_dir = sr_matrix_.inverse()*(-grad_);
      //getchar();
      // update variables
      vparms_ += search_tstep * search_dir;
      
      //vparms_ += search_tstep * (-grad_);
      // box constraint and max_norm (of components not hitting boundary) 
      vparms_ = lbound_.cwiseMax(vparms_.cwiseMin(ubound_));
      /*double gnorm = std::abs(grad_[0]);
      for (unsigned i=0; i<num_parms_; ++i) {
        double x = vparms_(i);
        double lb = lbound_(i);
        double ub = ubound_(i);
        if (x < lb) vparms_(i) = lb;
        else if (x > ub) vparms_(i) = ub;
        else {
         double norm = std::abs(grad_[i]);
         if (max_norm < norm) max_norm = norm;
        }
      } */
      // add data to Mann-Kendall statistic
      double gnorm = grad_.squaredNorm();
      if (gnorm/grad_.size() < 0.50) mk_statistic_ << vparms_;
      double mk_trend = mk_statistic_.elem_max_trend();
      if (print_progress_) {
        std::ios  state(NULL);
        state.copyfmt(std::cout);
        //double gnorm = grad_.squaredNorm();
        std::cout << std::fixed << std::setprecision(4);
        std::cout << " iter = " << iter << "\n";
        std::cout << " grad = " << grad_.transpose() << "\n";
        //std::cout << " search_dir = " << std::setw(6) << del_x.transpose() << "\n";
        std::cout << " search_dir = " << std::setw(6) << search_dir.transpose() << "\n";
        std::cout << " varp =\n" << vparms_.transpose() << "\n";
        std::cout.copyfmt(state);
        std::cout << " energy = " << en << "\n";
        std::cout << " gnorm = " << gnorm << "\n";
        std::cout << " trend = " << mk_trend << "\n"; 
      }
      // convergence criteria
      if (mk_statistic_.is_full() && mk_trend<mk_thresold_) {
        // converged, add data point to store
        mk_statistic_.get_series_avg(vparms_);
        optimal_parms_ << vparms_;
        // print
        if (print_log_) {
          //double gnorm = grad_.squaredNorm();
          logfile_ << "converged!" << std::endl;
          logfile_ << "iteration = "<< iter << std::endl;
          logfile_ << "MK trend = " << mk_statistic_.elem_max_trend() << "\n";
          logfile_ << "gnorm = " << gnorm << std::endl;
          if (num_parms_<10) {
            logfile_ << "varp = " << vparms_.transpose() << std::endl;
          }
          logfile_ << "energy = " << en << std::endl;
          logfile_ << std::endl;
        }
        if (print_progress_) {
          double gnorm = grad_.squaredNorm();
          std::cout << "converged!" << std::endl;
          std::cout << "iteration = "<< iter << std::endl;
          std::cout << "MK trend = " << mk_statistic_.elem_max_trend() << "\n";
          std::cout << "gnorm = " << gnorm << std::endl;
          if (num_parms_<10) {
            std::cout << "varp = " << vparms_.transpose() << std::endl;
          }
          std::cout << "energy = " << en << std::endl;
          std::cout << std::endl;
        }
        break;
      }
      // refinement if not converged early
      if (iter % refinement_cycle_ == 0) {
        if (print_log_) {
          logfile_ << "next refinement cycle" << std::endl;
        }
        if (print_progress_) {
          std::cout << "next refinement cycle" << std::endl;
        }
        mc_samples *= 2;
        search_tstep *= 0.5;
      }
    }
    if (iter>=max_iter_ && print_log_) {
      logfile_ << "NOT converged):" << std::endl << std::endl;
    }
    if (iter>=max_iter_ && print_progress_) {
      std::cout << "NOT converged):" << std::endl << std::endl;
    }
    // next sample
  }
  if (print_log_) {
    if (optimal_parms_.num_samples()==num_opt_samples_)
      logfile_ << "All done!" << std::endl;
    else if (optimal_parms_.num_samples()>0)
      logfile_ << "Done. Some not converged." << std::endl;
    else
      logfile_ << "None converged):" << std::endl;
  }
  if (print_progress_) {
    if (optimal_parms_.num_samples()==num_opt_samples_)
      std::cout << "All done!" << std::endl;
    else if (optimal_parms_.num_samples()>0)
      std::cout << "Done. Some not converged." << std::endl;
    else
      std::cout << "None converged):" << std::endl;
  }
  logfile_.close();
  // print results
  if (optimal_parms_.num_samples() > 0) {
    std::vector<double> xv({vmc.hole_doping()});
    optimal_parms_.print_result(xv);
    vparms_ = optimal_parms_.mean_data(); 
    return true;
  }
  else {
    std::cout << " StochasticReconf::optimize: optimization failed.\n";
    return false;
  }
}



} // end namespace vmc