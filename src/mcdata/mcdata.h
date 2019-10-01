/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-24 08:44:27
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-07 22:25:28
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef NEW_MCDATA_H
#define NEW_MCDATA_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

namespace mcdata {

using data_t = Eigen::Array<double,Eigen::Dynamic,1>;
using scalardata_t = Eigen::Array<double,1,1>;

/*----------------------DataBin class------------------*/
class DataBin 
{
public:
  DataBin(const int& size=1, const bool& no_error_bar=false);
  ~DataBin() {}
  void clear(void);
  void resize(const int& size);
  bool add_sample(const data_t& sample);
  bool has_samples(void) const { return (num_samples_ > 0); }
  bool has_carry_over(void) const { return !waiting_sample_exist_; }
  const int& size(void) const { return size_; }
  const int& num_samples(void) const { return num_samples_; }
  const data_t& carry(void) const { return carry_; } 
  bool have_new_samples(void) const { return num_samples_!=num_samples_last_; }
  void finalize(void) const;
  const data_t& mean(void) const { finalize(); return mean_; }
  const data_t& stddev(void) const { finalize(); return stddev_; }
private:
  int size_;
  int num_samples_{0};
  mutable int num_samples_last_{0};
  data_t ssum_;
  data_t sumsq_;
  data_t carry_;
  data_t waiting_sample_;
  bool waiting_sample_exist_{false};
  bool no_error_bar_{false};
  data_t MinusOne_;
  mutable data_t mean_;
  mutable data_t stddev_;
};

/*----------------------mcdata class------------------*/
class MC_Data : private std::vector<DataBin>
{
public:
  MC_Data() {}
  MC_Data(const std::string& name, const int& size=1, 
    const bool& no_error_bar=false) { init(name,size,no_error_bar); }
  ~MC_Data() {}
  virtual void init(const std::string& name, const int& size=1, 
    const bool& no_error_bar=false);
  virtual void resize(const int& size);
  void error_bar_off(void);
  void clear(void);
  void add_sample(const data_t& sample);
  void add_sample(const double& sample);
  void operator<<(const data_t& sample);
  void operator<<(const double& sample);
  const int& num_samples(void) const { return top_bin->num_samples(); }
  void finalize(void) const;
  const std::string& name(void) const { return name_; }
  unsigned size(void) const { return mean_.size(); }
  const data_t& mean_data(void) const; 
  double mean(void) const; 
  const double& mean(const int& n) const; 
  const data_t& stddev_data(void) const; 
  double stddev(void) const; // { return top_bin->stddev(); } 
  //const data_t& tau(void) const { finalize(); return tau_; } 
  const double& stddev(const int& n) const;
  const double& tau(void) const;
  std::string result_str(const int& n=0) const; 
  std::string conv_str(const int& n=0) const; 
  const MC_Data& with_statistic(void) const { show_statistic_=true; return *this; }
  void show_statistic(std::ostream& os=std::cout) const;
  friend std::ostream& operator<<(std::ostream& os, const MC_Data& obs);
private:
  std::string name_;
  int max_binlevel_default_ = 20;
  static const int good_sample_size_ = 30;
  std::vector<DataBin>::iterator top_bin;
  std::vector<DataBin>::iterator end_bin;
  mutable int dcorr_level_;
  mutable data_t mean_;
  mutable data_t stddev_;
  mutable double tau_;
  mutable bool show_statistic_;
  mutable std::string error_converged_;
  mutable std::string convergence_str_;

  void find_conv_and_tau(const unsigned& n=0) const;
  void check_convergence(const std::vector<double>& xv, const std::vector<double>& yv) const;
};



} // end namespace mcdata

#endif
