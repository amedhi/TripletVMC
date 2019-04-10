/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-13 11:23:45
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-17 10:45:06
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
//#define BOOST_FILESYSTEM_NO_DEPRECATED
//#include <boost/filesystem.hpp>
#include <cstdlib>
#include "./disorder.h"

namespace vmc {

const bool& SiteDisorder::check(const input::Parameters& inputs)
{
  int nowarn;
  exists_ = inputs.set_value("put_disorder",false,nowarn);
  return exists_;
}

int SiteDisorder::init(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model, const SysConfig& config,
    RandomGenerator& rng)
{
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  num_configs_ = inputs.set_value("disorder_configs",1);
  bandwidth_ = inputs.set_value("disorder_bandwidth",1.0);
  disorder_pot_.resize(num_configs_);
  for (auto& v : disorder_pot_) v.resize(num_sites_);
  num_opt_parms_ = graph.num_sites()+2*graph.num_bonds();
  optimal_parms_.resize(num_opt_parms_);
  // init out-of-range current config
  current_config_ = -1;
  current_config_it_ = disorder_pot_.end();
  // disorder potential file path
  datafile_prefix_ = inputs.set_value("disorder_inputs","./disorder");
  if (datafile_prefix_.back()!='/') datafile_prefix_ += "/";
  datafile_prefix_ += model.signature_str();
  optparm_path_ = config.signature_str();
  // load disorder potential from file, if available, otherwise 
  // generate the potentials
  std::string path = datafile_prefix_+"/cnf";
  for (unsigned i=0; i<num_configs_; ++i) {
    std::ostringstream cnf_id;
    cnf_id.width(3);
    cnf_id.fill('0');
    cnf_id << i;
    std::string config_path(path+cnf_id.str()+"/");
    std::string fname(config_path+"site_disorder.txt");
    std::ifstream fin(fname);
    if (fin.is_open()) {
      double x;
      unsigned n = 0;
      while (fin>>x && n<num_sites_) {
        disorder_pot_[i][n++] = x;
        //std::cout << i << "  " << x << "\n"; getchar();
      }
      fin.close();
      if (n != num_sites_) {
        throw std::range_error("SiteDisorder::init: disorder data size mismatch");
      }
    }
    else {
      //boost::filesystem::path dir(config_path);
      //boost::filesystem::create_directory(dir);
      std::string cmd("mkdir -p "+config_path);
      std::system(cmd.c_str());
      std::ofstream fout(fname);
      //std::cout << fname << "\n";
      if (!fout.is_open()) 
        throw std::runtime_error("SiteDisorder::init: file open failed");
      // random site potential
      fout<<std::scientific<<std::setprecision(12)<<std::uppercase<<std::right;
      for (unsigned n=0; n<num_sites_; ++n) {
        double V = bandwidth_ * (rng.random_real()-0.50);
        //fout << std::setw(6) << n << std::setw(22) << V << std::endl;
        // write to file
        fout << std::setw(22) << V << std::endl;
        // also keep the data
        disorder_pot_[i][n] = V;
        //std::cout << i << "  " << V << "\n";
      }
      fout.close();
    }

    // load optimized parameters if exists
    /*
    fname = config_path+"/"+optparm_path_+"/optimized_parms.txt";
    std::ifstream fopts(fname);
    if (fopts.is_open()) {
      double x;
      fopts.close();
      //if (n != num_sites_) {
      //  throw std::range_error("SiteDisorder::init: disorder data size mismatch");
      //}
    }*/
  }
  //std::cout << "SiteDisorder\n"; getchar();

  return 0;
}
  
void SiteDisorder::set_current_config(const unsigned& current_config)
{
  if (current_config >= num_configs_)
    throw std::range_error("SiteDisorder::set_current_config: out-of-range value");
  current_config_ = current_config;
  current_config_it_ = disorder_pot_.begin()+current_config_;
}

void SiteDisorder::save_optimal_parms(const var::parm_vector& optimal_parms)
{
  if (current_config_<0 || current_config_>= num_configs_)
    throw std::range_error("SiteDisorder::save_optimal_parms: out-of-range current config");
  std::string path = datafile_prefix_+"/cnf";
  std::ostringstream cnf_id;
  cnf_id.width(3);
  cnf_id.fill('0');
  cnf_id << current_config_;
  path += cnf_id.str()+"/";
  path += optparm_path_;
  std::string cmd("mkdir -p "+path);
  std::system(cmd.c_str());
  std::string fname = path+"/optimized_parms.txt";
  std::ofstream fout(fname);
  //std::cout << fname << "\n";
  if (!fout.is_open()) 
    throw std::runtime_error("SiteDisorder::save_optimal_parms: file open failed");
  fout<<std::scientific<<std::setprecision(6)<<std::uppercase<<std::right;
  unsigned t_start = num_sites_;
  unsigned delta_start = num_sites_+num_bonds_;
  fout<<"#"<<std::setw(5)<<"i"<<std::setw(20)<<"t"<<std::setw(20)<<"delta";
  fout<<std::setw(20)<<"mu\n";
  for (unsigned i=0; i<num_bonds_; ++i) {
    fout << std::setw(6) << i;
    // t-parameters
    fout << std::setw(20) << optimal_parms[t_start+i];
    // delta-parameters
    fout << std::setw(20) << optimal_parms[delta_start+i];
    // mu-parameters
    if (i < num_sites_)
      fout << std::setw(20) << optimal_parms[i] << "\n";
    else 
      fout << std::setw(20) << 0.0 << "\n";
  }
  fout.close();
}

const var::parm_vector& SiteDisorder::get_optimal_parms(void) 
{ 
  if (current_config_<0 || current_config_>= num_configs_)
    throw std::range_error("SiteDisorder::save_optimal_parms: out-of-range current config");

  // optimal parameter file name
  std::string path = datafile_prefix_+"/cnf";
  std::ostringstream cnf_id;
  cnf_id.width(3);
  cnf_id.fill('0');
  cnf_id << current_config_;
  path += cnf_id.str()+"/";
  path += optparm_path_;
  std::string fname = path+"/optimized_parms.txt";
  std::ifstream fin(fname);
  if (!fin.is_open()) 
    throw std::runtime_error("SiteDisorder::get_optimal_parms: file open failed");
  // discard the first line
  std::string line;
  std::getline(fin, line);
  int i;
  double mu, t, delta;
  //while (fin >> i >> mu >> t >> delta)
  // other lines
  unsigned m = num_sites_+num_bonds_;
  unsigned num_read = 0;
  while (std::getline(fin, line)) {
    std::istringstream ss(line);
    ss >> i;
    if (i < num_sites_) {
      ss >> t >> delta >> mu;
      optimal_parms_[i] = mu;
      optimal_parms_[num_sites_+i] = t;
      optimal_parms_[m+i] = delta;
      num_read += 3;
    }
    else {
      ss >> t >> delta;
      optimal_parms_[num_sites_+i] = t;
      optimal_parms_[m+i] = delta;
      num_read += 2;
    }
  }
  fin.close();
  if (num_read != optimal_parms_.size())
    throw std::range_error("SiteDisorder::get_optimal_parms: parameter size mismatch");
  /*for (int i=0; i<optimal_parms_.size(); ++i) {
    std::cout << optimal_parms_[i] << "\n"; 
  }*/
  return optimal_parms_;
}

bool SiteDisorder::optimal_parms_exists(const unsigned& config)
{
  // optimal parameter file name
  std::string path = datafile_prefix_+"/cnf";
  std::ostringstream cnf_id;
  cnf_id.width(3);
  cnf_id.fill('0');
  cnf_id << config;
  path += cnf_id.str()+"/";
  path += optparm_path_;
  std::string fname = path+"/optimized_parms.txt";
  std::ifstream fin(fname);
  if (fin.is_open()) {
    fin.close();
    return true; 
  }
  else return false;
}





} // end namespace vmc

