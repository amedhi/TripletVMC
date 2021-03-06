/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-17 07:24:49
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <boost/filesystem.hpp>
#include "vmc.h"
#include <Eigen/SVD>

namespace vmc {

VMC::VMC(const input::Parameters& inputs) 
  : graph(inputs) 
  , model(inputs, graph.lattice()) 
  , config(inputs, graph, model)
  , num_varparms_{config.num_varparms()}
{
  // seed random generator
  config.rng().seed(inputs.set_value("rng_seed", 1));
  // mc parameters
  num_sites_ = graph.num_sites();
  num_measure_steps_ = inputs.set_value("measure_steps", 0);
  num_warmup_steps_ = inputs.set_value("warmup_steps", 0);
  min_interval_ = inputs.set_value("min_interval", 0);
  max_interval_ = inputs.set_value("max_interval", 0);
  check_interval_ = std::max(1,num_measure_steps_/10);

  // output directory
  std::string outdir = inputs.set_value("prefix", "");
  prefix_ = "./"+outdir+"/";
  boost::filesystem::path prefix_dir(prefix_);
  boost::filesystem::create_directory(prefix_dir);

  // x-variables (assuming max 1)
  xvar_names_.resize(1);
  xvar_values_.resize(1);
  xvar_names_[0]=inputs.set_value("as_func_of","hole_doping");
  if (xvar_names_[0] == "hole_doping")
    xvar_values_[0] = config.hole_doping();
  else {
    xvar_values_[0] = inputs.set_value(xvar_names_[0],0.0);
  }

  // observables
  make_info_str(inputs);
  observables.headstream() << info_str_.str();
  observables.init(inputs,graph,model,config,prefix_);
  observables.as_functions_of(xvar_names_);
}

int VMC::start(const input::Parameters& inputs, const run_mode& mode, 
  const bool& silent)
{
  run_mode_ = mode;
  bool with_psi_grad = false;
  switch (mode) {
    case run_mode::normal:
      if (observables.energy_grad()) with_psi_grad = true;
      break;
    case run_mode::energy_function:
      observables.switch_off();
      observables.energy().setup(graph,model);
      observables.energy_grad().setup(config);
      break;
    case run_mode::sr_function:
      observables.switch_off();
      observables.energy().setup(graph,model);
      observables.energy_grad().setup(config);
      observables.sr_matrix().setup(graph,config);
      break;
  }
  silent_mode_ = silent;
  if (inputs.have_option_quiet()) silent_mode_ = true;
  // build config
  int stat = config.build(graph, inputs, with_psi_grad);
  // xvariables
  if (xvar_names_[0] == "hole_doping")
    xvar_values_[0] = config.hole_doping();
  else {
    xvar_values_[0] = inputs.set_value(xvar_names_[0],0.0);
  }
  return stat;
}

double VMC::energy_function(const Eigen::VectorXd& varp, Eigen::VectorXd& grad)
{
  // observables to calculate
  bool with_psi_grad;
  if (grad.size()>0) {
    with_psi_grad = true;
    observables.energy_grad().switch_on();
  }
  else {
    with_psi_grad = false;
    observables.energy_grad().switch_off();
  }
  // build the config from the variational parameters
  config.build(graph, varp, with_psi_grad);
  run_simulation();
  if (grad.size()>0) grad = observables.energy_grad().mean_data();
  //for (unsigned i=0; i<num_varparms_; ++i)
  //std::cout << " grad = " << grad.transpose() << "\n";
  //std::cout << " varp = " << x[0] << " " << x[1] << "\n";
  //std::cout << " energy = " << en << "\n";
  return observables.energy().mean_data().sum();
}

int VMC::sr_function(const Eigen::VectorXd& varp, double& en_mean, 
  double& en_stddev, Eigen::VectorXd& grad, Eigen::MatrixXd& sr_matrix, 
  const int& sample_size)
{
  // build the config from the variational parameters
  bool with_psi_grad = true;
  config.build(graph, varp, with_psi_grad);
  // run the simulation
  run_simulation(sample_size);
  // energy
  en_mean = observables.energy().mean();
  en_stddev = observables.energy().stddev();
  // gradient
  grad = observables.energy_grad().mean_data();
  // sr matrix
  observables.sr_matrix().get_matrix(sr_matrix);
  return 0;
}

// simulation after optimization
int VMC::run_simulation(const Eigen::VectorXd& varp)
{
  observables.switch_off();
  observables.energy().switch_on();
  config.build(graph, varp);
  run_simulation();
  return 0;
}

int VMC::run_simulation(const int& sample_size)
{
  // warmup run
  if (!silent_mode_) std::cout << "warming up... " << std::flush;
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
  if (!silent_mode_) std::cout << "done\n";
  // measuring run
  int num_measure_steps = num_measure_steps_;
  if (sample_size>0) num_measure_steps = sample_size;
  int measurement_count = 0;
  int skip_count = min_interval_;
  observables.reset();
  while (measurement_count < num_measure_steps) {
    config.update_state();
    if (skip_count >= min_interval_) {
      if (config.accept_ratio()>0.5 || skip_count==max_interval_) {
        skip_count = 0;
        config.reset_accept_ratio();
        observables.do_measurement(graph,model,config);
        ++measurement_count;
        if (!silent_mode_) print_progress(measurement_count, num_measure_steps);
      }
    }
    skip_count++;
  }
  observables.finalize();
  if (!silent_mode_) {
    std::cout << " simulation done\n";
    config.print_stats();
  }
  return 0;
}


void VMC::print_results(void) 
{
  //if (observables.energy_grad()) //finalize_energy_grad();
  // observables.print(config.wavefunc().hole_doping());
  // observables.print_results(config.vparm_vector());
  // observables.print_results(config.hole_doping());
  observables.print_results(xvar_values_);
  //std::cout << observables.energy().with_statistic() << "\n";
}

void VMC::print_progress(const int& num_measurement, const int& num_measure_steps) const
{
  if (num_measurement%check_interval_==0)
  std::cout<<" measurement = "<< double(100.0*num_measurement)/num_measure_steps<<" %\n";
}

void VMC::make_info_str(const input::Parameters& inputs)
{
  info_str_.clear();
  copyright_msg(info_str_);
  info_str_ << "# "<< inputs.job_id() <<"\n"; 
  info_str_ << model.info_str(); 
  info_str_ << config.info_str(); 
  info_str_ << "# Samples = " << num_measure_steps_;
  info_str_ << ", warmup = " << num_warmup_steps_;
  info_str_ << ", min_interval = " << min_interval_;
  info_str_ << ", max_interval = " << max_interval_ << "\n";
}

void VMC::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc





