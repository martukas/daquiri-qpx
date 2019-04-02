#pragma once

#include "gtest_color_print.h"
#include <core/util/clever_hist.h>

#include <core/fitting/data_model/data_model.h>
#include <core/fitting/parameter/abstract_param.h>

#include <core/fitting/optimizers/optlib_adapter.h>

#include <random>

struct ValueToVary
{
  ValueToVary() = default;
  ValueToVary(std::string var_name, DAQuiri::AbstractValue* var,
              double minimum, double maximum, double eps);

  DAQuiri::AbstractValue* variable;
  double min, max;
  double epsilon;
  std::uniform_real_distribution<double> distribution;
  std::string name;
  double goal;

  double max_delta{0};
  std::vector<double> deltas;

  std::string name_var() const;

  std::string declare() const;

  void vary(std::mt19937& rng);

  void record_delta();

  double get_delta() const;

  std::string print_delta();

  std::string summary() const;

  CleverHist deltas_hist() const;
};

class FunctionTest : public TestBase
{
 public:
  FunctionTest() = default;
  ~FunctionTest() = default;

 protected:
  DAQuiri::OptlibOptimizer optimizer;

  std::vector<double> val_proxy;
  std::vector<double> val_val;
  std::vector<double> chi_sq_norm;
  std::vector<double> gradient;
  std::vector<double> finite_gradient;
  std::vector<double> gradient_delta;
  std::vector<bool> gradient_ok;

  bool verbose {false};
  bool print_outside_tolerance {false};
  bool print_unconverged {true};

  size_t unconverged {0};
  size_t not_sane {0};
  size_t converged_finite {0};
  size_t converged_perturbed {0};
  size_t max_iterations_to_converge{0};
  size_t max_perturbations_to_converge{0};

  DAQuiri::WeightedData generate_data(
      const DAQuiri::FittableRegion* fittable, size_t bins) const;

  void visualize_data(const DAQuiri::WeightedData& data) const;

  void survey_grad(DAQuiri::FittableRegion* fittable,
                   DAQuiri::AbstractValue* variable,
                   double step_size = 0.1, double xmin = -M_PI_2, double xmax = M_PI_2);

  double check_chi_sq(bool print) const;

  double check_gradients(bool print) const;

  double check_gradient_deltas(bool print) const;

  void deterministic_test(size_t attempts,
                          DAQuiri::FittableRegion* fittable,
                          DAQuiri::AbstractValue* variable,
                          double wrong_value);

  void test_fit(size_t attempts,
                DAQuiri::FittableRegion* fittable,
                DAQuiri::AbstractValue* variable,
                double wrong_value,
                double epsilon);

  void test_fit_random(size_t attempts,
                       DAQuiri::FittableRegion* fittable,
                       ValueToVary var);

  void test_fit_random(size_t attempts,
                       DAQuiri::FittableRegion* fittable,
                       std::vector<ValueToVary> vals);
};
