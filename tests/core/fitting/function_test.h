#pragma once

#include "gtest_color_print.h"
#include "clever_hist.h"

#include <core/fitting/fittable_region.h>
#include <core/fitting/optimizers/abstract_optimizer.h>
#include <core/fitting/hypermet/Value.h>

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
 protected:
  std::vector<double> val_proxy;
  std::vector<double> val_val;
  std::vector<double> chi_sq_norm;
  std::vector<double> gradient;

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

  void survey_grad(const DAQuiri::FittableRegion* fittable,
                   DAQuiri::AbstractValue* variable,
                   double step_size = 0.1, double xmin = -M_PI_2, double xmax = M_PI_2);

  double check_chi_sq(bool print) const;

  double check_gradients(bool print) const;

  void deterministic_test(size_t attempts,
                          DAQuiri::AbstractOptimizer* optimizer,
                          DAQuiri::FittableRegion* fittable,
                          DAQuiri::AbstractValue* variable,
                          double wrong_value);

  void test_fit(size_t attempts,
                DAQuiri::AbstractOptimizer* optimizer,
                DAQuiri::FittableRegion* fittable,
                DAQuiri::AbstractValue* variable,
                double wrong_value,
                double epsilon);

  void test_fit_random(size_t attempts,
                       DAQuiri::AbstractOptimizer* optimizer,
                       DAQuiri::FittableRegion* fittable,
                       ValueToVary var);

  void test_fit_random(size_t attempts,
                       DAQuiri::AbstractOptimizer* optimizer,
                       DAQuiri::FittableRegion* fittable,
                       std::vector<ValueToVary> vals);
};
