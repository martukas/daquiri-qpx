#include "function_test.h"

#include <core/util/visualize_vector.h>

ValueToVary::ValueToVary(std::string var_name, DAQuiri::AbstractValue* var,
                         double minimum, double maximum, double eps)
{
  name = var_name;
  variable = var;
  min = minimum;
  max = maximum;
  distribution = std::uniform_real_distribution<double>(min, max);
  epsilon = eps;
  goal = variable->val();
}

std::string ValueToVary::name_var() const
{
  return fmt::format("{:<15} {}",
                     name, variable->to_string());
}

std::string ValueToVary::declare() const
{
  auto minmax = fmt::format("dist({},{})", min, max);
  return fmt::format("{:<55} {:<25}  epsilon={:>10}",
                     name_var(), minmax, epsilon);
}

void ValueToVary::vary(std::mt19937& rng)
{
  variable->val(distribution(rng));
}

void ValueToVary::record_delta()
{
  auto delta = get_delta();
  max_delta = std::max(delta, max_delta);
  deltas.push_back(delta);
}

double ValueToVary::get_delta() const
{
  return std::abs(goal - variable->val());
}

std::string ValueToVary::print_delta()
{
  return fmt::format("{:<55} delta={:>10}",
                     name_var(), get_delta());
}

std::string ValueToVary::summary() const
{
  return fmt::format("{:<15} max_delta={:>10}",
                     name, max_delta);
}

CleverHist ValueToVary::deltas_hist() const
{
  return CleverHist::make_linear(deltas, std::ceil(std::sqrt(deltas.size())));
//  return CleverHist::make_linear(deltas, 100);
}

DAQuiri::WeightedData FunctionTest::generate_data(
    const DAQuiri::FittableRegion* fittable, size_t bins) const
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < bins; ++i)
  {
    channels.push_back(i);
    y.push_back(fittable->eval(i));
  }
  return DAQuiri::WeightedData(channels, y);
}

void FunctionTest::visualize_data(const DAQuiri::WeightedData& data) const
{
  std::vector<double> channels;
  std::vector<double> counts;
  for (const auto& p : data.data)
  {
    channels.push_back(p.channel);
    counts.push_back(p.count);
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, counts, 100) << "\n";
}

void FunctionTest::survey_grad(const DAQuiri::FittableRegion* fittable,
                               DAQuiri::AbstractValue* variable,
                               double step_size, double xmin, double xmax)
{
  size_t chosen_var_idx = variable->index();

  val_proxy.clear();
  val_val.clear();
  chi_sq_norm.clear();
  gradient.clear();

  Eigen::VectorXd variables = fittable->variables();
  Eigen::VectorXd gradients;
  for (double proxy = xmin; proxy < xmax; proxy += step_size)
  {
    variables[chosen_var_idx] = proxy;

    val_proxy.push_back(proxy);
    val_val.push_back(variable->val_at(proxy));
    chi_sq_norm.push_back(fittable->chi_sq_gradient(variables, gradients));
    gradient.push_back(gradients[chosen_var_idx]);
  }
}

double FunctionTest::check_chi_sq(bool print) const
{
  auto min_chi = std::min_element(chi_sq_norm.begin(), chi_sq_norm.end());
  auto min_chi_i = std::distance(chi_sq_norm.begin(), min_chi);

  if (print)
  {
//      MESSAGE() << "chi_sq(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "chi_sq(val):\n" << visualize_all(val_val, chi_sq_norm, 100) << "\n";
  }
  MESSAGE() << "min(chi_sq)=" << chi_sq_norm[min_chi_i] << " at val=" << val_val[min_chi_i] << "\n";

  return val_val[min_chi_i];
}

double FunctionTest::check_gradients(bool print) const
{
  double min_abs = std::numeric_limits<double>::max();
  size_t grad_i = 0;
  for (size_t i = 0; i < gradient.size(); ++i)
  {
    if (std::abs(gradient[i]) < min_abs)
    {
      grad_i = i;
      min_abs = std::abs(gradient[i]);
    }
  }

  if (print)
  {
//      MESSAGE() << "gradient(proxy):\n" << visualize(val_proxy, gradient, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize_all(val_val, gradient, 100) << "\n";
  }
  MESSAGE() << "min(abs(grad))=" << gradient[grad_i] << " at val=" << val_val[grad_i] << "\n";

  return val_val[grad_i];
}

void FunctionTest::test_fit(size_t attempts,
                            DAQuiri::AbstractOptimizer* optimizer,
                            DAQuiri::FittableRegion* fittable,
                            DAQuiri::AbstractValue* variable,
                            double wrong_value,
                            double epsilon)
{
  fittable->update_indices();
  double goal_val = variable->val();

  MESSAGE() << "Will fit " << wrong_value << " --> " << variable->to_string()
            << " in " << attempts << " attempts with epsilon=" << epsilon << "\n";

  for (size_t i = 0; i < attempts; ++i)
  {
    variable->val(wrong_value);
    auto result = optimizer->minimize(fittable);
    if (optimizer->verbose)
      MESSAGE() << result.to_string() << "\n";
    fittable->save_fit(result);
    MESSAGE() << "Attempt[" << i << "] " << variable->to_string()
              << "  delta=" << (goal_val - variable->val()) << "\n";
    EXPECT_NEAR(variable->val(), goal_val, epsilon);
  }
}

void FunctionTest::test_fit_random(size_t attempts,
                                   DAQuiri::AbstractOptimizer* optimizer,
                                   DAQuiri::FittableRegion* fittable,
                                   DAQuiri::AbstractValue* variable,
                                   double wrong_min, double wrong_max,
                                   double epsilon)
{
  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());
  std::uniform_real_distribution<double> distribution(wrong_min, wrong_max);

  fittable->update_indices();
  double goal_val = variable->val();

  MESSAGE() << "Will fit random(" << wrong_min << "," << wrong_max << ") --> "
            << variable->to_string() << " in " << attempts
            << " attempts with epsilon=" << epsilon << "\n";

  for (size_t i = 0; i < attempts; ++i)
  {
    double wrong = distribution(random_generator);
    variable->val(wrong);
    auto result = optimizer->minimize(fittable);
    fittable->save_fit(result);
    MESSAGE() << "Attempt[" << i << "] " << wrong << "->" << variable->to_string()
              << "  delta=" << (goal_val - variable->val()) << "\n";
    if (optimizer->verbose)
      MESSAGE() << "     " << result.to_string(true) << "\n";
    EXPECT_NEAR(variable->val(), goal_val, epsilon);
  }
}

void FunctionTest::test_fit_random(size_t attempts,
                     DAQuiri::AbstractOptimizer* optimizer,
                     DAQuiri::FittableRegion* fittable,
                     ValueToVary var, bool verbose)
{
  std::vector<ValueToVary> vals;
  vals.push_back(var);
  test_fit_random(attempts, optimizer, fittable, vals, verbose);
}

void FunctionTest::test_fit_random(size_t attempts,
                                   DAQuiri::AbstractOptimizer* optimizer,
                                   DAQuiri::FittableRegion* fittable,
                                   std::vector<ValueToVary> vals,
                                   bool verbose)
{
  size_t unconverged {0};
  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());

  fittable->update_indices();

  MESSAGE() << "Will fit " << fittable->variable_count << " variables in "
            << attempts << " attempts\n";
  for (const auto& v : vals)
    MESSAGE() << " " << v.declare() << "\n";

  for (size_t i = 0; i < attempts; ++i)
  {
    for (auto& v : vals)
      v.vary(random_generator);

    if (verbose)
    {
      MESSAGE() << "   Attempt[" << i << "]\n";
      for (auto& v : vals)
        MESSAGE() << "     " << v.name_var() << "\n";
    }

    auto result = optimizer->minimize(fittable);
    if (optimizer->verbose)
      MESSAGE() << "        " << result.to_string(optimizer->verbose) << "\n";

    if (result.converged)
    {
      fittable->save_fit(result);
      for (auto& v : vals)
        v.record_delta();
    }
    else
    {
      MESSAGE() << "        " << result.to_string() << "\n";
      unconverged++;
    }

    if (verbose)
    {
      MESSAGE() << "     ----->\n";
      for (auto& v : vals)
        MESSAGE() << "          " << v.print_delta() << "\n";
    }
  }

  MESSAGE() << "Summary:\n";
  MESSAGE() << " Unconverged:" << unconverged << "\n";
  for (const auto& v : vals)
  {
    auto deltas = v.deltas_hist();
    MESSAGE() << " " << v.summary() << "\n"
              << visualize_all(deltas.bins, deltas.counts, 100) << "\n";
  }

  for (auto& v : vals)
    EXPECT_LE(v.max_delta, v.epsilon);

}
