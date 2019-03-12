#include "function_test.h"

#include <core/util/visualize_vector.h>
#include <core/util/UTF_extensions.h>

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
  return fmt::format("{:<96} {:<25}  epsilon={:>10}",
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
  return fmt::format("{:<55} \u0394={:>10}",
                     name_var(), get_delta());
}

std::string ValueToVary::summary() const
{
  return fmt::format("{:<15} max(\u0394)={:>10}",
                     name, max_delta);
}

CleverHist ValueToVary::deltas_hist() const
{
  auto ret = CleverHist::make_linear(deltas, std::ceil(std::sqrt(deltas.size())));
  while (!ret.counts.empty() && (ret.counts.back() == 0.0))
  {
    ret.counts.pop_back();
    ret.bins.pop_back();
  }

  return ret;
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

void FunctionTest::survey_grad(DAQuiri::FittableRegion* fittable,
                               DAQuiri::AbstractValue* variable,
                               double step_size, double xmin, double xmax)
{
  size_t chosen_var_idx = variable->index();

  val_proxy.clear();
  val_val.clear();
  chi_sq_norm.clear();
  gradient.clear();
  finite_gradient.clear();
  gradient_delta.clear();
  gradient_ok.clear();

  Eigen::VectorXd variables = fittable->variables();
  Eigen::VectorXd gradients;
  for (double proxy = xmin; proxy < xmax; proxy += step_size)
  {
    variables[chosen_var_idx] = proxy;

    val_proxy.push_back(proxy);
    val_val.push_back(variable->val_at(proxy));

    chi_sq_norm.push_back(fittable->chi_sq_gradient(variables, gradients));
    double analytical_grad = gradients[chosen_var_idx];
    gradient.push_back(analytical_grad);

    optimizer.finite_gradient(fittable, variables, gradients);
    double finite_grad = gradients[chosen_var_idx];
    finite_gradient.push_back(finite_grad);

    gradient_delta.push_back(analytical_grad - finite_grad);
    gradient_ok.push_back(optimizer.check_gradient(fittable, variables));
  }
}

double FunctionTest::check_chi_sq(bool print) const
{
  auto min_chi = std::min_element(chi_sq_norm.begin(), chi_sq_norm.end());
  auto min_chi_i = std::distance(chi_sq_norm.begin(), min_chi);

  if (print)
  {
//      MESSAGE() << "chi" << UTF_superscript(2) << "(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "chi" << UTF_superscript(2) << "(val):\n" << visualize_all(val_val, chi_sq_norm, 100) << "\n";
  }
  MESSAGE() << "min(chi" << UTF_superscript(2) << ")=" << chi_sq_norm[min_chi_i]
            << " at val=" << val_val[min_chi_i]
            << " x=" << val_proxy[min_chi_i] << "\n";

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
    MESSAGE() << "grad(val):\n" << visualize_all(val_val, gradient, 100) << "\n";
  }
  MESSAGE() << "min(|grad|)=" << gradient[grad_i]
            << " at val=" << val_val[grad_i]
            << " x=" << val_proxy[grad_i] << "\n";

  return val_val[grad_i];
}

double FunctionTest::check_gradient_deltas(bool print) const
{
  auto max_gd = std::max_element(gradient_delta.begin(), gradient_delta.end());
  auto max_gd_i = std::distance(gradient_delta.begin(), max_gd);

  if (print)
  {
//      MESSAGE() << "chi" << UTF_superscript(2) << "(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "finite_grad(val):\n" << visualize_all(val_val, finite_gradient, 100) << "\n";
    MESSAGE() << "\u0394grad(val):\n" << visualize_all(val_val, gradient_delta, 100) << "\n";
  }

  for (size_t i = 0; i < gradient_ok.size(); ++i)
    EXPECT_TRUE(gradient_ok[i]) << "Bad grad at " << val_val[i] << "\n";

  MESSAGE() << "max(\u0394grad)=" << gradient_delta[max_gd_i]
            << " at val=" << val_val[max_gd_i]
            << " x=" << val_proxy[max_gd_i] << "\n";

  return gradient_delta[max_gd_i];
}

void FunctionTest::test_fit(size_t attempts,
                            DAQuiri::FittableRegion* fittable,
                            DAQuiri::AbstractValue* variable,
                            double wrong_value,
                            double epsilon)
{
  double goal_val = variable->val();

  fittable->update_indices();

  if (verbose)
    MESSAGE() << "Will fit " << wrong_value << " --> " << variable->to_string()
              << " in " << attempts << " attempts with epsilon=" << epsilon << "\n";

  for (size_t i = 0; i < attempts; ++i)
  {
    variable->val(wrong_value);
    auto result = optimizer.minimize(fittable);
    if (verbose)
      MESSAGE() << result.to_string() << "\n";
    fittable->save_fit(result);
    EXPECT_NEAR(variable->val(), goal_val, epsilon)
              << "Attempt[" << i << "] " << variable->to_string()
              << "  \u0394=" << (goal_val - variable->val()) << "\n";
  }

  variable->val(goal_val);
}

void FunctionTest::deterministic_test(size_t attempts,
                                      DAQuiri::FittableRegion* fittable,
                                      DAQuiri::AbstractValue* variable,
                                      double wrong_value)
{
  double goal_val = variable->val();

  fittable->update_indices();

  std::vector<double> vals;
  std::vector<double> uncerts;
  for (size_t i = 0; i < attempts; ++i)
  {
    variable->val(wrong_value);
    fittable->save_fit(optimizer.minimize(fittable));
    vals.push_back(variable->val());
    uncerts.push_back(variable->uncert());
  }

  double val1 = vals.front();
  bool vals_ok = true;
  for (const auto& v : vals)
  {
    if (v != val1)
      vals_ok = false;
  }

  EXPECT_TRUE(vals_ok);
  if (!vals_ok || verbose)
    MESSAGE() << "Vals: " << print_vector(vals) << "\n";

  double uncert1 = uncerts.front();
  bool uncerts_ok = true;
  for (const auto& v : uncerts)
  {
    if (v != uncert1)
      uncerts_ok = false;
  }

  EXPECT_TRUE(uncerts_ok);
  if (!uncerts_ok || verbose)
    MESSAGE() << "Uncerts: " << print_vector(uncerts) << "\n";

  variable->val(goal_val);
}

void FunctionTest::test_fit_random(size_t attempts,
                                   DAQuiri::FittableRegion* fittable,
                                   ValueToVary var)
{
  test_fit_random(attempts, fittable, std::vector<ValueToVary>{var});
}

void FunctionTest::test_fit_random(size_t attempts,
                                   DAQuiri::FittableRegion* fittable,
                                   std::vector<ValueToVary> vals)
{
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

    auto result = optimizer.minimize(fittable);
    if (optimizer.verbosity)
      MESSAGE() << "        " << result.to_string(optimizer.verbosity) << "\n";

    fittable->save_fit(result);

    if (!fittable->sane())
      not_sane++;
    else if (result.converged)
    {
      max_iterations_to_converge =
          std::max(max_iterations_to_converge, result.iterations);
      max_perturbations_to_converge =
          std::max(max_perturbations_to_converge, result.total_perturbations);
      if (result.total_perturbations > 0)
        converged_perturbed++;
      else if (result.used_finite_grads)
        converged_finite++;
      for (auto& v : vals)
      {
        v.record_delta();
        if (print_outside_tolerance && (v.get_delta() > v.epsilon))
          MESSAGE() << "        Outside tolerance " << result.to_string() << "\n"
                    << "                          F=" <<
                    fittable->to_string("                         ") << "\n";
      }
    }
    else
    {
      if (print_unconverged)
        MESSAGE() << "        attempt:" << i << "  "
                  << result.to_string() << "\n"
                  << "                       F=" <<
                  fittable->to_string("                         ") << "\n";
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
  MESSAGE() << " Unconverged:" << unconverged << " ("
            << std::to_string(double(unconverged) / double(attempts) * 100.0) << "%)\n";
  MESSAGE() << " Not sane:" << not_sane << " ("
            << std::to_string(double(not_sane) / double(attempts) * 100.0) << "%)\n";
  MESSAGE() << " Converged as finite:" << converged_finite << " ("
            << std::to_string(double(converged_finite) / double(attempts) * 100.0) << "%)\n";
  MESSAGE() << " Converged perturbed:" << converged_perturbed << " ("
            << std::to_string(double(converged_perturbed) / double(attempts) * 100.0) << "%)\n";
  MESSAGE() << " Max iterations to converge:" << max_iterations_to_converge << "\n";
  MESSAGE() << " Max perturbations to converge:" << max_perturbations_to_converge << "\n";
  for (const auto& v : vals)
  {
    auto deltas = v.deltas_hist();
    MESSAGE() << " " << v.summary() << "\n"
              << visualize_all(deltas.bins, deltas.counts, 100) << "\n";
  }

  for (auto& v : vals)
    EXPECT_LE(v.max_delta, v.epsilon);

  for (auto& v : vals)
    v.variable->val(v.goal);
}
