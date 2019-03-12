#include "../function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Peak.h>

#include <core/fitting/optimizers/optlib_adapter.h>

class FittablePeak : public DAQuiri::FittableRegion
{
  std::uniform_real_distribution<double> x_dist {-M_PI_2, M_PI_2};

 public:
  DAQuiri::Peak peak;

  void update_indices() override
  {
    variable_count = 0;
    peak.update_indices(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    peak.put(ret);
    return ret;
  }

  bool perturb(std::mt19937& rng) override
  {
    if (peak.width.valid_index())
      peak.width.x(x_dist(rng));
    if (peak.position.valid_index())
      peak.position.x(x_dist(rng));
    if (peak.amplitude.valid_index())
      peak.amplitude.x(peak.amplitude.x() + x_dist(rng));

    if (peak.short_tail.amplitude.valid_index())
      peak.short_tail.amplitude.x(x_dist(rng));
    if (peak.short_tail.slope.valid_index())
      peak.short_tail.slope.x(x_dist(rng));

    if (peak.right_tail.amplitude.valid_index())
      peak.right_tail.amplitude.x(x_dist(rng));
    if (peak.right_tail.slope.valid_index())
      peak.right_tail.slope.x(x_dist(rng));

    if (peak.long_tail.amplitude.valid_index())
      peak.long_tail.amplitude.x(x_dist(rng));
    if (peak.long_tail.slope.valid_index())
      peak.long_tail.slope.x(x_dist(rng));

    if (peak.step.amplitude.valid_index())
      peak.step.amplitude.x(x_dist(rng));

    return true;
  }

  bool sane() const override
  {
    return peak.sane();
  }

  double eval(double chan) const override
  {
    return peak.eval(chan).all();
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return peak.eval_at(chan, fit).all();
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    return peak.eval_grad_at(chan, fit, grads).all();
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    peak.get(result.variables);

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();

    auto chi = chi_sq();
    peak.get_uncerts(diags, chi);
  }

  std::string to_string(std::string prepend = "") const override
  {
    Eigen::VectorXd grads;
    auto x = variables();
    chi_sq_gradient(x, grads);
    std::stringstream ss;

    ss << peak.to_string() << "\n";
    ss << prepend << "  chi2=" + std::to_string(chi_sq());
    ss << "  grads=" << grads.transpose() << "\n";
    return ss.str();
  }
};

class Peak : public FunctionTest
{
 protected:
  FittablePeak fp;
  size_t region_size{100};
  size_t random_samples{500};
  
  void SetUp() override
  {
//    optimizer.verbosity = 5;
    optimizer.maximum_iterations = 1000;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;
    optimizer.use_epsilon_check = false;
    optimizer.min_g_norm = 1e-7;
    optimizer.perform_sanity_checks = false;
    optimizer.maximum_perturbations = 0;

    fp.peak = fp.peak.gaussian_only();
    //fp.peak.amplitude.bound(0, 500);
    fp.peak.amplitude.val(40000);
    fp.peak.position.bound(44, 68);
    fp.peak.position.val(51);
    fp.peak.width_override = true;
    fp.peak.width.bound(0.8, 5.0);
    fp.peak.width.val(3.2);
  }
};

TEST_F(Peak, CheckSetup)
{
  MESSAGE() << "Peak:\n" << fp.peak.to_string() << "\n";
}

TEST_F(Peak, Visualize)
{
  auto data = generate_data(&fp, region_size);
  visualize_data(data);
}

TEST_F(Peak, WithinBounds)
{
  auto data = generate_data(&fp, region_size);
  EXPECT_NEAR(data.count_min, 0.0, 1e-14);
  EXPECT_NEAR(data.count_max, 40000.0, 1e-7);
}

TEST_F(Peak, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(fp.peak.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(fp.peak.update_indices(i));
}

TEST_F(Peak, UpdateIndex)
{
  int32_t i = 0;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 3);
  EXPECT_EQ(fp.peak.amplitude.index(), 4);
  EXPECT_EQ(fp.peak.width.index(), 5);
  EXPECT_EQ(i, 6);

  i = 42;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 42);
  EXPECT_EQ(fp.peak.amplitude.index(), 43);
  EXPECT_EQ(fp.peak.width.index(), 44);
  EXPECT_EQ(i, 45);
}

TEST_F(Peak, UpdateIndexInvalidates)
{
  int32_t i = 0;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fp.peak.position.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), -1);
  EXPECT_EQ(fp.peak.amplitude.index(), 3);
  EXPECT_EQ(fp.peak.width.index(), 4);
  EXPECT_EQ(i, 5);

  fp.peak.position.to_fit = true;
  fp.peak.amplitude.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 5);
  EXPECT_EQ(fp.peak.amplitude.index(), -1);
  EXPECT_EQ(fp.peak.width.index(), 6);
  EXPECT_EQ(i, 7);

  fp.peak.position.to_fit = true;
  fp.peak.amplitude.to_fit = true;
  fp.peak.width.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 7);
  EXPECT_EQ(fp.peak.amplitude.index(), 8);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 9);

  fp.peak.position.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.peak.width.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), -1);
  EXPECT_EQ(fp.peak.amplitude.index(), -1);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 9);
}

TEST_F(Peak, UpdateIndexDisabled)
{
  int32_t i = 0;

  fp.peak.width_override = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 2);

  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 2);
  EXPECT_EQ(fp.peak.amplitude.index(), 3);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 4);

  // \todo test resetting of indices
}

TEST_F(Peak, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 1.0);

  fp.peak.put(fit);
  EXPECT_EQ(fit[0], 1.0);
  EXPECT_NE(fit[0], fp.peak.position.x());
  EXPECT_EQ(fit[1], 1.0);
  EXPECT_NE(fit[1], fp.peak.amplitude.x());
  EXPECT_EQ(fit[2], 1.0);
  EXPECT_NE(fit[2], fp.peak.width.x());

  fp.peak.update_indices(fp.variable_count);
  fp.peak.put(fit);
  EXPECT_NE(fit[0], 1.0);
  EXPECT_EQ(fit[0], fp.peak.position.x());
  EXPECT_NE(fit[1], 1.0);
  EXPECT_EQ(fit[1], fp.peak.amplitude.x());
  EXPECT_NE(fit[2], 1.0);
  EXPECT_EQ(fit[2], fp.peak.width.x());
}

TEST_F(Peak, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);
  fit[0] = 0.5;
  fit[1] = 0.03;
  fit[2] = 0.01;

  fp.peak.get(fit);
  EXPECT_NEAR(fp.peak.position.val(), 51, 0.00001);
  EXPECT_NE(fp.peak.position.val(), fp.peak.position.val_at(0.5));
  EXPECT_NEAR(fp.peak.amplitude.val(), 40000, 0.00001);
  EXPECT_NE(fp.peak.amplitude.val(), fp.peak.amplitude.val_at(0.03));
  EXPECT_NEAR(fp.peak.width.val(), 3.2, 0.00001);
  EXPECT_NE(fp.peak.width.val(), fp.peak.width.val_at(0.01));

  fp.peak.update_indices(fp.variable_count);

  fp.peak.get(fit);
  EXPECT_EQ(fp.peak.position.val(), fp.peak.position.val_at(0.5));
  EXPECT_EQ(fp.peak.amplitude.val(), fp.peak.amplitude.val_at(0.03));
  EXPECT_EQ(fp.peak.width.val(), fp.peak.width.val_at(0.01));
}

TEST_F(Peak, EvalAt)
{
  auto goal = fp.peak.eval(20);

  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd fit;
  fit.setConstant(fp.variable_count, 0.0);
  fp.peak.put(fit);

  fp.peak.position.val(0.000001);
  fp.peak.amplitude.val(0.000001);
  fp.peak.width.val(0.000001);

  EXPECT_NE(fp.peak.eval(10).all(), goal.all());
  EXPECT_EQ(fp.peak.eval_at(20, fit).all(), goal.all());
}

TEST_F(Peak, EvalGrad)
{
  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(fp.variable_count, 0.0);

  auto result = fp.peak.eval_grad(20, grad);

  EXPECT_EQ(result.all(), fp.peak.eval(20).all());
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Peak, EvalGradAt)
{
  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(fp.variable_count, 0.0);
  fp.peak.eval_grad(20, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(fp.variable_count, 0.0);
  grad.setConstant(fp.variable_count, 0.0);

  fp.peak.put(fit);
  fp.peak.position.val(0.000001);
  fp.peak.amplitude.val(0.000001);
  fp.peak.width.val(0.000001);

  auto result = fp.peak.eval_grad_at(20, fit, grad);

  EXPECT_EQ(result.all(), fp.peak.eval_at(20, fit).all());
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
}

TEST_F(Peak, SurveyGradients)
{
  fp.data = generate_data(&fp, region_size);
  fp.update_indices();

  double goal_val = fp.peak.position.val();
  survey_grad(&fp, &fp.peak.position, 0.0001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);

  goal_val = fp.peak.width.val();
  survey_grad(&fp, &fp.peak.width, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.01);

  goal_val = fp.peak.amplitude.val();
  survey_grad(&fp, &fp.peak.amplitude, 0.001, std::sqrt(30000), std::sqrt(50000));
  EXPECT_NEAR(check_chi_sq(false), goal_val, 1);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.05);
}

TEST_F(Peak, FitPosition)
{
  fp.data = generate_data(&fp, region_size);
  fp.peak.width.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();

  SetUp();
  test_fit_random(random_samples, &fp,
                  {"position", &fp.peak.position,
                   fp.peak.position.min(), fp.peak.position.max(), 1e-10});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 15u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}


TEST_F(Peak, FitWidth)
{
  fp.data = generate_data(&fp, region_size);
  fp.peak.position.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();


  SetUp();
  test_fit_random(random_samples, &fp,
                  {"width", &fp.peak.width,
                   fp.peak.width.min(), fp.peak.width.max(), 1e-12});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 12u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Peak, FitAmplitude)
{
  fp.data = generate_data(&fp, region_size);
  fp.peak.width.to_fit = false;
  fp.peak.position.to_fit = false;
  fp.update_indices();

  SetUp();
  test_fit_random(random_samples, &fp,
                  {"amplitude", &fp.peak.amplitude,
                   30000, 40000, 1e-6});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 6u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Peak, FitAllThree)
{
  fp.data = generate_data(&fp, region_size);
  fp.update_indices();

  print_outside_tolerance = true;

  std::vector<ValueToVary> vals;
  vals.push_back({"width", &fp.peak.width,
                  fp.peak.width.min(), fp.peak.width.max(), 1e-11});
  vals.push_back({"position", &fp.peak.position,
                  fp.peak.position.min(), fp.peak.position.max(), 1e-12});
  vals.push_back({"amplitude", &fp.peak.amplitude,
                  30000, 50000, 1e-7});
  test_fit_random(random_samples, &fp, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 86u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Peak, FitWithSkews)
{
  fp.peak.short_tail.enabled = true;
  fp.peak.right_tail.enabled = true;
  fp.data = generate_data(&fp, region_size);
  //visualize_data(fp.data);

  fp.update_indices();

  //MESSAGE() << fp.peak.to_string() << "\n";

//  print_outside_tolerance = true;
//  verbose = true;
//  optimizer.verbose = true;

  std::vector<ValueToVary> vals;
  vals.push_back({"width", &fp.peak.width,
                  fp.peak.width.min(), fp.peak.width.max(), 1e-8});
  vals.push_back({"position", &fp.peak.position,
                  fp.peak.position.min(), fp.peak.position.max(), 1e-8});
  vals.push_back({"amplitude", &fp.peak.amplitude,
                  30000, 50000, 2e-4});
  vals.push_back({"ls_amp", &fp.peak.short_tail.amplitude,
                  fp.peak.short_tail.amplitude.min(),
                  fp.peak.short_tail.amplitude.max(), 1e-7});
  vals.push_back({"ls_slope", &fp.peak.short_tail.slope,
                  fp.peak.short_tail.slope.min(),
                  fp.peak.short_tail.slope.max(), 1e-6});
  vals.push_back({"rs_amp", &fp.peak.right_tail.amplitude,
                  fp.peak.right_tail.amplitude.min(),
                  fp.peak.right_tail.amplitude.max(), 1e-9});
  vals.push_back({"rs_slope", &fp.peak.right_tail.slope,
                  fp.peak.right_tail.slope.min(),
                  fp.peak.right_tail.slope.max(), 1e-6});
  test_fit_random(random_samples, &fp, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 100u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Peak, FitWithEverything)
{
  fp.peak.short_tail.enabled = true;
  fp.peak.right_tail.enabled = true;
  fp.peak.long_tail.enabled = true;
  fp.peak.step.enabled = true;
  fp.data = generate_data(&fp, region_size);
  //visualize_data(fp.data);

  fp.update_indices();

  //MESSAGE() << fp.peak.to_string() << "\n";

//  print_outside_tolerance = true;
//  verbose = true;0
//  optimizer.verbose = true;

  optimizer.maximum_iterations = 500;

  std::vector<ValueToVary> vals;
  vals.push_back({"width", &fp.peak.width,
                  fp.peak.width.min(), fp.peak.width.max(), 1e-7});
  vals.push_back({"position", &fp.peak.position,
                  fp.peak.position.min(), fp.peak.position.max(), 1e-7});
  vals.push_back({"amplitude", &fp.peak.amplitude,
                  30000, 50000, 1e-2});
  vals.push_back({"ls_amp", &fp.peak.short_tail.amplitude,
                  fp.peak.short_tail.amplitude.min(),
                  fp.peak.short_tail.amplitude.max(), 1e-4});
  vals.push_back({"ls_slope", &fp.peak.short_tail.slope,
                  fp.peak.short_tail.slope.min(),
                  fp.peak.short_tail.slope.max(), 1e-5});
  vals.push_back({"rs_amp", &fp.peak.right_tail.amplitude,
                  fp.peak.right_tail.amplitude.min(),
                  fp.peak.right_tail.amplitude.max(), 1e-7});
  vals.push_back({"rs_slope", &fp.peak.right_tail.slope,
                  fp.peak.right_tail.slope.min(),
                  fp.peak.right_tail.slope.max(), 1e-8});
  vals.push_back({"lt_amp", &fp.peak.long_tail.amplitude,
                  fp.peak.long_tail.amplitude.min(),
                  fp.peak.long_tail.amplitude.max(), 1e-7});
  vals.push_back({"lt_slope", &fp.peak.long_tail.slope,
                  fp.peak.long_tail.slope.min(),
                  fp.peak.long_tail.slope.max(), 1e-6});
  vals.push_back({"step_amp", &fp.peak.step.amplitude,
                  fp.peak.step.amplitude.min(),
                  fp.peak.step.amplitude.max(), 1e-8});
  test_fit_random(random_samples, &fp, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
//  EXPECT_LE(converged_finite, 0.97 * random_samples);
  EXPECT_LE(converged_perturbed, 0.40 * random_samples);
  EXPECT_LE(max_iterations_to_converge, 270u);
  EXPECT_LE(max_perturbations_to_converge, 3u);
}
