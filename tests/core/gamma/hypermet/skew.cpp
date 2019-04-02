#include "../../fitting/function_test.h"

#include <core/gamma/hypermet/skew.h>
#include <core/fitting/parameter/positive_param.h>
#include <core/fitting/parameter/sine_bounded_param.h>

#include <core/util/UTF_extensions.h>

class FittableTail : public DAQuiri::FittableRegion
{
  std::uniform_real_distribution<double> x_dist {-M_PI_2, M_PI_2};
 public:

  DAQuiri::SineBoundedValue position;
  DAQuiri::PositiveValue amplitude;
  DAQuiri::SineBoundedValue width;

  DAQuiri::Tail tail;

  void update_indices() override
  {
    variable_count = 0;
    amplitude.update_index(variable_count);
    width.update_index(variable_count);
    position.update_index(variable_count);
    tail.update_indices(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    position.put(ret);
    amplitude.put(ret);
    width.put(ret);
    tail.put(ret);
    return ret;
  }

  bool perturb(std::mt19937& rng) override
  {
    if (tail.amplitude.valid_index())
      tail.amplitude.x(x_dist(rng));
    if (tail.slope.valid_index())
      tail.slope.x(x_dist(rng));
    if (width.valid_index())
      width.x(x_dist(rng));
    if (position.valid_index())
      position.x(x_dist(rng));
    if (amplitude.valid_index())
      amplitude.x(amplitude.x() + x_dist(rng));
    return true;
  }

  bool sane() const override
  {
    if (!tail.sane(1e-5, 1e-4, 1e-3))
      return false;
    if (width.to_fit && width.at_extremum(1e-3, 1e-3))
      return false;
    if (position.to_fit && position.at_extremum(1e-2, 1e-2))
      return false;
    return true;
  }

  DAQuiri::PrecalcVals precalc(double chan) const
  {
    DAQuiri::PrecalcVals ret;
    ret.ampl = amplitude.val();
    ret.half_ampl = 0.5 * ret.ampl;
    ret.width = width.val();
    ret.spread = (chan - position.val()) / ret.width;

    ret.amp_grad = amplitude.grad();
    ret.width_grad = width.grad();
    ret.pos_grad = position.grad();

    ret.i_amp = amplitude.index();
    ret.i_width = width.index();
    ret.i_pos = position.index();
    return ret;
  }

  DAQuiri::PrecalcVals precalc_at(double chan, const Eigen::VectorXd& fit) const
  {
    DAQuiri::PrecalcVals ret;
    ret.ampl = amplitude.val_from(fit);
    ret.half_ampl = 0.5 * ret.ampl;
    ret.width = width.val_from(fit);
    ret.spread = (chan - position.val_from(fit)) / ret.width;

    ret.amp_grad = amplitude.grad_from(fit);
    ret.width_grad = width.grad_from(fit);
    ret.pos_grad = position.grad_from(fit);

    ret.i_amp = amplitude.index();
    ret.i_width = width.index();
    ret.i_pos = position.index();
    return ret;
  }

  double eval(double chan) const override
  {
    return tail.eval(precalc(chan));
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return tail.eval_at(precalc_at(chan, fit), fit);
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    return tail.eval_grad_at(precalc_at(chan, fit), fit, grads);
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    amplitude.get(result.variables);
    width.get(result.variables);
    position.get(result.variables);
    tail.get(result.variables);

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();

    auto chi = chi_sq();
    amplitude.get_uncert(diags, chi);
    width.get_uncert(diags, chi);
    position.get_uncert(diags, chi);
    tail.get_uncerts(diags, chi);
  }

  std::string to_string(std::string prepend = "") const override
  {
    Eigen::VectorXd grads;
    auto x = variables();
    chi_sq_gradient(x, grads);
    std::stringstream ss;

    ss << "parent_pos = " << position.to_string() << "\n";
    ss << prepend << "parent_amp = " << amplitude.to_string() << "\n";
    ss << prepend << "parent_width = "  << width.to_string() << "\n";
    ss << prepend << "tail = " << tail.to_string() << "\n";
    ss << prepend << "  chi" + UTF_superscript(2) + "=" + std::to_string(chi_sq());
    ss << "  grads=" << grads.transpose() << "\n";
    return ss.str();
  }
};

class Tail : public FunctionTest
{
 protected:
  FittableTail ft;
  size_t region_size{100};
  size_t random_samples{100};

  void SetUp() override
  {
//    print_outside_tolerance = true;

    //optimizer.verbosity = 5;
    optimizer.maximum_iterations = 1000;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;
//    optimizer.use_epsilon_check = false;
//    optimizer.min_g_norm = 1e-7;
    optimizer.perform_sanity_checks = false;
    optimizer.maximum_perturbations = 0;

    ft.tail.amplitude.bound(0.0001, 1.5);
    ft.tail.amplitude.val(0.05);
    ft.tail.slope.bound(0.2, 50);
    ft.tail.slope.val(30);

//    ft.amplitude.bound(0, 50000);
    ft.amplitude.val(40000);
    ft.amplitude.update_index(ft.variable_count);

    ft.width.bound(0.8, 5.0);
    ft.width.val(3.2);
    ft.width.update_index(ft.variable_count);

    ft.position.bound(44, 68);
    ft.position.val(51);
    ft.position.update_index(ft.variable_count);
  }
};

TEST_F(Tail, CheckSetup)
{
  MESSAGE() << "Gaussian amp: " << ft.amplitude.to_string() << "\n";
  MESSAGE() << "Gaussian pos: " << ft.position.to_string() << "\n";
  MESSAGE() << "Gaussian width: " << ft.width.to_string() << "\n";
  MESSAGE() << "Tail: " << ft.tail.to_string() << "\n";
}

TEST_F(Tail, Visualize)
{
  auto data = generate_data(&ft, region_size);
  visualize_data(data);
}


TEST_F(Tail, WithinBounds)
{
  auto data = generate_data(&ft, region_size);
  EXPECT_NEAR(data.count_min(), 0.0, 1.0);
  EXPECT_NEAR(data.count_max(), 1871, 1.0);
}

TEST_F(Tail, LeftOriented)
{
  ft.tail.side = DAQuiri::Side::left;
  auto data = generate_data(&ft, region_size);
  EXPECT_NEAR(data.count.front(), 1175, 1.0);
  EXPECT_NEAR(data.count.back(), 0.0, 1.0);
}

TEST_F(Tail, RightOriented)
{
  ft.tail.side = DAQuiri::Side::right;
  auto data = generate_data(&ft, region_size);
  EXPECT_NEAR(data.count.front(), 0.0, 1.0);
  EXPECT_NEAR(data.count.back(), 1213, 1.0);
}

TEST_F(Tail, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(ft.tail.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(ft.tail.update_indices(i));
}

TEST_F(Tail, UpdateIndex)
{
  int32_t i;

  i = 0;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 0);
  EXPECT_EQ(ft.tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 2);
  EXPECT_EQ(ft.tail.slope.index(), 3);
  EXPECT_EQ(i, 4);

  i = 42;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 42);
  EXPECT_EQ(ft.tail.slope.index(), 43);
  EXPECT_EQ(i, 44);
}

TEST_F(Tail, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 0);
  EXPECT_EQ(ft.tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  ft.tail.amplitude.to_fit = false;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), -1);
  EXPECT_EQ(ft.tail.slope.index(), 2);
  EXPECT_EQ(i, 3);

  ft.tail.amplitude.to_fit = true;
  ft.tail.slope.to_fit = false;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 3);
  EXPECT_EQ(ft.tail.slope.index(), -1);
  EXPECT_EQ(i, 4);

  ft.tail.slope.to_fit = true;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 4);
  EXPECT_EQ(ft.tail.slope.index(), 5);
  EXPECT_EQ(i, 6);
}

TEST_F(Tail, UpdateIndexDisabled)
{
  ft.tail.enabled = false;
  int32_t i;

  i = 0;
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), -1);
  EXPECT_EQ(ft.tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), -1);
  EXPECT_EQ(ft.tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  // \todo test resetting of indices
}

TEST_F(Tail, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);

  ft.tail.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], ft.tail.amplitude.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], ft.tail.slope.x());

  int32_t i{0};
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 0);
  EXPECT_EQ(ft.tail.slope.index(), 1);

  ft.tail.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], ft.tail.amplitude.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], ft.tail.slope.x());
}

TEST_F(Tail, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  fit[0] = 0.005;
  fit[1] = 0.03;

  ft.tail.get(fit);
  EXPECT_NEAR(ft.tail.amplitude.val(), 0.05, 0.00001);
  EXPECT_NEAR(ft.tail.slope.val(), 30, 0.00001);
  EXPECT_NE(ft.tail.amplitude.val(),ft.tail.amplitude.val_at(0.005));
  EXPECT_NE(ft.tail.slope.val(), ft.tail.slope.val_at(0.03));

  int32_t i{0};
  ft.tail.update_indices(i);
  EXPECT_EQ(ft.tail.amplitude.index(), 0);
  EXPECT_EQ(ft.tail.slope.index(), 1);

  ft.tail.get(fit);
  EXPECT_NE(ft.tail.amplitude.val(), 0.5);
  EXPECT_EQ(ft.tail.amplitude.val(), ft.tail.amplitude.val_at(0.005));
  EXPECT_NE(ft.tail.slope.val(), 30);
  EXPECT_EQ(ft.tail.slope.val(), ft.tail.slope.val_at(0.03));
}

TEST_F(Tail, EvalAt)
{
  auto pre = ft.precalc(20);

  auto goal = ft.tail.eval(pre);

  int32_t i{0};
  ft.tail.update_indices(i);

  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  ft.tail.put(fit);

  ft.tail.amplitude.val(0.000001);
  ft.tail.slope.val(0.000001);

  EXPECT_NE(ft.tail.eval(pre), goal);
  EXPECT_EQ(ft.tail.eval_at(pre, fit), goal);
}

TEST_F(Tail, EvalGrad)
{
  auto pre = ft.precalc(10);

  ft.tail.update_indices(ft.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(ft.variable_count, 0.0);

  auto result = ft.tail.eval_grad(pre, grad);

  EXPECT_EQ(result, ft.tail.eval(pre));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);
  EXPECT_NE(grad[3], 0.0);
  EXPECT_NE(grad[4], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Tail, EvalGradAt)
{
  auto pre = ft.precalc(10);

  ft.tail.update_indices(ft.variable_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(ft.variable_count, 0.0);
  ft.tail.eval_grad(pre, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(ft.variable_count, 0.0);
  grad.setConstant(ft.variable_count, 0.0);

  ft.tail.put(fit);
  ft.tail.amplitude.val(0.000001);
  ft.tail.slope.val(0.000001);

  auto result = ft.tail.eval_grad_at(pre, fit, grad);

  EXPECT_EQ(result, ft.tail.eval_at(pre, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
  EXPECT_EQ(grad[3], grad_goal[3]);
  EXPECT_EQ(grad[4], grad_goal[4]);
}

TEST_F(Tail, SurveyGradients)
{
  ft.data = generate_data(&ft, region_size);
  ft.update_indices();

  EXPECT_TRUE(optimizer.check_gradient(&ft));

  double goal_val = ft.tail.amplitude.val();
  survey_grad(&ft, &ft.tail.amplitude, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
  //EXPECT_NEAR(check_gradients(false), goal_val, 0.01);
  survey_grad(&ft, &ft.tail.amplitude, 0.05);
  check_gradients(true);
  check_gradient_deltas(true);

  goal_val = ft.tail.slope.val();
  survey_grad(&ft, &ft.tail.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
  survey_grad(&ft, &ft.tail.slope, 0.05);
  check_gradients(true);
  check_gradient_deltas(true);

  goal_val = ft.width.val();
  survey_grad(&ft, &ft.width, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.01);
  survey_grad(&ft, &ft.width, 0.05);
  check_gradients(true);
  check_gradient_deltas(true);

  goal_val = ft.position.val();
  survey_grad(&ft, &ft.position, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
//  EXPECT_NEAR(check_gradients(false), goal_val, 0.01);
  survey_grad(&ft, &ft.position, 0.05);
  check_gradients(true);
  check_gradient_deltas(true);

  goal_val = ft.amplitude.val();
  survey_grad(&ft, &ft.amplitude, 0.1, std::sqrt(30000), std::sqrt(40000));
  EXPECT_NEAR(check_chi_sq(false), goal_val, 50);
  EXPECT_NEAR(check_gradients(false), goal_val, 50);
  survey_grad(&ft, &ft.amplitude, 0.5, std::sqrt(30000), std::sqrt(40000));
  check_gradients(true);
  check_gradient_deltas(true);

}

TEST_F(Tail, FitAmplitude)
{
  ft.data = generate_data(&ft, region_size);
  ft.tail.slope.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  SetUp();
//  
  test_fit_random(random_samples, &ft,
                  {"amplitude", &ft.tail.amplitude,
                   ft.tail.amplitude.min(), ft.tail.amplitude.max(), 1e-7});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 13u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitSlope)
{
  ft.data = generate_data(&ft, region_size);
  ft.tail.amplitude.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  SetUp();
  test_fit_random(random_samples, &ft,
                  {"slope", &ft.tail.slope,
                   ft.tail.slope.min(), ft.tail.slope.max(), 1e-9});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 12u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitBoth)
{
  ft.data = generate_data(&ft, region_size);
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  
//  optimizer.verbosity=5;

  std::vector<ValueToVary> vals;
  vals.push_back({"amplitude", &ft.tail.amplitude,
                  ft.tail.amplitude.min(), ft.tail.amplitude.max(), 1e-8});
  vals.push_back({"slope", &ft.tail.slope,
                  ft.tail.slope.min(), ft.tail.slope.max(), 1e-5});
  test_fit_random(random_samples, &ft, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 48u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitParentWidth)
{
  ft.data = generate_data(&ft, region_size);
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();

  SetUp();
  test_fit_random(random_samples, &ft,
                  {"parent_width", &ft.width,
                   ft.width.min(), ft.width.max(), 1e-9});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 10u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitParentPosition)
{
  ft.data = generate_data(&ft, region_size);
  ft.width.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();

  SetUp();
  test_fit_random(random_samples, &ft,
                  {"parent_position", &ft.position,
                   ft.position.min(), ft.position.max(), 1e-7});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 13u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitParentAmplitude)
{
  ft.data = generate_data(&ft, region_size);
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();

  SetUp();
  test_fit_random(random_samples, &ft,
                  {"parent_amplitude", &ft.amplitude,
//                   ft.amplitude.min(), ft.amplitude.max(), 1e-3});
                   30000, 50000, 1e-15});

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 6u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

// the 2 amplitudes are degenerate. can't deconvolute without peak
// so we must test their fitting separately in these tests
TEST_F(Tail, FitFourA)
{
  ft.data = generate_data(&ft, region_size);
  ft.amplitude.to_fit = false;
  ft.update_indices();

  
//  optimizer.verbosity=5;

  std::vector<ValueToVary> vals;
  vals.push_back({"amplitude", &ft.tail.amplitude,
                  ft.tail.amplitude.min(), ft.tail.amplitude.max(), 1e-6});
  vals.push_back({"slope", &ft.tail.slope,
                  ft.tail.slope.min(), ft.tail.slope.max(), 1e-3});
  vals.push_back({"parent_width", &ft.width, ft.width.min(), ft.width.max(), 1e-4});
  vals.push_back({"parent_position", &ft.position,
                  ft.position.min(), ft.position.max(), 1e-4});
  test_fit_random(random_samples, &ft, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_LE(not_sane, 0u);
  EXPECT_LE(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 87u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Tail, FitFourB)
{
  ft.data = generate_data(&ft, region_size);
  ft.tail.amplitude.to_fit = false;
  ft.update_indices();

  print_outside_tolerance = true;

//  optimizer.verbosity=5;

  std::vector<ValueToVary> vals;
  vals.push_back({"slope", &ft.tail.slope,
                  ft.tail.slope.min(), ft.tail.slope.max(), 1e-3});
  vals.push_back({"parent_width", &ft.width, ft.width.min(), ft.width.max(), 1e-4});
  vals.push_back({"parent_position", &ft.position,
                  ft.position.min(), ft.position.max(), 1e-4});
  vals.push_back({"parent_amplitude", &ft.amplitude, 30000, 50000, 1e-5});
  test_fit_random(random_samples, &ft, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_LE(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 69u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}
