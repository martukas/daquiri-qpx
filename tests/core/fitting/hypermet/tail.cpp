#include "../function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Tail.h>

#include <core/fitting/optimizers/dlib_adapter.h>


class FittableTail : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::Value position;
  DAQuiri::ValuePositive amplitude;
  DAQuiri::Value width;

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

    // \todo uncerts
  }
};

class Tail : public FunctionTest
{
 protected:
  FittableTail ft;

  virtual void SetUp()
  {
    //ft.amplitude.bound(0, 10000);
    ft.amplitude.val(40000);
    ft.amplitude.update_index(ft.variable_count);

    ft.width.bound(0.8, 5.0);
    ft.width.val(3.2);
    ft.width.update_index(ft.variable_count);

    ft.position.bound(14, 28);
    ft.position.val(21);
    ft.position.update_index(ft.variable_count);

    ft.tail.amplitude.bound(0.0001, 1.5);
    ft.tail.amplitude.val(0.05);
    ft.tail.slope.bound(0.2, 50);
    ft.tail.slope.val(30);
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
  auto data = generate_data(&ft, 40);
  visualize_data(data);
}


TEST_F(Tail, WithinBounds)
{
  auto data = generate_data(&ft, 40);
  EXPECT_NEAR(data.count_min, 0.0, 1.0);
  EXPECT_NEAR(data.count_max, 1871, 1.0);
}

TEST_F(Tail, LeftOriented)
{
  ft.tail.side = DAQuiri::Side::left;
  auto data = generate_data(&ft, 40);
  EXPECT_NEAR(data.data.front().count, 1607, 1.0);
  EXPECT_NEAR(data.data.back().count, 0.0, 1.0);
}

TEST_F(Tail, RightOriented)
{
  ft.tail.side = DAQuiri::Side::right;
  auto data = generate_data(&ft, 40);
  EXPECT_NEAR(data.data.front().count, 0.0, 1.0);
  EXPECT_NEAR(data.data.back().count, 1658, 1.0);
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

TEST_F(Tail, GradAmplitude)
{
  double goal_val = ft.tail.amplitude.val();
  ft.data = generate_data(&ft, 40);

  ft.tail.slope.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();
  survey_grad(&ft, &ft.tail.amplitude, 0.000005);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.00001);
}

TEST_F(Tail, FitAmplitudeOnly)
{
  ft.data = generate_data(&ft, 40);

  ft.tail.slope.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &ft, &ft.tail.amplitude, 1.0, 1e-5);

  ft.tail.amplitude.val(0.05);
  test_fit_random(20, &optimizer, &ft, &ft.tail.amplitude,
                  ft.tail.amplitude.min(), ft.tail.amplitude.max(), 1e-4);
}

TEST_F(Tail, GradSlope)
{
  double goal_val = ft.tail.slope.val();
  ft.data = generate_data(&ft, 40);

  ft.tail.amplitude.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();
  survey_grad(&ft, &ft.tail.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
}

TEST_F(Tail, FitSlopeOnly)
{
  ft.data = generate_data(&ft, 40);

  ft.tail.amplitude.to_fit = false;
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &ft, &ft.tail.slope, 1.0, 1e-4);

  ft.tail.slope.val(30);
  test_fit_random(20, &optimizer, &ft, &ft.tail.slope,
                  ft.tail.slope.min(), ft.tail.slope.max(), 1e-3);
}

TEST_F(Tail, FitAmplitudeAndSlope)
{
  double goal_amplitude = ft.tail.amplitude.val();
  double goal_slope = ft.tail.slope.val();

  ft.data = generate_data(&ft, 40);
  ft.width.to_fit = false;
  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.update_indices();

  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());
  std::uniform_real_distribution<double> amplitude_dist(ft.tail.amplitude.min(),
                                                        ft.tail.amplitude.max());
  std::uniform_real_distribution<double> slope_dist(ft.tail.slope.min(),
                                                    ft.tail.slope.max());

  DAQuiri::DLibOptimizer optimizer;

  for (size_t i = 0; i < 50; ++i)
  {
    ft.tail.amplitude.val(amplitude_dist(random_generator));
    ft.tail.slope.val(slope_dist(random_generator));
    MESSAGE() << "Attempt[" << i << "] " << ft.tail.to_string() << "\n";
    ft.save_fit(optimizer.minimize(&ft));
    MESSAGE() << "Result:               " << ft.tail.to_string() << "\n";
    MESSAGE() << "      amplitude delta="
              << (goal_amplitude - ft.tail.amplitude.val()) << "\n";
    MESSAGE() << "          slope delta="
              << (goal_slope - ft.tail.slope.val()) << "\n";

    EXPECT_NEAR(ft.tail.amplitude.val(), goal_amplitude, 0.01);
    EXPECT_NEAR(ft.tail.slope.val(), goal_slope, 0.1);
  }
}

TEST_F(Tail, GradParentWidth)
{
  double goal_val = ft.width.val();
  ft.data = generate_data(&ft, 40);

  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();
  survey_grad(&ft, &ft.width, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.005);
}

TEST_F(Tail, FitParentWidthOnly)
{
  ft.data = generate_data(&ft, 40);

  ft.position.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &ft, &ft.width, 1.0, 0.05);

  ft.width.val(3.2);
  test_fit_random(20, &optimizer, &ft, &ft.width,
                  ft.width.min(), ft.width.max(), 0.05);
}

TEST_F(Tail, GradParentPosition)
{
  double goal_val = ft.position.val();
  ft.data = generate_data(&ft, 40);

  ft.width.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();
  survey_grad(&ft, &ft.position, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.005);
}

TEST_F(Tail, FitParentPositionOnly)
{
  ft.data = generate_data(&ft, 40);

  ft.width.to_fit = false;
  ft.amplitude.to_fit = false;
  ft.tail.amplitude.to_fit = false;
  ft.tail.slope.to_fit = false;
  ft.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &ft, &ft.position, 16, 0.05);

  ft.position.val(21);
  test_fit_random(20, &optimizer, &ft, &ft.position,
                  ft.position.min(), ft.position.max(), 0.05);
}

// \todo Nowhere close. Might be implemented wrong

//TEST_F(Tail, GradParentAmplitude)
//{
//  double goal_val = ft.amplitude.val();
//  ft.data = generate_data(&ft, 40);
//
//  ft.width.to_fit = false;
//  ft.position.to_fit = false;
//  ft.tail.amplitude.to_fit = false;
//  ft.tail.slope.to_fit = false;
//  ft.update_indices();
//  survey_grad(&ft, &ft.amplitude, 0.1, std::sqrt(30000), std::sqrt(40000));
//  EXPECT_NEAR(check_chi_sq(false), goal_val, 40);
//  EXPECT_NEAR(check_gradients(false), goal_val, 40);
//}
//
//TEST_F(Tail, FitParentAmplitudeOnly)
//{
//  ft.data = generate_data(&ft, 40);
//
//  ft.width.to_fit = false;
//  ft.position.to_fit = false;
//  ft.tail.amplitude.to_fit = false;
//  ft.tail.slope.to_fit = false;
//  ft.update_indices();
//
//  DAQuiri::DLibOptimizer optimizer;
//  optimizer.tolerance = 1e-12;
//  test_fit(5, &optimizer, &ft, &ft.amplitude, 30000, 0.5);
//
//  ft.width.val(40000);
//  test_fit_random(20, &optimizer, &ft, &ft.amplitude,
//                  30000, 50000, 260);
//}
