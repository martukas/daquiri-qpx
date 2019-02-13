#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/PolyBackground.h>
#include <core/fitting/weighted_data.h>

class PolyBackground : public FunctionTest
{
 protected:

  DAQuiri::PolyBackground background;

  virtual void SetUp()
  {
    background.x_offset = 0;
    background.base.bound(0, 7792);
    background.base.val(70);
    background.slope.bound(0, 198);
    background.slope.val(3);
    background.curve.bound(0, 15);
    background.curve.val(5);

    // \todo make these more permissive
  }

  double eval(double chan) const override
  {
    return background.eval(chan);
  }

  double eval_grad(double chan, Eigen::VectorXd& chan_gradients) const override
  {
    return background.eval_grad(chan, chan_gradients);
  }

};

TEST_F(PolyBackground, CheckSetup)
{
  MESSAGE() << "PolyBackground: " << background.to_string() << "\n";
}

TEST_F(PolyBackground, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(background.eval(i));
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";
}


TEST_F(PolyBackground, WithinBounds)
{
  auto data = generate_data(40);
  auto min = std::numeric_limits<double>::max();
  auto max = std::numeric_limits<double>::min();
  for (const auto& d : data.data)
  {
    min = std::min(min, d.y);
    max = std::max(max, d.y);
  }

  EXPECT_NEAR(max, 70 + 3 * 39 + 5 * square(39) , 1e-10);
  EXPECT_NEAR(min, 70, 1e-10);
}

TEST_F(PolyBackground, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(background.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(background.update_indices(i));
}

TEST_F(PolyBackground, UpdateIndex)
{
  int32_t i;

  i = 0;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 0);
  EXPECT_EQ(background.slope.index(), 1);
  EXPECT_EQ(background.curve.index(), 2);
  EXPECT_EQ(i, 3);

  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 3);
  EXPECT_EQ(background.slope.index(), 4);
  EXPECT_EQ(background.curve.index(), 5);
  EXPECT_EQ(i, 6);

  i = 42;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 42);
  EXPECT_EQ(background.slope.index(), 43);
  EXPECT_EQ(background.curve.index(), 44);
  EXPECT_EQ(i, 45);
}

TEST_F(PolyBackground, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 0);
  EXPECT_EQ(background.slope.index(), 1);
  EXPECT_EQ(background.curve.index(), 2);
  EXPECT_EQ(i, 3);

  background.base.to_fit = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), -1);
  EXPECT_EQ(background.slope.index(), 3);
  EXPECT_EQ(background.curve.index(), 4);
  EXPECT_EQ(i, 5);

  background.base.to_fit = true;
  background.slope.to_fit = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 5);
  EXPECT_EQ(background.slope.index(), -1);
  EXPECT_EQ(background.curve.index(), 6);
  EXPECT_EQ(i, 7);

  background.base.to_fit = true;
  background.slope.to_fit = true;
  background.curve.to_fit = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 7);
  EXPECT_EQ(background.slope.index(), 8);
  EXPECT_EQ(background.curve.index(), -1);
  EXPECT_EQ(i, 9);

  background.base.to_fit = false;
  background.slope.to_fit = false;
  background.curve.to_fit = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), -1);
  EXPECT_EQ(background.slope.index(), -1);
  EXPECT_EQ(background.curve.index(), -1);
  EXPECT_EQ(i, 9);
}

TEST_F(PolyBackground, UpdateIndexDisabled)
{
  int32_t i = 0;

  background.slope_enabled = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 0);
  EXPECT_EQ(background.slope.index(), -1);
  EXPECT_EQ(background.curve.index(), 1);
  EXPECT_EQ(i, 2);

  background.slope_enabled = true;
  background.curve_enabled = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 2);
  EXPECT_EQ(background.slope.index(), 3);
  EXPECT_EQ(background.curve.index(), -1);
  EXPECT_EQ(i, 4);

  background.slope_enabled = false;
  background.curve_enabled = false;
  background.update_indices(i);
  EXPECT_EQ(background.base.index(), 4);
  EXPECT_EQ(background.slope.index(), -1);
  EXPECT_EQ(background.curve.index(), -1);
  EXPECT_EQ(i, 5);
}

TEST_F(PolyBackground, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);

  background.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], background.base.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], background.slope.x());
  EXPECT_EQ(fit[2], 0.0);
  EXPECT_NE(fit[2], background.curve.x());

  background.update_indices(var_count);
  background.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], background.base.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], background.slope.x());
  EXPECT_NE(fit[2], 0.0);
  EXPECT_EQ(fit[2], background.curve.x());
}

TEST_F(PolyBackground, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);
  fit[0] = 10;
  fit[1] = 0.03;
  fit[2] = 0.05;

  background.get(fit);
  EXPECT_NEAR(background.base.val(), 70, 0.00001);
  EXPECT_NEAR(background.slope.val(), 3, 0.00001);
  EXPECT_NEAR(background.curve.val(), 5, 0.00001);
  EXPECT_NE(background.base.val(),background.base.val_at(10));
  EXPECT_NE(background.slope.val(), background.slope.val_at(0.03));
  EXPECT_NE(background.curve.val(), background.curve.val_at(0.05));

  background.update_indices(var_count);
  background.get(fit);
  EXPECT_NE(background.base.val(), 70);
  EXPECT_EQ(background.base.val(), background.base.val_at(10));
  EXPECT_NE(background.slope.val(), 3);
  EXPECT_EQ(background.slope.val(), background.slope.val_at(0.03));
  EXPECT_NE(background.curve.val(), 5);
  EXPECT_EQ(background.curve.val(), background.curve.val_at(0.05));
}

TEST_F(PolyBackground, EvalAt)
{
  auto goal = background.eval(20);

  background.update_indices(var_count);

  Eigen::VectorXd fit;
  fit.setConstant(var_count, 0.0);
  background.put(fit);

  background.base.val(10);
  background.slope.val(0.000001);
  background.curve.val(0.000001);

  EXPECT_NE(background.eval(20), goal);
  EXPECT_EQ(background.eval_at(20, fit), goal);
}

TEST_F(PolyBackground, EvalGrad)
{
  background.update_indices(var_count);

  Eigen::VectorXd grad;
  grad.setConstant(var_count, 0.0);

  auto result = background.eval_grad(10, grad);

  EXPECT_EQ(result, background.eval(10));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(PolyBackground, EvalGradAt)
{
  background.update_indices(var_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(var_count, 0.0);
  background.eval_grad(10, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(var_count, 0.0);
  grad.setConstant(var_count, 0.0);

  background.put(fit);
  background.base.val(0.000001);
  background.slope.val(0.000001);
  background.curve.val(0.000001);

  auto result = background.eval_grad_at(10, fit, grad);

  EXPECT_EQ(result, background.eval_at(10, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
}

TEST_F(PolyBackground, GradPolyBackgroundBaseOnePoint)
{
  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = background.base.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, background.base, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
  // false inflection point, but maybe ok, only one data point is not reliable anyways
}

TEST_F(PolyBackground, GradPolyBackgroundBase)
{
  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = background.base.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, background.base, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.5);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
  // \todo false gradient inflection point here
  // gradient is only good if EVEN SMALLER step sizes are used for more granularity
}

TEST_F(PolyBackground, GradPolyBackgroundSlopeOnePoint)
{
  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = background.slope.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, background.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
}

TEST_F(PolyBackground, GradPolyBackgroundSlope)
{
  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = background.slope.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, background.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
}

TEST_F(PolyBackground, GradPolyBackgroundCurveOnePoint)
{
  double goal_val = background.curve.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, background.curve);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.02);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.02);
}

TEST_F(PolyBackground, GradPolyBackgroundCurve)
{
  double goal_val = background.curve.val();
  background.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, background.curve);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.02);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.02);
}

