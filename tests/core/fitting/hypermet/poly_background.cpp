#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/PolyBackground.h>
#include <core/fitting/weighted_data.h>

class FittableBackground : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::PolyBackground background;

  double eval(double chan) const override
  {
    return background.eval(chan);
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return background.eval_at(chan, fit);
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    return background.eval_grad_at(chan, fit, grads);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    background.put(ret);
    return ret;
  }
};


class PolyBackground : public FunctionTest
{
 protected:

  FittableBackground fb;

  virtual void SetUp()
  {
    fb.background.x_offset = 0;
    fb.background.base.bound(0, 7792);
    fb.background.base.val(70);
    fb.background.slope.bound(0, 198);
    fb.background.slope.val(3);
    fb.background.curve.bound(0, 15);
    fb.background.curve.val(5);

    // \todo make these more permissive
  }

};

TEST_F(PolyBackground, CheckSetup)
{
  MESSAGE() << "PolyBackground: " << fb.background.to_string() << "\n";
}

TEST_F(PolyBackground, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(fb.background.eval(i));
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";
}


TEST_F(PolyBackground, WithinBounds)
{
  auto data = generate_data(&fb, 40);
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
  EXPECT_ANY_THROW(fb.background.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(fb.background.update_indices(i));
}

TEST_F(PolyBackground, UpdateIndex)
{
  int32_t i;

  i = 0;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 0);
  EXPECT_EQ(fb.background.slope.index(), 1);
  EXPECT_EQ(fb.background.curve.index(), 2);
  EXPECT_EQ(i, 3);

  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 3);
  EXPECT_EQ(fb.background.slope.index(), 4);
  EXPECT_EQ(fb.background.curve.index(), 5);
  EXPECT_EQ(i, 6);

  i = 42;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 42);
  EXPECT_EQ(fb.background.slope.index(), 43);
  EXPECT_EQ(fb.background.curve.index(), 44);
  EXPECT_EQ(i, 45);
}

TEST_F(PolyBackground, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 0);
  EXPECT_EQ(fb.background.slope.index(), 1);
  EXPECT_EQ(fb.background.curve.index(), 2);
  EXPECT_EQ(i, 3);

  fb.background.base.to_fit = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), -1);
  EXPECT_EQ(fb.background.slope.index(), 3);
  EXPECT_EQ(fb.background.curve.index(), 4);
  EXPECT_EQ(i, 5);

  fb.background.base.to_fit = true;
  fb.background.slope.to_fit = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 5);
  EXPECT_EQ(fb.background.slope.index(), -1);
  EXPECT_EQ(fb.background.curve.index(), 6);
  EXPECT_EQ(i, 7);

  fb.background.base.to_fit = true;
  fb.background.slope.to_fit = true;
  fb.background.curve.to_fit = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 7);
  EXPECT_EQ(fb.background.slope.index(), 8);
  EXPECT_EQ(fb.background.curve.index(), -1);
  EXPECT_EQ(i, 9);

  fb.background.base.to_fit = false;
  fb.background.slope.to_fit = false;
  fb.background.curve.to_fit = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), -1);
  EXPECT_EQ(fb.background.slope.index(), -1);
  EXPECT_EQ(fb.background.curve.index(), -1);
  EXPECT_EQ(i, 9);
}

TEST_F(PolyBackground, UpdateIndexDisabled)
{
  int32_t i = 0;

  fb.background.slope_enabled = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 0);
  EXPECT_EQ(fb.background.slope.index(), -1);
  EXPECT_EQ(fb.background.curve.index(), 1);
  EXPECT_EQ(i, 2);

  fb.background.slope_enabled = true;
  fb.background.curve_enabled = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 2);
  EXPECT_EQ(fb.background.slope.index(), 3);
  EXPECT_EQ(fb.background.curve.index(), -1);
  EXPECT_EQ(i, 4);

  fb.background.slope_enabled = false;
  fb.background.curve_enabled = false;
  fb.background.update_indices(i);
  EXPECT_EQ(fb.background.base.index(), 4);
  EXPECT_EQ(fb.background.slope.index(), -1);
  EXPECT_EQ(fb.background.curve.index(), -1);
  EXPECT_EQ(i, 5);
}

TEST_F(PolyBackground, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);

  fb.background.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], fb.background.base.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], fb.background.slope.x());
  EXPECT_EQ(fit[2], 0.0);
  EXPECT_NE(fit[2], fb.background.curve.x());

  fb.background.update_indices(fb.variable_count);
  fb.background.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], fb.background.base.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], fb.background.slope.x());
  EXPECT_NE(fit[2], 0.0);
  EXPECT_EQ(fit[2], fb.background.curve.x());
}

TEST_F(PolyBackground, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);
  fit[0] = 10;
  fit[1] = 0.03;
  fit[2] = 0.05;

  fb.background.get(fit);
  EXPECT_NEAR(fb.background.base.val(), 70, 0.00001);
  EXPECT_NEAR(fb.background.slope.val(), 3, 0.00001);
  EXPECT_NEAR(fb.background.curve.val(), 5, 0.00001);
  EXPECT_NE(fb.background.base.val(),fb.background.base.val_at(10));
  EXPECT_NE(fb.background.slope.val(), fb.background.slope.val_at(0.03));
  EXPECT_NE(fb.background.curve.val(), fb.background.curve.val_at(0.05));

  fb.background.update_indices(fb.variable_count);
  fb.background.get(fit);
  EXPECT_NE(fb.background.base.val(), 70);
  EXPECT_EQ(fb.background.base.val(), fb.background.base.val_at(10));
  EXPECT_NE(fb.background.slope.val(), 3);
  EXPECT_EQ(fb.background.slope.val(), fb.background.slope.val_at(0.03));
  EXPECT_NE(fb.background.curve.val(), 5);
  EXPECT_EQ(fb.background.curve.val(), fb.background.curve.val_at(0.05));
}

TEST_F(PolyBackground, EvalAt)
{
  auto goal = fb.background.eval(20);

  fb.background.update_indices(fb.variable_count);

  Eigen::VectorXd fit;
  fit.setConstant(fb.variable_count, 0.0);
  fb.background.put(fit);

  fb.background.base.val(10);
  fb.background.slope.val(0.000001);
  fb.background.curve.val(0.000001);

  EXPECT_NE(fb.background.eval(20), goal);
  EXPECT_EQ(fb.background.eval_at(20, fit), goal);
}

TEST_F(PolyBackground, EvalGrad)
{
  fb.background.update_indices(fb.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(fb.variable_count, 0.0);

  auto result = fb.background.eval_grad(10, grad);

  EXPECT_EQ(result, fb.background.eval(10));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(PolyBackground, EvalGradAt)
{
  fb.background.update_indices(fb.variable_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(fb.variable_count, 0.0);
  fb.background.eval_grad(10, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(fb.variable_count, 0.0);
  grad.setConstant(fb.variable_count, 0.0);

  fb.background.put(fit);
  fb.background.base.val(0.000001);
  fb.background.slope.val(0.000001);
  fb.background.curve.val(0.000001);

  auto result = fb.background.eval_grad_at(10, fit, grad);

  EXPECT_EQ(result, fb.background.eval_at(10, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
}

TEST_F(PolyBackground, GradPolyBackgroundBase)
{
  fb.data = generate_data(&fb, 40);

  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = fb.background.base.val();
  fb.background.update_indices(fb.variable_count);
  survey_grad(&fb, fb.background.base, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.5);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
  // \todo false gradient inflection point here
  // gradient is only good if EVEN SMALLER step sizes are used for more granularity
}

TEST_F(PolyBackground, GradPolyBackgroundSlope)
{
  fb.data = generate_data(&fb, 40);

  // chi-sq is only good if smaller step sizes are used for more granularity
  double goal_val = fb.background.slope.val();
  fb.background.update_indices(fb.variable_count);
  survey_grad(&fb, fb.background.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
}

TEST_F(PolyBackground, GradPolyBackgroundCurve)
{
  fb.data = generate_data(&fb, 40);

  double goal_val = fb.background.curve.val();
  fb.background.update_indices(fb.variable_count);
  survey_grad(&fb, fb.background.curve);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.02);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.02);
}

