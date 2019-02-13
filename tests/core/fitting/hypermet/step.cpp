#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/weighted_data.h>

class Step : public FunctionTest
{
 protected:

  DAQuiri::Step step;

  DAQuiri::Value position;
  DAQuiri::Value amplitude;
  DAQuiri::Value width;

  virtual void SetUp()
  {
    amplitude.bound(0, 1000);
    amplitude.val(40);
    amplitude.update_index(var_count);

    width.bound(0.8, 5.0);
    width.val(2);
    width.update_index(var_count);

    position.bound(0, 40);
    position.val(20);
    position.update_index(var_count);

    step.amplitude.bound(0.000001, 0.05);
    step.amplitude.val(0.05);
  }

  DAQuiri::PrecalcVals precalc_spoof(double chan) const
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
    ret.i_pos = position.index(); // should not matter?
    return ret;
  }

  double eval(double chan) const override
  {
    return step.eval(precalc_spoof(chan));
  }

  double eval_grad(double chan, Eigen::VectorXd& chan_gradients) const override
  {
    return step.eval_grad(precalc_spoof(chan), chan_gradients);
  }

};

TEST_F(Step, CheckSetup)
{
  MESSAGE() << "Gaussian amp: " << amplitude.to_string() << "\n";
  MESSAGE() << "Gaussian pos: " << position.to_string() << "\n";
  MESSAGE() << "Gaussian width: " << width.to_string() << "\n";
  MESSAGE() << "Step: " << step.to_string() << "\n";
}

TEST_F(Step, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(step.eval(precalc_spoof(i)));
  }
  MESSAGE() << "peak(channel):\n" << visualize(channels, y, 100) << "\n";
}

TEST_F(Step, WithinBounds)
{
  auto data = generate_data(40);
  auto min = std::numeric_limits<double>::max();
  auto max = std::numeric_limits<double>::min();
  for (const auto& d : data.data)
  {
    min = std::min(min, d.y);
    max = std::max(max, d.y);
  }

  EXPECT_NEAR(max, 2.0, 1e-15);
  EXPECT_NEAR(min, 0.0, 1e-40);
}

TEST_F(Step, LeftOriented)
{
  step.side = DAQuiri::Side::left;
  auto data = generate_data(40);
  EXPECT_NEAR(data.data.front().y, 2.0, 1e-15);
  EXPECT_NEAR(data.data.back().y, 0.0, 1e-40);
}

TEST_F(Step, RightOriented)
{
  step.side = DAQuiri::Side::right;
  auto data = generate_data(40);
  EXPECT_NEAR(data.data.front().y, 0.0, 1e-40);
  EXPECT_NEAR(data.data.back().y, 2.0, 1e-15);
}

TEST_F(Step, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(step.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(step.update_indices(i));
}

TEST_F(Step, UpdateIndex)
{
  int32_t i;

  i = 0;
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 0);
  EXPECT_EQ(i, 1);

  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 1);
  EXPECT_EQ(i, 2);

  i = 42;
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 42);
  EXPECT_EQ(i, 43);
}

TEST_F(Step, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 0);
  EXPECT_EQ(i, 1);

  step.amplitude.to_fit = false;
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), -1);
  EXPECT_EQ(i, 1);

  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), -1);
  EXPECT_EQ(i, 1);
}

TEST_F(Step, UpdateIndexDisabled)
{
  step.enabled = false;
  int32_t i;

  i = 0;
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), -1);
  EXPECT_EQ(i, 0);

  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), -1);
  EXPECT_EQ(i, 0);
}

TEST_F(Step, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(1, 0.0);

  step.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], step.amplitude.x());

  int32_t i{0};
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 0);

  step.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], step.amplitude.x());
}

TEST_F(Step, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(1, 0.005);

  step.get(fit);
  EXPECT_EQ(step.amplitude.val(), 0.05);

  int32_t i{0};
  step.update_indices(i);
  EXPECT_EQ(step.amplitude.index(), 0);

  step.get(fit);
  EXPECT_NE(step.amplitude.val(), 0.05);
  EXPECT_EQ(step.amplitude.val(), step.amplitude.val_at(0.005));
}

TEST_F(Step, EvalAt)
{
  auto pre = precalc_spoof(20);

  auto goal = step.eval(pre);

  int32_t i{0};
  step.update_indices(i);

  Eigen::VectorXd fit;
  fit.setConstant(1, 0.0);
  step.put(fit);

  step.amplitude.val(0.000001);

  EXPECT_NE(step.eval(pre), goal);
  EXPECT_EQ(step.eval_at(pre, fit), goal);
}

TEST_F(Step, EvalGrad)
{
  auto pre = precalc_spoof(10);

  step.update_indices(var_count);

  Eigen::VectorXd grad;
  grad.setConstant(var_count, 0.0);

  auto result = step.eval_grad(pre, grad);

  EXPECT_EQ(result, step.eval(pre));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_EQ(grad[2], 0.0); // pos gradient should be unaffected?
  EXPECT_NE(grad[3], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Step, EvalGradAt)
{
  auto pre = precalc_spoof(10);

  int32_t i{3};
  step.update_indices(i);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(i, 0.0);
  step.eval_grad(pre, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(i, 0.0);
  grad.setConstant(i, 0.0);

  step.put(fit);
  step.amplitude.val(0.000001);

  auto result = step.eval_grad_at(pre, fit, grad);

  EXPECT_EQ(result, step.eval_at(pre, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
  EXPECT_EQ(grad[3], grad_goal[3]);
}

TEST_F(Step, GradStepAmpOnePoint)
{
  double goal_val = step.amplitude.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, step.amplitude);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.00001);
}

TEST_F(Step, GradStepAmp)
{
  double goal_val = step.amplitude.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, step.amplitude);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.00001);
}

TEST_F(Step, GradWidthOnePoint)
{
  double goal_val = width.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, width);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.005);
  check_gradients(true);
  // false inflection point, but maybe ok, only one data point is not reliable anyways
}

TEST_F(Step, GradWidth)
{
  double goal_val = width.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, width);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.005);
}

TEST_F(Step, GradAmpOnePoint)
{
  double goal_val = amplitude.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 4);
  EXPECT_NEAR(check_gradients(true), goal_val, 4);
}

TEST_F(Step, GradAmp)
{
  double goal_val = amplitude.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 4);
  EXPECT_NEAR(check_gradients(true), goal_val, 4);
}

TEST_F(Step, GradPosOnePoint)
{
  double goal_val = position.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, position);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  check_gradients(true);
  // \todo gradient not affected!
}

TEST_F(Step, GradPos)
{
  double goal_val = position.val();
  step.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, position);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  check_gradients(true);
  // \todo gradient not affected!
}
