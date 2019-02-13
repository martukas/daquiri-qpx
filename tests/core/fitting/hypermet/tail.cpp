#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Tail.h>
#include <core/fitting/weighted_data.h>

class Tail : public FunctionTest
{
 protected:

  DAQuiri::Tail tail;

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

    tail.amplitude.bound(0.0001, 1.5);
    tail.amplitude.val(0.5);
    tail.slope.bound(0.2, 50);
    tail.slope.val(30);
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
    return tail.eval(precalc_spoof(chan));
  }

  double eval_grad(double chan, Eigen::VectorXd& chan_gradients) const override
  {
    return tail.eval_grad(precalc_spoof(chan), chan_gradients);
  }

};

TEST_F(Tail, CheckSetup)
{
  MESSAGE() << "Gaussian amp: " << amplitude.to_string() << "\n";
  MESSAGE() << "Gaussian pos: " << position.to_string() << "\n";
  MESSAGE() << "Gaussian width: " << width.to_string() << "\n";
  MESSAGE() << "Tail: " << tail.to_string() << "\n";
}

TEST_F(Tail, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(tail.eval(precalc_spoof(i)));
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";
}


TEST_F(Tail, WithinBounds)
{
  auto data = generate_data(40);
  auto min = std::numeric_limits<double>::max();
  auto max = std::numeric_limits<double>::min();
  for (const auto& d : data.data)
  {
    min = std::min(min, d.y);
    max = std::max(max, d.y);
  }

  EXPECT_LE(max, 20.0);
  EXPECT_NEAR(min, 0.0, 1e-39);
}

TEST_F(Tail, LeftOriented)
{
  tail.side = DAQuiri::Side::left;
  auto data = generate_data(40);
  EXPECT_NEAR(data.data.front().y, 14.0, 1.0);
  EXPECT_NEAR(data.data.back().y, 0.0, 1e-39);
}

TEST_F(Tail, RightOriented)
{
  tail.side = DAQuiri::Side::right;
  auto data = generate_data(40);
  EXPECT_NEAR(data.data.front().y, 0.0, 1e-39);
  EXPECT_NEAR(data.data.back().y, 14.0, 1.0);
}

TEST_F(Tail, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(tail.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(tail.update_indices(i));
}

TEST_F(Tail, UpdateIndex)
{
  int32_t i;

  i = 0;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 0);
  EXPECT_EQ(tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 2);
  EXPECT_EQ(tail.slope.index(), 3);
  EXPECT_EQ(i, 4);

  i = 42;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 42);
  EXPECT_EQ(tail.slope.index(), 43);
  EXPECT_EQ(i, 44);
}

TEST_F(Tail, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 0);
  EXPECT_EQ(tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  tail.amplitude.to_fit = false;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), -1);
  EXPECT_EQ(tail.slope.index(), 2);
  EXPECT_EQ(i, 3);

  tail.amplitude.to_fit = true;
  tail.slope.to_fit = false;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 3);
  EXPECT_EQ(tail.slope.index(), -1);
  EXPECT_EQ(i, 4);

  tail.slope.to_fit = true;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 4);
  EXPECT_EQ(tail.slope.index(), 5);
  EXPECT_EQ(i, 6);
}

TEST_F(Tail, UpdateIndexDisabled)
{
  tail.enabled = false;
  int32_t i;

  i = 0;
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), -1);
  EXPECT_EQ(tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), -1);
  EXPECT_EQ(tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  // \todo test resetting of indices
}

TEST_F(Tail, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);

  tail.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], tail.amplitude.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], tail.slope.x());

  int32_t i{0};
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 0);
  EXPECT_EQ(tail.slope.index(), 1);

  tail.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], tail.amplitude.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], tail.slope.x());
}

TEST_F(Tail, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  fit[0] = 0.005;
  fit[1] = 0.03;

  tail.get(fit);
  EXPECT_NEAR(tail.amplitude.val(), 0.5, 0.00001);
  EXPECT_NEAR(tail.slope.val(), 30, 0.00001);
  EXPECT_NE(tail.amplitude.val(),tail.amplitude.val_at(0.005));
  EXPECT_NE(tail.slope.val(), tail.slope.val_at(0.03));

  int32_t i{0};
  tail.update_indices(i);
  EXPECT_EQ(tail.amplitude.index(), 0);
  EXPECT_EQ(tail.slope.index(), 1);

  tail.get(fit);
  EXPECT_NE(tail.amplitude.val(), 0.5);
  EXPECT_EQ(tail.amplitude.val(), tail.amplitude.val_at(0.005));
  EXPECT_NE(tail.slope.val(), 30);
  EXPECT_EQ(tail.slope.val(), tail.slope.val_at(0.03));
}

TEST_F(Tail, EvalAt)
{
  auto pre = precalc_spoof(20);

  auto goal = tail.eval(pre);

  int32_t i{0};
  tail.update_indices(i);

  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  tail.put(fit);

  tail.amplitude.val(0.000001);
  tail.slope.val(0.000001);

  EXPECT_NE(tail.eval(pre), goal);
  EXPECT_EQ(tail.eval_at(pre, fit), goal);
}

TEST_F(Tail, EvalGrad)
{
  auto pre = precalc_spoof(10);

  tail.update_indices(var_count);

  Eigen::VectorXd grad;
  grad.setConstant(var_count, 0.0);

  auto result = tail.eval_grad(pre, grad);

  EXPECT_EQ(result, tail.eval(pre));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);
  EXPECT_NE(grad[3], 0.0);
  EXPECT_NE(grad[4], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Tail, EvalGradAt)
{
  auto pre = precalc_spoof(10);

  tail.update_indices(var_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(var_count, 0.0);
  tail.eval_grad(pre, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(var_count, 0.0);
  grad.setConstant(var_count, 0.0);

  tail.put(fit);
  tail.amplitude.val(0.000001);
  tail.slope.val(0.000001);

  auto result = tail.eval_grad_at(pre, fit, grad);

  EXPECT_EQ(result, tail.eval_at(pre, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
  EXPECT_EQ(grad[3], grad_goal[3]);
  EXPECT_EQ(grad[4], grad_goal[4]);
}

TEST_F(Tail, GradTailAmpOnePoint)
{
  double goal_val = tail.amplitude.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, tail.amplitude);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.002);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.002);
}

TEST_F(Tail, GradTailAmp)
{
  double goal_val = tail.amplitude.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, tail.amplitude);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.002);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.002);
}

TEST_F(Tail, GradTailSlopeOnePoint)
{
  double goal_val = tail.slope.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, tail.slope);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 1.0);
  check_gradients(true);
  // false inflection point, but maybe ok, only one data point is not reliable anyways
}

TEST_F(Tail, GradTailSlope)
{
  double goal_val = tail.slope.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, tail.slope, 0.05);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 1.0);
  check_gradients(true);
  // \todo false gradient inflection point here, even with multiple data points
}

TEST_F(Tail, GradWidthOnePoint)
{
  double goal_val = width.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, width);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.005);
  check_gradients(true);
}

TEST_F(Tail, GradWidth)
{
  double goal_val = width.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, width);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.005);
}

TEST_F(Tail, GradAmpOnePoint)
{
  double goal_val = amplitude.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 4);
  check_gradients(true);
  // false inflection point, but maybe ok, only one data point is not reliable anyways
}

TEST_F(Tail, GradAmp)
{
  double goal_val = amplitude.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 4);
  check_gradients(true);
  // \todo false gradient inflection point here, even with multiple data points
}

TEST_F(Tail, GradPosOnePoint)
{
  double goal_val = position.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 10, 10, position);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.00001);
}

TEST_F(Tail, GradPos)
{
  double goal_val = position.val();
  tail.update_indices(var_count);
  survey_grad(generate_data(40), 0, 40, position);
  EXPECT_NEAR(check_chi_sq(true), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(true), goal_val, 0.00001);
}
