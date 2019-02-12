#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/weighted_data.h>

class Step : public TestBase
{
 protected:

  virtual void SetUp()
  {
    step.amplitude.bound(0.000001, 0.05);
    step.amplitude.val(0.05);
  }

  DAQuiri::PrecalcVals precalc_spoof(double chan)
  {
    DAQuiri::Value position;
    position.bound(0, 40);
    position.val(20);
    EXPECT_DOUBLE_EQ(position.val(), 20.0);

    DAQuiri::Value amplitude;
    amplitude.bound(0, 1000);
    amplitude.val(40);
    EXPECT_DOUBLE_EQ(amplitude.val(), 40.0);

    DAQuiri::Value width_;
    width_.bound(0.8, 5.0);
    width_.val(2);
    EXPECT_DOUBLE_EQ(width_.val(), 2.0);

    DAQuiri::PrecalcVals ret;
    ret.ampl = amplitude.val();
    ret.half_ampl = 0.5 * ret.ampl;
    ret.width = width_.val();
    ret.spread = (chan - position.val()) / ret.width;

    ret.amp_grad = amplitude.grad();
    ret.width_grad = width_.grad();
    ret.pos_grad = position.grad();

    ret.i_amp = 0;
    ret.i_width = 1;
    ret.i_pos = 2; // should not matter?
    return ret;
  }

  DAQuiri::WeightedData generate_data()
  {
    std::vector<double> channels;
    std::vector<double> y;
    for (size_t i = 0; i < 40; ++i)
    {
      channels.push_back(i);
      auto pre = precalc_spoof(i);
      y.push_back(step.eval(pre));
    }
    return DAQuiri::WeightedData(channels, y);
  }

  void grad_at_point(DAQuiri::Value& variable, size_t channel_idx)
  {
    auto wdata = generate_data();

    int32_t idx {3};
    step.update_indices(idx);
    EXPECT_EQ(idx, 4);

    double degrees_freedom = wdata.data.size() - idx;
    EXPECT_EQ(degrees_freedom, 36);

    auto chosen_point = wdata.data.at(channel_idx);
    auto chosen_pre = precalc_spoof(chosen_point.x);
    EXPECT_EQ(chosen_point.x, channel_idx);
    EXPECT_EQ(chosen_point.y, step.eval(chosen_pre));

    size_t chosen_var_idx = variable.index();

    std::vector<double> val_proxy;
    std::vector<double> val_val;
    std::vector<double> chi_sq_norm;
    std::vector<double> gradient;

    Eigen::VectorXd chan_gradients;
    for (double proxy = -4; proxy < 4; proxy += 0.2)
    {
      val_proxy.push_back(proxy);
      val_val.push_back(variable.val_at(proxy));

      variable.x(proxy);
      chan_gradients.setConstant(idx, 0.0);
      auto pre = precalc_spoof(chosen_point.x);
      double FTotal = step.eval_grad(pre, chan_gradients);
      double t3 = -2.0 * (chosen_point.y - FTotal) / square(chosen_point.weight_phillips_marlow);

      gradient.push_back(chan_gradients[chosen_var_idx] * t3);
      chi_sq_norm.push_back(
          square((chosen_point.y - FTotal) / chosen_point.weight_phillips_marlow)
              / degrees_freedom);
    }

    MESSAGE() << "chi_sq_norm(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(proxy):\n" << visualize(val_proxy, gradient, 100) << "\n";
    MESSAGE() << "sq_norm(val):\n" << visualize(val_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize(val_val, gradient, 100) << "\n";
  }

  void grad(DAQuiri::Value& variable)
  {
    auto wdata = generate_data();

    int32_t idx {3};
    step.update_indices(idx);
    ASSERT_EQ(idx, 4);

    double degrees_freedom = wdata.data.size() - idx;
    ASSERT_EQ(degrees_freedom, 36);

    size_t chosen_var_idx = variable.index();

    std::vector<double> val_proxy;
    std::vector<double> val_val;
    std::vector<double> chi_sq_norm;
    std::vector<double> gradient;

    Eigen::VectorXd gradients;
    Eigen::VectorXd chan_gradients;
    for (double proxy = -4; proxy < 4; proxy += 0.2)
    {
      val_proxy.push_back(proxy);
      val_val.push_back(variable.val_at(proxy));
      variable.x(proxy);

      gradients.setConstant(idx, 0.0);
      Eigen::VectorXd chan_gradients;
      double Chisq{0.0};
      for (const auto& data : wdata.data)
      {
        chan_gradients.setConstant(idx, 0.0);

        auto pre = precalc_spoof(data.x);

        double FTotal = step.eval_grad(pre, chan_gradients);

        double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);

        for (size_t var = 0; var < static_cast<size_t>(idx); ++var)
          gradients[var] += chan_gradients[var] * t3;
        Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
      }

      gradient.push_back(gradients[chosen_var_idx]);
      chi_sq_norm.push_back(Chisq / degrees_freedom);
    }


    MESSAGE() << "chi_sq_norm(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(proxy):\n" << visualize(val_proxy, gradient, 100) << "\n";
    MESSAGE() << "chi_sq_norm(val):\n" << visualize(val_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize(val_val, gradient, 100) << "\n";
  }

  DAQuiri::Step step;
};

TEST_F(Step, CheckSetup)
{
  MESSAGE() << "Step:\n" << step.to_string() << "\n";
}

TEST_F(Step, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    auto pre = precalc_spoof(i);
    y.push_back(step.eval(pre));
  }
  MESSAGE() << "peak(channel):\n" << visualize(channels, y, 100) << "\n";
}

TEST_F(Step, WithinBounds)
{
  auto data = generate_data();
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
  auto data = generate_data();
  EXPECT_NEAR(data.data.front().y, 2.0, 1e-15);
  EXPECT_NEAR(data.data.back().y, 0.0, 1e-40);
}

TEST_F(Step, RightOriented)
{
  step.side = DAQuiri::Side::right;
  auto data = generate_data();
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

  i=42;
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

  int32_t i{3};
  step.update_indices(i);

  Eigen::VectorXd grad;
  grad.setConstant(i, 0.0);

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

TEST_F(Step, GradAmplitudeAtOnePoint)
{
  grad_at_point(step.amplitude, 10);
}

TEST_F(Step, GradAmplitudeOnly)
{
  grad(step.amplitude);
}

//TEST_F(Step, GradWidthAtOnePoint)
//{
//  grad_at_point(peak.width_, 5);
//}
//
//TEST_F(Step, GradWidthOnly)
//{
//  grad(peak.width_);
//}
//
//TEST_F(Step, GradPositionAtOnePoint)
//{
//  grad_at_point(peak.position, 10);
//}
//
//TEST_F(Step, GradPositionOnly)
//{
//  grad(peak.position);
//}
