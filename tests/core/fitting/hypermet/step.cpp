#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/weighted_data.h>

class Step : public TestBase
{
 protected:

  DAQuiri::Step step;

  DAQuiri::Value position;
  DAQuiri::Value amplitude;
  DAQuiri::Value width;

  virtual void SetUp()
  {
    int32_t i {0};

    amplitude.bound(0, 1000);
    amplitude.val(40);
    amplitude.update_index(i);

    width.bound(0.8, 5.0);
    width.val(2);
    width.update_index(i);

    position.bound(0, 40);
    position.val(20);
    position.update_index(i);

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

  double eval_grad(double chan, Eigen::VectorXd& chan_gradients) const
  {
    return step.eval_grad(precalc_spoof(chan), chan_gradients);
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

  void grad(const DAQuiri::WeightedData& wwdata,
      double from, double to,
      DAQuiri::Value& variable)
  {
    int32_t idx {3};
    step.update_indices(idx);
    ASSERT_EQ(idx, 4);

    double degrees_freedom = wwdata.data.size() - idx;
    ASSERT_EQ(degrees_freedom, 36);

    auto wdata = wwdata.subset(from, to);

    size_t chosen_var_idx = variable.index();

    std::vector<double> val_proxy;
    std::vector<double> val_val;
    std::vector<double> chi_sq_norm;
    std::vector<double> gradient;

    Eigen::VectorXd gradients;
    Eigen::VectorXd chan_gradients;
    for (double proxy = -M_PI; proxy < M_PI; proxy += 0.1)
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

        double FTotal = eval_grad(data.x, chan_gradients);
        double t3 = -2.0 * (data.y - FTotal) / square(data.weight_phillips_marlow);

        for (size_t var = 0; var < static_cast<size_t>(idx); ++var)
          gradients[var] += chan_gradients[var] * t3;
        Chisq += square((data.y - FTotal) / data.weight_phillips_marlow);
      }

      gradient.push_back(gradients[chosen_var_idx]);
      chi_sq_norm.push_back(Chisq / degrees_freedom);
    }


    auto min_chi = std::min_element(chi_sq_norm.begin(), chi_sq_norm.end());
    auto min_chi_i = std::distance(chi_sq_norm.begin(), min_chi);

    double min_abs = std::numeric_limits<double>::max();
    size_t grad_i = 0;
    for (size_t i=0; i < gradient.size(); ++i)
    {
      if (std::abs(gradient[i]) < min_abs)
      {
        grad_i = i;
        min_abs = std::abs(gradient[i]);
      }
    }


//    MESSAGE() << "chi_sq(proxy):\n" << visualize(val_proxy, chi_sq_norm, 100) << "\n";
    MESSAGE() << "chi_sq(val):\n" << visualize_all(val_val, chi_sq_norm, 100) << "\n";
    MESSAGE() << "min(chi_sq)=" << chi_sq_norm[min_chi_i] << " at val=" << val_val[min_chi_i] << "\n";
//    MESSAGE() << "gradient(proxy):\n" << visualize(val_proxy, gradient, 100) << "\n";
    MESSAGE() << "gradient(val):\n" << visualize_all(val_val, gradient, 100) << "\n";
    MESSAGE() << "min(abs(grad))=" << gradient[grad_i] << " at val=" << val_val[grad_i] << "\n";

  }
};

TEST_F(Step, CheckSetup)
{
  MESSAGE() << "Gauss-amp:\n" << amplitude.to_string() << "\n";
  MESSAGE() << "Gauss-pos:\n" << position.to_string() << "\n";
  MESSAGE() << "Gauss-width:\n" << width.to_string() << "\n";
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

TEST_F(Step, GradStepAmpOnePoint)
{
  auto wdata = generate_data();
  grad(wdata, 10, 10, step.amplitude);
}

TEST_F(Step, GradStepAmp)
{
  auto wdata = generate_data();
  grad(wdata, 0, 40, step.amplitude);
}

TEST_F(Step, GradWidthOnePoint)
{
  auto wdata = generate_data();
  grad(wdata, 10, 10, width);
}

TEST_F(Step, GradWidth)
{
  auto wdata = generate_data();
  grad(wdata, 0, 40, width);
}

TEST_F(Step, GradAmpOnePoint)
{
  auto wdata = generate_data();
  grad(wdata, 10, 10, amplitude);
}

TEST_F(Step, GradAmp)
{
  auto wdata = generate_data();
  grad(wdata, 0, 40, amplitude);
}

TEST_F(Step, GradPosOnePoint)
{
  auto wdata = generate_data();
  grad(wdata, 10, 10, position);
}

TEST_F(Step, GradPos)
{
  auto wdata = generate_data();
  grad(wdata, 0, 40, position);
}
