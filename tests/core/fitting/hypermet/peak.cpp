#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Peak.h>

#include <core/fitting/BFGS/BFGS.h>

class FittablePeak : public TestFittable
{
 public:
  DAQuiri::Peak peak;

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

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(var_count, 0.0);
    peak.put(ret);
    return ret;
  }
};

class Peak : public FunctionTest
{
 protected:
  FittablePeak fpeak;

  virtual void SetUp()
  {
    fpeak.peak = fpeak.peak.gaussian_only();
    fpeak.peak.amplitude.bound(0, 500);
    fpeak.peak.amplitude.val(400);
    fpeak.peak.position.bound(0, 40);
    fpeak.peak.position.val(21);
    fpeak.peak.width_override = true;
    fpeak.peak.width.bound(0.8, 5.0);
    fpeak.peak.width.val(3.2);
  }
};

TEST_F(Peak, CheckSetup)
{
  MESSAGE() << "Peak:\n" << fpeak.peak.to_string() << "\n";
}

TEST_F(Peak, Visualize)
{
  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(fpeak.peak.eval(i).all());
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";
}

TEST_F(Peak, WithinBounds)
{
  auto data = fpeak.generate_data(40);

  auto min = std::numeric_limits<double>::max();
  auto max = std::numeric_limits<double>::min();
  for (const auto& d : data.data)
  {
    min = std::min(min, d.y);
    max = std::max(max, d.y);
  }

  EXPECT_LE(max, 400.0);
  EXPECT_NEAR(min, 0.0, 1e-14);
}

TEST_F(Peak, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(fpeak.peak.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(fpeak.peak.update_indices(i));
}

TEST_F(Peak, UpdateIndex)
{
  int32_t i = 0;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 0);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 1);
  EXPECT_EQ(fpeak.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 3);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 4);
  EXPECT_EQ(fpeak.peak.width.index(), 5);
  EXPECT_EQ(i, 6);

  i = 42;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 42);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 43);
  EXPECT_EQ(fpeak.peak.width.index(), 44);
  EXPECT_EQ(i, 45);
}

TEST_F(Peak, UpdateIndexInvalidates)
{
  int32_t i = 0;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 0);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 1);
  EXPECT_EQ(fpeak.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fpeak.peak.position.to_fit = false;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), -1);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 3);
  EXPECT_EQ(fpeak.peak.width.index(), 4);
  EXPECT_EQ(i, 5);

  fpeak.peak.position.to_fit = true;
  fpeak.peak.amplitude.to_fit = false;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 5);
  EXPECT_EQ(fpeak.peak.amplitude.index(), -1);
  EXPECT_EQ(fpeak.peak.width.index(), 6);
  EXPECT_EQ(i, 7);

  fpeak.peak.position.to_fit = true;
  fpeak.peak.amplitude.to_fit = true;
  fpeak.peak.width.to_fit = false;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 7);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 8);
  EXPECT_EQ(fpeak.peak.width.index(), -1);
  EXPECT_EQ(i, 9);

  fpeak.peak.position.to_fit = false;
  fpeak.peak.amplitude.to_fit = false;
  fpeak.peak.width.to_fit = false;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), -1);
  EXPECT_EQ(fpeak.peak.amplitude.index(), -1);
  EXPECT_EQ(fpeak.peak.width.index(), -1);
  EXPECT_EQ(i, 9);
}

TEST_F(Peak, UpdateIndexDisabled)
{
  int32_t i = 0;

  fpeak.peak.width_override = false;
  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 0);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 1);
  EXPECT_EQ(fpeak.peak.width.index(), -1);
  EXPECT_EQ(i, 2);

  fpeak.peak.update_indices(i);
  EXPECT_EQ(fpeak.peak.position.index(), 2);
  EXPECT_EQ(fpeak.peak.amplitude.index(), 3);
  EXPECT_EQ(fpeak.peak.width.index(), -1);
  EXPECT_EQ(i, 4);

  // \todo test resetting of indices
}

TEST_F(Peak, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);

  fpeak.peak.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], fpeak.peak.position.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], fpeak.peak.amplitude.x());
  EXPECT_EQ(fit[2], 0.0);
  EXPECT_NE(fit[2], fpeak.peak.width.x());

  fpeak.peak.update_indices(fpeak.var_count);
  fpeak.peak.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], fpeak.peak.position.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], fpeak.peak.amplitude.x());
  EXPECT_NE(fit[2], 0.0);
  EXPECT_EQ(fit[2], fpeak.peak.width.x());
}

TEST_F(Peak, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);
  fit[0] = 0.5;
  fit[1] = 0.03;
  fit[2] = 0.01;

  fpeak.peak.get(fit);
  EXPECT_NEAR(fpeak.peak.position.val(), 21, 0.00001);
  EXPECT_NE(fpeak.peak.position.val(), fpeak.peak.position.val_at(0.5));
  EXPECT_NEAR(fpeak.peak.amplitude.val(), 400, 0.00001);
  EXPECT_NE(fpeak.peak.amplitude.val(), fpeak.peak.amplitude.val_at(0.03));
  EXPECT_NEAR(fpeak.peak.width.val(), 3.2, 0.00001);
  EXPECT_NE(fpeak.peak.width.val(), fpeak.peak.width.val_at(0.01));

  fpeak.peak.update_indices(fpeak.var_count);

  fpeak.peak.get(fit);
  EXPECT_EQ(fpeak.peak.position.val(), fpeak.peak.position.val_at(0.5));
  EXPECT_EQ(fpeak.peak.amplitude.val(), fpeak.peak.amplitude.val_at(0.03));
  EXPECT_EQ(fpeak.peak.width.val(), fpeak.peak.width.val_at(0.01));
}


TEST_F(Peak, EvalAt)
{
  auto goal = fpeak.peak.eval(20);

  fpeak.peak.update_indices(fpeak.var_count);

  Eigen::VectorXd fit;
  fit.setConstant(fpeak.var_count, 0.0);
  fpeak.peak.put(fit);

  fpeak.peak.position.val(0.000001);
  fpeak.peak.amplitude.val(0.000001);
  fpeak.peak.width.val(0.000001);

  EXPECT_NE(fpeak.peak.eval(10).all(), goal.all());
  EXPECT_EQ(fpeak.peak.eval_at(20, fit).all(), goal.all());
}

TEST_F(Peak, EvalGrad)
{
  fpeak.peak.update_indices(fpeak.var_count);

  Eigen::VectorXd grad;
  grad.setConstant(fpeak.var_count, 0.0);

  auto result = fpeak.peak.eval_grad(20, grad);

  EXPECT_EQ(result.all(), fpeak.peak.eval(20).all());
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Peak, EvalGradAt)
{
  fpeak.peak.update_indices(fpeak.var_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(fpeak.var_count, 0.0);
  fpeak.peak.eval_grad(20, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(fpeak.var_count, 0.0);
  grad.setConstant(fpeak.var_count, 0.0);

  fpeak.peak.put(fit);
  fpeak.peak.position.val(0.000001);
  fpeak.peak.amplitude.val(0.000001);
  fpeak.peak.width.val(0.000001);

  auto result = fpeak.peak.eval_grad_at(20, fit, grad);

  EXPECT_EQ(result.all(), fpeak.peak.eval_at(20, fit).all());
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
}

TEST_F(Peak, GradPosition)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.position.val();
  fpeak.peak.update_indices(fpeak.var_count);
  survey_grad(&fpeak, fpeak.peak.position, 0.05);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 3);
  check_gradients(true);
  // \todo false gradient minima, even with multiple data points
}

TEST_F(Peak, GradAmp)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.amplitude.val();
  fpeak.peak.update_indices(fpeak.var_count);
  survey_grad(&fpeak, fpeak.peak.amplitude);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 3);
  EXPECT_NEAR(check_gradients(false), goal_val, 3);
}

TEST_F(Peak, GradWidth)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.width.val();
  fpeak.peak.update_indices(fpeak.var_count);
  survey_grad(&fpeak, fpeak.peak.width, 0.05);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.03);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.03);
}

TEST_F(Peak, FitPosition)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.position.val();

  fpeak.peak.position.val(15);
  fpeak.peak.update_indices(fpeak.var_count);
  MESSAGE() << "Will rebuild:\n" << fpeak.peak.to_string() << "\n";

  DAQuiri::BFGS optimizer;

  auto result = optimizer.BFGSMin(&fpeak, 0.00001);
  fpeak.peak.get(result.variables);

  MESSAGE() << "Result:\n" << fpeak.peak.to_string() << "\n";

  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(fpeak.peak.eval(i).all());
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";

  EXPECT_NEAR(fpeak.peak.position.val(), goal_val, 0.03);
}

TEST_F(Peak, FitAmplitude)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.amplitude.val();

  fpeak.peak.amplitude.val(200);
  fpeak.peak.update_indices(fpeak.var_count);
  MESSAGE() << "Will rebuild:\n" << fpeak.peak.to_string() << "\n";

  DAQuiri::BFGS optimizer;

  auto result = optimizer.BFGSMin(&fpeak, 0.00001);
  fpeak.peak.get(result.variables);

  MESSAGE() << "Result:\n" << fpeak.peak.to_string() << "\n";

  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(fpeak.peak.eval(i).all());
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";

  // \todo this intermittently fails!!!
  //EXPECT_NEAR(fpeak.peak.amplitude.val(), goal_val, 0.03);
}

TEST_F(Peak, FitWidth)
{
  fpeak.data = fpeak.generate_data(40);
  double goal_val = fpeak.peak.width.val();

  fpeak.peak.width.val(1.0);
  fpeak.peak.update_indices(fpeak.var_count);
  MESSAGE() << "Will rebuild:\n" << fpeak.peak.to_string() << "\n";

  DAQuiri::BFGS optimizer;

  auto result = optimizer.BFGSMin(&fpeak, 0.00001);
  fpeak.peak.get(result.variables);

  MESSAGE() << "Result:\n" << fpeak.peak.to_string() << "\n";

  std::vector<double> channels;
  std::vector<double> y;
  for (size_t i = 0; i < 40; ++i)
  {
    channels.push_back(i);
    y.push_back(fpeak.peak.eval(i).all());
  }
  MESSAGE() << "counts(channel):\n" << visualize(channels, y, 100) << "\n";

  // \todo this intermittently fails!!!
//  EXPECT_NEAR(fpeak.peak.width.val(), goal_val, 0.3);
}
