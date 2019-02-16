#include "../function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Tail.h>
#include <core/fitting/weighted_data.h>

class FittableTail : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::Value position;
  DAQuiri::Value amplitude;
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
  FittableTail ftail;

  virtual void SetUp()
  {
    ftail.amplitude.bound(0, 1000);
    ftail.amplitude.val(40);
    ftail.amplitude.update_index(ftail.variable_count);

    ftail.width.bound(0.8, 5.0);
    ftail.width.val(2);
    ftail.width.update_index(ftail.variable_count);

    ftail.position.bound(0, 40);
    ftail.position.val(20);
    ftail.position.update_index(ftail.variable_count);

    ftail.tail.amplitude.bound(0.0001, 1.5);
    ftail.tail.amplitude.val(0.5);
    ftail.tail.slope.bound(0.2, 50);
    ftail.tail.slope.val(30);
  }
};

TEST_F(Tail, CheckSetup)
{
  MESSAGE() << "Gaussian amp: " << ftail.amplitude.to_string() << "\n";
  MESSAGE() << "Gaussian pos: " << ftail.position.to_string() << "\n";
  MESSAGE() << "Gaussian width: " << ftail.width.to_string() << "\n";
  MESSAGE() << "Tail: " << ftail.tail.to_string() << "\n";
}

TEST_F(Tail, Visualize)
{
  auto data = generate_data(&ftail, 40);
  visualize_data(data);
}


TEST_F(Tail, WithinBounds)
{
  auto data = generate_data(&ftail, 40);
  EXPECT_NEAR(data.count_min, 0.0, 1e-39);
  EXPECT_LE(data.count_max, 20.0);
}

TEST_F(Tail, LeftOriented)
{
  ftail.tail.side = DAQuiri::Side::left;
  auto data = generate_data(&ftail, 40);
  EXPECT_NEAR(data.data.front().count, 14.0, 1.0);
  EXPECT_NEAR(data.data.back().count, 0.0, 1e-39);
}

TEST_F(Tail, RightOriented)
{
  ftail.tail.side = DAQuiri::Side::right;
  auto data = generate_data(&ftail, 40);
  EXPECT_NEAR(data.data.front().count, 0.0, 1e-39);
  EXPECT_NEAR(data.data.back().count, 14.0, 1.0);
}

TEST_F(Tail, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(ftail.tail.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(ftail.tail.update_indices(i));
}

TEST_F(Tail, UpdateIndex)
{
  int32_t i;

  i = 0;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 0);
  EXPECT_EQ(ftail.tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 2);
  EXPECT_EQ(ftail.tail.slope.index(), 3);
  EXPECT_EQ(i, 4);

  i = 42;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 42);
  EXPECT_EQ(ftail.tail.slope.index(), 43);
  EXPECT_EQ(i, 44);
}

TEST_F(Tail, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 0);
  EXPECT_EQ(ftail.tail.slope.index(), 1);
  EXPECT_EQ(i, 2);

  ftail.tail.amplitude.to_fit = false;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), -1);
  EXPECT_EQ(ftail.tail.slope.index(), 2);
  EXPECT_EQ(i, 3);

  ftail.tail.amplitude.to_fit = true;
  ftail.tail.slope.to_fit = false;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 3);
  EXPECT_EQ(ftail.tail.slope.index(), -1);
  EXPECT_EQ(i, 4);

  ftail.tail.slope.to_fit = true;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 4);
  EXPECT_EQ(ftail.tail.slope.index(), 5);
  EXPECT_EQ(i, 6);
}

TEST_F(Tail, UpdateIndexDisabled)
{
  ftail.tail.enabled = false;
  int32_t i;

  i = 0;
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), -1);
  EXPECT_EQ(ftail.tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), -1);
  EXPECT_EQ(ftail.tail.slope.index(), -1);
  EXPECT_EQ(i, 0);

  // \todo test resetting of indices
}

TEST_F(Tail, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);

  ftail.tail.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], ftail.tail.amplitude.x());
  EXPECT_EQ(fit[1], 0.0);
  EXPECT_NE(fit[1], ftail.tail.slope.x());

  int32_t i{0};
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 0);
  EXPECT_EQ(ftail.tail.slope.index(), 1);

  ftail.tail.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], ftail.tail.amplitude.x());
  EXPECT_NE(fit[1], 0.0);
  EXPECT_EQ(fit[1], ftail.tail.slope.x());
}

TEST_F(Tail, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  fit[0] = 0.005;
  fit[1] = 0.03;

  ftail.tail.get(fit);
  EXPECT_NEAR(ftail.tail.amplitude.val(), 0.5, 0.00001);
  EXPECT_NEAR(ftail.tail.slope.val(), 30, 0.00001);
  EXPECT_NE(ftail.tail.amplitude.val(),ftail.tail.amplitude.val_at(0.005));
  EXPECT_NE(ftail.tail.slope.val(), ftail.tail.slope.val_at(0.03));

  int32_t i{0};
  ftail.tail.update_indices(i);
  EXPECT_EQ(ftail.tail.amplitude.index(), 0);
  EXPECT_EQ(ftail.tail.slope.index(), 1);

  ftail.tail.get(fit);
  EXPECT_NE(ftail.tail.amplitude.val(), 0.5);
  EXPECT_EQ(ftail.tail.amplitude.val(), ftail.tail.amplitude.val_at(0.005));
  EXPECT_NE(ftail.tail.slope.val(), 30);
  EXPECT_EQ(ftail.tail.slope.val(), ftail.tail.slope.val_at(0.03));
}

TEST_F(Tail, EvalAt)
{
  auto pre = ftail.precalc(20);

  auto goal = ftail.tail.eval(pre);

  int32_t i{0};
  ftail.tail.update_indices(i);

  Eigen::VectorXd fit;
  fit.setConstant(2, 0.0);
  ftail.tail.put(fit);

  ftail.tail.amplitude.val(0.000001);
  ftail.tail.slope.val(0.000001);

  EXPECT_NE(ftail.tail.eval(pre), goal);
  EXPECT_EQ(ftail.tail.eval_at(pre, fit), goal);
}

TEST_F(Tail, EvalGrad)
{
  auto pre = ftail.precalc(10);

  ftail.tail.update_indices(ftail.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(ftail.variable_count, 0.0);

  auto result = ftail.tail.eval_grad(pre, grad);

  EXPECT_EQ(result, ftail.tail.eval(pre));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);
  EXPECT_NE(grad[3], 0.0);
  EXPECT_NE(grad[4], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Tail, EvalGradAt)
{
  auto pre = ftail.precalc(10);

  ftail.tail.update_indices(ftail.variable_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(ftail.variable_count, 0.0);
  ftail.tail.eval_grad(pre, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(ftail.variable_count, 0.0);
  grad.setConstant(ftail.variable_count, 0.0);

  ftail.tail.put(fit);
  ftail.tail.amplitude.val(0.000001);
  ftail.tail.slope.val(0.000001);

  auto result = ftail.tail.eval_grad_at(pre, fit, grad);

  EXPECT_EQ(result, ftail.tail.eval_at(pre, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
  EXPECT_EQ(grad[3], grad_goal[3]);
  EXPECT_EQ(grad[4], grad_goal[4]);
}


TEST_F(Tail, GradTailAmp)
{
  ftail.data = generate_data(&ftail, 40);

  double goal_val = ftail.tail.amplitude.val();
  ftail.tail.update_indices(ftail.variable_count);
  survey_grad(&ftail, ftail.tail.amplitude);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.002);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.002);
}

TEST_F(Tail, GradTailSlope)
{
  ftail.data = generate_data(&ftail, 40);

  double goal_val = ftail.tail.slope.val();
  ftail.tail.update_indices(ftail.variable_count);
  survey_grad(&ftail, ftail.tail.slope, 0.01);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.01);
}

TEST_F(Tail, GradWidth)
{
  ftail.data = generate_data(&ftail, 40);

  double goal_val = ftail.width.val();
  ftail.tail.update_indices(ftail.variable_count);
  survey_grad(&ftail, ftail.width);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.005);
}

TEST_F(Tail, GradAmp)
{
  ftail.data = generate_data(&ftail, 40);

  double goal_val = ftail.amplitude.val();
  ftail.tail.update_indices(ftail.variable_count);
  survey_grad(&ftail, ftail.amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 4);
  check_gradients(true);
  // \todo false gradient inflection point here, even with multiple data points
}

TEST_F(Tail, GradPos)
{
  ftail.data = generate_data(&ftail, 40);

  double goal_val = ftail.position.val();
  ftail.tail.update_indices(ftail.variable_count);
  survey_grad(&ftail, ftail.position);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.00001);
}
