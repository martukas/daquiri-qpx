#include "../function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/weighted_data.h>

class FittableStep : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::Value position;
  DAQuiri::Value amplitude;
  DAQuiri::Value width;

  DAQuiri::Step step;

  void update_indices() override
  {
    variable_count = 0;
    amplitude.update_index(variable_count);
    width.update_index(variable_count);
    position.update_index(variable_count);
    step.update_indices(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    position.put(ret);
    amplitude.put(ret);
    width.put(ret);
    step.put(ret);
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
    ret.i_pos = position.index(); // should not matter?
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
    return step.eval(precalc(chan));
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return step.eval_at(precalc_at(chan, fit), fit);
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    return step.eval_grad_at(precalc_at(chan, fit), fit, grads);
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    amplitude.get(result.variables);
    width.get(result.variables);
    position.get(result.variables);
    step.get(result.variables);

    // \todo uncerts
  }
};

class Step : public FunctionTest
{
 protected:
  FittableStep fstep;

  virtual void SetUp()
  {
    fstep.amplitude.bound(0, 1000);
    fstep.amplitude.val(40);
    fstep.amplitude.update_index(fstep.variable_count);

    fstep.width.bound(0.8, 5.0);
    fstep.width.val(2);
    fstep.width.update_index(fstep.variable_count);

    fstep.position.bound(0, 40);
    fstep.position.val(20);
    fstep.position.update_index(fstep.variable_count);

    fstep.step.amplitude.bound(0.000001, 0.05);
    fstep.step.amplitude.val(0.05);
  }
};

TEST_F(Step, CheckSetup)
{
  MESSAGE() << "Gaussian amp: " << fstep.amplitude.to_string() << "\n";
  MESSAGE() << "Gaussian pos: " << fstep.position.to_string() << "\n";
  MESSAGE() << "Gaussian width: " << fstep.width.to_string() << "\n";
  MESSAGE() << "Step: " << fstep.step.to_string() << "\n";
}

TEST_F(Step, Visualize)
{
  auto data = generate_data(&fstep, 40);
  visualize_data(data);
}

TEST_F(Step, WithinBounds)
{
  auto data = generate_data(&fstep, 40);
  EXPECT_NEAR(data.count_min, 0.0, 1e-40);
  EXPECT_NEAR(data.count_max, 2.0, 1e-15);
}

TEST_F(Step, LeftOriented)
{
  fstep.step.side = DAQuiri::Side::left;
  auto data = generate_data(&fstep, 40);
  EXPECT_NEAR(data.data.front().count, 2.0, 1e-15);
  EXPECT_NEAR(data.data.back().count, 0.0, 1e-40);
}

TEST_F(Step, RightOriented)
{
  fstep.step.side = DAQuiri::Side::right;
  auto data = generate_data(&fstep, 40);
  EXPECT_NEAR(data.data.front().count, 0.0, 1e-40);
  EXPECT_NEAR(data.data.back().count, 2.0, 1e-15);
}

TEST_F(Step, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(fstep.step.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(fstep.step.update_indices(i));
}

TEST_F(Step, UpdateIndex)
{
  int32_t i;

  i = 0;
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 0);
  EXPECT_EQ(i, 1);

  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 1);
  EXPECT_EQ(i, 2);

  i = 42;
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 42);
  EXPECT_EQ(i, 43);
}

TEST_F(Step, UpdateIndexInvalidates)
{
  int32_t i;

  i = 0;
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 0);
  EXPECT_EQ(i, 1);

  fstep.step.amplitude.to_fit = false;
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), -1);
  EXPECT_EQ(i, 1);

  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), -1);
  EXPECT_EQ(i, 1);
}

TEST_F(Step, UpdateIndexDisabled)
{
  fstep.step.enabled = false;
  int32_t i;

  i = 0;
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), -1);
  EXPECT_EQ(i, 0);

  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), -1);
  EXPECT_EQ(i, 0);

  // \todo test resetting of indices
}

TEST_F(Step, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(1, 0.0);

  fstep.step.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], fstep.step.amplitude.x());

  int32_t i{0};
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 0);

  fstep.step.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], fstep.step.amplitude.x());
}

TEST_F(Step, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(1, 0.005);

  fstep.step.get(fit);
  EXPECT_EQ(fstep.step.amplitude.val(), 0.05);

  int32_t i{0};
  fstep.step.update_indices(i);
  EXPECT_EQ(fstep.step.amplitude.index(), 0);

  fstep.step.get(fit);
  EXPECT_NE(fstep.step.amplitude.val(), 0.05);
  EXPECT_EQ(fstep.step.amplitude.val(), fstep.step.amplitude.val_at(0.005));
}

TEST_F(Step, EvalAt)
{
  auto pre = fstep.precalc(20);

  auto goal = fstep.step.eval(pre);

  int32_t i{0};
  fstep.step.update_indices(i);

  Eigen::VectorXd fit;
  fit.setConstant(1, 0.0);
  fstep.step.put(fit);

  fstep.step.amplitude.val(0.000001);

  EXPECT_NE(fstep.step.eval(pre), goal);
  EXPECT_EQ(fstep.step.eval_at(pre, fit), goal);
}

TEST_F(Step, EvalGrad)
{
  auto pre = fstep.precalc(10);

  fstep.step.update_indices(fstep.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(fstep.variable_count, 0.0);

  auto result = fstep.step.eval_grad(pre, grad);

  EXPECT_EQ(result, fstep.step.eval(pre));
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_EQ(grad[2], 0.0); // pos gradient should be unaffected?
  EXPECT_NE(grad[3], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Step, EvalGradAt)
{
  auto pre = fstep.precalc(10);

  int32_t i{3};
  fstep.step.update_indices(i);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(i, 0.0);
  fstep.step.eval_grad(pre, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(i, 0.0);
  grad.setConstant(i, 0.0);

  fstep.step.put(fit);
  fstep.step.amplitude.val(0.000001);

  auto result = fstep.step.eval_grad_at(pre, fit, grad);

  EXPECT_EQ(result, fstep.step.eval_at(pre, fit));
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
  EXPECT_EQ(grad[3], grad_goal[3]);
}

TEST_F(Step, GradStepAmp)
{
  fstep.data = generate_data(&fstep, 40);

  double goal_val = fstep.step.amplitude.val();
  fstep.step.update_indices(fstep.variable_count);
  survey_grad(&fstep, &fstep.step.amplitude);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.00001);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.00001);
}

TEST_F(Step, GradWidth)
{
  fstep.data = generate_data(&fstep, 40);

  double goal_val = fstep.width.val();
  fstep.step.update_indices(fstep.variable_count);
  survey_grad(&fstep, &fstep.width);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.005);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.005);
}

TEST_F(Step, GradAmp)
{
  fstep.data = generate_data(&fstep, 40);

  double goal_val = fstep.amplitude.val();
  fstep.step.update_indices(fstep.variable_count);
  survey_grad(&fstep, &fstep.amplitude, 0.05);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 4);
  EXPECT_NEAR(check_gradients(false), goal_val, 4);
}

TEST_F(Step, GradPos)
{
  fstep.data = generate_data(&fstep, 40);

  double goal_val = fstep.position.val();
  fstep.step.update_indices(fstep.variable_count);
  survey_grad(&fstep, &fstep.position);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.00001);
  check_gradients(true);
  // \todo gradient not affected!
}
