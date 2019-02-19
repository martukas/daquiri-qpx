#include "../function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Peak.h>

#include <core/fitting/optimizers/dlib_adapter.h>

class FittablePeak : public DAQuiri::FittableRegion
{
 public:
  DAQuiri::Peak peak;

  void update_indices() override
  {
    variable_count = 0;
    peak.update_indices(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    peak.put(ret);
    return ret;
  }

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

  void save_fit(const DAQuiri::FitResult& result) override
  {
    peak.get(result.variables);
    // \todo uncerts
  }

};

class Peak : public FunctionTest
{
 protected:
  FittablePeak fp;

  virtual void SetUp()
  {
    fp.peak = fp.peak.gaussian_only();
    //fp.peak.amplitude.bound(0, 500);
    fp.peak.amplitude.val(4000);
    fp.peak.position.bound(14, 28);
    fp.peak.position.val(21);
    fp.peak.width_override = true;
    fp.peak.width.bound(0.8, 5.0);
    fp.peak.width.val(3.2);
  }
};

TEST_F(Peak, CheckSetup)
{
  MESSAGE() << "Peak:\n" << fp.peak.to_string() << "\n";
}

TEST_F(Peak, Visualize)
{
  auto data = generate_data(&fp, 40);
  visualize_data(data);
}

TEST_F(Peak, WithinBounds)
{
  auto data = generate_data(&fp, 40);
  EXPECT_NEAR(data.count_min, 0.0, 1e-14);
  EXPECT_NEAR(data.count_max, 4000.0, 1e-7);
}

TEST_F(Peak, UpdateIndexInvalidThrows)
{
  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(fp.peak.update_indices(i));

  i = -42;
  EXPECT_ANY_THROW(fp.peak.update_indices(i));
}

TEST_F(Peak, UpdateIndex)
{
  int32_t i = 0;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 3);
  EXPECT_EQ(fp.peak.amplitude.index(), 4);
  EXPECT_EQ(fp.peak.width.index(), 5);
  EXPECT_EQ(i, 6);

  i = 42;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 42);
  EXPECT_EQ(fp.peak.amplitude.index(), 43);
  EXPECT_EQ(fp.peak.width.index(), 44);
  EXPECT_EQ(i, 45);
}

TEST_F(Peak, UpdateIndexInvalidates)
{
  int32_t i = 0;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), 2);
  EXPECT_EQ(i, 3);

  fp.peak.position.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), -1);
  EXPECT_EQ(fp.peak.amplitude.index(), 3);
  EXPECT_EQ(fp.peak.width.index(), 4);
  EXPECT_EQ(i, 5);

  fp.peak.position.to_fit = true;
  fp.peak.amplitude.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 5);
  EXPECT_EQ(fp.peak.amplitude.index(), -1);
  EXPECT_EQ(fp.peak.width.index(), 6);
  EXPECT_EQ(i, 7);

  fp.peak.position.to_fit = true;
  fp.peak.amplitude.to_fit = true;
  fp.peak.width.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 7);
  EXPECT_EQ(fp.peak.amplitude.index(), 8);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 9);

  fp.peak.position.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.peak.width.to_fit = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), -1);
  EXPECT_EQ(fp.peak.amplitude.index(), -1);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 9);
}

TEST_F(Peak, UpdateIndexDisabled)
{
  int32_t i = 0;

  fp.peak.width_override = false;
  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 0);
  EXPECT_EQ(fp.peak.amplitude.index(), 1);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 2);

  fp.peak.update_indices(i);
  EXPECT_EQ(fp.peak.position.index(), 2);
  EXPECT_EQ(fp.peak.amplitude.index(), 3);
  EXPECT_EQ(fp.peak.width.index(), -1);
  EXPECT_EQ(i, 4);

  // \todo test resetting of indices
}

TEST_F(Peak, Put)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 1.0);

  fp.peak.put(fit);
  EXPECT_EQ(fit[0], 1.0);
  EXPECT_NE(fit[0], fp.peak.position.x());
  EXPECT_EQ(fit[1], 1.0);
  EXPECT_NE(fit[1], fp.peak.amplitude.x());
  EXPECT_EQ(fit[2], 1.0);
  EXPECT_NE(fit[2], fp.peak.width.x());

  fp.peak.update_indices(fp.variable_count);
  fp.peak.put(fit);
  EXPECT_NE(fit[0], 1.0);
  EXPECT_EQ(fit[0], fp.peak.position.x());
  EXPECT_NE(fit[1], 1.0);
  EXPECT_EQ(fit[1], fp.peak.amplitude.x());
  EXPECT_NE(fit[2], 1.0);
  EXPECT_EQ(fit[2], fp.peak.width.x());
}

TEST_F(Peak, Get)
{
  Eigen::VectorXd fit;
  fit.setConstant(3, 0.0);
  fit[0] = 0.5;
  fit[1] = 0.03;
  fit[2] = 0.01;

  fp.peak.get(fit);
  EXPECT_NEAR(fp.peak.position.val(), 21, 0.00001);
  EXPECT_NE(fp.peak.position.val(), fp.peak.position.val_at(0.5));
  EXPECT_NEAR(fp.peak.amplitude.val(), 4000, 0.00001);
  EXPECT_NE(fp.peak.amplitude.val(), fp.peak.amplitude.val_at(0.03));
  EXPECT_NEAR(fp.peak.width.val(), 3.2, 0.00001);
  EXPECT_NE(fp.peak.width.val(), fp.peak.width.val_at(0.01));

  fp.peak.update_indices(fp.variable_count);

  fp.peak.get(fit);
  EXPECT_EQ(fp.peak.position.val(), fp.peak.position.val_at(0.5));
  EXPECT_EQ(fp.peak.amplitude.val(), fp.peak.amplitude.val_at(0.03));
  EXPECT_EQ(fp.peak.width.val(), fp.peak.width.val_at(0.01));
}

TEST_F(Peak, EvalAt)
{
  auto goal = fp.peak.eval(20);

  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd fit;
  fit.setConstant(fp.variable_count, 0.0);
  fp.peak.put(fit);

  fp.peak.position.val(0.000001);
  fp.peak.amplitude.val(0.000001);
  fp.peak.width.val(0.000001);

  EXPECT_NE(fp.peak.eval(10).all(), goal.all());
  EXPECT_EQ(fp.peak.eval_at(20, fit).all(), goal.all());
}

TEST_F(Peak, EvalGrad)
{
  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd grad;
  grad.setConstant(fp.variable_count, 0.0);

  auto result = fp.peak.eval_grad(20, grad);

  EXPECT_EQ(result.all(), fp.peak.eval(20).all());
  EXPECT_NE(grad[0], 0.0);
  EXPECT_NE(grad[1], 0.0);
  EXPECT_NE(grad[2], 0.0);

  // \todo confirm that gradient is meaningful?
}

TEST_F(Peak, EvalGradAt)
{
  fp.peak.update_indices(fp.variable_count);

  Eigen::VectorXd grad_goal;
  grad_goal.setConstant(fp.variable_count, 0.0);
  fp.peak.eval_grad(20, grad_goal);

  Eigen::VectorXd fit, grad;
  fit.setConstant(fp.variable_count, 0.0);
  grad.setConstant(fp.variable_count, 0.0);

  fp.peak.put(fit);
  fp.peak.position.val(0.000001);
  fp.peak.amplitude.val(0.000001);
  fp.peak.width.val(0.000001);

  auto result = fp.peak.eval_grad_at(20, fit, grad);

  EXPECT_EQ(result.all(), fp.peak.eval_at(20, fit).all());
  EXPECT_EQ(grad[0], grad_goal[0]);
  EXPECT_EQ(grad[1], grad_goal[1]);
  EXPECT_EQ(grad[2], grad_goal[2]);
}

TEST_F(Peak, GradPosition)
{
  double goal_val = fp.peak.position.val();
  fp.data = generate_data(&fp, 40);

  fp.peak.width.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();
  survey_grad(&fp, &fp.peak.position, 0.0001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.1);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.1);
}

TEST_F(Peak, FitPositionOnly)
{
  fp.data = generate_data(&fp, 40);

  fp.peak.width.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.position, 17, 0.5);

  fp.peak.position.val(21);
  test_fit_random(20, &optimizer, &fp, &fp.peak.position, 14, 28, 0.5);
}

TEST_F(Peak, FitPositionRelaxed)
{
  fp.data = generate_data(&fp, 40);
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.position, 17, 0.5);

  fp.peak.position.val(21);
  test_fit_random(20, &optimizer, &fp, &fp.peak.position, 14, 28, 0.5);
}

TEST_F(Peak, GradWidth)
{
  double goal_val = fp.peak.width.val();
  fp.data = generate_data(&fp, 40);

  fp.peak.position.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();
  survey_grad(&fp, &fp.peak.width, 0.001);
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.01);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.01);
}

TEST_F(Peak, FitWidthOnly)
{
  fp.data = generate_data(&fp, 40);

  fp.peak.position.to_fit = false;
  fp.peak.amplitude.to_fit = false;
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.width, 1.0, 0.5);

  fp.peak.width.val(3.2);
  test_fit_random(20, &optimizer, &fp, &fp.peak.width,
                  fp.peak.width.min(), fp.peak.width.max(), 0.5);
}

TEST_F(Peak, FitWidthRelaxed)
{
  fp.data = generate_data(&fp, 40);
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.width, 1.0, 0.5);

  fp.peak.width.val(3.2);
  test_fit_random(20, &optimizer, &fp, &fp.peak.width,
                  fp.peak.width.min(), fp.peak.width.max(), 0.5);
}

TEST_F(Peak, GradAmplitude)
{
  double goal_val = fp.peak.amplitude.val();
  fp.data = generate_data(&fp, 40);

  fp.peak.width.to_fit = false;
  fp.peak.position.to_fit = false;
  fp.update_indices();
  survey_grad(&fp, &fp.peak.amplitude, 0.001, std::sqrt(3000), std::sqrt(5000));
  EXPECT_NEAR(check_chi_sq(false), goal_val, 0.05);
  EXPECT_NEAR(check_gradients(false), goal_val, 0.05);
}

TEST_F(Peak, FitAmplitudeOnly)
{
  fp.data = generate_data(&fp, 40);

  fp.peak.width.to_fit = false;
  fp.peak.position.to_fit = false;
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.amplitude, 3000, 0.5);

  fp.peak.amplitude.val(4000);
  test_fit_random(20, &optimizer, &fp, &fp.peak.amplitude,
                  3000, 4000, 0.5);
}

TEST_F(Peak, FitAmplitudeRelaxed)
{
  fp.data = generate_data(&fp, 40);
  fp.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fp, &fp.peak.amplitude, 3000, 0.5);

  fp.peak.amplitude.val(4000);
  test_fit_random(20, &optimizer, &fp, &fp.peak.amplitude,
                  3000, 4000, 0.5);
}

TEST_F(Peak, FitAllThree)
{
  double goal_pos = fp.peak.position.val();
  double goal_width = fp.peak.width.val();
  double goal_amplitude = fp.peak.amplitude.val();

  fp.data = generate_data(&fp, 40);
  fp.update_indices();

  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());
  std::uniform_real_distribution<double> pos_dist(fp.peak.position.min(),
                                                  fp.peak.position.max());
  std::uniform_real_distribution<double> width_dist(fp.peak.width.min(),
                                                    fp.peak.width.max());
  std::uniform_real_distribution<double> amplitude_dist(3000, 5000);

  DAQuiri::DLibOptimizer optimizer;

  for (size_t i = 0; i < 50; ++i)
  {
    fp.peak.position.val(pos_dist(random_generator));
    fp.peak.width.val(width_dist(random_generator));
    fp.peak.amplitude.val(amplitude_dist(random_generator));
    MESSAGE() << "Attempt[" << i << "]\n" << fp.peak.to_string();
    fp.save_fit(optimizer.minimize(&fp));
    MESSAGE() << "Result:               \n"
              << fp.peak.to_string("                      ");
    MESSAGE() << "       position delta="
              << (goal_pos - fp.peak.position.val()) << "\n";
    MESSAGE() << "          width delta="
              << (goal_width - fp.peak.width.val()) << "\n";
    MESSAGE() << "      amplitude delta="
              << (goal_amplitude - fp.peak.amplitude.val()) << "\n";

    EXPECT_NEAR(fp.peak.position.val(), goal_pos, 0.001);
    EXPECT_NEAR(fp.peak.width.val(), goal_width, 0.001);
    EXPECT_NEAR(fp.peak.amplitude.val(), goal_amplitude, 0.01);
  }
}

TEST_F(Peak, WithSkews)
{
  fp.peak.short_tail.enabled = true;
  fp.peak.right_tail.enabled = true;
  fp.data = generate_data(&fp, 40);
  fp.update_indices();

  MESSAGE() << "Peak:\n" << fp.peak.to_string() << "\n";
  visualize_data(fp.data);

  double goal_pos = fp.peak.position.val();
  double goal_width = fp.peak.width.val();
  double goal_amplitude = fp.peak.amplitude.val();
  double goal_ls_amp = fp.peak.short_tail.amplitude.val();
  double goal_ls_slope = fp.peak.short_tail.slope.val();
  double goal_rs_amp = fp.peak.right_tail.amplitude.val();
  double goal_rs_slope = fp.peak.right_tail.slope.val();

  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());
  std::uniform_real_distribution<double> pos_dist(fp.peak.position.min(),
                                                  fp.peak.position.max());
  std::uniform_real_distribution<double> width_dist(fp.peak.width.min(),
                                                    fp.peak.width.max());
  std::uniform_real_distribution<double> amplitude_dist(3000, 5000);
  std::uniform_real_distribution<double> ls_amplitude_dist(fp.peak.short_tail.amplitude.min(),
                                                           fp.peak.short_tail.amplitude.max());
  std::uniform_real_distribution<double> ls_slope_dist(fp.peak.short_tail.slope.min(),
                                                       fp.peak.short_tail.slope.max());
  std::uniform_real_distribution<double> rs_amplitude_dist(fp.peak.right_tail.amplitude.min(),
                                                           fp.peak.right_tail.amplitude.max());
  std::uniform_real_distribution<double> rs_slope_dist(fp.peak.right_tail.slope.min(),
                                                       fp.peak.right_tail.slope.max());

  DAQuiri::DLibOptimizer optimizer;

  for (size_t i = 0; i < 10; ++i)
  {
    fp.peak.position.val(pos_dist(random_generator));
    fp.peak.width.val(width_dist(random_generator));
    fp.peak.amplitude.val(amplitude_dist(random_generator));
    fp.peak.short_tail.amplitude.val(ls_amplitude_dist(random_generator));
    fp.peak.short_tail.slope.val(ls_slope_dist(random_generator));
    fp.peak.right_tail.amplitude.val(rs_amplitude_dist(random_generator));
    fp.peak.right_tail.slope.val(rs_slope_dist(random_generator));

    MESSAGE() << "Attempt[" << i << "]\n" << fp.peak.to_string();
    fp.save_fit(optimizer.minimize(&fp));
    MESSAGE() << "Result:               \n"
              << fp.peak.to_string("                      ");
    MESSAGE() << "       position delta="
              << (goal_pos - fp.peak.position.val()) << "\n";
    MESSAGE() << "          width delta="
              << (goal_width - fp.peak.width.val()) << "\n";
    MESSAGE() << "      amplitude delta="
              << (goal_amplitude - fp.peak.amplitude.val()) << "\n";
    MESSAGE() << "   ls_amplitude delta="
              << (goal_ls_amp - fp.peak.short_tail.amplitude.val()) << "\n";
    MESSAGE() << "       ls_slope delta="
              << (goal_ls_slope - fp.peak.short_tail.slope.val()) << "\n";
    MESSAGE() << "   rs_amplitude delta="
              << (goal_rs_amp - fp.peak.right_tail.amplitude.val()) << "\n";
    MESSAGE() << "       rs_slope delta="
              << (goal_rs_slope - fp.peak.right_tail.slope.val()) << "\n";

    EXPECT_NEAR(fp.peak.position.val(), goal_pos, 7);
    EXPECT_NEAR(fp.peak.width.val(), goal_width, 2.4);
    EXPECT_NEAR(fp.peak.amplitude.val(), goal_amplitude, 1500);
    EXPECT_NEAR(fp.peak.short_tail.amplitude.val(), goal_ls_amp, 0.75);
    EXPECT_NEAR(fp.peak.short_tail.slope.val(), goal_ls_slope, 0.15);
    EXPECT_NEAR(fp.peak.right_tail.amplitude.val(), goal_rs_amp, 0.45);
    EXPECT_NEAR(fp.peak.right_tail.slope.val(), goal_rs_slope, 0.57);
  }

}