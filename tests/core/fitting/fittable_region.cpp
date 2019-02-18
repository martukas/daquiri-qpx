#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/fittable_region.h>
#include <core/fitting/hypermet/Value.h>

#include <core/fitting/optimizers/dlib_adapter.h>

template <typename T>
class ConstFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    return val.val();
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit);
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.index() >= 0)
      grads[val.index()] = val.grad_from(fit);
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);
    // \todo uncerts
  }
};

template <typename T>
class LinearFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    return val.val() * chan;
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit) * chan;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.index() >= 0)
      grads[val.index()] = val.grad_from(fit) * chan;
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);
    // \todo uncerts
  }
};

template <typename T>
class QuadraticFunction : public DAQuiri::FittableRegion
{
 public:
  T val;

  void update_indices() override
  {
    variable_count = 0;
    val.update_index(variable_count);
  }

  Eigen::VectorXd variables() const override
  {
    Eigen::VectorXd ret;
    ret.setConstant(variable_count, 0.0);
    val.put(ret);
    return ret;
  }

  double eval(double chan) const override
  {
    return val.val() * chan * chan;
  }

  double eval_at(double chan, const Eigen::VectorXd& fit) const override
  {
    return val.val_from(fit) * chan * chan;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    double ret = eval_at(chan, fit);
    if (val.index() >= 0)
      grads[val.index()] = val.grad_from(fit) * chan * chan;
    return ret;
  }

  void save_fit(const DAQuiri::FitResult& result) override
  {
    val.get(result.variables);
    // \todo uncerts
  }
};

class FittableRegion : public FunctionTest
{
 protected:
  virtual void SetUp()
  {

  }
};

TEST_F(FittableRegion, UnboundedConstFunctionSurvey)
{
  ConstFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  EXPECT_EQ(fl.variable_count, 0);
  fl.update_indices();
  EXPECT_EQ(fl.variable_count, 1);

  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-5);
  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
}

TEST_F(FittableRegion, UnboundedConstFunctionFit)
{
  ConstFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-5);
  test_fit_random(20, &optimizer, &fl, &fl.val, 0, 40, 1e-4);
}

TEST_F(FittableRegion, PositiveConstFunctionSurvey)
{
  ConstFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  EXPECT_EQ(fl.variable_count, 0);
  fl.update_indices();
  EXPECT_EQ(fl.variable_count, 1);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 2e-3);
  //EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
}

TEST_F(FittableRegion, PositiveConstFunctionFit)
{
  ConstFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-5);
  test_fit_random(20, &optimizer, &fl, &fl.val, 0, 40, 1e-4);
}

TEST_F(FittableRegion, BoundedConstFunctionSurvey)
{
  ConstFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  EXPECT_EQ(fl.variable_count, 0);
  fl.update_indices();
  EXPECT_EQ(fl.variable_count, 1);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-3);
  //EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
}

TEST_F(FittableRegion, BoundedConstFunctionFit)
{
  ConstFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-5);
  test_fit_random(20, &optimizer, &fl, &fl.val, 0, 40, 1e-3);
}

TEST_F(FittableRegion, BoundedLinearFunctionSurvey)
{
  LinearFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-10);

  EXPECT_EQ(fl.variable_count, 0);
  fl.update_indices();
  EXPECT_EQ(fl.variable_count, 1);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
}

TEST_F(FittableRegion, BoundedLinearFunctionFit)
{
  LinearFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-5);
  test_fit_random(20, &optimizer, &fl, &fl.val, 0, 40, 1e-4);
}

TEST_F(FittableRegion, BoundedQuadraticFunctionSurvey)
{
  QuadraticFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  EXPECT_EQ(fl.variable_count, 0);
  fl.update_indices();
  EXPECT_EQ(fl.variable_count, 1);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
}

TEST_F(FittableRegion, BoundedQuadraticFunctionFit)
{
  QuadraticFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  DAQuiri::DLibOptimizer optimizer;
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-5);
  test_fit_random(20, &optimizer, &fl, &fl.val, 0, 40, 1e-4);
}
