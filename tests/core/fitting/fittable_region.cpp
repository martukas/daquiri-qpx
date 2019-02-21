#include "function_test.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/fittable_region.h>
#include <core/fitting/hypermet/Value.h>

#include <core/fitting/optimizers/optlib_adapter.h>

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

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
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

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
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

    if (!result.inv_hessian.innerSize() || !result.inv_hessian.outerSize())
      return;

    Eigen::VectorXd diags = result.inv_hessian.diagonal();
    diags *= degrees_of_freedom();
    val.get_uncert(diags, chi_sq());
  }
};

class FittableRegion : public FunctionTest
{
 protected:
  DAQuiri::OptlibOptimizer optimizer;
  size_t random_samples{100};

  virtual void SetUp()
  {
//    optimizer.verbose = true;
  }
};

TEST_F(FittableRegion, UnboundedConst)
{
  ConstFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-20);
  EXPECT_NEAR(check_gradients(false), 10.0, 1e-20);
  
  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(10);

  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-20);
  fl.val.val(10);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-20});
}

TEST_F(FittableRegion, PositiveConst)
{
  ConstFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 2e-3);
//  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(10);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-9);
  fl.val.val(10);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-8});
}

TEST_F(FittableRegion, BoundedConst)
{
  ConstFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(10);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-3);
//  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(10);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-13);
  fl.val.val(10);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-9});
}

TEST_F(FittableRegion, UnboundedLinear)
{
  LinearFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-20);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-20);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-14);
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-11});
}

TEST_F(FittableRegion, PositiveLinear)
{
  LinearFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-20);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-13);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  //EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-12);
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-10});
}

TEST_F(FittableRegion, BoundedLinear)
{
  LinearFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-14);
  // \todo why is this so bad?!?!?!
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 22});
}

TEST_F(FittableRegion, UnboundedQuadratic)
{
  QuadraticFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-12);
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-11});
}

TEST_F(FittableRegion, PositiveQuadratic)
{
  QuadraticFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
//  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  test_fit(5, &optimizer, &fl, &fl.val, 30, 1e-12);
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 1e-11});
}

TEST_F(FittableRegion, BoundedQuadratic)
{
  QuadraticFunction<DAQuiri::Value> fl;
  fl.val.bound(0,40);
  fl.val.val(5);
  fl.data = generate_data(&fl, 40);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);

  EXPECT_TRUE(optimizer.check_gradient(&fl));

  fl.val.val(5);
  // \todo This is terrible...
  test_fit(5, &optimizer, &fl, &fl.val, 30, 25);
  fl.val.val(5);
  test_fit_random(random_samples, &optimizer, &fl, {"val", &fl.val, 0, 40, 34});
}

