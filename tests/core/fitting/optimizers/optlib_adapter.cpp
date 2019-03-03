#include "../function_test.h"

#include "rosenbrock.h"

#include <core/fitting/optimizers/optlib_adapter.h>

#include <random>

class OptlibOptimizer : public TestBase
{
 protected:
  DAQuiri::OptlibOptimizer optimizer;
  Rosenbrock rb {10};

  virtual void SetUp()
  {
//    optimizer.verbose = true;
//    optimizer.maximum_iterations = 200;
  }
};

TEST_F(OptlibOptimizer, CheckSetup)
{
  MESSAGE() << "Rosenbrock: " << rb.variables().transpose() << "\n";
}

TEST_F(OptlibOptimizer, CheckGradient)
{
  EXPECT_TRUE(optimizer.check_gradient(&rb));
}

TEST_F(OptlibOptimizer, Fit10)
{
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 25);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10nonzero)
{
  for (size_t i=0; i < rb.vals_.size(); ++i)
    rb.vals_[i] = 3;

  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 67);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10nonzero2)
{
  for (size_t i=0; i < rb.vals_.size(); ++i)
    rb.vals_[i] = i;

  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 210);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10Random)
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(0, 10);
  for (size_t i=0; i < rb.vals_.size(); ++i)
    rb.vals_[i] = dist6(rng);

  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 300);
  EXPECT_FALSE(result.used_finite_grads);
}

