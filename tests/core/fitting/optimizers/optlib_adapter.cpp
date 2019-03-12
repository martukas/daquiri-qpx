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
    optimizer.verbosity = 4;
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

TEST_F(OptlibOptimizer, Fit10zerosEpsilon)
{
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 28u);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10zerosMinDX)
{
  optimizer.use_epsilon_check = false;
  optimizer.min_x_delta = 100 * std::numeric_limits<double>::min();
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 28u);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10zerosMinDF)
{
  optimizer.use_epsilon_check = false;
  optimizer.min_f_delta = 100 * std::numeric_limits<double>::min();
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 28u);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10zerosMinGNorm)
{
  optimizer.use_epsilon_check = false;
  optimizer.min_g_norm = 100 * std::numeric_limits<double>::min();
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 27u);
  EXPECT_FALSE(result.used_finite_grads);
}

TEST_F(OptlibOptimizer, Fit10Random)
{
  std::mt19937 rng;
  rng.seed(std::random_device()());
  std::uniform_real_distribution<double> dist6 (-10.0, 10.0);
  for (size_t i=0; i < (size_t)rb.vals_.size(); ++i)
    rb.vals_[i] = dist6(rng);

  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";
  auto result = optimizer.minimize(&rb);
  MESSAGE() << "Result: " << result.to_string(true) << "\n";
  EXPECT_TRUE(result.converged);
  EXPECT_LE(result.iterations, 300u);
  EXPECT_FALSE(result.used_finite_grads);
}

