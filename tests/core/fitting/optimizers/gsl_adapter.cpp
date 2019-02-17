#include "../function_test.h"

#include "rosenbrock.h"

#include <core/fitting/optimizers/gsl_adapter.h>

#include <random>

class GSLAdapter : public TestBase
{
 protected:

  virtual void SetUp()
  {
  }
};

TEST_F(GSLAdapter, CheckSetup)
{
  Rosenbrock rb(10);
  MESSAGE() << "Rosenbrock: " << rb.variables().transpose() << "\n";
}

TEST_F(GSLAdapter, Fit10)
{
  Rosenbrock rb(10);
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";

  DAQuiri::GSLAdapter optimizer;
  auto result = optimizer.minimize(&rb, 0.00001);

  MESSAGE() << "Final: " << result.variables.transpose() << "\n";
}

TEST_F(GSLAdapter, Fit10nonzero)
{
  Rosenbrock rb(10);
  for (size_t i=0; i < rb.vals_.size(); ++i)
    rb.vals_[i] = 3;
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";

  DAQuiri::GSLAdapter optimizer;
  auto result = optimizer.minimize(&rb, 0.00001);

  MESSAGE() << "Final: " << result.variables.transpose() << "\n";
}

TEST_F(GSLAdapter, Fit10Random)
{
  Rosenbrock rb(10);

  std::mt19937 rng;
  rng.seed(std::random_device()());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(0, 10);
  for (size_t i=0; i < rb.vals_.size(); ++i)
    rb.vals_[i] = dist6(rng);
  MESSAGE() << "Starting: " << rb.variables().transpose() << "\n";

  DAQuiri::GSLAdapter optimizer;
  auto result = optimizer.minimize(&rb, 0.00001);

  MESSAGE() << "Final: " << result.variables.transpose() << "\n";
}

