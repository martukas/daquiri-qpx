#include "../function_test.h"
#include "../simple_functions.h"

class DataModel : public FunctionTest
{
 protected:
  size_t region_size{40};
  size_t random_samples{1000};
  bool perform_fitting_tests {false};

  void SetUp() override
  {
    perform_fitting_tests = true;
    //  mpfr::mpreal::set_default_prec(mpfr::digits2bits(100));

    optimizer.maximum_iterations = 200;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;

    optimizer.use_epsilon_check = false;
    optimizer.min_g_norm = 1e-7;

//    optimizer.verbosity = 5;
//    optimizer.tolerance = 1e-14;
  }
};

// \todo tests...