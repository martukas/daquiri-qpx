#include "../function_test.h"
#include "../simple_functions.h"

#include <core/fitting/parameter/positive_param.h>


class PositiveParam : public FunctionTest
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

//////////////
/// Positive
//////////////

TEST_F(PositiveParam, PositiveConst)
{
  ConstFunction<DAQuiri::PositiveParam> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min(), 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max(), 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0.001, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 2e-3);
  EXPECT_NEAR(check_gradients(false), 10.0, 2e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
//  check_gradients(true);
  check_gradient_deltas(false);

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(10);
    test_fit(1, &fl, &fl.val, 30, 1e-9);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-8});
    EXPECT_EQ(unconverged, 0u);
    EXPECT_EQ(not_sane, 0u);
    EXPECT_EQ(converged_finite, 0u);
    EXPECT_LE(max_iterations_to_converge,13u);
    EXPECT_LE(max_perturbations_to_converge, 0u);
  }
}

TEST_F(PositiveParam, PositiveLinear)
{
  LinearFunction<DAQuiri::PositiveParam> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min(), 0, 1e-20);
  EXPECT_NEAR(fl.data.count_max(), 39 * 5, 1e-13);

  survey_grad(&fl, &fl.val, 0.001, 0.001, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
//  check_gradients(true);
  check_gradient_deltas(false);

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-14);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-10});
    EXPECT_EQ(unconverged, 0u);
    EXPECT_EQ(not_sane, 0u);
    EXPECT_EQ(converged_finite, 0u);
    EXPECT_LE(max_iterations_to_converge, 14u);
    EXPECT_LE(max_perturbations_to_converge, 0u);
  }
}

TEST_F(PositiveParam, PositiveQuadratic)
{
  QuadraticFunction<DAQuiri::PositiveParam> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min(), 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max(), 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0.001, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
//  check_gradients(true);
  check_gradient_deltas(false);

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-13);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-11});
    EXPECT_EQ(unconverged, 0u);
    EXPECT_EQ(not_sane, 0u);
    EXPECT_EQ(converged_finite, 0u);
    EXPECT_LE(max_iterations_to_converge, 13u);
    EXPECT_LE(max_perturbations_to_converge, 0u);
  }
}

