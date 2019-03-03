#include "function_test.h"
#include "simple_functions.h"

class FittableRegion : public FunctionTest
{
 protected:
  size_t region_size{40};
  size_t random_samples{10000};
  bool perform_fitting_tests {false};

  virtual void SetUp()
  {
    //  mpfr::mpreal::set_default_prec(mpfr::digits2bits(100));

    optimizer.maximum_iterations = 20;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;
//    optimizer.verbose = true;
//    optimizer.tolerance = 1e-14;
  }
};

//////////////
/// Unbounded
//////////////

TEST_F(FittableRegion, UnboundedConst)
{
  ConstFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-20);
  EXPECT_NEAR(check_gradients(false), 10.0, 1e-20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(10);
    test_fit(1, &fl, &fl.val, 30, 1e-90);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-90});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 2);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, UnboundedLinear)
{
  LinearFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-20);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-20);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-14);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-14});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 2);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, UnboundedQuadratic)
{
  QuadraticFunction<DAQuiri::ValueSimple> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-12);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-11});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 2);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

//////////////
/// Positive
//////////////

TEST_F(FittableRegion, PositiveConst)
{
  ConstFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(10);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 2e-3);
//  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(10);
    test_fit(1, &fl, &fl.val, 30, 1e-14);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-8});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 5);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, PositiveLinear)
{
  LinearFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-20);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-13);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  //EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-14);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-10});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 6);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, PositiveQuadratic)
{
  QuadraticFunction<DAQuiri::ValuePositive> fl;
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, 0, 20);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
//  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.5, 0, 20);
  check_gradients(true);
  check_gradient_deltas(true);

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-90);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-11});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 6);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}


//////////////
/// Bounded
//////////////

TEST_F(FittableRegion, BoundedConst)
{
  ConstFunction<DAQuiri::Value> fl;
  fl.val.bound(0, region_size);
  fl.val.val(10);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, -M_PI, M_PI);
  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-3);
//  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
  survey_grad(&fl, &fl.val, 0.05, -M_PI, M_PI);
  check_gradients(true);
  check_gradient_deltas(true);

  optimizer.gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::FiniteAlways;
//  optimizer.gradient_selection =
//      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(10);
    test_fit(1, &fl, &fl.val, 30, 1e-9);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-9});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_EQ(converged_finite, 0);
    EXPECT_LE(max_iterations_to_converge, 5);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, BoundedLinear)
{
  LinearFunction<DAQuiri::Value> fl;
  fl.val.bound(0, region_size);
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, -M_PI, M_PI);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.05, -M_PI, M_PI);
  check_gradients(true);
  check_gradient_deltas(true);

  optimizer.gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::FiniteAlways;
//  optimizer.gradient_selection =
//      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-11);
    deterministic_test(10, &fl, &fl.val, 35);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-9});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_LE(converged_finite, 0.03 * random_samples);
    EXPECT_LE(max_iterations_to_converge, 6);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, BoundedQuadratic)
{
  QuadraticFunction<DAQuiri::Value> fl;
  fl.val.bound(0, region_size);
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, -M_PI, M_PI);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.05, -M_PI, M_PI);
  check_gradients(true);
  check_gradient_deltas(true);

  optimizer.gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::FiniteAlways;
//  optimizer.gradient_selection =
//      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-13);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-10});
    EXPECT_LE(unconverged, 0.05 * random_samples);
    EXPECT_EQ(not_sane, 0);
    EXPECT_LE(converged_finite, 0.45 * random_samples);
    EXPECT_EQ(converged_perturbed, 0);
    EXPECT_LE(max_iterations_to_converge, 7);
    EXPECT_LE(max_perturbations_to_converge, 0);
  }
}

TEST_F(FittableRegion, BoundedQuadraticAutoperturb)
{
  class BoundedQuad : public QuadraticFunction<DAQuiri::Value>
  {
    std::uniform_real_distribution<double> x_dist{-M_PI_2, M_PI_2};

    bool perturb(std::mt19937& rng) override
    {
      val.x(x_dist(rng));
      return true;
    }
  };
  BoundedQuad fl;
  fl.val.bound(0, region_size);
  fl.val.val(5);
  fl.data = generate_data(&fl, region_size);
  fl.update_indices();

  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
  EXPECT_NEAR(fl.data.count_max, 39 * 39 * 5, 1e-10);

  survey_grad(&fl, &fl.val, 0.001, -M_PI, M_PI);
  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
  survey_grad(&fl, &fl.val, 0.05, -M_PI, M_PI);
  check_gradients(true);
  check_gradient_deltas(true);

//  mpfr::mpreal::set_default_prec(mpfr::digits2bits(100));
  optimizer.gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::FiniteAlways;
//  optimizer.gradient_selection =
//      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;

  print_outside_tolerance = true;

  if (perform_fitting_tests)
  {
    fl.val.val(5);
    test_fit(1, &fl, &fl.val, 30, 1e-13);
    deterministic_test(10, &fl, &fl.val, 30);
    test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-10});
    EXPECT_EQ(unconverged, 0);
    EXPECT_EQ(not_sane, 0);
    EXPECT_LE(converged_finite, 0.45 * random_samples);
    EXPECT_LE(converged_perturbed, 0.05 * random_samples);
    EXPECT_LE(max_iterations_to_converge, 7);
    EXPECT_LE(max_perturbations_to_converge, 2);
  }
}

//////////////
/// Bounded2
//////////////
//
//TEST_F(FittableRegion, Bounded2Const)
//{
//  ConstFunction<DAQuiri::Value2> fl;
//  fl.val.bound(0, region_size);
//  fl.val.val(10);
//  fl.data = generate_data(&fl, region_size);
//  fl.update_indices();
//
//  EXPECT_NEAR(fl.data.count_min, 10, 1e-10);
//  EXPECT_NEAR(fl.data.count_max, 10, 1e-10);
//
//  survey_grad(&fl, &fl.val, 0.001, 0, 20);
//  EXPECT_NEAR(check_chi_sq(false), 10.0, 1e-3);
////  EXPECT_NEAR(check_gradients(false), 10.0, 1e-5);
//
//  EXPECT_TRUE(optimizer.check_gradient(&fl));
//
//  fl.val.val(10);
//  test_fit(5, &fl, &fl.val, 30, 1e-13);
//  fl.val.val(10);
//  test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-9});
//  EXPECT_EQ(unconverged, 0);
//  EXPECT_EQ(converged_finite, 0);
//  EXPECT_LE(max_iterations_to_converge, 5);
//  EXPECT_LE(max_perturbations_to_converge, 0);
//}
//
//TEST_F(FittableRegion, Bounded2Linear)
//{
//  LinearFunction<DAQuiri::Value2> fl;
//  fl.val.bound(0, region_size);
//  fl.val.val(5);
//  fl.data = generate_data(&fl, region_size);
//  fl.update_indices();
//
//  EXPECT_NEAR(fl.data.count_min, 0, 1e-10);
//  EXPECT_NEAR(fl.data.count_max, 39 * 5, 1e-10);
//
//  survey_grad(&fl, &fl.val, 0.001, 0, 20);
//  EXPECT_NEAR(check_chi_sq(false), 5.0, 1e-3);
//  EXPECT_NEAR(check_gradients(false), 5.0, 1e-3);
//
//  EXPECT_TRUE(optimizer.check_gradient(&fl));
//
////  optimizer.gradient_selection =
////      DAQuiri::OptlibOptimizer::GradientSelection::FiniteAlways;
//  optimizer.gradient_selection =
//      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;
//
//  fl.val.val(5);
//  test_fit(5, &fl, &fl.val, 30, 1e-10);
//  fl.val.val(5);
//  test_fit_random(random_samples, &fl, {"val", &fl.val, 0, 40, 1e-9});
//  EXPECT_EQ(unconverged, 0);
//  EXPECT_LE(converged_finite, 0.60 * random_samples);
//  EXPECT_LE(converged_perturbed, 0);
//  EXPECT_LE(max_iterations_to_converge, 6);
//  EXPECT_LE(max_perturbations_to_converge, 0);
//}
