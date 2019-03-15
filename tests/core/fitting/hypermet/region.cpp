#include "../function_test.h"
#include <core/util/visualize_vector.h>
#include <random>

#include <core/fitting/hypermet/Region.h>

#include <core/fitting/optimizers/optlib_adapter.h>

class Region : public FunctionTest
{
 protected:
  DAQuiri::Region region;
  size_t region_size{100};
  size_t random_samples{100};


  void SetUp() override
  {
//    optimizer.verbosity = 5;
    optimizer.maximum_iterations = 1000;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;
    optimizer.epsilon = 1e-10;
    optimizer.tolerance = 1e-4;
    optimizer.perform_sanity_checks = false;
    optimizer.maximum_perturbations = 0;

    region.background.x_offset = 0;
//    fb.background.base.bound(0, 7792);
    region.background.base.val(70);
//    fb.background.slope.bound(0, 198);
    region.background.slope.val(0.7);
//    fb.background.curve.bound(0, 15);
    region.background.curve.val(-0.01);

    //region.default_peak_ = region.default_peak_.gaussian_only();
    //fp.peak.amplitude.bound(0, 500);
    //region.default_peak_.amplitude.val(40000);
    //region.default_peak_.position.bound(44, 68);
    //region.default_peak_.position.val(51);
    //region.default_peak_.width.val(3.2);
  }

  DAQuiri::Peak make_peak(double pos, double amp)
  {
    DAQuiri::Peak p = region.default_peak_;
    //fp.peak.amplitude.bound(0, 500);
    p.amplitude.val(amp);
    p.position.bound(pos-5, pos+5);
    p.position.val(pos);
    return p;
  }
};

TEST_F(Region, IndexBackgroundOnly)
{
  DAQuiri::Region region;
  EXPECT_EQ(region.variable_count, 0);
  region.update_indices();
  EXPECT_EQ(region.variable_count, 3);
}

TEST_F(Region, VisualizeBackgroundOnly)
{
  auto data = generate_data(&region, region_size);
  visualize_data(data);
}

TEST_F(Region, FitBackgroundOnly)
{
  region.data = generate_data(&region, region_size);
  region.update_indices();

  //auto_bound();
  std::vector<ValueToVary> vals;
  vals.push_back({"base", &region.background.base, 50, 7792, 1e-6});
  vals.push_back({"slope", &region.background.slope, -460, 460, 1e-7});
  vals.push_back({"curve", &region.background.curve, -10, 10, 1e-8});
  test_fit_random(random_samples, &region, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 20u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Region, VisualizeOnePeak)
{
  region.default_peak_ = region.default_peak_.gaussian_only();
  auto p = make_peak(51, 40000);
  region.peaks_[p.id()] = p;
  MESSAGE() << region.to_string() << "\n";
  auto data = generate_data(&region, region_size);
  visualize_data(data);
}

TEST_F(Region, FitOnePeak)
{
  region.default_peak_ = region.default_peak_.gaussian_only();
  auto p = make_peak(51, 40000);
  region.peaks_[p.id()] = p;
  auto data = generate_data(&region, region_size);

  region.peaks_.clear();
  region.data = data;
  region.add_peak(41, 61, 40000);
  region.update_indices();

  auto& pp = region.peaks_.begin()->second;

  MESSAGE() << region.to_string() << "\n";

  std::vector<ValueToVary> vals;
  vals.push_back({"base", &region.background.base, 50, 7792, 1e-6});
  vals.push_back({"slope", &region.background.slope, -460, 460, 1e-7});
  vals.push_back({"curve", &region.background.curve, -10, 10, 1e-8});
  vals.push_back({"width", &region.default_peak_.width,
                  region.default_peak_.width.min(),
                  region.default_peak_.width.max(), 1e-8});
  vals.push_back({"position", &pp.position,
                  pp.position.min(), pp.position.max(), 1e-8});
  vals.push_back({"amplitude", &pp.amplitude,
                  30000, 50000, 2e-4});
  test_fit_random(random_samples, &region, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_EQ(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 50u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Region, FitTwoPeaks)
{
  region.default_peak_ = region.default_peak_.gaussian_only();
  auto p = make_peak(35, 40000);
  region.peaks_[p.id()] = p;
  auto p2 = make_peak(51, 35000);
  region.peaks_[p2.id()] = p2;
  auto data = generate_data(&region, region_size);
  visualize_data(data);

  region.peaks_.clear();
  region.data = data;
  region.add_peak(25, 45, 40000);
  region.add_peak(41, 61, 35000);
  region.update_indices();

  auto& pp = region.peaks_.begin()->second;
  auto& pp2 = region.peaks_.rbegin()->second;

  MESSAGE() << region.to_string() << "\n";

  std::vector<ValueToVary> vals;
  vals.push_back({"base", &region.background.base, 50, 7792, 1e-6});
  vals.push_back({"slope", &region.background.slope, -460, 460, 1e-7});
  vals.push_back({"curve", &region.background.curve, -10, 10, 1e-8});
  vals.push_back({"width", &region.default_peak_.width,
                  region.default_peak_.width.min(),
                  region.default_peak_.width.max(), 1e-8});
  vals.push_back({"position1", &pp.position,
                  pp.position.min(), pp.position.max(), 1e-8});
  vals.push_back({"amplitude1", &pp.amplitude,
                  30000, 50000, 2e-4});
  vals.push_back({"position2", &pp2.position,
                  pp2.position.min(), pp2.position.max(), 1e-8});
  vals.push_back({"amplitude2", &pp2.amplitude,
                  25000, 45000, 2e-4});
  test_fit_random(random_samples, &region, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 120u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}

TEST_F(Region, FitTwoPeaksWithEverything)
{
  auto p = make_peak(35, 40000);
  region.peaks_[p.id()] = p;
  auto p2 = make_peak(51, 35000);
  region.peaks_[p2.id()] = p2;
  auto data = generate_data(&region, region_size);
  visualize_data(data);

  region.peaks_.clear();
  region.data = data;
  region.add_peak(25, 45, 40000);
  region.add_peak(41, 61, 35000);
  region.update_indices();

  auto& pp = region.peaks_.begin()->second;
  auto& pp2 = region.peaks_.rbegin()->second;

  MESSAGE() << region.to_string() << "\n";

  optimizer.maximum_iterations = 600;
  optimizer.perform_sanity_checks = true;

  std::vector<ValueToVary> vals;
  vals.push_back({"base", &region.background.base, 50, 7792, 1e-5});
  vals.push_back({"slope", &region.background.slope, -460, 460, 1e-7});
  vals.push_back({"curve", &region.background.curve, -10, 10, 1e-8});
  vals.push_back({"width", &region.default_peak_.width,
                  region.default_peak_.width.min(),
                  region.default_peak_.width.max(), 1e-8});
  vals.push_back({"position1", &pp.position,
                  pp.position.min(), pp.position.max(), 1e-8});
  vals.push_back({"amplitude1", &pp.amplitude,
                  30000, 50000, 2e-4});
  vals.push_back({"position2", &pp2.position,
                  pp2.position.min(), pp2.position.max(), 1e-8});
  vals.push_back({"amplitude2", &pp2.amplitude,
                  25000, 45000, 2e-4});
  test_fit_random(random_samples, &region, vals);

  EXPECT_EQ(unconverged, 0u);
  EXPECT_EQ(not_sane, 0u);
  EXPECT_EQ(converged_finite, 0u);
  EXPECT_LE(converged_perturbed, 0u);
  EXPECT_LE(max_iterations_to_converge, 500u);
  EXPECT_LE(max_perturbations_to_converge, 0u);
}
