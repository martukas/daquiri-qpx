#include "../function_test.h"
#include <core/util/visualize_vector.h>
#include <random>

#include <core/fitting/hypermet/Region.h>
#include <core/fitting/region_manager.h>

#include <core/fitting/optimizers/dlib_adapter.h>

class Region : public FunctionTest
{
 protected:
  DAQuiri::Region region;
  std::mt19937 rng;

  virtual void SetUp()
  {
    rng.seed(std::random_device()());
  }
};

TEST_F(Region, IndexBackgroundOnly)
{
  DAQuiri::Region region;
  EXPECT_EQ(region.variable_count, 0);
  region.update_indices();
  EXPECT_EQ(region.variable_count, 3);
}

TEST_F(Region, FitBackgroundOnly)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(7000);
  region.background.slope.val(1.2);
  region.background.curve.val(-0.02);

  MESSAGE() << "Start region:\n" << region.to_string(" ") << "\n";
  auto data = generate_data(&region, 100);
  visualize_data(data);

  double goal_base = region.background.base.val();
  double goal_slope = region.background.slope.val();
  double goal_curve = region.background.curve.val();

  region.data = data;
  region.update_indices();

  survey_grad(&region, &region.background.base, 0.01, std::sqrt(6900), std::sqrt(7100));
  EXPECT_NEAR(check_chi_sq(false), goal_base, 0.5);
  EXPECT_NEAR(check_gradients(false), goal_base, 0.1);

  survey_grad(&region, &region.background.slope, 0.01, 0, 5);
  EXPECT_NEAR(check_chi_sq(false), goal_slope, 0.5);
  EXPECT_NEAR(check_gradients(false), goal_slope, 0.1);

  survey_grad(&region, &region.background.curve, 0.01, -5, 5);
  EXPECT_NEAR(check_chi_sq(false), goal_curve, 0.5);
  EXPECT_NEAR(check_gradients(false), goal_curve, 0.1);

  region.replace_data(data, 3, 3);
  region.update_indices();
  MESSAGE() << "Auto region:\n" << region.to_string(" ") << "\n";

  std::uniform_real_distribution<double> base_dist(6900, 7100);
  std::uniform_real_distribution<double> slope_dist(0, 5);
  std::uniform_real_distribution<double> curve_dist(-5, 5);

  DAQuiri::DLibOptimizer optimizer;

  for (size_t i = 0; i < 10; ++i)
  {
    region.background.base.val(base_dist(rng));
    region.background.slope.val(slope_dist(rng));
    region.background.curve.val(curve_dist(rng));
    MESSAGE() << "Attempt[" << i << "]\n" << region.background.to_string();

    region.save_fit(optimizer.minimize(&region));
    MESSAGE() << "Result:               \n"
              << region.background.to_string("                      ");
    MESSAGE() << "           base delta="
              << (goal_base - region.background.base.val()) << "\n";
    MESSAGE() << "          slope delta="
              << (goal_slope - region.background.slope.val()) << "\n";
    MESSAGE() << "          curve delta="
              << (goal_curve - region.background.curve.val()) << "\n";

    EXPECT_NEAR(region.background.base.val(), goal_base, 1.0);
    EXPECT_NEAR(region.background.slope.val(), goal_slope, 0.02);
    EXPECT_NEAR(region.background.curve.val(), goal_curve, 0.001);
  }
}

TEST_F(Region, FitGaussianOnly)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.slope_enabled = false;
  region.background.curve_enabled = false;

  auto pk = DAQuiri::Peak().gaussian_only();
  pk.position.bound(44, 58);
  pk.position.val(51);
  pk.amplitude.val(4000);
  pk.width.val(3.2);

  double goal_pos = pk.position.val();
  double goal_width = pk.width.val();
  double goal_amplitude = pk.amplitude.val();

  region.default_peak_ = pk;

  region.peaks_[pk.id()] = pk;

  MESSAGE() << "Start region:\n" << region.to_string(" ") << "\n";
  auto data = generate_data(&region, 100);
  visualize_data(data);

  region = DAQuiri::Region();
  region.default_peak_ = DAQuiri::Peak().gaussian_only();

  region.replace_data(data, 3, 3);
  region.add_peak(44.5, 58.7);
  region.update_indices();

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";


  std::uniform_real_distribution<double> pos_dist(pk.position.min(),
                                                  pk.position.max());
  std::uniform_real_distribution<double> width_dist(pk.width.min(),
                                                    pk.width.max());
  std::uniform_real_distribution<double> amplitude_dist(3000, 5000);

  DAQuiri::DLibOptimizer optimizer;

  for (size_t i = 0; i < 10; ++i)
  {
    auto& pk = region.peaks_.begin()->second;

    pk.position.val(pos_dist(rng));
    pk.width.val(width_dist(rng));
    pk.amplitude.val(amplitude_dist(rng));
    MESSAGE() << "Attempt[" << i << "]\n" << pk.to_string();
    region.save_fit(optimizer.minimize(&region));
    MESSAGE() << "Result:               \n"
              << region.to_string("                      ");
    MESSAGE() << "       position delta="
              << (goal_pos - pk.position.val()) << "\n";
    MESSAGE() << "          width delta="
              << (goal_width - pk.width.val()) << "\n";
    MESSAGE() << "      amplitude delta="
              << (goal_amplitude - pk.amplitude.val()) << "\n";

    EXPECT_NEAR(pk.position.val(), goal_pos, 0.001);
    EXPECT_NEAR(pk.width.val(), goal_width, 0.001);
    EXPECT_NEAR(pk.amplitude.val(), goal_amplitude, 0.2);
  }
}

TEST_F(Region, FitGaussianBackground)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(7000);
  region.background.slope.val(1.2);
  region.background.curve.val(-0.02);

  auto pk = DAQuiri::Peak().gaussian_only();
  pk.position.bound(44, 58);
  pk.position.val(51);
  pk.amplitude.val(4000);
  pk.width.val(3.2);

  double goal_base = region.background.base.val();
  double goal_slope = region.background.slope.val();
  double goal_curve = region.background.curve.val();
  double goal_pos = pk.position.val();
  double goal_width = pk.width.val();
  double goal_amplitude = pk.amplitude.val();

  region.default_peak_ = pk;

  region.peaks_[pk.id()] = pk;

  MESSAGE() << "Start region:\n" << region.to_string(" ") << "\n";
  auto data = generate_data(&region, 100);
  visualize_data(data);

  region = DAQuiri::Region();
  region.default_peak_ = DAQuiri::Peak().gaussian_only();
  region.replace_data(data, 3, 3);
  region.add_peak(44.5, 58.7);
  region.update_indices();

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  std::uniform_real_distribution<double> base_dist(6900, 7100);
  std::uniform_real_distribution<double> slope_dist(0, 5);
  std::uniform_real_distribution<double> curve_dist(-5, 5);
  std::uniform_real_distribution<double> pos_dist(pk.position.min(),
                                                  pk.position.max());
  std::uniform_real_distribution<double> width_dist(pk.width.min(),
                                                    pk.width.max());
  std::uniform_real_distribution<double> amplitude_dist(3000, 5000);

  DAQuiri::DLibOptimizer optimizer;
  optimizer.verbose = true;

  for (size_t i = 0; i < 10; ++i)
  {
    auto& pk = region.peaks_.begin()->second;

    region.background.base.val(base_dist(rng));
    region.background.slope.val(slope_dist(rng));
    region.background.curve.val(curve_dist(rng));
    pk.position.val(pos_dist(rng));
    pk.width.val(width_dist(rng));
    pk.amplitude.val(amplitude_dist(rng));
    MESSAGE() << "Attempt[" << i << "]\n" << region.to_string();
    region.save_fit(optimizer.minimize(&region));
    MESSAGE() << "Result:               \n"
              << region.to_string("                      ");
    MESSAGE() << "           base delta="
              << (goal_base - region.background.base.val()) << "\n";
    MESSAGE() << "          slope delta="
              << (goal_slope - region.background.slope.val()) << "\n";
    MESSAGE() << "          curve delta="
              << (goal_curve - region.background.curve.val()) << "\n";
    MESSAGE() << "       position delta="
              << (goal_pos - pk.position.val()) << "\n";
    MESSAGE() << "          width delta="
              << (goal_width - pk.width.val()) << "\n";
    MESSAGE() << "      amplitude delta="
              << (goal_amplitude - pk.amplitude.val()) << "\n";

    EXPECT_NEAR(region.background.base.val(), goal_base, 1.0);
    EXPECT_NEAR(region.background.slope.val(), goal_slope, 0.02);
    EXPECT_NEAR(region.background.curve.val(), goal_curve, 0.001);
    EXPECT_NEAR(pk.position.val(), goal_pos, 0.001);
    EXPECT_NEAR(pk.width.val(), goal_width, 0.001);
    EXPECT_NEAR(pk.amplitude.val(), goal_amplitude, 0.4);
  }
}

