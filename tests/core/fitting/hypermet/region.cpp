#include "gtest_color_print.h"
#include "visualize_vector.h"
#include <random>

#include <core/fitting/hypermet/Region.h>
#include <core/fitting/region_manager.h>

class Region : public TestBase
{
  virtual void SetUp()
  {
    rng.seed(std::random_device()());
  }
 protected:
  std::mt19937 rng;
};

TEST_F(Region, IndexBackgroundOnly)
{
  DAQuiri::Region region;
  EXPECT_EQ(region.variable_count(), 0);
  region.map_fit();
  EXPECT_EQ(region.variable_count(), 3);
}

TEST_F(Region, EvalBackground)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(70);
  region.background.slope.val(3);
  region.background.curve.val(5);

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();

  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";
  auto y2 = region.background.eval(x);

  MESSAGE() << "Region:\n" << region.to_string(" ") << "\n";

  MESSAGE() << "Estmate val:\n" << visualize(x, y2, 80) << "\n";

  size_t fits = 100;
  double delta_slope = (3.0 - region.background.slope.val()) / static_cast<double>(fits);
  double delta_curve = (5.0 - region.background.curve.val()) / static_cast<double>(fits);

  std::vector<double> chi, base, slope, curve;
  Eigen::VectorXd grad(region.variable_count());
  for (size_t i=0; i < fits; ++i)
  {
    region.background.slope.val(region.background.slope.val() + delta_slope);
    region.background.curve.val(region.background.curve.val() + delta_curve);
    Eigen::VectorXd vars = region.variables();
    chi.push_back(region(vars, grad));
    base.push_back(grad[0]);
    slope.push_back(grad[1]);
    curve.push_back(grad[2]);
  }

  MESSAGE() << "Chi:\n" << visualize(chi, 80) << "\n";
  MESSAGE() << "Grad(base):\n" << visualize(base, 80) << "\n";
  MESSAGE() << "Grad(slope):\n" << visualize(slope, 80) << "\n";
  MESSAGE() << "Grad(curve):\n" << visualize(curve, 80) << "\n";
}

TEST_F(Region, FitBackground)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(70);
  region.background.slope.val(3);
  region.background.curve.val(5);

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);
  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  MESSAGE() << "Estimate val:\n" << visualize(x, y2, 80) << "\n";


  DAQuiri::fitter_vector starting_point = dlib::mat(region.variables());

  MESSAGE() << "Starting point:\n" <<
            fmt::format(" size={}, 0={}, 1={}, 2={}\n", starting_point.size(),
                starting_point(0), starting_point(1), starting_point(2));

  const auto& fe = std::bind(&DAQuiri::Region::eval, region, std::placeholders::_1);
  const auto& fd = std::bind(&DAQuiri::Region::derivative, region, std::placeholders::_1);
  dlib::find_min(dlib::bfgs_search_strategy(),
                 dlib::objective_delta_stop_strategy(1e-10).be_verbose(),
                 fe, fd, starting_point, -1);

  MESSAGE() << "End point:\n" <<
            fmt::format(" size={}, 0={}, 1={}, 2={}\n", starting_point.size(),
                        starting_point(0), starting_point(1), starting_point(2));

  Eigen::VectorXd v;
  v.setConstant(starting_point.size(), 0.0);
  for (long i = 0; i < starting_point.size(); ++i)
    v[i] = starting_point(i);
  region.save_fit(v);

  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
}


TEST_F(Region, FitBackgroundWithNoise)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(70);
  region.background.slope.val(3);
  region.background.curve.val(5);

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);
  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";

  std::normal_distribution<double> normal_dist(0.0, 5.0);
  for (auto& p : y)
    p += normal_dist(rng);

  MESSAGE() << "Noisy val:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  DAQuiri::fitter_vector starting_point = dlib::mat(region.variables());

  const auto& fe = std::bind(&DAQuiri::Region::eval, region, std::placeholders::_1);
  const auto& fd = std::bind(&DAQuiri::Region::derivative, region, std::placeholders::_1);
  dlib::find_min(dlib::bfgs_search_strategy(),
                 dlib::objective_delta_stop_strategy(1e-10).be_verbose(),
                 fe, fd, starting_point, -1);

  Eigen::VectorXd v;
  v.setConstant(starting_point.size(), 0.0);
  for (long i = 0; i < starting_point.size(); ++i)
    v[i] = starting_point(i);
  region.save_fit(v);

  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
}

TEST_F(Region, FitBackgroundOurImplementation)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(70);
  region.background.slope.val(3);
  region.background.curve.val(5);

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);
  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  DAQuiri::OptimizerType optimizer;
  auto result = optimizer.BFGSMin(&region, 1e-10);
  region.save_fit_uncerts(result);

  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
}


TEST_F(Region, FitOnePeakGaussianOnly)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(70);
  region.background.slope.val(1);
  region.background.curve.val(0);

  DAQuiri::Peak pk;
  pk.position.bound(3, 27);
  pk.position.val(15);
  pk.amplitude.val(500);
  pk.width_.val(3);
  pk.step.enabled = false;
  pk.short_tail.enabled = false;
  pk.right_tail.enabled = false;
  pk.long_tail.enabled = false;

  region.default_peak_ = pk;

//  pk.step.override = true;
//  pk.short_tail.override = true;
//  pk.right_tail.override = true;
//  pk.long_tail.override = true;
  region.peaks_[pk.id()] = pk;

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);
  for (size_t i=0; i < 30; ++i)
    y[i] += pk.eval(x[i]).all();

  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";

  region.remove_peak(pk.id());

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.add_peak(3, 27, 300);
  region.map_fit();
  auto y2 = region.background.eval(x);
  auto pk2 = region.peaks_.begin()->second;
  for (size_t i=0; i < 30; ++i)
    y2[i] += pk2.eval(x[i]).all();

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  MESSAGE() << "Estimate val:\n" << visualize(x, y2, 80) << "\n";

  DAQuiri::fitter_vector starting_point = dlib::mat(region.variables());

  MESSAGE() << "Starting point:\n" <<
            fmt::format(" size={}, 0={}, 1={}, 2={}\n", starting_point.size(),
                        starting_point(0), starting_point(1), starting_point(2));

  const auto& fe = std::bind(&DAQuiri::Region::eval, region, std::placeholders::_1);
  const auto& fd = std::bind(&DAQuiri::Region::derivative, region, std::placeholders::_1);
  dlib::find_min(dlib::bfgs_search_strategy(),
                 dlib::objective_delta_stop_strategy(1e-10).be_verbose(),
                 fe, fd, starting_point, -1);

  MESSAGE() << "End point:\n" <<
            fmt::format(" size={}, 0={}, 1={}, 2={}\n", starting_point.size(),
                        starting_point(0), starting_point(1), starting_point(2));

  Eigen::VectorXd v;
  v.setConstant(starting_point.size(), 0.0);
  for (long i = 0; i < starting_point.size(); ++i)
    v[i] = starting_point(i);
  region.save_fit(v);

  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
}
