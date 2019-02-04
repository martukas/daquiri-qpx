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

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";

  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);

  MESSAGE() << "Goal eval:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();

  MESSAGE() << "Auto region:\n" << region.to_string(" ") << "\n";

  auto y2 = region.background.eval(x);
  MESSAGE() << "Estimate eval:\n" << visualize(x, y2, 80) << "\n";

  size_t fits = 20;
  double delta_slope = (3.0 - region.background.slope.val()) / static_cast<double>(fits);
  double delta_curve = (5.0 - region.background.curve.val()) / static_cast<double>(fits);

  std::vector<double> chi, base, slope, curve;
  Eigen::VectorXd grad(region.variable_count());
  for (size_t i=0; i < 2*fits; ++i)
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

  MESSAGE() << "Goal eval:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Auto region:\n" << region.to_string(" ") << "\n";

  MESSAGE() << "Estimate eval:\n" << visualize(x, y2, 80) << "\n";


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


//TEST_F(Region, FitBackgroundWithNoise)
//{
//  DAQuiri::Region region;
//  region.background.x_offset = 0;
//  region.background.base.val(70);
//  region.background.slope.val(3);
//  region.background.curve.val(5);
//
//  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";
//
//  std::vector<double> x;
//  for (size_t i=0; i < 30; ++i)
//    x.push_back(i);
//  auto y = region.background.eval(x);
//  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";
//
//  std::normal_distribution<double> normal_dist(0.0, 5.0);
//  for (auto& p : y)
//    p += normal_dist(rng);
//
//  MESSAGE() << "Noisy val:\n" << visualize(x, y, 80) << "\n";
//
//  DAQuiri::WeightedData wd(x, y);
//  region.replace_data(wd);
//  region.map_fit();
//  auto y2 = region.background.eval(x);
//
//  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";
//
//  DAQuiri::fitter_vector starting_point = dlib::mat(region.variables());
//
//  const auto& fe = std::bind(&DAQuiri::Region::eval, region, std::placeholders::_1);
//  const auto& fd = std::bind(&DAQuiri::Region::derivative, region, std::placeholders::_1);
//  dlib::find_min(dlib::bfgs_search_strategy(),
//                 dlib::objective_delta_stop_strategy(1e-10).be_verbose(),
//                 fe, fd, starting_point, -1);
//
//  Eigen::VectorXd v;
//  v.setConstant(starting_point.size(), 0.0);
//  for (long i = 0; i < starting_point.size(); ++i)
//    v[i] = starting_point(i);
//  region.save_fit(v);
//
//  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";
//
//  auto y3 = region.background.eval(x);
//  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
//}

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
  MESSAGE() << "Goal eval:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Auto region:\n" << region.to_string(" ") << "\n";

  DAQuiri::OptimizerType optimizer;
  auto result = optimizer.BFGSMin(&region, 1e-10);
  region.save_fit_uncerts(result);

  MESSAGE() << "Fit region:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Fit eval:\n" << visualize(x, y3, 80) << "\n";
}

TEST_F(Region, FitNoisyBackgroundOurImplementation)
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
  MESSAGE() << "Clean eval:\n" << visualize(x, y, 80) << "\n";

  std::normal_distribution<double> normal_dist(0.0, 5.0);
  for (auto& p : y)
    p += normal_dist(rng);

  MESSAGE() << "Noisy eval:\n" << visualize(x, y, 80) << "\n";

  DAQuiri::WeightedData wd(x, y);
  region.replace_data(wd, 3, 3);
  region.map_fit();
  auto y2 = region.background.eval(x);

  MESSAGE() << "Auto region:\n" << region.to_string(" ") << "\n";

  DAQuiri::OptimizerType optimizer;
  auto result = optimizer.BFGSMin(&region, 1e-10);
  region.save_fit_uncerts(result);

  MESSAGE() << "Fit region:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Fit eval:\n" << visualize(x, y3, 80) << "\n";
}


TEST_F(Region, EvalOnePeakGaussianOnly)
{
  DAQuiri::Region region;
  region.background.x_offset = 0;
  region.background.base.val(0);
  region.background.slope.val(0);
  region.background.curve.val(0);

  auto pk = DAQuiri::Peak().gaussian_only();
  pk.position.bound(3, 27);
  pk.position.val(15);
  pk.amplitude.val(500);
  pk.width_.val(3);

  region.default_peak_ = pk;

  region.peaks_[pk.id()] = pk;

  MESSAGE() << "Goal region:\n" << region.to_string(" ") << "\n";


  std::vector<double> x;
  for (size_t i=0; i < 30; ++i)
    x.push_back(i);
  auto y = region.background.eval(x);
  for (size_t i=0; i < 30; ++i)
    y[i] += pk.eval(x[i]).all();

  MESSAGE() << "Original val:\n" << visualize(x, y, 80) << "\n";

  region = DAQuiri::Region();
  region.default_peak_ = DAQuiri::Peak().gaussian_only();

  region.replace_data(DAQuiri::WeightedData(x, y));
  region.add_peak(5, 27, 300);
  region.map_fit();

  MESSAGE() << "Region with data:\n" << region.to_string(" ") << "\n";

  auto y2 = region.background.eval(x);
  auto pk2 = region.peaks_.begin()->second;
  for (size_t i=0; i < 30; ++i)
    y2[i] += pk2.eval(x[i]).all();

  MESSAGE() << "Estimate val:\n" << visualize(x, y2, 80) << "\n";

  size_t fits = 20;
  double delta_width = (3.0 - pk2.width_.val()) / static_cast<double>(fits);
  double delta_pos = (5.0 - pk2.position.val()) / static_cast<double>(fits);
  double delta_amp = (500.0 - pk2.amplitude.val()) / static_cast<double>(fits);

  std::vector<double> chi, width, pos, amp;
  Eigen::VectorXd grad(region.variable_count());
  for (size_t i=0; i < 2*fits; ++i)
  {
    auto& pkk = region.peaks_.begin()->second;
    pkk.width_.val(pkk.width_.val() + delta_width);
    pkk.position.val(pkk.position.val() + delta_pos);
    pkk.amplitude.val(pkk.amplitude.val() + delta_amp);
    Eigen::VectorXd vars = region.variables();
    chi.push_back(region(vars, grad));
    width.push_back(grad[3]);
    pos.push_back(grad[4]);
    amp.push_back(grad[5]);
  }

  MESSAGE() << "Chi:\n" << visualize(chi, 80) << "\n";
  MESSAGE() << "Grad(with):\n" << visualize(width, 80) << "\n";
  MESSAGE() << "Grad(pos):\n" << visualize(pos, 80) << "\n";
  MESSAGE() << "Grad(amp):\n" << visualize(amp, 80) << "\n";

  region = DAQuiri::Region();
  region.default_peak_ = DAQuiri::Peak().gaussian_only();
  region.replace_data(DAQuiri::WeightedData(x, y));
  region.add_peak(5, 27, 300);
  region.map_fit();

  DAQuiri::OptimizerType optimizer;
  auto result = optimizer.BFGSMin(&region, 1e-10);
  region.save_fit_uncerts(result);

  MESSAGE() << "Region after fit:\n" << region.to_string(" ") << "\n";

  auto y3 = region.background.eval(x);
  MESSAGE() << "Final fit:\n" << visualize(x, y3, 80) << "\n";
}
