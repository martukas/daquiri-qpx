#include "../function_test.h"
#include <core/util/visualize_vector.h>
#include <random>

#include <core/fitting/hypermet/Region.h>

#include <core/fitting/optimizers/optlib_adapter.h>

class Region : public FunctionTest
{
 protected:
  DAQuiri::Region region;
  DAQuiri::OptlibOptimizer optimizer;
  size_t region_size{100};
  size_t random_samples{200};


  void SetUp() override
  {
    optimizer.tolerance = 1e-7;
    optimizer.maximum_perturbations = 10;
    optimizer.maximum_iterations = 200;
    optimizer.gradient_selection =
        DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;

    region.default_peak_ = region.default_peak_.gaussian_only();
    //fp.peak.amplitude.bound(0, 500);
    //region.default_peak_.amplitude.val(40000);
    //region.default_peak_.position.bound(44, 68);
    //region.default_peak_.position.val(51);
    region.default_peak_.width.val(3.2);
  }
};

TEST_F(Region, IndexBackgroundOnly)
{
  DAQuiri::Region region;
  EXPECT_EQ(region.variable_count, 0);
  region.update_indices();
  EXPECT_EQ(region.variable_count, 3);
}
