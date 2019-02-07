#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/weighted_data.h>

class WeightedData : public TestBase
{
};


TEST_F(WeightedData, DefaultConstructed)
{
  DAQuiri::WeightedData wd;
  EXPECT_TRUE(wd.empty());
}

