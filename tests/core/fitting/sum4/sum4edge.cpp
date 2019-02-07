#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/sum4/sum4edge.h>

class SUM4Edge : public TestBase
{
};

TEST_F(SUM4Edge, DefaultConstructed)
{
  DAQuiri::SUM4Edge e;
  EXPECT_EQ(e.width(), 0);
  EXPECT_EQ(e.sum(), UncertainDouble(0, 0));
}
