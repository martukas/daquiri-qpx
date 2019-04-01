#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/gamma/sum4/peak.h>

class SUM4 : public TestBase
{
};


TEST_F(SUM4, DefaultConstructed)
{
  DAQuiri::SUM4 e;
  EXPECT_EQ(e.peak_width(), 0u);
  EXPECT_EQ(e.peak_area(), UncertainDouble(0, 0));
}
