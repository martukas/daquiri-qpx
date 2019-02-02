#include "gtest_color_print.h"

#include <core/fitting/hypermet/Value.h>

class Value : public TestBase
{
};

TEST_F(Value, GamValueTrue)
{
  DAQuiri::ValueGam vg;
  vg.val(4.0);
  EXPECT_DOUBLE_EQ(vg.val(), 4.0);
  EXPECT_DOUBLE_EQ(vg.x(), 2.0);
  EXPECT_DOUBLE_EQ(vg.val_at(9.0), 81.0);
}

TEST_F(Value, GamGradient)
{
  DAQuiri::ValueGam vg;
  vg.val(9.0);
  EXPECT_DOUBLE_EQ(vg.grad(), 6.0);
  EXPECT_DOUBLE_EQ(vg.grad_at(5.0), 10.0);
}
