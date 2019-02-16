#include "gtest_color_print.h"

#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/GaussianPrecalc.h>

TEST(Side, Flip)
{
  EXPECT_DOUBLE_EQ(DAQuiri::flip(DAQuiri::Side::left, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(DAQuiri::flip(DAQuiri::Side::left, -1.0), -1.0);
  EXPECT_DOUBLE_EQ(DAQuiri::flip(DAQuiri::Side::right, 1.0), -1.0);
  EXPECT_DOUBLE_EQ(DAQuiri::flip(DAQuiri::Side::right, -1.0), 1.0);
}

TEST(Side, ToFromString)
{
  EXPECT_EQ(DAQuiri::to_string(DAQuiri::Side::left), "left");
  EXPECT_EQ(DAQuiri::to_string(DAQuiri::Side::right), "right");

  EXPECT_EQ(DAQuiri::to_side("left"), DAQuiri::Side::left);
  EXPECT_EQ(DAQuiri::to_side("right"), DAQuiri::Side::right);

  EXPECT_EQ(DAQuiri::to_side(DAQuiri::to_string(DAQuiri::Side::left)),
      DAQuiri::Side::left);
  EXPECT_EQ(DAQuiri::to_side(DAQuiri::to_string(DAQuiri::Side::right)),
      DAQuiri::Side::right);

  EXPECT_EQ(DAQuiri::to_string(DAQuiri::to_side("left")), "left");
  EXPECT_EQ(DAQuiri::to_string(DAQuiri::to_side("right")), "right");
}
