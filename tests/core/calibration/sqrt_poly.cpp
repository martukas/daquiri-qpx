#include "gtest_color_print.h"
#include <core/calibration/sqrt_poly.h>

class SqrtPoly : public TestBase
{
};

TEST_F(SqrtPoly, DefaultConstr)
{
  DAQuiri::SqrtPoly cf;
  EXPECT_FALSE(cf.valid());
}

TEST_F(SqrtPoly, SetCoef)
{
  DAQuiri::SqrtPoly cf({2});
  EXPECT_TRUE(cf.valid());
}

TEST_F(SqrtPoly, InitCoefs)
{
  DAQuiri::SqrtPoly cf({5.0, 3.0});
  EXPECT_TRUE(cf.valid());
}

TEST_F(SqrtPoly, Eval)
{
  DAQuiri::SqrtPoly cf{{5.0, 2.0, 1.0}};
  EXPECT_DOUBLE_EQ(sqrt(5 + 2 + 1), cf.eval(1.0));
  EXPECT_DOUBLE_EQ(sqrt(5 + 2*2 + 1*4), cf.eval(2.0));
  EXPECT_DOUBLE_EQ(sqrt(5 + 2*3 + 1*9), cf.eval(3.0));
}

TEST_F(SqrtPoly, Debug)
{
  DAQuiri::SqrtPoly cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_string() << "\n";
}

TEST_F(SqrtPoly, UTF8)
{
  DAQuiri::SqrtPoly cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_UTF8(3) << "\n";
}

TEST_F(SqrtPoly, Markup)
{
  DAQuiri::SqrtPoly cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_markup(3) << "\n";
}
