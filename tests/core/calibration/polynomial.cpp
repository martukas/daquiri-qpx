#include "gtest_color_print.h"
#include <core/calibration/polynomial.h>

class Polynomial : public TestBase
{
};

TEST_F(Polynomial, DefaultConstr)
{
  DAQuiri::Polynomial cf;
  EXPECT_FALSE(cf.valid());
}

TEST_F(Polynomial, SetCoef)
{
  DAQuiri::Polynomial cf({2});
  EXPECT_TRUE(cf.valid());
}

TEST_F(Polynomial, InitCoefs)
{
  DAQuiri::Polynomial cf({5.0, 3.0});
  EXPECT_TRUE(cf.valid());
}

TEST_F(Polynomial, Eval)
{
  DAQuiri::Polynomial cf{{5.0, 2.0, 1.0}};
  EXPECT_DOUBLE_EQ(8.0, cf.eval(1.0));
  EXPECT_DOUBLE_EQ(13.0, cf.eval(2.0));
  EXPECT_DOUBLE_EQ(20.0, cf.eval(3.0));
}

TEST_F(Polynomial, Debug)
{
  DAQuiri::Polynomial cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_string() << "\n";
}

TEST_F(Polynomial, UTF8)
{
  DAQuiri::Polynomial cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_UTF8(3) << "\n";
}

TEST_F(Polynomial, Markup)
{
  DAQuiri::Polynomial cf{{5.0, 2.0, 1.0}};
  MESSAGE() << cf.to_markup(3) << "\n";
}
