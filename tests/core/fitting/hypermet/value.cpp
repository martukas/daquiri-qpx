#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Value.h>

class ValueGam : public TestBase
{
};

TEST_F(ValueGam, GamValueTrue)
{
  DAQuiri::ValueGam vg;
  vg.val(4.0);
  EXPECT_DOUBLE_EQ(vg.val(), 4.0);
  EXPECT_DOUBLE_EQ(vg.x(), 2.0);
  EXPECT_DOUBLE_EQ(vg.val_at(9.0), 81.0);
}

TEST_F(ValueGam, GamGradient)
{
  DAQuiri::ValueGam vg;
  vg.val(9.0);
  EXPECT_DOUBLE_EQ(vg.grad(), 6.0);
  EXPECT_DOUBLE_EQ(vg.grad_at(5.0), 10.0);
}


class Value : public TestBase
{
};


TEST_F(Value, ValBounds)
{
  DAQuiri::Value v;
  v.bound(10, 20);

  v.val(10);
  EXPECT_EQ(v.val(), 10);
  v.val(20);
  EXPECT_EQ(v.val(), 20);
  v.val(5);
  EXPECT_EQ(v.val(), 10);
  v.val(30);
  EXPECT_EQ(v.val(), 20);
}

TEST_F(Value, ValAt)
{
  DAQuiri::Value val;
  val.bound(10, 20);

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> grad;
  for (double v = -7.0; v < 7.0; v+=0.1)
  {
    x.push_back(v);
    y.push_back(val.val_at(v));
    grad.push_back(val.grad_at(v));
  }

  MESSAGE() << "Value::val_at(x):\n" << visualize(x, y, 100) << "\n";
  MESSAGE() << "grad(val):\n" << visualize(y, grad, 100) << "\n";
}

