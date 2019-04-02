#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/parameter/abstract_bounded_param.h>

class DummyBounded : public DAQuiri::BoundedValue
{
 public:
  DummyBounded() = default;

  using AbstractValue::x;
  using AbstractValue::val;
  using AbstractValue::grad;
  using BoundedValue::min;
  using BoundedValue::max;
  using BoundedValue::bound;

  void val(double new_val) override
  {
    double t = (min() + max() - 2.0 * new_val) / (min() - max());
    if (std::abs(t) <= 1)
      x(std::asin((min() + max() - 2.0 * new_val) / (min() - max())));
    else if (signum(t) < 0)
      x(std::asin(-1));
    else
      x(std::asin(1));
  }

  double val_at(double at_x) const override
  {
    return (1.0 + std::sin(at_x)) * (max() - min()) / 2.0 + min();
  }

  double grad_at(double at_x) const override
  {
    return std::cos(at_x) * (max() - min()) / 2.0;
  }
};

class AbstractBoundedParam : public TestBase
{
 protected:

  void SetUp() override
  {
  }
};


TEST_F(AbstractBoundedParam, ValMapsOneToOne)
{
  DummyBounded val;

  for (double v = val.min(); v < val.max(); v += 0.0001)
  {
    val.val(v);
    EXPECT_NEAR(v, val.val(), 0.00000000000000012);
  }
}

TEST_F(AbstractBoundedParam, BoundsEnforcedOnAssignment)
{
  DummyBounded v;
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

TEST_F(AbstractBoundedParam, MinEnforced)
{
  DummyBounded v;
  v.bound(10, 20);
  v.val(10);
  EXPECT_EQ(v.val(), 10);

  v.min(11);
  EXPECT_EQ(v.val(), 11);
  v.min(15);
  EXPECT_EQ(v.val(), 15);
  v.min(17);
  EXPECT_EQ(v.val(), 17);
}

TEST_F(AbstractBoundedParam, MaxEnforced)
{
  DummyBounded v;
  v.bound(10, 20);
  v.val(20);
  EXPECT_EQ(v.val(), 20);

  v.max(19);
  EXPECT_EQ(v.val(), 19);
  v.max(18);
  EXPECT_EQ(v.val(), 18);
  v.max(12);
  EXPECT_EQ(v.val(), 12);
}

TEST_F(AbstractBoundedParam, MinMaxEnforced)
{
  DummyBounded v;

  v.bound(10, 20);
  EXPECT_EQ(v.min(), 10);
  EXPECT_EQ(v.max(), 20);

  v.bound(-10, -20);
  EXPECT_EQ(v.min(), -20);
  EXPECT_EQ(v.max(), -10);
}

TEST_F(AbstractBoundedParam, Print)
{
  DummyBounded val;
  val.bound(-87654321e19, 87654321e19);
  val.val(1235678e5);
  MESSAGE() << "|" << val.to_string() << "|\n";
}
