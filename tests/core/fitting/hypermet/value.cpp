#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/hypermet/Value.h>

class Value : public TestBase
{
};

TEST_F(Value, CheckErfc)
{
  std::vector<double> x;
  std::vector<double> y;
  for (double v = -10; v < 10; v+=0.1)
  {
    x.push_back(v);
    y.push_back(std::erfc(v));
  }

  MESSAGE() << "Value::val_at(x):\n" << visualize(x, y, 100) << "\n";
}

TEST_F(Value, DefaultConstruct)
{
  DAQuiri::SineBoundedValue v;
  EXPECT_EQ(v.index(), -1);
  EXPECT_EQ(v.uncert(), 0.0);
  EXPECT_EQ(v.x(), 0.0);
  EXPECT_EQ(v.val(), 0.5);
}

TEST_F(Value, UpdateIndexInvalidThrows)
{
  DAQuiri::SineBoundedValue v;
  EXPECT_EQ(v.index(), -1);

  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(v.update_index(i));

  i = -42;
  EXPECT_ANY_THROW(v.update_index(i));
}

TEST_F(Value, UpdateIndex)
{
  DAQuiri::SineBoundedValue v;
  EXPECT_EQ(v.index(), -1);

  int32_t i;

  i = 0;
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);
  EXPECT_EQ(i, 1);

  v.update_index(i);
  EXPECT_EQ(v.index(), 1);
  EXPECT_EQ(i, 2);

  i=42;
  v.update_index(i);
  EXPECT_EQ(v.index(), 42);
  EXPECT_EQ(i, 43);
}

TEST_F(Value, UpdateIndexInvalidates)
{
  DAQuiri::SineBoundedValue v;
  EXPECT_EQ(v.index(), -1);

  int32_t i;

  i = 0;
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);
  EXPECT_EQ(i, 1);

  v.to_fit = false;
  v.update_index(i);
  EXPECT_EQ(v.index(), -1);
  EXPECT_EQ(i, 1);

  v.update_index(i);
  EXPECT_EQ(v.index(), -1);
  EXPECT_EQ(i, 1);
}

TEST_F(Value, AcceptAllX)
{
  DAQuiri::SineBoundedValue val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    EXPECT_EQ(val.x(), x);
  }
}

TEST_F(Value, ValAtWithinBoundsForAllX)
{
  DAQuiri::SineBoundedValue val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    auto vv = val.val_at(x);
    EXPECT_GE(vv, 0.0);
    EXPECT_LE(vv, 1.0);
  }
}

TEST_F(Value, ValWithinBoundsForAllX)
{
  DAQuiri::SineBoundedValue val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    auto vv = val.val();
    EXPECT_GE(vv, 0.0);
    EXPECT_LE(vv, 1.0);
  }
}

TEST_F(Value, GradAtWithinBoundsForAllX)
{
  DAQuiri::SineBoundedValue val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    auto g = val.grad_at(x);
    EXPECT_GE(g, -0.5);
    EXPECT_LE(g, 0.5);
  }
}

TEST_F(Value, GradWithinBoundsForAllX)
{
  DAQuiri::SineBoundedValue val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    auto g = val.grad();
    EXPECT_GE(g, -0.5);
    EXPECT_LE(g, 0.5);
  }
}

TEST_F(Value, ValMapsOneToOne)
{
  DAQuiri::SineBoundedValue val;

  for (double v = val.min(); v < val.max(); v += 0.0001)
  {
    val.val(v);
    EXPECT_NEAR(v, val.val(), 0.00000000000000012);
  }
}

TEST_F(Value, BoundsEnforcedOnAssignment)
{
  DAQuiri::SineBoundedValue v;
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

TEST_F(Value, MinEnforced)
{
  DAQuiri::SineBoundedValue v;
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

TEST_F(Value, MaxEnforced)
{
  DAQuiri::SineBoundedValue v;
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

TEST_F(Value, MinMaxEnforced)
{
  DAQuiri::SineBoundedValue v;

  v.bound(10, 20);
  EXPECT_EQ(v.min(), 10);
  EXPECT_EQ(v.max(), 20);

  v.bound(-10, -20);
  EXPECT_EQ(v.min(), -20);
  EXPECT_EQ(v.max(), -10);
}

TEST_F(Value, Put)
{
  DAQuiri::SineBoundedValue v;
  v.val(0.4);

  Eigen::VectorXd fit;
  fit.setConstant(1, 0.0);

  v.put(fit);
  EXPECT_EQ(fit[0], 0.0);
  EXPECT_NE(fit[0], v.x());

  int32_t i{0};
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);

  v.put(fit);
  EXPECT_NE(fit[0], 0.0);
  EXPECT_EQ(fit[0], v.x());
}

TEST_F(Value, Get)
{
  DAQuiri::SineBoundedValue v;
  v.val(0.4);

  Eigen::VectorXd fit;
  fit.setConstant(1, 4.0);

  v.get(fit);
  EXPECT_EQ(v.val(), 0.4);

  int32_t i{0};
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);

  v.get(fit);
  EXPECT_NE(v.val(), 0.4);
  EXPECT_EQ(v.val(), v.val_at(4.0));
}

TEST_F(Value, ValFrom)
{
  DAQuiri::SineBoundedValue v;
  v.val(0.4);

  Eigen::VectorXd fit;
  fit.setConstant(1, 4.0);

  EXPECT_EQ(v.val_from(fit), 0.4);

  int32_t i{0};
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);

  EXPECT_NE(v.val_from(fit), 0.4);
  EXPECT_EQ(v.val_from(fit), v.val_at(4.0));
}

TEST_F(Value, GradFrom)
{
  DAQuiri::SineBoundedValue v;
  v.val(0.4);

  Eigen::VectorXd fit;
  fit.setConstant(1, 4.0);

  EXPECT_EQ(v.grad_from(fit), v.grad());

  auto old_grad = v.grad();

  int32_t i{0};
  v.update_index(i);
  EXPECT_EQ(v.index(), 0);

  EXPECT_NE(v.grad_from(fit), old_grad);
  EXPECT_EQ(v.grad_from(fit), v.grad_at(4.0));
}

TEST_F(Value, Visualize)
{
  DAQuiri::SineBoundedValue val;
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

//  MESSAGE() << "Value::val_at(x):\n" << visualize(x, y, 100) << "\n";
  MESSAGE() << "grad(val):\n" << visualize(y, grad, 100) << "\n";
}

TEST_F(Value, Print)
{
  DAQuiri::UnboundedValue sval;
  sval.val(1235678e99);
  MESSAGE() << "|" << sval.to_string() << "|\n";

  DAQuiri::SineBoundedValue val;
  val.bound(-87654321e19, 87654321e19);
  val.val(1235678e5);
  MESSAGE() << "|" << val.to_string() << "|\n";
}

TEST_F(Value, Val2MapsOneToOne)
{
  DAQuiri::Value2 val;
  val.slope_ = 0.001;

  for (double v = val.min(); v < val.max(); v += 0.0001)
  {
    val.val(v);
    EXPECT_NEAR(v, val.val(), 0.00000000000000012);
  }
}

TEST_F(Value, Visualize2)
{
  DAQuiri::Value2 val;
  val.slope_ = 0.001;
  val.bound(10, 20);

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> grad;
  for (double v = -5000; v < 5000; v+=100)
  {
    x.push_back(v);
    y.push_back(val.val_at(v));
    grad.push_back(val.grad_at(v));
  }

  MESSAGE() << "Value::val_at(x):\n" << visualize(x, y, 100) << "\n";
  MESSAGE() << "grad(val):\n" << visualize(y, grad, 100) << "\n";
}