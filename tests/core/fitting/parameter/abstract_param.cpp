#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/parameter/abstract_param.h>


class DummyParam : public DAQuiri::AbstractParam
{
 public:
  DummyParam() = default;

  using DAQuiri::AbstractParam::x;
  using DAQuiri::AbstractParam::val;
  using DAQuiri::AbstractParam::grad;

  void val(double new_val) override
  {
    x(new_val / 2.0);
  }

  double val_at(double at_x) const override
  {
    return 2.0 * at_x;
  }

  double grad_at(double at_x) const override
  {
    (void) at_x;
    return 2.0;
  }
};

class AbstractParam : public TestBase
{
 protected:

  void SetUp() override
  {
  }
};


TEST_F(AbstractParam, DefaultConstruct)
{
  DummyParam v;
  EXPECT_EQ(v.index(), -1);
  EXPECT_EQ(v.uncert(), 0.0);
  EXPECT_EQ(v.x(), 0.0);
  EXPECT_EQ(v.val(), 0.5);
}

TEST_F(AbstractParam, UpdateIndexInvalidThrows)
{
  DummyParam v;
  EXPECT_EQ(v.index(), -1);

  int32_t i;

  i = -1;
  EXPECT_ANY_THROW(v.update_index(i));

  i = -42;
  EXPECT_ANY_THROW(v.update_index(i));
}

TEST_F(AbstractParam, UpdateIndex)
{
  DummyParam v;
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

TEST_F(AbstractParam, UpdateIndexInvalidates)
{
  DummyParam v;
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

TEST_F(AbstractParam, AcceptAllX)
{
  DummyParam val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    EXPECT_EQ(val.x(), x);
  }
}

TEST_F(AbstractParam, ValAtWithinBoundsForAllX)
{
  DummyParam val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    auto vv = val.val_at(x);
    EXPECT_GE(vv, 0.0);
    EXPECT_LE(vv, 1.0);
  }
}

TEST_F(AbstractParam, ValWithinBoundsForAllX)
{
  DummyParam val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    auto vv = val.val();
    EXPECT_GE(vv, 0.0);
    EXPECT_LE(vv, 1.0);
  }
}

TEST_F(AbstractParam, GradAtWithinBoundsForAllX)
{
  DummyParam val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    auto g = val.grad_at(x);
    EXPECT_GE(g, -0.5);
    EXPECT_LE(g, 0.5);
  }
}

TEST_F(AbstractParam, GradWithinBoundsForAllX)
{
  DummyParam val;

  double max = 4242;
  for (double x = -max; x < max; x += 0.01)
  {
    val.x(x);
    auto g = val.grad();
    EXPECT_GE(g, -0.5);
    EXPECT_LE(g, 0.5);
  }
}

TEST_F(AbstractParam, ValMapsOneToOne)
{
  DummyParam val;

  for (double v = -1e5; v < 1e5; v += 1e-4)
  {
    val.val(v);
    EXPECT_NEAR(v, val.val(), 1e-15);
  }
}


TEST_F(AbstractParam, Put)
{
  DummyParam v;
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

TEST_F(AbstractParam, Get)
{
  DummyParam v;
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

TEST_F(AbstractParam, ValFrom)
{
  DummyParam v;
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

TEST_F(AbstractParam, GradFrom)
{
  DummyParam v;
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

TEST_F(AbstractParam, Visualize)
{
  DummyParam val;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> grad;
  for (double v = -7.0; v < 7.0; v+=0.1)
  {
    x.push_back(v);
    y.push_back(val.val_at(v));
    grad.push_back(val.grad_at(v));
  }

//  MESSAGE() << "AbstractParam::val_at(x):\n" << visualize(x, y, 100) << "\n";
  MESSAGE() << "grad(val):\n" << visualize(y, grad, 100) << "\n";
}

TEST_F(AbstractParam, Print)
{
  DummyParam sval;
  sval.val(1235678e99);
  MESSAGE() << "|" << sval.to_string() << "|\n";
}
