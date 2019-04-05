#include "gtest_color_print.h"
#include <core/calibration/function_factory.h>

class CalibFunctionFactory : public TestBase
{
  virtual void TearDown() {
    DAQuiri::CalibFunctionFactory::singleton().clear();
  }
};

class FakeFunction : public DAQuiri::CalibFunction
{
 public:
  // Inherit constructors
  using DAQuiri::CalibFunction::CalibFunction;

  bool valid() const override { return true; }

  bool is_equal(CalibFunction* other) const override { return  true; }

  double d_dx(double) const override { return 1; }

  void update_indices() override {}

  Eigen::VectorXd variables() const override
  {
    return {};
  }

  double eval(double x) const override
  {
    return 1 + x;
  }

  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override
  {
    return 1 + chan;
  }

  void save_fit(const DAQuiri::FitResult& result) override {}

  std::string to_string(std::string prepend) const override { return ""; }
  std::string to_UTF8(int precision) const override { return ""; }
  std::string to_markup(int precision) const override { return ""; }

  void from_json(const nlohmann::json&) override {};
};

class Function1 : public FakeFunction
{
 public:
  Function1* clone() const override { return new Function1(*this); }
  std::string type() const override { return "Function1"; }
};

class Function2 : public FakeFunction
{
 public:
  Function2* clone() const override { return new Function2(*this); }
  std::string type() const override { return "Function2"; }
};

class Function3 : public FakeFunction
{
 public:
  Function3* clone() const override { return new Function3(*this); }
  std::string type() const override { return "Function3"; }
};


TEST_F(CalibFunctionFactory, Singleton)
{
  auto& a = DAQuiri::CalibFunctionFactory::singleton();
  auto& b = DAQuiri::CalibFunctionFactory::singleton();
  EXPECT_EQ(&a, &b);
}


TEST_F(CalibFunctionFactory, types)
{
  auto& cf = DAQuiri::CalibFunctionFactory::singleton();

  EXPECT_TRUE(cf.types().empty());
  EXPECT_EQ(cf.types().size(), 0UL);

  DAQUIRI_REGISTER_COEF_FUNCTION(Function1);
  EXPECT_EQ(cf.types().size(), 1UL);

  DAQUIRI_REGISTER_COEF_FUNCTION(Function2);
  EXPECT_EQ(cf.types().size(), 2UL);

  DAQUIRI_REGISTER_COEF_FUNCTION(Function3);
  EXPECT_EQ(cf.types().size(), 3UL);

  cf.clear();
  EXPECT_EQ(cf.types().size(), 0UL);
  EXPECT_TRUE(cf.types().empty());
}

TEST_F(CalibFunctionFactory, create_type)
{
  auto& cf = DAQuiri::CalibFunctionFactory::singleton();
  DAQUIRI_REGISTER_COEF_FUNCTION(Function1);
  DAQUIRI_REGISTER_COEF_FUNCTION(Function2);

  EXPECT_FALSE(cf.create_type("bad_id"));
  EXPECT_EQ(cf.create_type("Function1")->type(), "Function1");
  EXPECT_EQ(cf.create_type("Function2")->type(), "Function2");
}

TEST_F(CalibFunctionFactory, create_copy)
{
  auto& cf = DAQuiri::CalibFunctionFactory::singleton();
  DAQUIRI_REGISTER_COEF_FUNCTION(Function1);

  auto c1 = cf.create_type("Function1");
  auto c2 = cf.create_copy(c1);

  EXPECT_EQ(c2->type(), "Function1");
  EXPECT_NE(c1.get(), c2.get());

  EXPECT_FALSE(cf.create_copy(nullptr));
}

TEST_F(CalibFunctionFactory, create_from_json)
{
  auto& cf = DAQuiri::CalibFunctionFactory::singleton();
  DAQUIRI_REGISTER_COEF_FUNCTION(Function1);

  auto c1 = cf.create_type("Function1");

  nlohmann::json j = c1->to_json();

  auto c2 = cf.create_from_json(j);

  EXPECT_EQ(c2->type(), "Function1");
  EXPECT_NE(c1.get(), c2.get());
}
