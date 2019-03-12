#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/sum4/sum4edge.h>

class SUM4Edge : public TestBase
{
 protected:
  void SetUp() override
  {
    std::vector<double> x;
    std::vector<double> y;
    for (size_t i=0; i < 10; ++i)
    {
      x.push_back(i + 10);
      y.push_back(2.5);
    }
    wd = DAQuiri::WeightedData(x,y);
  }

  DAQuiri::WeightedData wd;
};

TEST_F(SUM4Edge, DefaultConstructed)
{
  DAQuiri::SUM4Edge e;
  EXPECT_EQ(e.width(), 0);
  EXPECT_EQ(e.sum(), UncertainDouble(0, 0));
}

TEST_F(SUM4Edge, Construct)
{
  DAQuiri::SUM4Edge e(wd);
  EXPECT_DOUBLE_EQ(e.width(), 10);
  EXPECT_EQ(e.left(), 10);
  EXPECT_EQ(e.right(), 19);
  //EXPECT_EQ(e.sum(), UncertainDouble(25.0, 5.0)) << e.sum().debug();
  auto sum = e.sum();
  EXPECT_DOUBLE_EQ(sum.value(), 25);
  EXPECT_DOUBLE_EQ(sum.sigma(), 5);
  //EXPECT_EQ(e.average(), UncertainDouble(2.5, 0.5)) << e.average().debug();
  auto av = e.average();
  EXPECT_DOUBLE_EQ(av.value(), 2.5);
  EXPECT_DOUBLE_EQ(av.sigma(), 0.5);

  EXPECT_DOUBLE_EQ(e.variance(), 0.25);

// \todo more tests for confirming math; talk to Dick
}

TEST_F(SUM4Edge, GeneratePolynomial)
{
  DAQuiri::SUM4Edge e1(wd);

  std::vector<double> x;
  std::vector<double> y;
  for (size_t i=19; i <= 29; ++i)
  {
    x.push_back(i + 10);
    y.push_back(12.5);
  }
  wd = DAQuiri::WeightedData(x,y);

  DAQuiri::SUM4Edge e2(wd);

  auto poly = DAQuiri::SUM4Edge::sum4_background(e1, e2);
  auto coeffs = poly.coeffs();
  EXPECT_EQ(coeffs.size(), 2u);
  EXPECT_DOUBLE_EQ(coeffs[0].value(), 2.5);
  EXPECT_DOUBLE_EQ(coeffs[1].value(), 1.0);
}

TEST_F(SUM4Edge, JsonEmpty)
{
  DAQuiri::SUM4Edge ee;
  nlohmann::json j = ee;
  DAQuiri::SUM4Edge ee2 = j;
  EXPECT_EQ(ee2.width(), 0);
}

TEST_F(SUM4Edge, JsonNonEmpty)
{
  DAQuiri::SUM4Edge ee(wd);
  nlohmann::json j = ee;
  DAQuiri::SUM4Edge ee2 = j;
  EXPECT_DOUBLE_EQ(ee2.width(), 10);
}
