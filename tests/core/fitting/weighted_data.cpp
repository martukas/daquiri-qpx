#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/weighted_data.h>

class WeightedData : public TestBase
{
};


TEST_F(WeightedData, DefaultConstructed)
{
  DAQuiri::WeightedData wd;
  EXPECT_TRUE(wd.empty());
}

TEST_F(WeightedData, ConstructionFailsIfEmpty)
{
  std::vector<double> x;
  std::vector<double> y;
  EXPECT_ANY_THROW(DAQuiri::WeightedData(x,y));
}

TEST_F(WeightedData, ConstructionFailsIfSizesDiffer)
{
  std::vector<double> x(3);
  std::vector<double> y(4);
  EXPECT_ANY_THROW(DAQuiri::WeightedData(x,y));
}

TEST_F(WeightedData, ConstructionSucceeds)
{
  std::vector<double> x(4);
  std::vector<double> y(4);
  DAQuiri::WeightedData wd;
  EXPECT_NO_THROW(wd = DAQuiri::WeightedData(x,y));
  EXPECT_EQ(wd.data.size(), 4u);
  EXPECT_FALSE(wd.empty());
  EXPECT_DOUBLE_EQ(wd.count_min(), 0.0);
  EXPECT_DOUBLE_EQ(wd.count_max(), 0.0);
}

TEST_F(WeightedData, SubsetEmptyOutOfRange)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.data.size(), 10u);
  EXPECT_DOUBLE_EQ(wd.count_min(), 10.0);
  EXPECT_DOUBLE_EQ(wd.count_max(), 19.0);

  auto indices_too_low = wd.subset(4,9);
  EXPECT_TRUE(indices_too_low.empty());

  auto indices_too_high = wd.subset(20,25);
  EXPECT_TRUE(indices_too_high.empty());
}

TEST_F(WeightedData, SubsetValidRange)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.data.size(), 10u);

  auto indices_ordered = wd.subset(11,15);
  EXPECT_FALSE(indices_ordered.empty());
  EXPECT_EQ(indices_ordered.data.front().channel, 11u);
  EXPECT_EQ(indices_ordered.data.back().channel, 15u);
  EXPECT_DOUBLE_EQ(indices_ordered.count_min(), 11.0);
  EXPECT_DOUBLE_EQ(indices_ordered.count_max(), 15.0);

  auto indices_reversed = wd.subset(15,11);
  EXPECT_FALSE(indices_reversed.empty());
  EXPECT_EQ(indices_reversed.data.front().channel, 11u);
  EXPECT_EQ(indices_reversed.data.back().channel, 15u);
  EXPECT_DOUBLE_EQ(indices_reversed.count_min(), 11.0);
  EXPECT_DOUBLE_EQ(indices_reversed.count_max(), 15.0);
}

TEST_F(WeightedData, SubsetPartial)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.data.size(), 10u);

  auto front_overlap = wd.subset(4,10);
  EXPECT_FALSE(front_overlap.empty());
  EXPECT_EQ(front_overlap.data.front().channel, 10u);
  EXPECT_EQ(front_overlap.data.back().channel, 10u);
  EXPECT_DOUBLE_EQ(front_overlap.count_min(), 10.0);
  EXPECT_DOUBLE_EQ(front_overlap.count_max(), 10.0);

  auto back_overlap = wd.subset(19,25);
  EXPECT_FALSE(back_overlap.empty());
  EXPECT_EQ(back_overlap.data.front().channel, 19u);
  EXPECT_EQ(back_overlap.data.back().channel, 19u);
  EXPECT_DOUBLE_EQ(back_overlap.count_min(), 19.0);
  EXPECT_DOUBLE_EQ(back_overlap.count_max(), 19.0);
}

TEST_F(WeightedData, LeftSubset)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.data.size(), 10u);

  auto left_subset = wd.left(3);
  EXPECT_FALSE(left_subset.empty());
  EXPECT_EQ(left_subset.data.size(), 3u);
  EXPECT_EQ(left_subset.data.front().channel, 10u);
  EXPECT_EQ(left_subset.data.back().channel, 12u);
  EXPECT_DOUBLE_EQ(left_subset.count_min(), 10.0);
  EXPECT_DOUBLE_EQ(left_subset.count_max(), 12.0);
}

TEST_F(WeightedData, RightSubset)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.data.size(), 10u);

  auto right_subset = wd.right(3);
  EXPECT_FALSE(right_subset.empty());
  EXPECT_EQ(right_subset.data.size(), 3u);
  EXPECT_EQ(right_subset.data.front().channel, 17u);
  EXPECT_EQ(right_subset.data.back().channel, 19u);
  EXPECT_DOUBLE_EQ(right_subset.count_min(), 17.0);
  EXPECT_DOUBLE_EQ(right_subset.count_max(), 19.0);
}
