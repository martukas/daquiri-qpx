#include "gtest_color_print.h"
#include <core/util/visualize_vector.h>

#include <core/fitting/data_model/weighted_data.h>

class WeightedData : public TestBase
{
};


TEST_F(WeightedData, DefaultConstructed)
{
  DAQuiri::WeightedData wd;
  EXPECT_FALSE(wd.valid());
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
  EXPECT_EQ(wd.chan.size(), 4u);
  EXPECT_FALSE(wd.valid());
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
  EXPECT_EQ(wd.chan.size(), 10u);
  EXPECT_DOUBLE_EQ(wd.count_min(), 10.0);
  EXPECT_DOUBLE_EQ(wd.count_max(), 19.0);

  auto indices_too_low = wd.subset(4,9);
  EXPECT_EQ(indices_too_low.chan.size(), 0u);

  auto indices_too_high = wd.subset(20,25);
  EXPECT_EQ(indices_too_high.chan.size(), 0u);
}

TEST_F(WeightedData, SubsetValidRange)
{
  std::vector<double> x;
  for (size_t i=0; i < 10; ++i)
    x.push_back(i+10);
  auto y = x;
  DAQuiri::WeightedData wd(x,y);
  EXPECT_EQ(wd.chan.size(), 10u);

  auto indices_ordered = wd.subset(11,15);
  EXPECT_FALSE(indices_ordered.chan.empty());
  EXPECT_EQ(indices_ordered.chan.front(), 11u);
  EXPECT_EQ(indices_ordered.chan.back(), 15u);
  EXPECT_EQ(indices_ordered.count.front(), 11u);
  EXPECT_EQ(indices_ordered.count.back(), 15u);
  EXPECT_DOUBLE_EQ(indices_ordered.count_min(), 11.0);
  EXPECT_DOUBLE_EQ(indices_ordered.count_max(), 15.0);

  auto indices_reversed = wd.subset(15,11);
  EXPECT_FALSE(indices_reversed.chan.empty());
  EXPECT_EQ(indices_reversed.chan.front(), 11u);
  EXPECT_EQ(indices_reversed.chan.back(), 15u);
  EXPECT_EQ(indices_reversed.count.front(), 11u);
  EXPECT_EQ(indices_reversed.count.back(), 15u);
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
  EXPECT_EQ(wd.chan.size(), 10u);

  auto front_overlap = wd.subset(4,10);
  EXPECT_FALSE(front_overlap.chan.empty());
  EXPECT_EQ(front_overlap.chan.front(), 10u);
  EXPECT_EQ(front_overlap.chan.back(), 10u);
  EXPECT_DOUBLE_EQ(front_overlap.count_min(), 10.0);
  EXPECT_DOUBLE_EQ(front_overlap.count_max(), 10.0);

  auto back_overlap = wd.subset(19,25);
  EXPECT_FALSE(back_overlap.chan.empty());
  EXPECT_EQ(back_overlap.chan.front(), 19u);
  EXPECT_EQ(back_overlap.chan.back(), 19u);
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
  EXPECT_EQ(wd.chan.size(), 10u);

  auto left_subset = wd.left(3);
  EXPECT_FALSE(left_subset.chan.empty());
  EXPECT_EQ(left_subset.chan.size(), 3u);
  EXPECT_EQ(left_subset.chan.front(), 10u);
  EXPECT_EQ(left_subset.chan.back(), 12u);
  EXPECT_EQ(left_subset.count.front(), 10u);
  EXPECT_EQ(left_subset.count.back(), 12u);
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
  EXPECT_EQ(wd.chan.size(), 10u);

  auto right_subset = wd.right(3);
  EXPECT_FALSE(right_subset.chan.empty());
  EXPECT_EQ(right_subset.chan.size(), 3u);
  EXPECT_EQ(right_subset.chan.front(), 17u);
  EXPECT_EQ(right_subset.chan.back(), 19u);
  EXPECT_EQ(right_subset.count.front(), 17u);
  EXPECT_EQ(right_subset.count.back(), 19u);
  EXPECT_DOUBLE_EQ(right_subset.count_min(), 17.0);
  EXPECT_DOUBLE_EQ(right_subset.count_max(), 19.0);
}
