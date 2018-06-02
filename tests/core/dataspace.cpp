#include <gtest/gtest.h>
#include "dataspace.h"

using namespace DAQuiri;

class MockDataspace : public Dataspace
{
  public:
    MockDataspace() {}

    MockDataspace* clone() const override { return new MockDataspace(*this); }

    bool empty() const override { return  true; }

    void clear() override {}
    void add(const Entry&) override {}
    void add_one(const Coords&) override {}
    PreciseFloat get(const Coords&) const override { return 0; }
    EntryList range(std::vector<Pair> list) const override { return EntryList(); }
    void recalc_axes() override {}

    void save(std::ostream& os) override {}

  protected:
    void data_save(hdf5::node::Group) const override {}
    void data_load(hdf5::node::Group) override {}
};

TEST(DataAxis, Init)
{
  DataAxis a;
  EXPECT_TRUE(a.label().empty());
}

TEST(Dataspace, Init)
{
  MockDataspace d;
  EXPECT_TRUE(d.empty());
}