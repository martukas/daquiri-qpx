#include "gtest_color_print.h"
#include <consumers/spectrum_time.h>

class TimeSpectrum : public TestBase
{
  protected:
    virtual void SetUp()
    {
      h.set_attribute(DAQuiri::Setting::text("stream_id", "stream"));
      h.set_attribute(DAQuiri::Setting::text("value_latch/value_id", "val"));
      h.set_attribute(DAQuiri::Setting::floating("time_resolution", 1));
      h.set_attribute(DAQuiri::Setting::integer("time_units", 0));

      s.event_model.add_value("val", 100);
      s.event_model.add_value("val2", 100);
      s.events.reserve(3, s.event_model);
      s.events.last().set_time(10);
      s.events.last().set_value(0, 0);
      s.events.last().set_value(1, 0);
      ++s.events;
      s.events.last().set_time(20);
      s.events.last().set_value(0, 1);
      s.events.last().set_value(1, 15);
      ++s.events;
      s.events.last().set_time(20);
      s.events.last().set_value(0, 1);
      s.events.last().set_value(1, 30);
      ++s.events;
      s.events.finalize();
    }

    DAQuiri::TimeSpectrum h;
    DAQuiri::Spill s{"stream", DAQuiri::Spill::Type::start};
};

TEST_F(TimeSpectrum, DefaultConstructed)
{
  DAQuiri::TimeSpectrum default_h;
  EXPECT_FALSE(default_h.changed());
  EXPECT_EQ(default_h.type(), "TimeSpectrum 2D");
  EXPECT_EQ(default_h.dimensions(), 2);
}

TEST_F(TimeSpectrum, HistogramsEvents)
{
  h.push_spill(s);

  EXPECT_EQ(h.metadata().get_attribute("total_count").get_number(), 3);

  auto data = h.data();
  EXPECT_EQ(data->get({10,0}), 1);
  EXPECT_EQ(data->get({20,1}), 2);
}

TEST_F(TimeSpectrum, HistogramDecimate)
{
  h.set_attribute(DAQuiri::Setting::floating("time_resolution", 10));

  h.push_spill(s);

  EXPECT_EQ(h.metadata().get_attribute("total_count").get_number(), 3);

  auto data = h.data();
  EXPECT_EQ(data->get({1,0}), 1);
  EXPECT_EQ(data->get({2,1}), 2);
}

TEST_F(TimeSpectrum, LatchesStream)
{
  h.set_attribute(DAQuiri::Setting::text("stream_id", "N/A"));

  h.push_spill(s);

  EXPECT_EQ(h.metadata().get_attribute("total_count").get_number(), 0);
}

TEST_F(TimeSpectrum, ZeroResolutionBinsNothing)
{
  h.set_attribute(DAQuiri::Setting::floating("time_resolution", 0));
  h.push_spill(s);

  EXPECT_EQ(h.metadata().get_attribute("total_count").get_number(), 0);
}

TEST_F(TimeSpectrum, FilterByValue)
{
  h.set_attribute(DAQuiri::Setting::integer("filter_count", 1));

  auto filters = h.metadata().get_attribute("filters");

  auto fe = DAQuiri::Setting::boolean("filter/enabled", true);
  fe.set_indices({0});
  h.set_attribute(fe);

  auto fn = DAQuiri::Setting::text("filter/value_id", "val2");
  fn.set_indices({0});
  h.set_attribute(fn);

  auto fmin = DAQuiri::Setting::integer("filter/min", 10);
  fmin.set_indices({0});
  h.set_attribute(fmin);

  auto fmax = DAQuiri::Setting::integer("filter/max", 20);
  fmax.set_indices({0});
  h.set_attribute(fmax);

  h.push_spill(s);

  EXPECT_EQ(h.metadata().get_attribute("total_count").get_number(), 1);
}

//TODO: test other time units
//TODO: test data axes

TEST_F(TimeSpectrum, Clone)
{
  h.push_spill(s);

  auto h_copy = std::shared_ptr<DAQuiri::TimeSpectrum>(h.clone());

  EXPECT_NE(h_copy.get(), &h);
  EXPECT_EQ(h_copy->metadata().get_attribute("total_count").get_number(), 3);

  auto data = h_copy->data();
  EXPECT_EQ(data->get({10,0}), 1);
  EXPECT_EQ(data->get({20,1}), 2);
}
