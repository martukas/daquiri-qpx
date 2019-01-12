#include "gtest_color_print.h"

#include <core/consumer_factory.h>
#include <consumers/histogram_1d.h>

#include <core/calibration/coef_function_factory.h>
#include <core/calibration/polynomial.h>

#include <importers/ImporterMCA.h>
#include <date/date.h>

class ImportMCA : public TestBase
{
 protected:
  virtual void SetUp()
  {
    DAQuiri::ConsumerRegistrar<DAQuiri::Histogram1D> h1d;
    DAQuiri::CoefFunctionRegistrar<DAQuiri::Polynomial> pol;
    p = std::make_shared<DAQuiri::Project>();
  }

  virtual void TearDown()
  {
    DAQuiri::ImporterFactory::singleton().clear();
    DAQuiri::CoefFunctionFactory::singleton().clear();
    DAQuiri::ConsumerFactory::singleton().clear();
  }

  DAQuiri::ProjectPtr p;
  ImporterMCA importer;
};


TEST_F(ImportMCA, ImportB79)
{
  importer.import(std::string(TEST_DATA_PATH) + "/B79hossz.mca", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(13912630)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(13972200)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2004);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 4u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 21u);
  EXPECT_EQ(tod.hours().count(), 12);
  EXPECT_EQ(tod.minutes().count(), 46);
  EXPECT_EQ(tod.seconds().count(), 24);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -0.043826, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.708278, 0.000001);
  EXPECT_NEAR(axis.domain.front(), -0.043826, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 11165.251511, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 15765u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 0);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 6);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 478);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 14);
}
