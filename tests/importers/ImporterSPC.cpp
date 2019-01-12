#include "gtest_color_print.h"

#include <core/consumer_factory.h>
#include <consumers/histogram_1d.h>

#include <core/calibration/coef_function_factory.h>
#include <core/calibration/polynomial.h>

#include <importers/ImporterSPC.h>
#include <date/date.h>

class ImportSPC : public TestBase
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
  ImporterSPC importer;
};


TEST_F(ImportSPC, ImportMarine)
{
  importer.import(std::string(TEST_DATA_PATH) + "/IAEA433_Marine_sediment_Q39.spc", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 2u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(10486700)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(11595799)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2014);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 2u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 11u);
  EXPECT_EQ(tod.hours().count(), 13);
  EXPECT_EQ(tod.minutes().count(), 14);
  EXPECT_EQ(tod.seconds().count(), 46);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), 0, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.212068, 0.000001);
  EXPECT_NEAR(axis.domain.front(), 0, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 3473.039937, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 16378u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 0);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 0);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 14606);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 3);



  auto c2 = cs.get(1);
  EXPECT_EQ(c2->type(), "Histogram 1D");

  auto md2 = c2->metadata();

  auto lt2 = md2.get_attribute("live_time");
  EXPECT_EQ(lt2.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(lt2.duration());
  EXPECT_EQ(lt_ms2, std::chrono::milliseconds(11595799)) << lt_ms2.count();

  auto rt2 = md2.get_attribute("real_time");
  EXPECT_EQ(rt2.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(rt2.duration());
  EXPECT_EQ(rt_ms2, std::chrono::milliseconds(11595799)) << rt_ms2.count();

  auto st2 = md2.get_attribute("start_time");
  EXPECT_EQ(st2.metadata().type(), DAQuiri::SettingType::time);
  auto converted2 = st2.time();
  EXPECT_EQ(converted, converted2);

  auto data2 = c2->data();

  auto axis2 = data2->axis(0);
  EXPECT_EQ(axis.calibration, axis2.calibration);

  auto list2 = data2->all_data();
  EXPECT_EQ(list2->size(), 16378u);

  auto p10b = list2->at(10);
  EXPECT_EQ(p10b.first[0], 10u);
  EXPECT_EQ(p10b.second, 0);

  auto p100b = list2->at(100);
  EXPECT_EQ(p100b.first[0], 100u);
  EXPECT_EQ(p100b.second, 0);

  auto p1000b = list2->at(1000);
  EXPECT_EQ(p1000b.first[0], 1000u);
  EXPECT_EQ(p1000b.second, 16111);

  auto p10000b = list2->at(10000);
  EXPECT_EQ(p10000b.first[0], 10000u);
  EXPECT_EQ(p10000b.second, 3);
}
