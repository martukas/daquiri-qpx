#include "gtest_color_print.h"

#include <core/consumer_factory.h>
#include <consumers/histogram_1d.h>

#include <core/calibration/function_factory.h>
#include <core/calibration/polynomial.h>

#include <importers/ImporterN42.h>
#include <date/date.h>

class ImportN42 : public TestBase
{
 protected:
  virtual void SetUp()
  {
    DAQuiri::ConsumerRegistrar<DAQuiri::Histogram1D> h1d;
    DAQuiri::CalibFunctionRegistrar<DAQuiri::Polynomial> pol;
    p = std::make_shared<DAQuiri::Project>();
  }

  virtual void TearDown()
  {
    DAQuiri::ImporterFactory::singleton().clear();
    DAQuiri::CalibFunctionFactory::singleton().clear();
    DAQuiri::ConsumerFactory::singleton().clear();
  }

  DAQuiri::ProjectPtr p;
  ImporterN42 importer;
};


TEST_F(ImportN42, ImportDaphne)
{
  importer.import(std::string(TEST_DATA_PATH) + "/Daphne_Eu152.n42", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(21284000)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(21600000)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2014);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 12u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 11u);
  EXPECT_EQ(tod.hours().count(), 12);
  EXPECT_EQ(tod.minutes().count(), 51);
  EXPECT_EQ(tod.seconds().count(), 28);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -0.422413, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.249778, 0.000001);
  EXPECT_NEAR(func->coeffs()[2].value(), 7.32964e-09, 1e-14);
  EXPECT_NEAR(axis.domain.front(), -0.422413, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 3855.391725, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 15431u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 155);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 31);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 31066);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 20);
}

TEST_F(ImportN42, ImportSophia)
{
  importer.import(std::string(TEST_DATA_PATH) + "/Sophia_Eu152.n42", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(21273000)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(21600000)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2014);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 12u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 11u);
  EXPECT_EQ(tod.hours().count(), 12);
  EXPECT_EQ(tod.minutes().count(), 51);
  EXPECT_EQ(tod.seconds().count(), 44);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -0.505959, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.250237, 0.000001);
  EXPECT_NEAR(func->coeffs()[2].value(), 1.05746e-08, 1e-13);
  EXPECT_NEAR(axis.domain.front(), -0.505959, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 3856.660671, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 15405u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 150);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 102);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 37456);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 14);
}
