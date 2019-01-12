#include "gtest_color_print.h"

#include <core/consumer_factory.h>
#include <consumers/histogram_1d.h>

#include <core/calibration/coef_function_factory.h>
#include <core/calibration/polynomial.h>

#include <importers/ImporterCNF.h>
#include <date/date.h>

class ImportCNF : public TestBase
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
  ImporterCNF importer;
};


TEST_F(ImportCNF, ImportDaphne)
{
  importer.import(std::string(TEST_DATA_PATH) + "/Daphne_Co60.cnf", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(1789000)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(1801000)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 1858);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 11u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 17u);
  EXPECT_EQ(tod.hours().count(), 0);
  EXPECT_EQ(tod.minutes().count(), 0);
  EXPECT_EQ(tod.seconds().count(), 0);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -0.760052, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.249888, 0.000001);
  EXPECT_NEAR(axis.domain.front(), -0.760052, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 3854.015321, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 15427u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 0);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 0);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 414);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 5);
}


TEST_F(ImportCNF, ImportPromptLithium)
{
  importer.import(std::string(TEST_DATA_PATH) + "/prompt lithium.cnf", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(40026530)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(40093990)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2016);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 2u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 14u);
  EXPECT_EQ(tod.hours().count(), 22);
  EXPECT_EQ(tod.minutes().count(), 35);
  EXPECT_EQ(tod.seconds().count(), 48);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -1.50776, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.713635, 0.000001);
  EXPECT_NEAR(axis.domain.front(), -1.50776, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 11684.978999, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 16377u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 327);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 7052);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 685);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 0);
}

TEST_F(ImportCNF, ImportPrompt0)
{
  importer.import(std::string(TEST_DATA_PATH) + "/prompt0.cnf", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1u);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();

  auto lt = md.get_attribute("live_time");
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(4812860)) << lt_ms.count();

  auto rt = md.get_attribute("real_time");
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(4819640)) << rt_ms.count();

  auto st = md.get_attribute("start_time");
  EXPECT_EQ(st.metadata().type(), DAQuiri::SettingType::time);
  auto converted = st.time();
  auto daypoint = date::floor<date::days>(converted);
  auto ymd = date::year_month_day(daypoint);   // calendar date
  auto tod = date::make_time(converted - daypoint); // Yields time_of_day type
  EXPECT_EQ(static_cast<int>(ymd.year()), 2015);
  EXPECT_EQ(static_cast<unsigned>(ymd.month()), 8u);
  EXPECT_EQ(static_cast<unsigned>(ymd.day()), 19u);
  EXPECT_EQ(tod.hours().count(), 14);
  EXPECT_EQ(tod.minutes().count(), 35);
  EXPECT_EQ(tod.seconds().count(), 26);

  auto data = c->data();

  auto axis = data->axis(0);
  EXPECT_EQ(axis.calibration.from(), DAQuiri::CalibID("energy", "unknown", ""));
  EXPECT_EQ(axis.calibration.to(), DAQuiri::CalibID("energy", "unknown", "keV"));
  auto func = axis.calibration.function();
  EXPECT_TRUE(func);
  EXPECT_EQ(func->type(), "Polynomial");
  EXPECT_NEAR(func->coeffs()[0].value(), -1.58029, 0.000001);
  EXPECT_NEAR(func->coeffs()[1].value(), 0.713815, 0.000001);
  EXPECT_NEAR(axis.domain.front(), -1.58029, 0.000001);
  EXPECT_NEAR(axis.domain.back(), 11644.313068, 0.000001);

  auto list = data->all_data();
  EXPECT_EQ(list->size(), 16316u);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10u);
  EXPECT_EQ(p10.second, 1);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100u);
  EXPECT_EQ(p100.second, 127);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000u);
  EXPECT_EQ(p1000.second, 21);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000u);
  EXPECT_EQ(p10000.second, 1);
}
