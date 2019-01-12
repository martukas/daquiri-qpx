#include "gtest_color_print.h"

#include <core/consumer_factory.h>
#include <consumers/histogram_1d.h>

#include <importers/ImporterTKA.h>

class ImportTKA : public TestBase
{
  virtual void SetUp()
  {
    using namespace DAQuiri;
    DAQUIRI_REGISTER_CONSUMER(Histogram1D)
  }

  virtual void TearDown()
  {
    DAQuiri::ImporterFactory::singleton().clear();
  }
};


TEST_F(ImportTKA, ImportPrompt1)
{
  ImporterTKA importer;
  auto p = std::make_shared<DAQuiri::Project>();
  importer.import(std::string(TEST_DATA_PATH) + "/prompt1.tka", p);

  auto cs = p->get_consumers();
  EXPECT_EQ(cs.size(), 1);

  auto c = cs.get(0);
  EXPECT_EQ(c->type(), "Histogram 1D");

  auto md = c->metadata();
  auto lt = md.get_attribute("live_time");
  auto lt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(lt.duration());
  auto rt = md.get_attribute("real_time");
  auto rt_ms = std::chrono::duration_cast<std::chrono::milliseconds>(rt.duration());
  EXPECT_EQ(lt.metadata().type(), DAQuiri::SettingType::duration);
  EXPECT_EQ(lt_ms, std::chrono::milliseconds(65694899)) << lt_ms.count();
  EXPECT_EQ(rt.metadata().type(), DAQuiri::SettingType::duration);
  EXPECT_EQ(rt_ms, std::chrono::milliseconds(65884300)) << rt_ms.count();

  auto data = c->data();
  auto list = data->all_data();
  EXPECT_EQ(list->size(), 16331);

  auto p10 = list->at(10);
  EXPECT_EQ(p10.first[0], 10);
  EXPECT_EQ(p10.second, 0);

  auto p100 = list->at(100);
  EXPECT_EQ(p100.first[0], 100);
  EXPECT_EQ(p100.second, 4144);

  auto p1000 = list->at(1000);
  EXPECT_EQ(p1000.first[0], 1000);
  EXPECT_EQ(p1000.second, 437);

  auto p10000 = list->at(10000);
  EXPECT_EQ(p10000.first[0], 10000);
  EXPECT_EQ(p10000.second, 5);
}
