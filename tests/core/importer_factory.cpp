#include "gtest_color_print.h"
#include <core/importer_factory.h>

class ImporterFactory : public TestBase
{
  virtual void TearDown()
  {
    DAQuiri::ImporterFactory::singleton().clear();
  }
};

class MockImporter : public DAQuiri::Importer
{
 public:
  MockImporter() : DAQuiri::Importer() {}
  virtual ~MockImporter() {}

  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override
  {
    (void) path;
    (void) project;
  }
};

class Importer1 : public MockImporter {
 public:
  Importer1() : MockImporter() {}
  ~Importer1() {}

  std::string ext() const override { return "ext1"; }
  std::string description() const override { return "descr1"; }
  bool validate(const boost::filesystem::path& path) const override
  {
    (void) path;
    return true;
  }
};

class Importer2a : public MockImporter {
 public:
  Importer2a() : MockImporter() {}
  ~Importer2a() {}

  std::string ext() const override { return "ext2"; }
  std::string description() const override { return "descr2"; }
  bool validate(const boost::filesystem::path& path) const override { return !path.empty(); }
};

class Importer2b : public MockImporter {
 public:
  Importer2b() : MockImporter() {}
  ~Importer2b() {}

  std::string ext() const override { return "ext2"; }
  std::string description() const override { return "descr3"; }
  bool validate(const boost::filesystem::path& path) const override { return path.empty(); }
};

TEST_F(ImporterFactory, Singleton)
{
  auto& a = DAQuiri::ImporterFactory::singleton();
  auto& b = DAQuiri::ImporterFactory::singleton();
  EXPECT_EQ(&a, &b);
}

TEST_F(ImporterFactory, empty)
{
  auto& cf = DAQuiri::ImporterFactory::singleton();
  EXPECT_TRUE(cf.extensions().empty());
  EXPECT_TRUE(cf.descriptions().empty());
}

TEST_F(ImporterFactory, registered)
{
  auto& cf = DAQuiri::ImporterFactory::singleton();

  DAQUIRI_REGISTER_IMPORTER(Importer1);
  EXPECT_EQ(cf.extensions().size(), 1u);
  EXPECT_TRUE(cf.extensions().count("ext1"));
  EXPECT_EQ(cf.descriptions().size(), 2u);
  EXPECT_EQ(cf.descriptions()[1], "descr1 (*.ext1)");

  DAQUIRI_REGISTER_IMPORTER(Importer2a);
  EXPECT_EQ(cf.extensions().size(), 2u);
  EXPECT_TRUE(cf.extensions().count("ext2"));
  EXPECT_EQ(cf.descriptions().size(), 3u);
  EXPECT_EQ(cf.descriptions()[2], "descr2 (*.ext2)");

  DAQUIRI_REGISTER_IMPORTER(Importer2b);
  EXPECT_EQ(cf.extensions().size(), 2u);
  EXPECT_TRUE(cf.extensions().count("ext2"));
  EXPECT_EQ(cf.descriptions().size(), 4u);
  EXPECT_EQ(cf.descriptions()[3], "descr3 (*.ext2)");
}

TEST_F(ImporterFactory, importer_filtering)
{
  auto& cf = DAQuiri::ImporterFactory::singleton();

  DAQUIRI_REGISTER_IMPORTER(Importer1);
  DAQUIRI_REGISTER_IMPORTER(Importer2a);
  DAQUIRI_REGISTER_IMPORTER(Importer2b);

  auto l1 = cf.attempt_import("whatever.ext1");
  EXPECT_EQ(l1.size(), 1u);
  EXPECT_EQ(l1[0]->description(), "descr1");

  auto l2 = cf.attempt_import("nonempty.ext2");
  EXPECT_EQ(l2.size(), 1u);
  EXPECT_EQ(l2[0]->description(), "descr2");

  auto l3 = cf.attempt_import(".ext2");
  EXPECT_EQ(l3.size(), 1u);
  EXPECT_EQ(l3[0]->description(), "descr3");
}
