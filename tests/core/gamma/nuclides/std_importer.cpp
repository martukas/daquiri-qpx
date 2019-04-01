#include "gtest_color_print.h"

#include <core/gamma/nuclides/std_importer.h>

#include <date/date.h>

class NuclidStdImporter : public TestBase
{
 protected:
  DAQuiri::NuclideStdImporter importer;
};

TEST_F(NuclidStdImporter, OpenAndScan)
{
  importer.open(std::string(TEST_DATA_PATH) + "/Nuclid.std");

  EXPECT_EQ(importer.stats_total_lines, 1173u);
  EXPECT_EQ(importer.nuclide_lines.size(), 849);
  EXPECT_EQ(importer.certificate_lines.size(), 125);
}

TEST_F(NuclidStdImporter, ReadNuclides)
{
  importer.open(std::string(TEST_DATA_PATH) + "/Nuclid.std");
  importer.parse_nuclides();

  EXPECT_EQ(importer.nuclides.size(), 67);

  for (const auto&n : importer.nuclides)
  {
    MESSAGE()
    << "  " << n.name
//    << "\n   comment:" << n.comment
    << "   HL:" << n.half_life.debug()
//    << "   hl_comment:" << n.half_life_comment
    << "\n";
    for (const auto& g: n.gammas)
    {
      MESSAGE()
          << "     E=" << g.energy.debug()
          << "     I=" << g.abundance.debug()
//          << "     C:" << g.comment
          << "\n";
    }
  }
}

TEST_F(NuclidStdImporter, ReadCertificates)
{
  importer.open(std::string(TEST_DATA_PATH) + "/Nuclid.std");
  importer.parse_certificates();

  EXPECT_EQ(importer.certificates.size(), 31);

  for (const auto&n : importer.certificates)
  {
    MESSAGE()
        << "  " << n.name
        << "   comment:" << n.comment
        << "   Nuclide:" << n.nuclide_id
        << "   Activity:" << n.activity.debug()
        << "   comment:" << n.activity_comment
        << "   Time:" << date::format("%F %T", n.reference_date)
        << "\n";
  }
}
