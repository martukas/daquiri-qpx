#pragma once

#include <core/gamma/nuclides/nuclide.h>
#include <boost/filesystem.hpp>
#include <core/util/time_extensions.h>

namespace DAQuiri
{

struct Certificate
{
  std::string name;
  std::string comment;
  std::string nuclide_id;
  UncertainDouble activity;
  std::string activity_comment;
  hr_time_t reference_date;
};

class NuclideStdImporter
{
 public:
  NuclideStdImporter() = default;
  void open(const boost::filesystem::path& path);

  void parse_nuclides();

  void parse_certificates();

  size_t stats_total_lines{0};
  std::vector<std::string> nuclide_lines;
  std::vector<std::string> certificate_lines;

  std::vector<Isotope> nuclides;
  std::vector<Certificate> certificates;

  static void split_comments(const std::string& line,
                             std::string& stuff, std::string& comments);

  static Isotope parse_nuclide_line(const std::string& line);
  static void parse_hl_line(const std::string& line, Isotope& nuclide);
  static void parse_comment_line(const std::string& line, Isotope& nuclide);
  static Radiation parse_gamma_line(const std::string& line);

  static Certificate parse_certificate_line(const std::string& line);
  static void parse_activity_line(const std::string& line, Certificate& certificate);
  static hr_time_t parse_reference_time(const std::string& line);

};

}
