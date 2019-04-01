#include <core/gamma/nuclides/std_importer.h>
#include <core/util/string_extensions.h>
#include <date/date.h>

namespace DAQuiri
{

void NuclideStdImporter::open(const boost::filesystem::path& path)
{
  std::ifstream myfile(path.string(), std::ios::in);

  bool preamble_finished{false};
  bool certificates_started{false};
  for (std::string line; getline(myfile, line);)
  {
    stats_total_lines++;
    trim(line);
    if (line.empty())
      continue;
    if (line[0] == '#')
      preamble_finished = true;
    if (line.substr(0, 14) == "#CERTIFICATES:")
    {
      certificates_started = true;
      continue;
    }
    if (certificates_started)
      certificate_lines.push_back(line);
    else if (preamble_finished)
      nuclide_lines.push_back(line);
  }

  // last non-empty line of file is just termination string
  certificate_lines.pop_back();
}

void NuclideStdImporter::split_comments(const std::string& line,
                                        std::string& stuff, std::string& comments)
{
  auto comment_start = line.find(";");
  if (comment_start != std::string::npos)
  {
    stuff = trim_copy(line.substr(0, comment_start));
    comments = trim_copy(line.substr(comment_start + 1,
                                     line.size() - comment_start - 1));
  }
  else
    stuff = trim_copy(line);
}

Isotope NuclideStdImporter::parse_nuclide_line(const std::string& line)
{
  Isotope ret;
  split_comments(line, ret.name, ret.comment);
  return ret;
}

void NuclideStdImporter::parse_hl_line(const std::string& line, Isotope& nuclide)
{
  std::string data;
  split_comments(line, data, nuclide.half_life_comment);
  std::stringstream ss(data);
  double hl, hl_sigma;
  ss >> hl >> hl_sigma;
  nuclide.half_life = {hl, hl_sigma};
}

void NuclideStdImporter::parse_comment_line(const std::string& line, Isotope& nuclide)
{
  auto truncated = trim_copy(line);
  //if (!nuclide.comment.empty())
  nuclide.comment += "\n";
  nuclide.comment += truncated;
}

Radiation NuclideStdImporter::parse_gamma_line(const std::string& line)
{
  Radiation ret;
  std::string data;
  split_comments(line, data, ret.comment);
  std::stringstream ss(data);
  double e, e_sigma, i, i_sigma;
  ss >> e >> e_sigma >> i >> i_sigma;
  ret.energy = {e, e_sigma};
  ret.abundance = {i, i_sigma};
  return ret;
}

void NuclideStdImporter::parse_nuclides()
{
  Isotope nuclide;
  for (auto line : nuclide_lines)
  {
    if (line[0] == '#')
    {
      if (nuclide.name.size())
        nuclides.push_back(nuclide);
      nuclide = parse_nuclide_line(line.substr(1, line.size() - 1));
    }
    else if ((line.size() > 4) && line.substr(0, 5) == "T1/2:")
      parse_hl_line(line.substr(5, line.size() - 5), nuclide);
    else if (line[0] == ';')
      parse_comment_line(line.substr(1, line.size() - 1), nuclide);
    else
      nuclide.gammas.add_a(parse_gamma_line(line));
  }

  if (nuclide.name.size())
    nuclides.push_back(nuclide);
}

Certificate NuclideStdImporter::parse_certificate_line(const std::string& line)
{
  Certificate ret;
  split_comments(line, ret.name, ret.comment);
  return ret;
}

void NuclideStdImporter::parse_activity_line(const std::string& line,
                                             Certificate& certificate)
{
  std::string data;
  split_comments(line, data, certificate.activity_comment);
  std::stringstream ss(trim_copy(data));
  double a, a_sigma;
  ss >> a >> data;
  std::stringstream ss2(data.substr(1, data.size() - 3));
  ss2 >> a_sigma;
  certificate.activity = {a, a_sigma};
}

hr_time_t NuclideStdImporter::parse_reference_time(const std::string& line)
{
  // mm-dd-yyyy hh:mm:ss
  hr_time_t ret;
  std::istringstream stream{line};
  stream >> date::parse("%m-%d-%Y %H:%M:%S", ret);
  return ret;
}

void NuclideStdImporter::parse_certificates()
{
  Certificate certificate;
  for (auto line : certificate_lines)
  {
    if (line[0] == '#')
    {
      if (certificate.name.size())
        certificates.push_back(certificate);
      certificate = parse_certificate_line(line);
    }
    else if ((line.size() > 7) && line.substr(0, 8) == "Nuclide:")
      certificate.nuclide_id = trim_copy(line.substr(8, line.size() - 8));
    else if ((line.size() > 8) && line.substr(0, 9) == "Activity:")
      parse_activity_line(trim_copy(line.substr(9, line.size() - 9)), certificate);
    else if ((line.size() > 7) && line.substr(0, 8) == "RefDate:")
      certificate.reference_date =
          parse_reference_time(trim_copy(line.substr(8, line.size() - 8)));
  }

  if (certificate.name.size())
    certificates.push_back(certificate);
}

}
