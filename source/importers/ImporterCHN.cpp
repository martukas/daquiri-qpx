#include "ImporterCHN.h"
#include <core/consumer_factory.h>
#include <date/date.h>
#include <iostream>

#include <core/util/custom_logger.h>

struct __attribute__ ((packed)) ChnHeader
{
  int16_t must_be_neg1;
  int16_t detector;
  int16_t segment;
  char start_seconds[2];
  int32_t real_time;
  int32_t live_time;
  char start_date[8];
  char start_time[4];
  int16_t channel_offset;
  int16_t num_channels;
};

bool ImporterCHN::validate(const boost::filesystem::path& path) const
{
  return true;
  std::ifstream file(path.string(), std::ios::binary);
  int16_t must_be_neg1;
  file.read(reinterpret_cast<char*>(&must_be_neg1), sizeof(int16_t));
  return (must_be_neg1 == -1);
}

void ImporterCHN::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::ifstream file(path.string(), std::ios::binary);

  ChnHeader header;
  file.read(reinterpret_cast<char*>(&header), sizeof(ChnHeader));

  if (header.must_be_neg1 != -1)
    throw std::runtime_error("<ImporterCHN> invalid CHN file, marker != -1");

  std::string start_day(&header.start_date[0], 2);
  std::string start_month(&header.start_date[2], 3);
  std::string start_year =
      ((header.start_date[7] == '1') ? "20" : "19")
          + std::string(&header.start_date[5], 2);

  std::string start_hours(&header.start_time[0], 2);
  std::string start_mins(&header.start_time[2], 2);
  std::string start_secs(&header.start_seconds[0], 2);

  date::sys_time<std::chrono::seconds> t;
  std::stringstream stream;
  stream << start_year << "-" << start_month << "-" << start_day
         << " " << start_hours << ":" << start_mins << ":" << start_secs;

  INFO("Detector = {}", header.detector);
  INFO("Segment = {}", header.segment);

  INFO("Start time = {}", stream.str());
  stream >> date::parse("%Y-%b-%d %H:%M:%S", t);
  if (stream.fail())
    throw std::runtime_error("<ImporterCHN> failed to parse date");
  INFO("parsed start time = {}", date::format("%FT%TZ", t));

  std::vector<int32_t> spectrum(static_cast<size_t>(header.num_channels));
  file.read(reinterpret_cast<char*>(spectrum.data()), spectrum.size() * sizeof(int32_t));

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("<ImporterCHN> could not get a valid Histogram 1D from factory");

  hist->set_attribute(DAQuiri::Setting("start_time", t));

  auto lt_ms = static_cast<int64_t>(header.live_time) * 20;
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));

  auto rt_ms = static_cast<int64_t>(header.real_time) * 20;
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));

  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    DAQuiri::Entry new_entry;
    new_entry.first.resize(1);
    new_entry.first[0] = i + header.channel_offset;
    new_entry.second = PreciseFloat(spectrum[i]);
    entry_list.push_back(new_entry);
  }

  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));

  hist->import(*this);

  project->add_consumer(hist);
}
