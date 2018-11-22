#include "ImporterCHN.h"
#include <core/consumer_factory.h>
#include <date/date.h>
#include <iostream>

#include <core/util/custom_logger.h>

struct __attribute__ ((packed)) ChnHeader
{
  int16_t must_be_neg1;
  int16_t detector_num;
  int16_t segment_num;
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
}

void ImporterCHN::import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project)
{
  std::ifstream file (path.string(), std::ios::binary);

  std::vector<char> buf(sizeof(ChnHeader));

  INFO("Read array size: {}", buf.size());
  file.read(buf.data(), buf.size());

  ChnHeader *header = reinterpret_cast<ChnHeader *>(buf.data());

  if (header->must_be_neg1 != -1)
    throw std::runtime_error("<ImporterCHN> invalid CHN file, marker != -1");

  INFO("start.date: {}", std::string(header->start_date));
  INFO("start.time: {}", std::string(header->start_time));
  INFO("start.secs: {}", std::string(header->start_seconds));

  std::vector<int32_t> spectrum(static_cast<size_t>(header->num_channels));
  file.read(reinterpret_cast<char*>(spectrum.data()), spectrum.size()*4);

  auto hist = DAQuiri::ConsumerFactory::singleton().create_type("Histogram 1D");
  if (!hist)
    throw std::runtime_error("ImporterCHN could not get a valid Histogram 1D from factory");

  auto lt_ms = static_cast<int64_t>(header->live_time * 20.0);
  hist->set_attribute(DAQuiri::Setting("live_time", std::chrono::milliseconds(lt_ms)));

  auto rt_ms = static_cast<int64_t>(header->real_time * 20.0);
  hist->set_attribute(DAQuiri::Setting("real_time", std::chrono::milliseconds(rt_ms)));

  for (size_t i=0; i < spectrum.size(); ++i)
  {
    DAQuiri::Entry new_entry;
    new_entry.first.resize(1);
    new_entry.first[0] = i + header->channel_offset;
    new_entry.second = PreciseFloat(spectrum[i]);
    entry_list.push_back(new_entry);
  }

  hist->set_attribute(DAQuiri::Setting::text("name", path.stem().string()));
  hist->set_attribute(DAQuiri::Setting::boolean("visible", true));

  hist->import(*this);

  project->add_consumer(hist);
}
