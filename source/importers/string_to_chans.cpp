#include "string_to_chans.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

DAQuiri::EntryList_t string_to_chans(std::istream& data_stream)
{
  DAQuiri::EntryList_t ret;

  int i = 0;
  std::string numero;
  while (data_stream.rdbuf()->in_avail())
  {
    data_stream >> numero;
    DAQuiri::Entry new_entry;
    new_entry.first.resize(1);
    new_entry.first[0] = i;
    PreciseFloat nr{0};
    try { nr = std::stold(numero); }
    catch (...) {}
    new_entry.second = nr;
    ret.push_back(new_entry);
    i++;
  }

  return ret;
}

DAQuiri::EntryList_t string_to_chans_zero_suppressed(std::istream& data_stream)
{
  DAQuiri::EntryList_t ret;

  int i = 0;

  std::string numero, numero_z;
  while (data_stream.rdbuf()->in_avail())
  {
    data_stream >> numero;
    if (numero == "0")
    {
      data_stream >> numero_z;
      i += boost::lexical_cast<uint16_t>(numero_z);
    }
    else
    {
      DAQuiri::Entry new_entry;
      new_entry.first.resize(1);
      new_entry.first[0] = i;
      PreciseFloat nr{0};
      try { nr = std::stold(numero); }
      catch (...) {}
      new_entry.second = nr;
      ret.push_back(new_entry);
      i++;
    }
  }

  return ret;
}
