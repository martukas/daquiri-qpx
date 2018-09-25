#pragma once

#include <core/project.h>

namespace DAQuiri
{

class Importer
{
 public:
  Importer() {}
  virtual ~Importer() {}

  virtual std::string ext() const = 0;
  virtual std::string description() const = 0;
  virtual bool validate(const boost::filesystem::path& path) const = 0;
  virtual void import(const boost::filesystem::path& path, ProjectPtr project) = 0;

  EntryList_t entry_list;
};

using ImporterPtr = std::shared_ptr<Importer>;

}