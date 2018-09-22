#include <core/importer_factory.h>
#include <core/util/custom_logger.h>

namespace DAQuiri {

void ImporterFactory::register_importer(const std::string& ext,
                                        const std::string& descr,
                                        std::function<Importer*(void)> constructor)
{
  if (ext.empty())
    WARN("<ImporterFactory> failed to register: no extension");
  else if (descr.empty())
    WARN("<ImporterFactory> failed to register: no description");
  else
  {
    constructors_.emplace(ext, constructor);
    DBG("<ImporterFactory> registered [{}]:'{}'", ext, descr);
  }
}

std::unordered_set<std::string> ImporterFactory::extensions() const
{
  std::unordered_set<std::string> ret;
  for (auto &q : constructors_)
    ret.insert(q.first);
  return ret;
}

std::vector<std::string> ImporterFactory::descriptions() const
{
  std::vector<std::string> ret;
  for (auto &q : constructors_)
  {
    // \ todo add ext wildcard
    ret.push_back(ImporterPtr(q.second())->description());
  }
  return ret;
}

std::vector<ImporterPtr> ImporterFactory::attempt_import(const boost::filesystem::path& path) const
{
  std::vector<ImporterPtr> ret;
  for (auto &q : constructors_)
  {
    auto i = ImporterPtr(q.second());
    if (i->validate(path))
      ret.push_back(i);
  }
  return ret;
}

}
