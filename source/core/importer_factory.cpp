#include <core/importer_factory.h>
#include <core/util/string_extensions.h>

#include <core/util/logger.h>

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
  std::string all_files;
  for (auto &q : constructors_)
  {
    if (!all_files.empty())
      all_files += " ";
    all_files += "*." + q.first;
  }

  std::vector<std::string> ret;
  if (!all_files.empty())
    ret.push_back("All available formats (" + all_files + ")");
  for (auto &q : constructors_)
  {
    ImporterPtr imp(q.second());
    ret.push_back(imp->description() + " (*." + imp->ext() + ")" );
  }

  return ret;
}

std::vector<ImporterPtr> ImporterFactory::attempt_import(const boost::filesystem::path& path) const
{
  std::vector<ImporterPtr> ret;
  for (auto &q : constructors_)
  {
    auto i = ImporterPtr(q.second());
    if (!iequals(("." + i->ext()), path.extension().string()))
      continue;
    auto path2 = path.stem();
    if (i->validate(path2))
      ret.push_back(i);
  }
  return ret;
}

void ImporterFactory::clear()
{
  constructors_.clear();
}

}
