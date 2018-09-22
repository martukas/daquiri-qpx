#pragma once

#include <core/project.h>
#include <core/util/unique_mangle.h>
#include <unordered_set>

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
  virtual void import(const boost::filesystem::path& path, ProjectPtr project) const = 0;
};

using ImporterPtr = std::shared_ptr<Importer>;

class ImporterFactory
{
 public:
  static ImporterFactory& singleton()
  {
    static ImporterFactory singleton_instance;
    return singleton_instance;
  }

  void register_importer(const std::string& ext, const std::string& descr,
                         std::function<Importer*(void)> constructor);

  std::unordered_set<std::string> extensions() const;
  std::vector<std::string> descriptions() const;
  std::vector<ImporterPtr> attempt_import(const boost::filesystem::path& path) const;

 private:
  std::multimap<std::string, std::function<Importer*(void)>> constructors_;

  //singleton assurance
  ImporterFactory() {}
  ImporterFactory(ImporterFactory const&);
  void operator=(ImporterFactory const&);
};

template<class T>
class ImporterRegistrar
{
 public:
  ImporterRegistrar()
  {
    ImporterFactory::singleton().register_importer(T().ext(), T().description(),
                                                   [](void) -> Importer* { return new T(); });
  }
};

#define DAQUIRI_REGISTER_IMPORTER(T) static DAQuiri::ImporterRegistrar< T >\
  UNIQUE_MANGLE(MangledDAQuiriImporterReg) ;

}
