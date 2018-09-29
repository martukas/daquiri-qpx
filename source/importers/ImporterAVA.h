#include <core/importer_factory.h>
#include <boost/lexical_cast.hpp>
#include <core/project.h>

class ImporterAVA : public DAQuiri::Importer
{
 public:
  ImporterAVA() : DAQuiri::Importer() {}
  ~ImporterAVA() {}

  std::string ext() const override { return "ava"; }
  std::string description() const override { return "The glorious AVA format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

