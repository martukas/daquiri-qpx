#include <core/importer_factory.h>

class ImporterCNF : public DAQuiri::Importer
{
 public:
  ImporterCNF() : DAQuiri::Importer() {}
  ~ImporterCNF() {}

  std::string ext() const override { return "cnf"; }
  std::string description() const override { return "The glorious CNF format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

