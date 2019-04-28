#include <core/importer_factory.h>

class ImporterN42 : public DAQuiri::Importer
{
 public:
  ImporterN42() : DAQuiri::Importer() {}
  ~ImporterN42() {}

  std::string ext() const override { return "n42"; }
  std::string description() const override { return "The glorious N42 format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

