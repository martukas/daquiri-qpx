#include <core/importer_factory.h>

class ImporterSPC : public DAQuiri::Importer
{
 public:
  ImporterSPC() : DAQuiri::Importer() {}
  ~ImporterSPC() {}

  std::string ext() const override { return "spc"; }
  std::string description() const override { return "The glorious SPC format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

