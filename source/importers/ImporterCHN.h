#include <core/importer_factory.h>

class ImporterCHN : public DAQuiri::Importer
{
 public:
  ImporterCHN() : DAQuiri::Importer() {}
  ~ImporterCHN() {}

  std::string ext() const override { return "chn"; }
  std::string description() const override { return "The glorious CHN format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

