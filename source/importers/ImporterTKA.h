#include <core/importer_factory.h>
#include <core/project.h>

class ImporterTKA : public DAQuiri::Importer
{
 public:
  ImporterTKA() : DAQuiri::Importer() {}
  ~ImporterTKA() {}

  std::string ext() const override { return "tka"; }
  std::string description() const override { return "The glorious TKA format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

