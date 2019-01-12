#include <core/importer_factory.h>
#include <boost/lexical_cast.hpp>
#include <core/project.h>

class ImporterMCA : public DAQuiri::Importer
{
 public:
  ImporterMCA() : DAQuiri::Importer() {}
  ~ImporterMCA() {}

  std::string ext() const override { return "mca"; }
  std::string description() const override { return "The glorious MCA format"; }
  bool validate(const boost::filesystem::path& path) const override;
  void import(const boost::filesystem::path& path, DAQuiri::ProjectPtr project) override;
};

