#include <core/gamma/nuclides/transition.h>

namespace DAQuiri {

void to_json(nlohmann::json& j, const Radiation& s)
{
  j["energy"] = s.energy;
  j["abundance"] = s.abundance;
}

void from_json(const nlohmann::json& j, Radiation& s)
{
  s.energy = j["energy"];
  s.abundance = j["abundance"];
}

}
