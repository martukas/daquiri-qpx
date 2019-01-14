#include <gui/analysis/isotope.h>
#include <core/util/string_extensions.h>

namespace RadTypes {

void AbstractRadiation::to_xml(pugi::xml_node &node) const {
  pugi::xml_node child = node.append_child();
  child.set_name(this->xml_element_name().c_str());

  child.append_child("energy");
  child.last_child().append_child(pugi::node_pcdata).set_value(dbl2str(energy).c_str());
  child.append_child("abundance");
  child.last_child().append_child(pugi::node_pcdata).set_value(dbl2str(abundance).c_str());
}

void AbstractRadiation::from_xml(const pugi::xml_node &node) {
  std::string str;

  str = std::string(node.child_value("energy"));
  if (!str.empty())
    energy = std::stod(trim_copy(str));

  str = std::string(node.child_value("abundance"));
  if (!str.empty())
    abundance = std::stod(trim_copy(str));
}


void Isotope::to_xml(pugi::xml_node &node) const {
  pugi::xml_node child = node.append_child();
  child.set_name(this->xml_element_name().c_str());

  child.append_child("name").append_child(pugi::node_pcdata).set_value(name.c_str());

  if (gamma_constant.size()) {
    child.append_child("gammaConstant");
    child.last_child().append_child(pugi::node_pcdata).set_value(gamma_constant.c_str());
  }

  child.append_child("halfLife");
  child.last_child().append_child(pugi::node_pcdata).set_value(dbl2str(half_life).c_str());

  if (!gammas.empty())
    gammas.to_xml(child);

  if (beta.energy > 0)
    beta.to_xml(child);
}

void Isotope::from_xml(const pugi::xml_node &node) {
  name = std::string(node.child_value("name"));
  gamma_constant = std::string(node.child_value("gammaConstant"));
  std::string hl(node.child_value("halfLife"));
  if (!hl.empty())
    half_life = std::stod(hl);
  gammas.from_xml(node.child(gammas.xml_element_name().c_str()));
  beta.from_xml(node.child(beta.xml_element_name().c_str()));
}

}
