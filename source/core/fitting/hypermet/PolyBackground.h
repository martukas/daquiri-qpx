#pragma once

#include <core/fitting/hypermet/Value.h>

namespace DAQuiri
{

struct PolyBackground
{
  PolyBackground() = default;
  // \todo construct with default vals

  // \todo why are these unbounded?
  double x_offset {0};
  ValueGam base;
  bool slope_enabled{true};
  ValueBkg slope;
  bool curve_enabled{true};
  ValueBkg curve;

  void update_indices(int32_t& i);
  void put(std::vector<double>& fit) const;
  void get(const std::vector<double>& fit);
  void get_uncerts(const std::vector<double>& diagonals, double chisq_norm);

  double eval(double bin) const;
  double eval_at(double bin, const std::vector<double>& fit) const;
  double eval_grad(double bin, std::vector<double>& grads) const;
  double eval_grad_at(double bin,
      const std::vector<double>& fit,
                      std::vector<double>& grads) const;

  double eval_add(const std::vector<double>& bins, std::vector<double>& vals) const;
  std::vector<double> eval(const std::vector<double>& bins) const;

  std::string to_string(std::string prepend = "") const;
};

void to_json(nlohmann::json& j, const PolyBackground& s);
void from_json(const nlohmann::json& j, PolyBackground& s);


}
