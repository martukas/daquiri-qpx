#pragma once

#include <core/fitting/hypermet/Value.h>

namespace Hypermet
{

struct PolyBackground
{
  PolyBackground() = default;

  double bin_offset {0};
  ValueGam background_base_;
  bool slope_enabled_{true};
  ValueBkg background_slope_;
  bool curve_enabled_{true};
  ValueBkg background_curve_;


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

  std::string to_string() const;
};

void to_json(nlohmann::json& j, const PolyBackground& s);
void from_json(const nlohmann::json& j, PolyBackground& s);


}
