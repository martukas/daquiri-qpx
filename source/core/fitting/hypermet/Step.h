#pragma once

#include <core/fitting/hypermet/Value.h>

namespace DAQuiri
{

struct Step
{
  Step() = default;
  Step(Side s) : side(s) {}

  bool override{false};
  bool enabled{true};
  Value amplitude;
  Side side {Side::left};

  inline double flip(double spread) const
  {
    if (side == Side::right)
      return -spread;
    return spread;
  }

  void update_indices(int32_t& i);
  void put(Eigen::VectorXd& fit) const;
  void get(const Eigen::VectorXd& fit);
  void get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm);

  double eval_with(const PrecalcVals& pre, double ampl) const;

  double eval(const PrecalcVals& pre) const;

  double eval_at(const PrecalcVals& pre, const Eigen::VectorXd& fit) const;

  double eval_grad(const PrecalcVals& pre,
                   Eigen::VectorXd& grads,
                   size_t i_width, size_t i_pos, size_t i_amp) const;

  double eval_grad_at(const PrecalcVals& pre, const Eigen::VectorXd& fit,
                   Eigen::VectorXd& grads,
                   size_t i_width, size_t i_pos, size_t i_amp) const;

  std::string to_string() const;
};

void to_json(nlohmann::json& j, const Step& s);
void from_json(const nlohmann::json& j, Step& s);

}
