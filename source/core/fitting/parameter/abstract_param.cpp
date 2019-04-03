#include <core/fitting/parameter/abstract_param.h>
#include <core/util/more_math.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

void AbstractParam::update_index(int32_t& idx)
{
  if (idx < 0)
    throw std::runtime_error("Value cannot save negative variable index");

  if (to_fit)
    index_ = idx++;
  else
    reset_index();
}

void AbstractParam::reset_index()
{
  index_ = InvalidIndex;
}

int32_t AbstractParam::index() const
{
  return index_;
}

bool AbstractParam::valid_index() const
{
  return index_ > InvalidIndex;
}

double AbstractParam::x() const
{
  return x_;
}

void AbstractParam::x(double new_x)
{
  x_ = new_x;
}

double AbstractParam::val() const
{
  return this->val_at(x_);
}

double AbstractParam::grad() const
{
  return this->grad_at(x_);
}

double AbstractParam::val_from(const Eigen::VectorXd& fit) const
{
  // \todo access without range checking once we have tests
  if (valid_index())
    return this->val_at(fit(static_cast<size_t>(index_)));
  return val();
}

double AbstractParam::grad_from(const Eigen::VectorXd& fit) const
{
  // \todo access without range checking once we have tests
  if (valid_index())
    return this->grad_at(fit(static_cast<size_t>(index_)));
  return grad();
}

double AbstractParam::uncert() const
{
  return val_uncert_;
}

void AbstractParam::uncert(double new_uncert)
{
  val_uncert_ = new_uncert;
}

void AbstractParam::put(Eigen::VectorXd& fit) const
{
  if (valid_index())
    fit[index_] = x();
}

void AbstractParam::get(const Eigen::VectorXd& fit)
{
  if (valid_index())
    x(fit[index_]);
}

void AbstractParam::get_uncert(const Eigen::VectorXd& diagonals, double chisq_norm)
{
  if (valid_index())
    uncert(std::sqrt(std::abs(diagonals[index_] * square(this->grad()) * chisq_norm)));
}

std::string AbstractParam::to_string() const
{
  auto val_part = fmt::format("{}\u00B1{}", val(), val_uncert_);
  auto x_part = fmt::format("(x={})", x_);
  auto i_part = fmt::format("{}[{}]",
      (to_fit ? "F" : ""),
      (valid_index() ? std::to_string(index_) : "-"));
  return fmt::format("{:>20} {:<17} {:>8}",
                     val_part, x_part, i_part);
}

void to_json(nlohmann::json& j, const AbstractParam& s)
{
  j["x"] = s.x_;
  j["to_fit"] = s.to_fit;
  j["uncert_value"] = s.val_uncert_;
}

void from_json(const nlohmann::json& j, AbstractParam& s)
{
  s.x_ = j["x"];
  s.to_fit = j["to_fit"];
  s.val_uncert_ = j["uncert_value"];
}


}
