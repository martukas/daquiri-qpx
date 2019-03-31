#include <core/fitting/fit_evaluation.h>
#include <cmath>
#include <range/v3/all.hpp>

namespace DAQuiri
{

FitEvaluation::FitEvaluation(const WeightedData& data)
{
  setNewData(data);
}

void FitEvaluation::cloneRange(const FitEvaluation& other, double l, double r)
{
  setNewData(other.weighted_data.subset(l, r));
}

void FitEvaluation::setNewData(const WeightedData& d)
{
  clear();
  weighted_data = d;
  x_ = d.chan;
  y_ = d.count;
  reset();
}

void FitEvaluation::clear()
{
  x_.clear();
  y_.clear();
  y_fit_.clear();
  y_background_.clear();
  y_resid_.clear();
  y_resid_on_background_.clear();

  weighted_data.clear();
}

void FitEvaluation::reset()
{
  y_resid_on_background_ = y_resid_ = y_;
  y_fit_.assign(x_.size(), 0.0);
  y_background_.assign(x_.size(), 0.0);
  y_resid_weighted_.assign(x_.size(), std::numeric_limits<double>::quiet_NaN());
}

bool FitEvaluation::empty() const
{
  return x_.empty();
}

void FitEvaluation::update_fit(const std::vector<double>& y_fit,
                               const std::vector<double>& y_background)
{
  if ((y_.size() != y_background.size())
      || (y_fit.size() != y_background.size())
      || (y_fit.empty()))
    return;

  for (size_t i = 0; i < y_fit.size(); ++i)
  {
    y_fit_[i] = y_fit[i];
    y_background_[i] = y_background[i];
    double resid = y_[i] - y_fit[i];
    y_resid_[i] = resid;
    y_resid_on_background_[i] = y_background[i] + resid;
    // \todo chose which type of weight you want to use for this?
    y_resid_weighted_[i] = y_resid_[i] / std::sqrt(y_[i] + 1.0);
  }
}

void FitEvaluation::merge_fit(const FitEvaluation& other)
{
  if (other.x_.empty() ||
      (other.x_.front() < x_.front()) ||
      (other.x_.back() > x_.back()))
    return;

  // find starting point
  size_t i = 0;
  while (x_[i] < other.x_.front())
    i++;

  for (size_t j = 0; j < other.x_.size(); ++j)
  {
    y_fit_[i + j] = other.y_fit_[j];
    y_background_[i + j] = other.y_background_[j];
    y_resid_[i + j] = other.y_resid_[j];
    y_resid_on_background_[i + j] = other.y_resid_on_background_[j];
    y_resid_weighted_[i + j] = other.y_resid_weighted_[j];
  }
}

}
