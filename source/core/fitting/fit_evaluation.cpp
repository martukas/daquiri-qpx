#include <core/fitting/fit_evaluation.h>

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
  for (const auto& p : d.data)
  {
    x_.push_back(p.x);
    y_.push_back(p.y);
  }
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
  }
}


}
