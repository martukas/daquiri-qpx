#pragma once

#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

class FitEvaluation
{
 public:
  FitEvaluation() = default;
  FitEvaluation(const WeightedData& data);

  void clear();
  void reset();

  bool empty() const;

  bool cloneRange(const FitEvaluation& other, double l, double r);
  void update_fit(const std::vector<double>& y_fit,
                  const std::vector<double>& y_background);

  //DATA

  std::vector<double> x_;
  std::vector<double> y_;

  WeightedData weighted_data;

  std::vector<double> y_fit_, y_background_, y_resid_, y_resid_on_background_;

 private:

  void setNewData(const WeightedData& d);
};

}
