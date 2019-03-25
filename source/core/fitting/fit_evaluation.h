#pragma once

#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

/// \class FitEvaluation fit_evaluation.h <core/fitting/fit_evaluation.h>
/// \brief stores histogram data for a region as well as an evaluation of the fit,
///         background and residuals. Provides interface for extracting subsets.
class FitEvaluation
{
 public:
  FitEvaluation() = default;

  /// \brief constructs fit evaluation from weighted data without fit.
  /// \param data weighted region data
  /// \post y_fit_ and y_background_ are zeroed out, y_resid_ and y_resid_on_background_
  ///         are made equal to y_
  FitEvaluation(const WeightedData& data);

  /// \brief clears all data
  void clear();

  /// \brief resets the evaluation assuming valid data but no fit
  /// \post y_fit_ and y_background_ are zeroed out, y_resid_ and y_resid_on_background_
  ///         are made equal to y_
  void reset();

  /// \returns true if region contains no data
  bool empty() const;

  // \todo replace this
  void cloneRange(const FitEvaluation& other, double l, double r);

  /// \brief updates evaluation with new background and total fit evaluation
  /// \param y_fit total fit evaluation for the region
  /// \param y_background background fit evaluation for the region
  /// \post y_resid_ is made equal to (count - fit) for each channel
  /// \post y_resid_on_background_ is made equal to (background + residual) for each channel
  void update_fit(const std::vector<double>& y_fit,
                  const std::vector<double>& y_background);

  void merge_fit(const FitEvaluation& other);

  std::vector<double> x_; /// < channel
  std::vector<double> y_; /// < counts

  WeightedData weighted_data; /// < weighted data for region

  std::vector<double> y_fit_;                  /// < evaluation of fit
  std::vector<double> y_background_;           /// < background portion of fit
  std::vector<double> y_resid_;                /// < y - fit
  std::vector<double> y_resid_on_background_;  /// < background + residual

  std::vector<double> y_resid_weighted_;       /// < (y - fit) / weight

 private:

  void setNewData(const WeightedData& d);
};

}
