#pragma once

#include <core/fitting/hypermet/Peak.h>
#include <core/fitting/hypermet/PolyBackground.h>
#include <core/fitting/weighted_data.h>

#include <core/fitting/BFGS/Fittable.h>

#include <set>

namespace DAQuiri
{

//enum class Type : uint16_t
//{
//  Normal = 0,
//  Annihilation = 1,
//  Boron = 2,
//  GeTriangle = 3,
//  IntPeakShape = 4
//};
//Type region_type{Type::Normal};


class Region : public Fittable
{
 public:
  PolyBackground background;
  Peak default_peak_;
  std::map<double, Peak> peaks_;
  SUM4Edge LB_, RB_;

 public:
  Region();
  Region(const WeightedData& data, uint16_t background_samples);
  void replace_data(const WeightedData& data, const SUM4Edge& lb, const SUM4Edge& rb);
  void replace_data(const WeightedData& data, uint16_t left_samples, uint16_t right_samples);
  void replace_data(const WeightedData& data);
  void adjust_LB(const SUM4Edge& lb);
  void adjust_RB(const SUM4Edge& rb);

  double left() const;
  double right() const;
  bool empty() const;
  bool dirty() const;

  bool add_peak(double l, double r, double amp_hint = 10);
  bool adjust_sum4(double peakID, double left, double right);
  bool auto_sum4();
  bool replace_hypermet(double peakID, Peak hyp);
  bool remove_peak(double peakID);
  bool remove_peaks(const std::set<double>& ids);

  double chi_sq_normalized() const;

  // Fitting related
  void map_fit();
  void save_fit_uncerts(const FitResult& result);
  // Fittable implementation
  Eigen::VectorXd variables() const override;
  double degrees_of_freedom() const  override;
  double chi_sq(const Eigen::VectorXd& fit) const override;
  double grad_chi_sq(const Eigen::VectorXd& fit,
                     Eigen::VectorXd& gradients,
                     Eigen::VectorXd& chan_gradients) const override;


  std::string to_string(std::string prepend = "") const;
  friend void to_json(nlohmann::json& j, const Region& s);
  friend void from_json(const nlohmann::json& j, Region& s);

 private:
  int32_t var_count_{0};
  WeightedData data_;
  bool dirty_{false};

  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  size_t fit_var_count() const;
  void save_fit(const Eigen::VectorXd& variables);
  double chi_sq() const;
  double grad_chi_sq(Eigen::VectorXd& gradients,
                     Eigen::VectorXd& chan_gradients) const;
  void init_background();
  void cull_peaks();
  void reindex_peaks();
};

}
