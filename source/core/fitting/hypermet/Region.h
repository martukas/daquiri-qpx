#pragma once

#include <core/fitting/hypermet/Peak.h>
#include <core/fitting/hypermet/PolyBackground.h>
#include <core/fitting/finder.h>

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

  bool dirty{false};

 public:
  Region();
  Region(const SpectrumData& data, uint16_t background_samples);
  void replace_data(const SpectrumData& data, const SUM4Edge& lb, const SUM4Edge& rb);
  void replace_data(const SpectrumData& data, uint16_t left_samples, uint16_t right_samples);
  void replace_data(const SpectrumData& data);
  void adjust_LB(const SUM4Edge& lb);
  void adjust_RB(const SUM4Edge& rb);

  double left() const;
  double right() const;

  bool adjust_sum4(double peakID, double left, double right);
  bool replace_hypermet(double peakID, Peak hyp);
  bool remove_peak(double peakID);
  bool remove_peaks(const std::set<double>& ids);

  // \todo add peak
  void reindex_peaks();

  double chi_sq_normalized() const;

  // Fitting related
  void map_fit();
  void save_fit_uncerts(const FitResult& result);
  // Fittable implementation
  std::vector<double> variables() const override;
  double degrees_of_freedom() const  override;
  double chi_sq(const std::vector<double>& fit) const override;
  double grad_chi_sq(const std::vector<double>& fit,
                     std::vector<double>& gradients) const override;

  friend void to_json(nlohmann::json& j, const Region& s);
  friend void from_json(const nlohmann::json& j, Region& s);

 private:
  int32_t var_count_{0};
  SpectrumData data_;

  //public: BoronPeak As CBoronPeak
  //public: AnnPeak As CAnnPeak

  size_t fit_var_count() const;
  void save_fit(const std::vector<double>& variables);
  double chi_sq() const;
  double grad_chi_sq(std::vector<double>& gradients) const;
  void init_background();
  void cull_peaks();
};

}
