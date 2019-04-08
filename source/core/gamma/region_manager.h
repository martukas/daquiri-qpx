#pragma once

#include <core/gamma/hypermet/gamma_region.h>
#include <core/gamma/fit_evaluation.h>
#include <core/gamma/fit_settings.h>

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri
{

struct FitDescription
{
  std::string description;
  size_t peaknum{0};
  double chi_sq_norm{0};
  double sum4aggregate{0};
};

class Fit
{
 public:
  Fit(const Region& r, std::string descr);

  FitDescription description;
  Region region;
};

class RegionManager
{
 public:
  RegionManager() = default;
  RegionManager(const Region& initial_region);
  RegionManager(const nlohmann::json& j, const WeightedData& super_region);

  double id() const;

  //access peaks
  size_t peak_count() const;
  bool contains(double peakID) const;
  Peak peak(double peakID) const;

  // current state
  const Region& region() const;
  void modify_region(const Region& new_region, std::string message = "");
  //access history
  size_t current_fit() const;
  std::vector<FitDescription> history() const;
  bool rollback(size_t i);

  //access other
  FitEvaluation eval() const;

  //manupulation, may invoke optimizer
  bool refit(AbstractOptimizer* optimizer);

  nlohmann::json to_json(const FitEvaluation& parent_finder) const;

 private:
//  FitSettings settings_;

  //history
  std::vector<Fit> fits_;
  size_t current_fit_{0};

  //intrinsic, these are saved
  Region region_;

//  bool add_from_resid(AbstractOptimizer* optimizer);
//  void iterative_fit(AbstractOptimizer* optimizer);

  void save_current_fit(std::string description);
};

}
