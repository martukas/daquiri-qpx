#pragma once

#include <QThread>
#include <QMutex>
#include <core/fitting/fitter.h>

enum FitterAction {kFit, kStop, kIdle, kAddPeak, kRemovePeaks, kRefit,
                  kAdjustLB, kAdjustRB, kOverrideSettingsROI, kMergeRegions};


enum class RefitPolicy { kAlways, kAsk, kNever };

class ThreadFitter : public QThread
{
  Q_OBJECT
public:
  explicit ThreadFitter(QObject *parent = 0);
  void terminate();

  void begin();
  void set_data(const DAQuiri::Fitter &data);

  void fit_peaks();
  void stop_work();
  void refit(double target_ROI);
  void add_peak(double L, double R);
  void merge_regions(double L, double R);
  void adjust_LB(double target_ROI, double L, double R);
  void adjust_RB(double target_ROI, double L, double R);
  void override_ROI_settings(double regionID, DAQuiri::FitSettings fs);
  void remove_peaks(std::set<double> chosen_peaks);

signals:
  void fit_updated(DAQuiri::Fitter data);
  void dirty(double region_id);
  void fitting_done();

protected:
  void run();

private:
  DAQuiri::Fitter fitter_;

  std::shared_ptr<DAQuiri::AbstractOptimizer> optimizer_;

  QMutex mutex_;
  FitterAction action_;

  RefitPolicy refit_policy_{RefitPolicy::kAsk};

  double LL, RR;
  double target_;
  DAQuiri::Peak hypermet_;
  DAQuiri::FitSettings settings_;
  std::set<double> chosen_peaks_;

  std::atomic<bool> running_;
  std::atomic<bool> terminating_;

  void conditional_refit(double region_id);
};
