#pragma once

#include <QThread>
#include <QMutex>
#include <core/gamma/fitter.h>

enum FitterAction {kFit, kStop, kIdle, kRefit,
                  kMergeRegions};


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
  void merge_regions(double L, double R);

signals:
  void fit_updated(DAQuiri::Fitter data);
  void fitting_done();

protected:
  void run();

private:
  DAQuiri::Fitter fitter_;

  std::shared_ptr<DAQuiri::AbstractOptimizer> optimizer_;

  QMutex mutex_;
  FitterAction action_;

  double LL, RR;
  double target_;
  DAQuiri::FitSettings settings_;
  std::set<double> chosen_peaks_;

  std::atomic<bool> running_;
  std::atomic<bool> terminating_;
};
