#include <gui/fitter/thread_fitter.h>
#include <core/util/timer.h>
#include <core/util/custom_logger.h>

#include <core/fitting/optimizers/optlib_adapter.h>

ThreadFitter::ThreadFitter(QObject *parent) :
  QThread(parent),
  terminating_(false),
  running_(false)
{
  auto opt = std::make_shared<DAQuiri::OptlibOptimizer>();
//  opt->verbose = true;
  opt->gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::DefaultToFinite;
  opt->maximum_iterations = 200;
  opt->maximum_perturbations = 10;
//  opt->perform_sanity_checks = true;

  optimizer_ = opt;

  action_ = kIdle;
  start(HighPriority);
}

void ThreadFitter::terminate() {
  terminating_.store(true);
  wait();
}

void ThreadFitter::begin() {
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::set_data(const DAQuiri::Fitter &data)
{
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  fitter_ = data;
  action_ = kIdle;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::fit_peaks() {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kFit;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::add_peak(double L, double R) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kAddPeak;
  LL = L;
  RR = R;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::refit(double target_ROI) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kRefit;
  target_ = target_ROI;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::override_ROI_settings(double regionID, DAQuiri::FitSettings fs)
 {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kOverrideSettingsROI;
  target_ = regionID;
  settings_ = fs;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::adjust_LB(double target_ROI, double L, double R) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kAdjustLB;
  target_ = target_ROI;
  LL = L;
  RR = R;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::adjust_RB(double target_ROI, double L, double R) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kAdjustRB;
  target_ = target_ROI;
  LL = L;
  RR = R;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::merge_regions(double L, double R) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kMergeRegions;
  LL = L;
  RR = R;
  if (!isRunning())
    start(HighPriority);
}


void ThreadFitter::remove_peaks(std::set<double> chosen_peaks) {
  if (running_.load()) {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kRemovePeaks;
  chosen_peaks_ = chosen_peaks;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::stop_work() {
  QMutexLocker locker(&mutex_);
  action_ = kStop; //not thread safe
  optimizer_->cancel.store(true);
}

void ThreadFitter::run() {

  while (!terminating_.load()) {
    if (action_ != kIdle) {
      running_.store(true);
      optimizer_->cancel.store(false);
    }

    if (action_ == kFit) {
      int current = 1;
      Timer total_timer(true);
      std::shared_ptr<Timer> timer(new Timer(true));
      for (auto &q : fitter_.regions())
      {
        fitter_.find_and_fit(q.first, optimizer_.get());
        current++;
        if (timer->s() > 2) {
          emit fit_updated(fitter_);
          timer = std::shared_ptr<Timer>(new Timer(true));
          DBG("<Fitter> {} of {} regions completed",
              current, fitter_.region_count());
        }
        if ((action_ == kStop) || terminating_.load())
          break;
      }
      DBG("<Fitter> Fitting spectrum was on average {} s/peak",
          total_timer.s() / double(fitter_.peaks().size()));
      emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kRefit) {
      if (fitter_.refit_region(target_, optimizer_.get()))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAddPeak) {
      fitter_.add_peak(LL, RR, optimizer_.get());
      emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAdjustLB) {
      if (fitter_.adj_LB(target_, LL, RR, optimizer_.get()))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAdjustRB) {
      if (fitter_.adj_RB(target_, LL, RR, optimizer_.get()))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kOverrideSettingsROI) {
      if (fitter_.override_ROI_settings(target_, settings_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kMergeRegions) {
      if (fitter_.merge_regions(LL, RR, optimizer_.get()))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kRemovePeaks) {
      if (fitter_.remove_peaks(chosen_peaks_, optimizer_.get()))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else {
      QThread::sleep(2);
    }
    running_.store(false);
  }
}


