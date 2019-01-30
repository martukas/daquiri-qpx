#include <gui/fitter/thread_fitter.h>
#include <core/util/timer.h>
#include <core/util/custom_logger.h>

ThreadFitter::ThreadFitter(QObject *parent) :
  QThread(parent),
  terminating_(false),
  running_(false)
{
  action_ = kIdle;
  optimizer_ = DAQuiri::OptimizerFactory::getInstance().create_any();
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
  interruptor_.store(true);
}

void ThreadFitter::run() {

  while (!terminating_.load()) {
    if (action_ != kIdle) {
      running_.store(true);
      interruptor_.store(false);
    }

    if (action_ == kFit) {
      int current = 1;
      Timer total_timer(true);
      std::shared_ptr<Timer> timer(new Timer(true));
      for (auto &q : fitter_.regions())
      {
        if (optimizer_)
          fitter_.find_and_fit(q.first, optimizer_, interruptor_);
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
      if (optimizer_ && fitter_.refit_region(target_, optimizer_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAddPeak) {
      if (optimizer_)
        fitter_.add_peak(LL, RR, optimizer_, interruptor_);
      emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAdjustLB) {
      if (optimizer_ && fitter_.adj_LB(target_, LL, RR, optimizer_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kAdjustRB) {
      if (optimizer_ && fitter_.adj_RB(target_, LL, RR, optimizer_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kOverrideSettingsROI) {
      if (fitter_.override_ROI_settings(target_, settings_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kMergeRegions) {
      if (optimizer_ && fitter_.merge_regions(LL, RR, optimizer_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else if (action_ == kRemovePeaks) {
      if (optimizer_ && fitter_.remove_peaks(chosen_peaks_, optimizer_, interruptor_))
        emit fit_updated(fitter_);
      emit fitting_done();
      action_ = kIdle;
    } else {
      QThread::sleep(2);
    }
    running_.store(false);
  }
}


