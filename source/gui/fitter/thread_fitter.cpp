#include <gui/fitter/thread_fitter.h>
#include <core/util/timer.h>
#include <core/util/logger.h>

#include <core/fitting/optimizers/optlib_adapter.h>

ThreadFitter::ThreadFitter(QObject* parent) :
    QThread(parent), terminating_(false), running_(false)
{
  auto opt = std::make_shared<DAQuiri::OptlibOptimizer>();
//  opt->verbose = true;
  opt->maximum_iterations = 1000;
  opt->gradient_selection =
      DAQuiri::OptlibOptimizer::GradientSelection::AnalyticalAlways;
//  opt->epsilon = 1e-10;
//  opt->tolerance = 1e-4;
  opt->use_epsilon_check = false;
  opt->min_g_norm = 1e-7;

  opt->perform_sanity_checks = false;
  opt->maximum_perturbations = 0;

  optimizer_ = opt;

  action_ = kIdle;
  start(HighPriority);
}

void ThreadFitter::terminate()
{
  terminating_.store(true);
  wait();
}

void ThreadFitter::begin()
{
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::set_data(const DAQuiri::Fitter& data)
{
  if (running_.load())
  {
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

void ThreadFitter::fit_peaks()
{
  if (running_.load())
  {
    WARN("Fitter busy");
    return;
  }
  QMutexLocker locker(&mutex_);
  terminating_.store(false);
  action_ = kFit;
  if (!isRunning())
    start(HighPriority);
}

void ThreadFitter::refit(double target_ROI)
{
  if (running_.load())
  {
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

void ThreadFitter::stop_work()
{
  QMutexLocker locker(&mutex_);
  action_ = kStop; //not thread safe
  optimizer_->cancel.store(true);
}

void ThreadFitter::run()
{

  while (!terminating_.load())
  {
    if (action_ != kIdle)
    {
      running_.store(true);
      optimizer_->cancel.store(false);
    }

    if (action_ == kFit)
    {
      int current = 1;
      Timer total_timer(true);
      std::shared_ptr<Timer> timer(new Timer(true));
      for (auto& q : fitter_.regions())
      {
        fitter_.find_and_fit(q.first, optimizer_.get());
        current++;
        if (timer->s() > 2)
        {
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
    }
    else if (action_ == kRefit)
    {
      auto r = fitter_.region(target_);
      if (r.region().peaks_.empty())
      {
        if (fitter_.find_and_fit(target_, optimizer_.get()))
            emit fit_updated(fitter_);
      }
      else
      {
        if (fitter_.refit_region(target_, optimizer_.get()))
            emit fit_updated(fitter_);
      }
      emit fitting_done();
      action_ = kIdle;
    }
    else
    {
      QThread::sleep(2);
    }
    running_.store(false);
  }
}


