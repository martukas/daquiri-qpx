#pragma once

#include <QWidget>
#include <gui/widgets/SettingDelegate.h>
#include <QItemSelection>
#include <QSettings>
#include <gui/fitter/thread_fitter.h>
#include <gui/fitter/qp_fitter.h>

//#include <QMediaPlayer>

namespace Ui {
class FormFitter;
}

class FormFitter : public QWidget
{
  Q_OBJECT

  enum class RefitPolicy { kAlways, kAsk, kNever };

 public:
  explicit FormFitter(QWidget *parent = 0);
  ~FormFitter();

  void clear();

  void setFit(DAQuiri::Fitter *fit);
  void update_spectrum();

  bool busy() { return busy_; }

  void clearSelection();
  std::set<double> get_selected_peaks();

  void make_range(double energy);

  void perform_fit();

  void loadSettings(QSettings &settings_);
  void saveSettings(QSettings &settings_);

public slots:
  void set_selected_peaks(std::set<double> selected_peaks);
  void updateData();

signals:

  void peak_selection_changed(std::set<double> selected_peaks);
  void data_changed();
  void fitting_done();
  void fitter_busy(bool);

  void range_selection_changed(double l, double r);

private slots:

  void selection_changed();
  void update_range_selection(double l, double r);

  void fit_updated(DAQuiri::Fitter);
  void fitting_complete();
  void dirty(double region_id);

  void add_peak(double l, double r);
  void delete_selected_peaks();
  void adjust_sum4(double peak_id, double l, double r);
  void adjust_background_L(double roi_id, double l, double r);
  void adjust_background_R(double roi_id, double l, double r);
  void peak_info(double peak);

  void rollback_ROI(double);
  void roi_settings(double region);
  void refit_ROI(double);
  void delete_ROI(double);

  void toggle_push(bool busy);

  void on_pushFindPeaks_clicked();
  void on_pushStopFitter_clicked();
  void on_pushSettings_clicked();

  void on_pushClearAll_clicked();

private:
  Ui::FormFitter *ui;

  DAQuiri::Fitter *fit_;

  bool busy_ {false};

  ThreadFitter thread_fitter_;
//  QMediaPlayer *player;

  void merge_regions(std::set<double> rois);

};
