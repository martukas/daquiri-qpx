#pragma once

#include <QWidget>
#include <QItemSelection>
#include <core/fitting/fitter.h>

#include <gui/widgets/SettingDelegate.h>
#include <gui/analysis/widget_plot_calib.h>

namespace Ui
{
class FormFwhmCalibration;
}

class FormFwhmCalibration : public QWidget
{
 Q_OBJECT

 public:
  explicit FormFwhmCalibration(DAQuiri::Detector&, DAQuiri::Fitter&, QWidget* parent = 0);
  ~FormFwhmCalibration();

  void newSpectrum();
  DAQuiri::Calibration get_new_calibration() { return new_calibration_; }

  void clear();
  bool save_close();

 public slots:
  void update_selection(std::set<double> selected_peaks);
  void update_data();

 signals:
  void selection_changed(std::set<double> selected_peaks);
  void detectorsChanged();
  void update_detector();
  void new_fit();

 private slots:
  void selection_changed_in_table();
  void selection_changed_in_plot();

  void toggle_push();
  void detectorsUpdated() { emit detectorsChanged(); }

  void on_pushApplyCalib_clicked();
  void on_pushFit_clicked();
  void on_pushFromDB_clicked();
  void on_pushDetDB_clicked();

  void on_doubleMaxFitErr_valueChanged(double arg1);
  void on_doubleMaxWidthErr_valueChanged(double arg1);

 private:
  Ui::FormFwhmCalibration* ui;

  //from parent
  QString data_directory_;

  DAQuiri::Detector& detector_;
  DAQuiri::Fitter& fit_data_;
  std::set<double> selected_peaks_;

  DAQuiri::Calibration new_calibration_;
  QPlot::Appearance style_pts, style_relevant, style_fit;

  void loadSettings();
  void saveSettings();
  void fit_calibration();

  void replot_calib();
  void rebuild_table();
  void select_in_table();
  void select_in_plot();

  void add_peak_to_table(const DAQuiri::Peak&, int, bool);
};

