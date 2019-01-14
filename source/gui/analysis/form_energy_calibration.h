#pragma once

#include <QWidget>
#include <gui/analysis/isotope.h>
#include <gui/analysis/widget_plot_calib.h>
#include <QItemSelection>
#include <core/fitting/fitter.h>

#include <gui/widgets/SettingDelegate.h>

namespace Ui {
class FormEnergyCalibration;
}

class FormEnergyCalibration : public QWidget
{
  Q_OBJECT

public:
  explicit FormEnergyCalibration(XMLableDB<DAQuiri::Detector>&, DAQuiri::Fitter&, QWidget *parent = 0);
  ~FormEnergyCalibration();

  DAQuiri::Calibration get_new_calibration() {return new_calibration_;}

  void newSpectrum();
  bool save_close();

  void clear();

public slots:
  void update_selection(std::set<double> selected_peaks);
  void update_data();

signals:
  void selection_changed(std::set<double> selected_peaks);
  void change_peaks();
  void detectorsChanged();
  void update_detector();
  void new_fit();

private slots:
  void selection_changed_in_table();
  void selection_changed_in_plot();

  void toggle_push();
  void isotope_energies_chosen();
  void on_pushApplyCalib_clicked();

  void detectorsUpdated() {emit detectorsChanged();}

  void on_pushFit_clicked();
  void on_pushFromDB_clicked();
  void on_pushDetDB_clicked();
  void on_pushPeaksToNuclide_clicked();
  void on_pushEnergiesToPeaks_clicked();

private:
  Ui::FormEnergyCalibration *ui;

  QString data_directory_;

  XMLableDB<DAQuiri::Detector> &detectors_;
  DAQuiri::Fitter &fit_data_;
  std::set<double> selected_peaks_;

  DAQuiri::Calibration new_calibration_;
  QPlot::Appearance style_fit, style_pts;

  void loadSettings();
  void saveSettings();

  void replot_calib();
  void rebuild_table();
  void select_in_table();
  void select_in_plot();

  void add_peak_to_table(const DAQuiri::Peak &, int, bool);
};

