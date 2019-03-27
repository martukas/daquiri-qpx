#pragma once

#include <QWidget>
#include <core/project.h>
#include <gui/analysis/form_energy_calibration.h>
#include <gui/widgets/SelectorWidget.h>

namespace Ui
{
class FormEfficiencyCalibration;
}

class FormEfficiencyCalibration : public QWidget
{
 Q_OBJECT

 public:
  explicit FormEfficiencyCalibration(QWidget* parent = 0);
  ~FormEfficiencyCalibration();

  void setDetector(DAQuiri::ProjectPtr newset, QString detector);
//  void clear();

 signals:
  void calibrationComplete();
  void detectorsChanged();

 public slots:
  void update_spectrum();
  void update_detector_calibs();
  void update_selection(std::set<double> selected_peaks);

 private slots:
  void update_data();

  void update_spectra();
  void detectorsUpdated() { emit detectorsChanged(); }

  void setSpectrum(int64_t idx);
  void spectrumLooksChanged(SelectorItem);
  void spectrumDetails(SelectorItem);

  void selection_changed_in_table();
  void selection_changed_in_calib_plot();

  void replot_calib();
  void rebuild_table(bool);
  void toggle_push();
  void isotope_chosen();
  void on_pushApplyCalib_clicked();

  void on_pushFit_clicked();
  void on_pushDetDB_clicked();

  void on_pushCullPeaks_clicked();

  void on_doubleEpsilonE_editingFinished();

  void on_doubleScaleFactor_editingFinished();

  void on_doubleScaleFactor_valueChanged(double arg1);

  void on_pushFit_2_clicked();

  void on_pushFitEffit_clicked();

 protected:
  void closeEvent(QCloseEvent*);

 private:
  Ui::FormEfficiencyCalibration* ui;

  std::map<int64_t, DAQuiri::Fitter> peak_sets_;
  std::set<double> selected_peaks_;
  std::set<double> falgged_peaks_;

  DAQuiri::Fitter fit_data_;

  QString mca_load_formats_;  //valid mca file formats that can be opened

  QString data_directory_;

  DAQuiri::ProjectPtr project_;
  int64_t current_spectrum_;
  std::string current_detector_;
//  DAQuiri::OptimizerPtr optimizer_;

  void loadSettings();
  void saveSettings();

  DAQuiri::Calibration new_calibration_;
  QPlot::Appearance style_fit, style_pts;
  void add_peak_to_table(const DAQuiri::Peak&, int, QColor);
  void select_in_table();
};
