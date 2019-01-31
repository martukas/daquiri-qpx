#pragma once

#include <QWidget>
#include <core/project.h>
#include <gui/analysis/form_energy_calibration.h>
//#include <gui/analysis/form_fwhm_calibration.h>
//#include <gui/analysis/form_fit_results.h>

namespace Ui {
class FormAnalysis1D;
}

class FormAnalysis1D : public QWidget
{
  Q_OBJECT

public:
  explicit FormAnalysis1D(QWidget *parent = 0);
  ~FormAnalysis1D();

  void setSpectrum(DAQuiri::ProjectPtr newset, int64_t idx);

  void clear();

signals:
  void calibrationComplete();
  void detectorsChanged();

public slots:
  void update_spectrum();
  void update_detector_calibs();
  void save_report();
  void update_fit();

private slots:

  void detectorsUpdated() {emit detectorsChanged();}
  void update_fits();

protected:
  void closeEvent(QCloseEvent*);

private:
  Ui::FormAnalysis1D *ui;

  FormEnergyCalibration *form_energy_calibration_;
//  FormFwhmCalibration *form_fwhm_calibration_;
//  FormFitResults *form_fit_results_;

  DAQuiri::Fitter fit_data_;

  DAQuiri::Calibration new_energy_calibration_;
  DAQuiri::Calibration new_fwhm_calibration_;

  //from parent
  QString data_directory_;
  DAQuiri::ProjectPtr spectra_;
  int64_t current_spectrum_;

  DAQuiri::Detector detector_;


  void loadSettings();
  void saveSettings();
};
