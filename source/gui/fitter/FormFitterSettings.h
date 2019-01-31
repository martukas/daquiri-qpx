#pragma once

#include <QDialog>
#include <QCloseEvent>
#include <core/fitting/fit_settings.h>

namespace Ui {
class FormFitterSettings;
}

class FormFitterSettings : public QDialog
{
  Q_OBJECT

public:
  explicit FormFitterSettings(DAQuiri::FitSettings &fs, QWidget *parent = 0);
  ~FormFitterSettings();

protected:
  void closeEvent(QCloseEvent*);


private slots:
  void on_buttonBox_accepted();

  void on_buttonBox_rejected();

  void on_doubleMinRskewSlope_valueChanged(double);
  void on_doubleMaxRskewSlope_valueChanged(double);
  void on_doubleMinRskewAmp_valueChanged(double);
  void on_doubleMaxRskewAmp_valueChanged(double);
  void on_doubleMinLskewSlope_valueChanged(double);
  void on_doubleMaxLskewSlope_valueChanged(double);
  void on_doubleMinLskewAmp_valueChanged(double);
  void on_doubleMaxLskewAmp_valueChanged(double);
  void on_doubleMinTailSlope_valueChanged(double);
  void on_doubleMaxTailSlope_valueChanged(double);
  void on_doubleMinTailAmp_valueChanged(double);
  void on_doubleMaxTailAmp_valueChanged(double);
  void on_doubleMinStep_valueChanged(double);
  void on_doubleMaxStep_valueChanged(double);

  void on_checkGaussOnly_clicked();

private:
  Ui::FormFitterSettings *ui;

  DAQuiri::FitSettings &fit_settings_;

  void enforce_bounds();
};
