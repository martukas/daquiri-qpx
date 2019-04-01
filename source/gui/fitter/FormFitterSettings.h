#pragma once

#include <QDialog>
#include <core/gamma/hypermet/peak.h>
#include <core/gamma/fit_settings.h>

class QCloseEvent;
class FitParameterWidget;
class TailWidget;
class UncertainDoubleWidget;

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
  void update();

  void on_buttonBox_accepted();

  void on_buttonBox_rejected();

//  void on_checkGaussOnly_clicked();

private:
  Ui::FormFitterSettings *ui;

  DAQuiri::FitSettings &fit_settings_;
  DAQuiri::FitSettings backup_settings_;

  FitParameterWidget* width_;
  UncertainDoubleWidget* fwhm_;
  UncertainDoubleWidget* fwhm_energy_;

  FitParameterWidget* step_amp_;

  TailWidget* left_skew_;
  TailWidget* right_skew_;

  TailWidget* tail_;
};
