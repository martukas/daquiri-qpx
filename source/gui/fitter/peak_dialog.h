#pragma once

#include <QDialog>
#include <core/gamma/hypermet/peak.h>
#include <core/gamma/fit_settings.h>

class QCloseEvent;
class BoundedParameterWidget;
class PositiveParameterWidget;
class SkewWidget;
class StepWidget;
class UncertainDoubleWidget;

namespace Ui
{
class PeakDialog;
}

class PeakDialog : public QDialog
{
 Q_OBJECT

 public:
  explicit PeakDialog(DAQuiri::Peak& peak,
                      const DAQuiri::FCalibration& calibration,
                      QWidget* parent = 0);
  ~PeakDialog();

 protected:
  void closeEvent(QCloseEvent*);

 private slots:
  void update();
  void on_buttonBox_accepted();
  void on_buttonBox_rejected();

 private:
  Ui::PeakDialog* ui;

  DAQuiri::FCalibration calibration_;

  DAQuiri::Peak& peak_;
  DAQuiri::Peak peak_backup_;

  BoundedParameterWidget* position_;
  UncertainDoubleWidget* energy_;

  PositiveParameterWidget* amplitude_;
  UncertainDoubleWidget* area_;

  BoundedParameterWidget* width_;
  UncertainDoubleWidget* fwhm_;
  UncertainDoubleWidget* fwhm_energy_;

  SkewWidget* left_skew_;
  SkewWidget* right_skew_;

  StepWidget* step_;
  SkewWidget* tail_;
};
