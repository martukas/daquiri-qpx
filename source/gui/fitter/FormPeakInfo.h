#pragma once

#include <QDialog>
#include <core/gamma/hypermet/peak.h>
#include <core/gamma/fit_settings.h>

class QCloseEvent;
class FitParameterWidget;
class TailWidget;
class UncertainDoubleWidget;

namespace Ui
{
class FormPeakInfo;
}

class FormPeakInfo : public QDialog
{
 Q_OBJECT

 public:
  explicit FormPeakInfo(DAQuiri::Peak& hm,
                        const DAQuiri::FCalibration& cal,
                        QWidget* parent = 0);
  ~FormPeakInfo();

 protected:
  void closeEvent(QCloseEvent*);

 private slots:
  void update();
  void on_buttonBox_accepted();
  void on_buttonBox_rejected();

 private:
  Ui::FormPeakInfo* ui;

  DAQuiri::FCalibration calib_;

  DAQuiri::Peak& peak_;
  DAQuiri::Peak peak_backup_;

  FitParameterWidget* position_;
  UncertainDoubleWidget* energy_;

  UncertainDoubleWidget* area_;

  FitParameterWidget* width_;
  UncertainDoubleWidget* fwhm_;
  UncertainDoubleWidget* fwhm_energy_;

  FitParameterWidget* step_amp_;

  TailWidget* left_skew_;
  TailWidget* right_skew_;

  TailWidget* tail_;

};
