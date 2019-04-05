#pragma once

#include <QDialog>
#include <core/gamma/hypermet/gamma_region.h>
#include <core/gamma/fit_settings.h>

class QCloseEvent;
class BoundedParameterWidget;
class PositiveParameterWidget;
class SkewWidget;
class StepWidget;
class UncertainDoubleWidget;

namespace Ui
{
class RegionDialog;
}

class RegionDialog : public QDialog
{
 Q_OBJECT

 public:
  explicit RegionDialog(DAQuiri::Region& region,
                      const DAQuiri::FCalibration& calibration,
                      QWidget* parent = 0);
  ~RegionDialog();

 protected:
  void closeEvent(QCloseEvent*);

 private slots:
  void update();
  void on_buttonBox_accepted();
  void on_buttonBox_rejected();

 private:
  Ui::RegionDialog* ui;

  DAQuiri::FCalibration calibration_;

  DAQuiri::Region& region_;
  DAQuiri::Region region_backup_;

  BoundedParameterWidget* background_base_;
  BoundedParameterWidget* background_slope_;
  BoundedParameterWidget* background_curve_;

  BoundedParameterWidget* width_;
  UncertainDoubleWidget* fwhm_;
  UncertainDoubleWidget* fwhm_energy_;

  SkewWidget* left_skew_;
  SkewWidget* right_skew_;

  StepWidget* step_;
  SkewWidget* tail_;
};
