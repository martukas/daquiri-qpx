#pragma once

#include <QDialog>
#include <QCloseEvent>
#include <core/fitting/hypermet.h>

namespace Ui {
class FormPeakInfo;
}

class FormPeakInfo : public QDialog
{
  Q_OBJECT

public:
  explicit FormPeakInfo(DAQuiri::Hypermet &hm, QWidget *parent = 0);
  ~FormPeakInfo();

protected:
  void closeEvent(QCloseEvent*);


private slots:
  void on_buttonBox_accepted();

  void on_buttonBox_rejected();

  void on_doubleMinRskewSlope_valueChanged(double arg1);
  void on_doubleMaxRskewSlope_valueChanged(double arg1);
  void on_doubleMinRskewAmp_valueChanged(double arg1);
  void on_doubleMaxRskewAmp_valueChanged(double arg1);
  void on_doubleMinLskewSlope_valueChanged(double arg1);
  void on_doubleMaxLskewSlope_valueChanged(double arg1);
  void on_doubleMinLskewAmp_valueChanged(double arg1);
  void on_doubleMaxLskewAmp_valueChanged(double arg1);
  void on_doubleMinTailSlope_valueChanged(double arg1);
  void on_doubleMaxTailSlope_valueChanged(double arg1);
  void on_doubleMinTailAmp_valueChanged(double arg1);
  void on_doubleMaxTailAmp_valueChanged(double arg1);
  void on_doubleMinStep_valueChanged(double arg1);
  void on_doubleMaxStep_valueChanged(double arg1);

private:
  Ui::FormPeakInfo *ui;

  DAQuiri::Hypermet &hm_;

  void enforce_bounds();
};
