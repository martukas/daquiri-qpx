#pragma once

// \todo generalize for other implementations of param
#include <core/fitting/parameter/sine_bounded_param.h>

#include <QWidget>

class QLabel;
class QCheckBox;
class QDoubleSpinBox;
class UncertainDoubleWidget;

class FitParameterWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit FitParameterWidget(const DAQuiri::SineBoundedParam& val,
      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::SineBoundedParam parameter() const { return parameter_; }

 signals:
  void updated();

 private slots:
  void to_fit_clicked();
  void valChanged(double);
  void minChanged(double);
  void maxChanged(double);

 private:
  DAQuiri::SineBoundedParam original_;
  DAQuiri::SineBoundedParam parameter_;

  QLabel* modified_ {nullptr};
  QCheckBox* to_fit_{nullptr};
  QDoubleSpinBox* val_{nullptr};
  QDoubleSpinBox* min_{nullptr};
  QDoubleSpinBox* max_{nullptr};
  UncertainDoubleWidget* unc_val_ {nullptr};

  void update();
};
