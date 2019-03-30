#pragma once

#include <core/fitting/hypermet/Value.h>

#include <QWidget>

class QLabel;
class QCheckBox;
class QDoubleSpinBox;
class UncertainDoubleWidget;

class FitParameterWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit FitParameterWidget(const DAQuiri::SineBoundedValue& val,
      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::SineBoundedValue parameter() const { return parameter_; }

 signals:
  void updated();

 private slots:
  void to_fit_clicked();
  void valChanged(double);
  void minChanged(double);
  void maxChanged(double);

 private:
  DAQuiri::SineBoundedValue original_;
  DAQuiri::SineBoundedValue parameter_;

  QLabel* modified_ {nullptr};
  QCheckBox* to_fit_{nullptr};
  QDoubleSpinBox* val_{nullptr};
  QDoubleSpinBox* min_{nullptr};
  QDoubleSpinBox* max_{nullptr};
  UncertainDoubleWidget* unc_val_ {nullptr};

  void update();
};
