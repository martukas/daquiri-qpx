#pragma once

#include <core/fitting/parameter/positive_param.h>

#include <QWidget>

class QLabel;
class QCheckBox;
class QDoubleSpinBox;
class UncertainDoubleWidget;

class PositiveParameterWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit PositiveParameterWidget(const DAQuiri::PositiveParam& val,
      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  void set_decimals(int decimals);
  void set_step(double step);
  DAQuiri::PositiveParam parameter() const { return current_; }

 signals:
  void updated();

 private slots:
  void to_fit_clicked();
  void valChanged(double);

 private:
  DAQuiri::PositiveParam original_;
  DAQuiri::PositiveParam current_;

  QLabel* modified_ {nullptr};
  QCheckBox* to_fit_{nullptr};
  QDoubleSpinBox* val_{nullptr};
  UncertainDoubleWidget* unc_val_ {nullptr};

  void update_values();
};
