#pragma once

#include <core/fitting/parameter/unbounded_param.h>

#include <QWidget>

class QLabel;
class QCheckBox;
class QDoubleSpinBox;
class UncertainDoubleWidget;

class UnboundedParameterWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit UnboundedParameterWidget(const DAQuiri::UnboundedParam& val,
      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  void set_decimals(int decimals);
  void set_step(double step);
  DAQuiri::UnboundedParam parameter() const { return current_; }

 signals:
  void updated();

 private slots:
  void to_fit_clicked();
  void valChanged(double);

 private:
  DAQuiri::UnboundedParam original_;
  DAQuiri::UnboundedParam current_;

  QLabel* modified_ {nullptr};
  QCheckBox* to_fit_{nullptr};
  QDoubleSpinBox* val_{nullptr};
  UncertainDoubleWidget* unc_val_ {nullptr};

  void update_values();
};
