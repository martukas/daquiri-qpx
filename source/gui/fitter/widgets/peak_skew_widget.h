#pragma once

#include <core/gamma/hypermet/skew.h>
#include <QWidget>

class QLabel;
class QCheckBox;
class BoundedParameterWidget;

class SkewWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit SkewWidget(QString name, const DAQuiri::Skew& t,
                      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::Skew skew() const { return current_; }

 signals:
  void updated();

 private slots:
  void enabled_clicked();
  void override_clicked();
  void param_changed();

 private:
  DAQuiri::Skew original_;
  DAQuiri::Skew current_;

  QLabel* modified_ {nullptr};
  QCheckBox* enabled_{nullptr};
  QCheckBox* override_{nullptr};
  BoundedParameterWidget* amplitude_ {nullptr};
  BoundedParameterWidget* slope_ {nullptr};

  void update_values();
};
