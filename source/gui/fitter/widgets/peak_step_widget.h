#pragma once

#include <core/gamma/hypermet/step.h>
#include <QWidget>

class QLabel;
class QCheckBox;
class BoundedParameterWidget;

class StepWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit StepWidget(QString name, const DAQuiri::Step& s,
                      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::Step step() const { return current_; }

 signals:
  void updated();

 private slots:
  void enabled_clicked();
  void override_clicked();
  void param_changed();

 private:
  DAQuiri::Step original_;
  DAQuiri::Step current_;

  QLabel* modified_ {nullptr};
  QCheckBox* enabled_{nullptr};
  QCheckBox* override_{nullptr};
  BoundedParameterWidget* amplitude_ {nullptr};

  void update_values();
};
