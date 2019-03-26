#pragma once

#include <core/fitting/hypermet/Tail.h>
#include <QWidget>

class QLabel;
class QCheckBox;
class FitParameterWidget;

class TailWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit TailWidget(QString name, const DAQuiri::Tail& t,
                      double spin_width, double label_width, QWidget* parent = 0);
  bool changed() const;
  DAQuiri::Tail tail() const { return tail_; }

 signals:
  void updated();

 private slots:
  void enabled_clicked();
  void override_clicked();
  void param_changed();

 private:
  DAQuiri::Tail original_;
  DAQuiri::Tail tail_;

  QLabel* modified_ {nullptr};
  QCheckBox* enabled_{nullptr};
  QCheckBox* override_{nullptr};
  FitParameterWidget* amplitude_ {nullptr};
  FitParameterWidget* slope_ {nullptr};

  void update();
};
