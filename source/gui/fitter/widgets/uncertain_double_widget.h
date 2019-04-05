#pragma once

#include <core/fitting/uncertain.h>

#include <QWidget>

class QLabel;

class UncertainDoubleWidget : public QWidget
{
 Q_OBJECT

 public:
  explicit UncertainDoubleWidget(double width, QWidget* parent = 0);
  void update(const UncertainDouble& val);

 private:
  QLabel* val_ {nullptr};
  QLabel* sigma_{nullptr};
  QLabel* percent_{nullptr};

  static QLabel* make_framed_label(double width);
};
