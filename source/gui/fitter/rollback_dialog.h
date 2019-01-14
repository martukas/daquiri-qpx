#pragma once

#include <QDialog>
#include <QRadioButton>
#include <core/fitting/roi.h>

class RollbackDialog : public QDialog {
  Q_OBJECT

public:
  explicit RollbackDialog(DAQuiri::ROI setting, QWidget *parent = 0);
  int get_choice();

private:
  DAQuiri::ROI      roi_;
  std::vector<QRadioButton*> radios_;
};
