#pragma once

#include <QDialog>
#include <QRadioButton>
#include <core/fitting/region_manager.h>

class RollbackDialog : public QDialog {
  Q_OBJECT

public:
  explicit RollbackDialog(DAQuiri::RegionManager setting, QWidget *parent = 0);
  int get_choice();

private:
  DAQuiri::RegionManager      roi_;
  std::vector<QRadioButton*> radios_;
};
