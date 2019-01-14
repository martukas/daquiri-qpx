#pragma once

#include <QWidget>
#include <QSettings>
#include <core/fitting/peak.h>

namespace Ui {
class WidgetFitterAutoselect;
}

class WidgetFitterAutoselect : public QWidget
{
  Q_OBJECT

public:

  explicit WidgetFitterAutoselect(QWidget *parent = 0);
  ~WidgetFitterAutoselect();

  std::set<double> autoselect(const std::map<double, DAQuiri::Peak> &peaks,
                              const std::set<double> &prevsel) const;

  std::set<double> reselect(const std::map<double, DAQuiri::Peak> &peaks,
                            const std::set<double> &prevsel) const;

  void loadSettings(QSettings &set, QString name = "");
  void saveSettings(QSettings &set, QString name = "");

private:
  Ui::WidgetFitterAutoselect *ui;
};
