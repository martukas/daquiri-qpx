#include <gui/fitter/WidgetFitterAutoselect.h>
#include "ui_WidgetFitterAutoselect.h"

#include <core/util/custom_logger.h>

WidgetFitterAutoselect::WidgetFitterAutoselect(QWidget *parent)
  : QWidget(parent)
  , ui(new Ui::WidgetFitterAutoselect)
{
  ui->setupUi(this);

  ui->comboAutoMethod->addItem("None");
  ui->comboAutoMethod->addItem("Max(area)");
  ui->comboAutoMethod->addItem("Max(height)");
  ui->comboAutoMethod->addItem("Nearest");
}

void WidgetFitterAutoselect::loadSettings(QSettings &set, QString name)
{
  if (name.isEmpty())
    name = "WidgetFitterAutoselect";
  set.beginGroup(name);
  ui->comboAutoMethod->setCurrentText(set.value("autosel_method", "None").toString());
  ui->checkReselect->setChecked(set.value("reselect", false).toBool());
  ui->doubleSlack->setValue(set.value("reselect_slack", 0.5).toDouble());
  set.endGroup();
}

void WidgetFitterAutoselect::saveSettings(QSettings &set, QString name)
{
  if (name.isEmpty())
    name = "WidgetFitterAutoselect";
  set.beginGroup(name);
  set.setValue("autosel_method", ui->comboAutoMethod->currentText());
  set.setValue("reselect", ui->checkReselect->isChecked());
  set.setValue("reselect_slack", ui->doubleSlack->value());
  set.endGroup();
}

WidgetFitterAutoselect::~WidgetFitterAutoselect()
{
  delete ui;
}

std::set<double> WidgetFitterAutoselect::autoselect(const std::map<double, DAQuiri::Peak> &peaks,
                                              const std::set<double> &prevsel) const
{
  double selection = peaks.begin()->first;

  if (ui->comboAutoMethod->currentText() == "Max(height)")
  {
    double maxheight = peaks.begin()->first;
    for (auto &p : peaks) {
      if (p.second.hypermet().height().value() > maxheight)
      {
        maxheight = p.second.hypermet().height().value();
        selection = p.first;
      }
    }
  }
  else if (ui->comboAutoMethod->currentText() == "Max(area)")
  {
    double maxarea = peaks.begin()->first;
    for (auto &p : peaks) {
      if (p.second.area_best() > maxarea) {
        maxarea = p.second.area_best();
        selection = p.first;
      }
    }
  }
  else if (ui->comboAutoMethod->currentText() == "Nearest")
  {
    if (prevsel.empty())
      return std::set<double>();

    //make this work for multiple peaks
    double prevpeak = *prevsel.begin();
    double delta = std::numeric_limits<double>::max();
    for (auto &p : peaks) {
      if (abs(p.second.center() - prevpeak) < delta) {
        delta = abs(p.second.center() - prevpeak);
        selection = p.first;
      }
    }
  }

  std::set<double> newselr;
  if (ui->comboAutoMethod->currentText() != "None")
    newselr.insert(selection);
  return newselr;
}

std::set<double> WidgetFitterAutoselect::reselect(const std::map<double, DAQuiri::Peak> &peaks,
                                              const std::set<double> &prevsel) const
{
  if (!ui->checkReselect->isChecked() || prevsel.empty() || peaks.empty())
    return prevsel;

  std::set<double> newselr;

  for (auto &prev : prevsel)
  {
    double selection = peaks.begin()->first;
    double delta = std::numeric_limits<double>::max();
    for (auto &p : peaks) {
      if (abs(p.second.center() - prev) < delta) {
        delta = abs(p.second.center() - prev);
        selection = p.first;
      }
    }
    if (abs(selection - prev) < ui->doubleSlack->value() * peaks.at(selection).fwhm())
      newselr.insert(selection);
  }

  if (newselr.size() == prevsel.size())
    return newselr;
  else
    return prevsel;
}

