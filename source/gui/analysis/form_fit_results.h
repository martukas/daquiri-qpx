#pragma once

#include <QWidget>
#include <QItemSelection>
#include <core/fitting/fitter.h>

#include <gui/widgets/SettingDelegate.h>

namespace Ui {
class FormFitResults;
}

class FormFitResults : public QWidget
{
  Q_OBJECT

public:
  explicit FormFitResults(DAQuiri::Fitter&, QWidget *parent = 0);
  ~FormFitResults();

  void clear();
  bool save_close();

public slots:
  void update_selection(std::set<double> selected_peaks);
  void update_data();

signals:
  void selection_changed(std::set<double> selected_peaks);
  void save_peaks_request();

private slots:
  void selection_changed_in_table();
  void toggle_push();

  void on_pushSaveReport_clicked();


private:
  Ui::FormFitResults *ui;

  //from parent
  QString data_directory_;

  DAQuiri::Fitter &fit_data_;
  std::set<double> selected_peaks_;

  void loadSettings();
  void saveSettings();
  void select_in_table();

  void add_peak_to_table(const DAQuiri::Peak &, int, bool);

};

