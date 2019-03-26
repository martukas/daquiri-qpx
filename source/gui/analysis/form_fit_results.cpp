#include <gui/analysis/form_fit_results.h>
//#include "widget_detectors.h"
#include "ui_form_fit_results.h"
#include <QSettings>
#include <gui/widgets/qt_util.h>

FormFitResults::FormFitResults(DAQuiri::Fitter& fit, QWidget* parent) :
    QWidget(parent), ui(new Ui::FormFitResults), fit_data_(fit)
{
  ui->setupUi(this);

  loadSettings();

  ui->tablePeaks->verticalHeader()->hide();
  ui->tablePeaks->setColumnCount(11);
  ui->tablePeaks->setHorizontalHeaderLabels({"energy", "\u03C3", "err",
                                             "area(hyp)", "\u03C3", "err",
                                             "\u03C7\u00B2-norm",
                                             "area(S4)", "\u03C3", "err",
                                             "Quality(S4)"});
  ui->tablePeaks->setSelectionBehavior(QAbstractItemView::SelectRows);
  ui->tablePeaks->setSelectionMode(QAbstractItemView::ExtendedSelection);
  ui->tablePeaks->setEditTriggers(QTableView::NoEditTriggers);
  ui->tablePeaks->horizontalHeader()->setStretchLastSection(true);
  ui->tablePeaks->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  ui->tablePeaks->show();

  connect(ui->tablePeaks, SIGNAL(itemSelectionChanged()), this, SLOT(selection_changed_in_table()));
}

FormFitResults::~FormFitResults()
{
  delete ui;
}

bool FormFitResults::save_close()
{
  saveSettings();
  return true;
}

void FormFitResults::loadSettings()
{
  QSettings settings_;

  settings_.beginGroup("Program");
  data_directory_ = settings_.value("save_directory", QDir::homePath() + "/qpx/data").toString();
  settings_.endGroup();

  settings_.beginGroup("peak_fitter");
  settings_.endGroup();

}

void FormFitResults::saveSettings()
{
  QSettings settings_;

  settings_.beginGroup("peak_fitter");
  settings_.endGroup();
}

void FormFitResults::clear()
{
  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(0);

  toggle_push();
}

void FormFitResults::update_data()
{
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);

  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(fit_data_.peaks().size());
  int i = 0;
  for (auto& q : fit_data_.peaks())
  {
    add_peak_to_table(q.second, i, false);
    ++i;
  }

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);

  select_in_table();

  toggle_push();
}

void FormFitResults::add_peak_to_table(const DAQuiri::Peak& p, int row, bool gray)
{
  double max_param_err = 10.0;
  double max_chi2 = 10.0;

  QBrush background(gray ? Qt::lightGray : Qt::white);

  QColor yellow;
  yellow.setNamedColor("#EBDD8D");
  QColor orange;
  orange.setNamedColor("#E4B372");
  QColor red;
  red.setNamedColor("#E46D59");

  auto calib = fit_data_.settings().calib;

  auto energy = p.peak_energy(calib.cali_nrg_);
  auto width = p.fwhm_energy(calib.cali_nrg_);

  int qual_energy = 1;
  if (energy.error_percent() > max_param_err)
    qual_energy = 3;
  else if (!std::isfinite(energy.sigma()) || !energy.sigma())
    qual_energy = 2;

  int qual_width = 1;
  if (width.error_percent() > max_param_err)
    qual_width = 3;
  else if (!std::isfinite(width.sigma()) || !width.sigma())
    qual_width = 2;

  int qual_chi = 1;
  if (p.chi_sq_norm > max_chi2)
    qual_chi = 3;
  else if (!std::isfinite(p.chi_sq_norm) || !p.chi_sq_norm)
    qual_chi = 2;

  bool good = ((p.sum4.quality() == 1)
      && (qual_chi == 1)
      && (qual_energy == 1)
      && (qual_width == 1));

  if (!good)
    background = yellow;

  QBrush s4back = background;
  if (p.sum4.quality() > 3)
    s4back = red;
  else if (p.sum4.quality() > 1)
    s4back = orange;

  QBrush hyback = background;
  if (qual_chi > 2)
    hyback = red;
  else if (qual_chi > 1)
    hyback = orange;

  QBrush eback = background;
  if (qual_energy > 2)
    eback = red;
  else if (qual_energy > 1)
    eback = orange;

  QBrush wback = background;
  if (qual_width > 2)
    wback = red;
  else if (qual_width > 1)
    wback = orange;

  auto area_hyp = p.area();
  auto area_s4 = p.sum4.peak_area();

  add_to_table(ui->tablePeaks, row, 0, QString::number(energy.value()), QVariant::fromValue(p.id()), eback);
  add_to_table(ui->tablePeaks, row, 1, QString::number(energy.sigma()), {}, eback);
  add_to_table(ui->tablePeaks, row, 2, QString::number(energy.error_percent()), {}, eback);
  add_to_table(ui->tablePeaks, row, 3, QString::number(area_hyp.value()), {}, hyback);
  add_to_table(ui->tablePeaks, row, 4, QString::number(area_hyp.sigma()), {}, hyback);
  add_to_table(ui->tablePeaks, row, 5, QString::number(area_hyp.error_percent()), {}, hyback);
  add_to_table(ui->tablePeaks, row, 6, QString::number(p.chi_sq_norm), {}, background);
  add_to_table(ui->tablePeaks, row, 7, QString::number(area_s4.value()), {}, s4back);
  add_to_table(ui->tablePeaks, row, 8, QString::number(area_s4.sigma()), {}, s4back);
  add_to_table(ui->tablePeaks, row, 9, QString::number(area_s4.error_percent()), {}, s4back);
  add_to_table(ui->tablePeaks, row, 10, QString::number(p.sum4.quality()), {}, s4back);
}

void FormFitResults::selection_changed_in_table()
{
  selected_peaks_.clear();
      foreach (QModelIndex i,
               ui->tablePeaks->selectionModel()->selectedRows())selected_peaks_.insert(ui->tablePeaks->item(i.row(),
                                                                                                            0)->data(Qt::UserRole).toDouble());
  if (isVisible())
      emit selection_changed(selected_peaks_);
  toggle_push();
}

void FormFitResults::toggle_push()
{
  ui->pushSaveReport->setEnabled(!fit_data_.peaks().empty());
}

void FormFitResults::update_selection(std::set<double> selected_peaks)
{
  bool changed = (selected_peaks_ != selected_peaks);
  selected_peaks_ = selected_peaks;

  if (changed)
    select_in_table();
}

void FormFitResults::select_in_table()
{
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);
  ui->tablePeaks->clearSelection();

  QItemSelectionModel* selectionModel = ui->tablePeaks->selectionModel();
  QItemSelection itemSelection = selectionModel->selection();

  for (int i = 0; i < ui->tablePeaks->rowCount(); ++i)
    if (selected_peaks_.count(ui->tablePeaks->item(i, 0)->data(Qt::UserRole).toDouble()))
    {
      ui->tablePeaks->selectRow(i);
      itemSelection.merge(selectionModel->selection(), QItemSelectionModel::Select);
    }

  selectionModel->clearSelection();
  selectionModel->select(itemSelection, QItemSelectionModel::Select);

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormFitResults::on_pushSaveReport_clicked()
{
  emit save_peaks_request();
}
