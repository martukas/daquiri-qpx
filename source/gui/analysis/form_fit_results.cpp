#include <gui/analysis/form_fit_results.h>
//#include "widget_detectors.h"
#include "ui_form_fit_results.h"
#include <QSettings>
#include <gui/widgets/qt_util.h>

FormFitResults::FormFitResults(DAQuiri::Fitter &fit, QWidget *parent) :
  QWidget(parent),
  ui(new Ui::FormFitResults),
  fit_data_(fit)
{
  ui->setupUi(this);

  loadSettings();

  ui->tablePeaks->verticalHeader()->hide();
  ui->tablePeaks->setColumnCount(7);
  ui->tablePeaks->setHorizontalHeaderLabels({"energy", "err(energy)", "fwhm", "err(fwhm)",
                                             "cps(hyp)", /*"err(hyp)",*/ "cps(S4)", "err(S4)"
                                             /*, "Quality"*/});
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

bool FormFitResults::save_close() {
  saveSettings();
  return true;
}

void FormFitResults::loadSettings() {
  QSettings settings_;

  settings_.beginGroup("Program");
  data_directory_ = settings_.value("save_directory", QDir::homePath() + "/qpx/data").toString();
  settings_.endGroup();

  settings_.beginGroup("peak_fitter");
  settings_.endGroup();

}

void FormFitResults::saveSettings() {
  QSettings settings_;

  settings_.beginGroup("peak_fitter");
  settings_.endGroup();
}

void FormFitResults::clear() {
  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(0);

  toggle_push();
}

void FormFitResults::update_data() {
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);

  ui->tablePeaks->clearContents();
  ui->tablePeaks->setRowCount(fit_data_.peaks().size());
  int i=0;
  for (auto &q : fit_data_.peaks()) {
    add_peak_to_table(q.second, i, false);
    ++i;
  }

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
  
  select_in_table();

  toggle_push();
}

void FormFitResults::add_peak_to_table(const DAQuiri::Peak &p, int row, bool gray) {
  QBrush background(gray ? Qt::lightGray : Qt::white);

  QColor yellow;
  yellow.setNamedColor("#EBDD8D");
  QColor orange;
  orange.setNamedColor("#E4B372");
  QColor red;
  red.setNamedColor("#E46D59");

  if (!p.good())
    background = yellow;

  QBrush eback = background;
  QBrush wback = background;
  QBrush s4back = background;

  if (p.sum4().quality() > 3)
    s4back = red;
  else if (p.sum4().quality() > 1)
    s4back = orange;

  if (p.quality_energy() > 2)
    eback = red;
  else if (p.quality_energy() > 1)
    eback = orange;

  if (p.quality_fwhm() > 2)
    wback = red;
  else if (p.quality_fwhm() > 1)
    wback = orange;

  // \todo reintroduce uncertainties
  add_to_table(ui->tablePeaks, row, 0, QString::number(p.energy()),
               QVariant::fromValue(p.center()), eback);
  //add_to_table(ui->tablePeaks, row, 1, p.energy().error_percent(), QVariant(), eback);
  add_to_table(ui->tablePeaks, row, 2, QString::number(p.fwhm()), QVariant(), wback);
//  add_to_table(ui->tablePeaks, row, 3, p.fwhm().error_percent(), QVariant(), wback);
  add_to_table(ui->tablePeaks, row, 4, QString::number(p.cps_hyp()), QVariant(), background);
  //  add_to_table(ui->tablePeaks, row, 5, p.cps_hyp.error_percent(), QVariant(), background);
  add_to_table(ui->tablePeaks, row, 5, QString::number(p.cps_sum4()), QVariant(), s4back);
//  add_to_table(ui->tablePeaks, row, 6, p.cps_sum4().error_percent(), QVariant(), s4back);
  //  add_to_table(ui->tablePeaks, row, 7, std::to_string(p.quality_),
  //               QVariant(), background);
}

void FormFitResults::selection_changed_in_table() {
  selected_peaks_.clear();
  foreach (QModelIndex i, ui->tablePeaks->selectionModel()->selectedRows())
    selected_peaks_.insert(ui->tablePeaks->item(i.row(), 0)->data(Qt::UserRole).toDouble());
  if (isVisible())
    emit selection_changed(selected_peaks_);
  toggle_push();
}

void FormFitResults::toggle_push() {
  ui->pushSaveReport->setEnabled(!fit_data_.peaks().empty());
}

void FormFitResults::update_selection(std::set<double> selected_peaks) {
  bool changed = (selected_peaks_ != selected_peaks);
  selected_peaks_ = selected_peaks;

  if (changed)
    select_in_table();
}

void FormFitResults::select_in_table() {
  ui->tablePeaks->blockSignals(true);
  this->blockSignals(true);
  ui->tablePeaks->clearSelection();

  QItemSelectionModel *selectionModel = ui->tablePeaks->selectionModel();
  QItemSelection itemSelection = selectionModel->selection();

  for (int i=0; i < ui->tablePeaks->rowCount(); ++i)
    if (selected_peaks_.count(ui->tablePeaks->item(i, 0)->data(Qt::UserRole).toDouble())) {
      ui->tablePeaks->selectRow(i);
      itemSelection.merge(selectionModel->selection(), QItemSelectionModel::Select);
    }

  selectionModel->clearSelection();
  selectionModel->select(itemSelection,QItemSelectionModel::Select);

  ui->tablePeaks->blockSignals(false);
  this->blockSignals(false);
}

void FormFitResults::on_pushSaveReport_clicked()
{
  emit save_peaks_request();
}
