#include <gui/fitter/FormFitter.h>
#include "ui_FormFitter.h"

//#include "widget_detectors.h"

#include <core/gamma/fitter.h>
#include <QPlot/QPlotButton.h>
#include <gui/fitter/FormFitterSettings.h>
#include <gui/fitter/peak_dialog.h>
#include <gui/fitter/region_dialog.h>
#include <gui/fitter/rollback_dialog.h>

FormFitter::FormFitter(QWidget* parent) :
    QWidget(parent), fit_(nullptr), ui(new Ui::FormFitter)
{
  ui->setupUi(this);
//  player = new QMediaPlayer(this);

  connect(ui->plot, SIGNAL(selectionChangedByUser()), this, SLOT(selection_changed()));
  connect(ui->plot, SIGNAL(range_selection_changed(double, double)),
          this, SLOT(update_range_selection(double, double)));

  connect(ui->plot, SIGNAL(add_peak(double, double)), this, SLOT(add_peak(double, double)));
  connect(ui->plot, SIGNAL(delete_selected_peaks()), this, SLOT(delete_selected_peaks()));
  connect(ui->plot, SIGNAL(adjust_sum4(double, double, double)),
          this, SLOT(adjust_sum4(double, double, double)));
  connect(ui->plot, SIGNAL(adjust_background_L(double, double, double)),
          this, SLOT(adjust_background_L(double, double, double)));
  connect(ui->plot, SIGNAL(adjust_background_R(double, double, double)),
          this, SLOT(adjust_background_R(double, double, double)));
  connect(ui->plot, SIGNAL(peak_info(double)), this, SLOT(peak_info(double)));

  connect(ui->plot, SIGNAL(rollback_ROI(double)), this, SLOT(rollback_ROI(double)));
  connect(ui->plot, SIGNAL(roi_settings(double)), this, SLOT(roi_settings(double)));
  connect(ui->plot, SIGNAL(refit_ROI(double)), this, SLOT(refit_ROI(double)));
  connect(ui->plot, SIGNAL(delete_ROI(double)), this, SLOT(delete_ROI(double)));

  QShortcut* shortcut = new QShortcut(QKeySequence(Qt::Key_Backspace), ui->plot);
  connect(shortcut, SIGNAL(activated()), ui->plot, SLOT(zoomOut()));

//  QShortcut *shortcut2 = new QShortcut(QKeySequence(Qt::Key_Insert), ui->plot);
//  connect(shortcut2, SIGNAL(activated()), this, SLOT(add_peak()));

  QShortcut* shortcut3 = new QShortcut(QKeySequence(QKeySequence::Delete), ui->plot);
  connect(shortcut3, SIGNAL(activated()), this, SLOT(delete_selected_peaks()));

  //  QShortcut* shortcut4 = new QShortcut(QKeySequence(QKeySequence(Qt::Key_L)), ui->plot);
  //  connect(shortcut4, SIGNAL(activated()), this, SLOT(switch_scale_type()));

  connect(&thread_fitter_, SIGNAL(fit_updated(DAQuiri::Fitter)), this, SLOT(fit_updated(DAQuiri::Fitter)));
  connect(&thread_fitter_, SIGNAL(fitting_done()), this, SLOT(fitting_complete()));

  QMovie* movie = new QMovie(":/icons/loader.gif");
  ui->labelMovie->setMovie(movie);
  ui->labelMovie->show();
  movie->start();
  ui->labelMovie->setVisible(false);

  thread_fitter_.start();
}

FormFitter::~FormFitter()
{
  thread_fitter_.stop_work();
  thread_fitter_.terminate(); //in thread itself?

  delete ui;
}

void FormFitter::setFit(DAQuiri::Fitter* fit)
{
  fit_ = fit;
  ui->plot->setFit(fit);
  update_spectrum();
  updateData();
}

void FormFitter::loadSettings(QSettings& settings_)
{
  settings_.beginGroup("Peaks");
  //  scale_log_ = settings_.value("scale_log", true).toBool();
  settings_.endGroup();
}

void FormFitter::saveSettings(QSettings& settings_)
{
  settings_.beginGroup("Peaks");
  //  settings_.setValue("scale_log", scale_log_);
  settings_.endGroup();
}

void FormFitter::clear()
{
  if (!fit_ || busy_)
    return;

  //DBG << "FormFitter::clear()";

  ui->plot->clearAll();

  clearSelection();
  ui->plot->replot();
  toggle_push(busy_);
}

void FormFitter::update_spectrum()
{
  ui->plot->update_spectrum();
  toggle_push(busy_);
}

void FormFitter::refit_ROI(double ROI_bin)
{
  if (!fit_ || busy_)
    return;

  toggle_push(true);

  thread_fitter_.set_data(*fit_);
  thread_fitter_.refit(ROI_bin);
}

void FormFitter::rollback_ROI(double ROI_bin)
{
  if (!fit_ || busy_)
    return;

  if (fit_->contains_region(ROI_bin))
  {
    RollbackDialog* editor = new RollbackDialog(fit_->region(ROI_bin), qobject_cast<QWidget*>(parent()));
    int ret = editor->exec();
    if (ret == QDialog::Accepted)
    {
      fit_->rollback_ROI(ROI_bin, editor->get_choice());
      toggle_push(false);
      updateData();

      emit data_changed();
      emit fitting_done();
    }
  }
}

void FormFitter::delete_ROI(double ROI_bin)
{
  if (!fit_ || busy_)
    return;

  fit_->delete_ROI(ROI_bin);
  toggle_push(false);
  updateData();

  emit data_changed();
  emit fitting_done();
}

void FormFitter::on_pushClearAll_clicked()
{
  if (!fit_ || busy_)
    return;

  fit_->clear_all_ROIs();
  toggle_push(false);
  updateData();

  emit data_changed();
  emit fitting_done();
}

void FormFitter::make_range(double energy)
{
  if (!fit_ || busy_)
    return;

  ui->plot->make_range_selection(energy);
  ui->plot->follow_selection();
}

void FormFitter::set_selected_peaks(std::set<double> selected_peaks)
{
  ui->plot->set_selected_peaks(selected_peaks);
  toggle_push(busy_);
}

std::set<double> FormFitter::get_selected_peaks()
{
  return ui->plot->get_selected_peaks();
}

void FormFitter::on_pushFindPeaks_clicked()
{
  perform_fit();
}

void FormFitter::perform_fit()
{
  if (!fit_ || busy_)
    return;

  fit_->find_regions();
  //  DBG << "number of peaks found " << fit_data_->peaks_.size();

  toggle_push(true);

  thread_fitter_.set_data(*fit_);
  thread_fitter_.fit_peaks();

}

void FormFitter::add_peak(double l, double r)
{
  if (!fit_ || busy_)
    return;

  ui->plot->clear_range_selection();

  double l_bin = fit_->settings().calib.nrg_to_bin(l);
  double r_bin = fit_->settings().calib.nrg_to_bin(r);

  auto relevant_regions = fit_->relevant_regions(l_bin, r_bin);

  if (relevant_regions.size() > 1)
  {
    // ambiguous which region
    merge_regions(relevant_regions);
    // \todo and also add to merged region
    return;
  }
  else if (relevant_regions.size() == 1)
  {
    // Add peak to one region
    auto id = fit_->add_peak(*relevant_regions.begin(), l_bin, r_bin);
    dirty(id);
    return;
  }
  else
  {
    // Create new region
    auto id = fit_->create_region(l_bin, r_bin);
    dirty(id);
    return;
  }
}

void FormFitter::adjust_sum4(double peak_id, double l, double r)
{
  if (!fit_ || busy_)
    return;

  if (fit_->adjust_sum4(peak_id,
                             fit_->settings().calib.nrg_to_bin(l),
                             fit_->settings().calib.nrg_to_bin(r)))
  {
    updateData();
    std::set<double> selected_peaks;
    selected_peaks.insert(peak_id);
    ui->plot->set_selected_peaks(selected_peaks);

    emit data_changed();
    emit fitting_done();
  }
  else
    ui->plot->clear_range_selection();
}

void FormFitter::adjust_background_L(double ROI_id, double l, double r)
{
  if (!fit_ || busy_)
    return;

  ui->plot->clear_range_selection();

  if (!fit_->contains_region(ROI_id))
  {
    WARN("No such region: id={}", ROI_id);
    return;
  }

  std::set<double> rois = fit_->relevant_regions(
      fit_->settings().calib.nrg_to_bin(l),
      fit_->region(ROI_id).right_bin());

  if (!rois.count(ROI_id))
  {
    QMessageBox::information(this, "Out of bounds", "Background sample bad. Very bad...");
    return;
  }

  bool merge = ((rois.size() > 1) &&
      (QMessageBox::question(this, "Merge?", "Regions overlap. Merge them?") == QMessageBox::Yes));

  thread_fitter_.set_data(*fit_);

  toggle_push(true);

  if (merge)
    merge_regions(rois);
  else
  {
    fit_->adj_LB(ROI_id,
                      fit_->settings().calib.nrg_to_bin(l),
                      fit_->settings().calib.nrg_to_bin(r));

    dirty(ROI_id);
  }
}

void FormFitter::adjust_background_R(double ROI_id, double l, double r)
{
  if (!fit_ || busy_)
    return;

  ui->plot->clear_range_selection();

  if (!fit_->contains_region(ROI_id))
  {
    WARN("No such region: id={}", ROI_id);
    return;
  }

  std::set<double> rois = fit_->relevant_regions(
      fit_->region(ROI_id).left_bin(),
      fit_->settings().calib.nrg_to_bin(r));

  if (!rois.count(ROI_id))
  {
    QMessageBox::information(this, "Out of bounds", "Background sample bad. Very bad...");
    return;
  }

  bool merge = ((rois.size() > 1) &&
      (QMessageBox::question(this, "Merge?", "Regions overlap. Merge them?") == QMessageBox::Yes));

  thread_fitter_.set_data(*fit_);

  toggle_push(true);

  if (merge)
    merge_regions(rois);
  else
  {
    fit_->adj_RB(ROI_id,
                      fit_->settings().calib.nrg_to_bin(l),
                      fit_->settings().calib.nrg_to_bin(r));

    dirty(ROI_id);
  }
}

void FormFitter::merge_regions(std::set<double> rois)
{
  while (rois.size() > 1)
  {
    auto id1 = *rois.begin();
    rois.erase(id1);
    auto id2 = *rois.begin();

    auto region1 = fit_->region(id1);
    auto region2 = fit_->region(id2);

    std::string r1str =
        fmt::format("[{}, {}]",
                    fit_->settings().calib.nrg_to_bin(region1.left_bin()),
                    fit_->settings().calib.nrg_to_bin(region1.right_bin()));

    std::string r2str =
        fmt::format("[{}, {}]",
                    fit_->settings().calib.nrg_to_bin(region2.left_bin()),
                    fit_->settings().calib.nrg_to_bin(region2.right_bin()));

    std::string message = fmt::format("Merge regions {} and {} ?", r1str, r2str);

    bool merge = QMessageBox::question(this, "Merge?",
        QString::fromStdString(message)) == QMessageBox::Yes;

    if (!merge)
      break;

    auto newr = fit_->merge_regions(id1, id2);
    if (newr == -1)
    {
      WARN("Merge failed");
      break;
    }

    rois.erase(id2);
    rois.insert(newr);
  }

  // \todo refit merged regions?
}

void FormFitter::delete_selected_peaks()
{
  if (!fit_ || busy_)
    return;

  std::set<double> chosen_peaks = this->get_selected_peaks();
  if (chosen_peaks.empty())
    return;

  clearSelection();

  fit_->remove_peaks(chosen_peaks);

  fitting_complete();

  // \todo refirt if some regions dirty
  // \todo delete some empty regions
}

void FormFitter::fit_updated(DAQuiri::Fitter fitter)
{
  //  while (player->state() == QMediaPlayer::PlayingState)
  //    player->stop();

  //  if (player->state() != QMediaPlayer::PlayingState) {
  //    player->setMedia(QUrl("qrc:/sounds/laser6.wav"));
  //    player->setVolume(100);
  //    player->setPosition(0);
  //    player->play();
  ////    while (player->state() == QMediaPlayer::PlayingState) {}
  //  }

  *fit_ = fitter;
  toggle_push(busy_);
  updateData();;
  emit data_changed();
}

void FormFitter::fitting_complete()
{
  //  while (player->state() == QMediaPlayer::PlayingState)
  //    player->stop();

  //  if (player->state() != QMediaPlayer::PlayingState) {
  //    player->setMedia(QUrl("qrc:/sounds/laser12.wav"));
  //    player->setVolume(100);
  //    player->setPosition(0);
  //    player->play();
  ////    while (player->state() == QMediaPlayer::PlayingState) {}
  //  }

  //  busy_= false;
  //  calc_visible();
  //  ui->plot->replot();
  toggle_push(false);
  emit data_changed();
  emit fitting_done();
}

void FormFitter::dirty(double region_id)
{
  updateData();
  fitting_complete();

  bool refit = (QMessageBox::question(this, "Refit?",
                                      "Region at bin=" + QString::number(region_id) + " modified. Refit?")
      == QMessageBox::Yes);

  thread_fitter_.set_data(*fit_);

  if (refit)
  {
    toggle_push(true);
    thread_fitter_.set_data(*fit_);
    thread_fitter_.refit(region_id);
  }
}

void FormFitter::toggle_push(bool busy)
{
  busy_ = busy;
  bool hasdata = (fit_ && fit_->region_count());

  ui->pushStopFitter->setEnabled(busy_);
  ui->pushFindPeaks->setEnabled(!busy_);
  ui->pushSettings->setEnabled(!busy_);
  ui->pushClearAll->setEnabled(!busy_ && hasdata);
  ui->labelMovie->setVisible(busy_);
  ui->plot->set_busy(busy);

  emit fitter_busy(busy_);
}

void FormFitter::on_pushStopFitter_clicked()
{
  ui->pushStopFitter->setEnabled(false);
  thread_fitter_.stop_work();
}

void FormFitter::clearSelection()
{
  ui->plot->clearSelection();
}

void FormFitter::selection_changed()
{
  toggle_push(busy_);
  emit peak_selection_changed(ui->plot->get_selected_peaks());
}

void FormFitter::update_range_selection(double l, double r)
{
  emit range_selection_changed(l, r);
}

void FormFitter::updateData()
{
  ui->plot->updateData();
}

void FormFitter::on_pushSettings_clicked()
{
  if (!fit_)
    return;

  DAQuiri::FitSettings fs = fit_->settings();
  FormFitterSettings* FitterSettings = new FormFitterSettings(fs, this);
  FitterSettings->setWindowTitle("Fitter settings");
  int ret = FitterSettings->exec();

  if (ret == QDialog::Accepted)
  {
    fit_->apply_settings(fs);
    updateData();
  }
}

void FormFitter::roi_settings(double roi)
{
  if (!fit_ || !fit_->contains_region(roi))
    return;

  auto region = fit_->region(roi).region();

  RegionDialog* regionDialog = new RegionDialog(region, fit_->settings().calib, this);
  regionDialog->setWindowTitle("Region settings");
  int ret = regionDialog->exec();

  if (ret == QDialog::Accepted)
  {
    ui->plot->clear_range_selection();

    toggle_push(true);
    fit_->override_region(roi, region);
    emit peak_selection_changed(ui->plot->get_selected_peaks());
    emit data_changed();
    emit fitting_done();
    dirty(roi);
  }
}

void FormFitter::peak_info(double bin)
{
  if (!fit_ || !fit_->contains_peak(bin))
    return;

  DAQuiri::Peak hm = fit_->peaks().at(bin);

  PeakDialog* peakInfo = new PeakDialog(hm, fit_->settings().calib, this);
  auto nrg = fit_->peaks().at(bin).peak_energy(fit_->settings().calib.cali_nrg_);
  peakInfo->setWindowTitle("Parameters for peak at " + QString::number(nrg.value()));
  int ret = peakInfo->exec();

  if ((ret == QDialog::Accepted) && fit_->replace_hypermet(bin, hm))
  {
    updateData();
    std::set<double> selected_peaks;
    selected_peaks.insert(bin);
    ui->plot->set_selected_peaks(selected_peaks);

    emit data_changed();
    emit fitting_done();
  }
}


