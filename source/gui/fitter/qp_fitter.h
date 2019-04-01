#pragma once

#include <gui/fitter/qp_range_selector.h>
#include <core/gamma/fitter.h>
//#include "QPlot1D.h"

class QpFitter : public QPlot::Multi1D
{
  Q_OBJECT

public:
  explicit QpFitter(QWidget *parent = 0);

  void clear();

  void setFit(DAQuiri::Fitter *fit);
  void update_spectrum();

  void set_busy(bool);

  void clearSelection();
  std::set<double> get_selected_peaks();

  void make_range_selection(double energy);
  void clear_range_selection();

  void follow_selection();

  QCPRange getRange(QCPRange domain = QCPRange()) Q_DECL_OVERRIDE;

public slots:
  void set_selected_peaks(std::set<double> selected_peaks);
  void updateData();

signals:

  void delete_selected_peaks();
  void peak_info(double peak_id);

  void rollback_ROI(double roi_id);
  void roi_settings(double roi_id);
  void refit_ROI(double roi_id);
  void delete_ROI(double roi_id);

  void add_peak(double l, double r);
  void adjust_sum4(double peak_id, double l, double r);
  void adjust_background_L(double roi_id, double l, double r);
  void adjust_background_R(double roi_id, double l, double r);

  void peak_selection_changed(std::set<double> selected_peaks);
  void range_selection_changed(double l, double r);

protected slots:

  void update_range_selection();

  void selection_changed();
  void changeROI(QAction*);

  void adjustY() Q_DECL_OVERRIDE;


protected:
  void executeButton(QPlot::Button *button) Q_DECL_OVERRIDE;
  void mouseClicked(double x, double y, QMouseEvent *event) Q_DECL_OVERRIDE;

  QCPAxisRect* residualsAxisRect {nullptr};

  DAQuiri::Fitter *fit_data_ {nullptr};
  bool busy_ {false};

  QMenu menuROI;

  std::set<double> selected_peaks_;
  double selected_roi_ {-1};
  bool hold_selection_ {false};

  RangeSelector* range_ {nullptr};

  void trim_log_lower(std::vector<double> &array);
  void calc_visible();

  QCPGraph *addGraph(const std::vector<double>& x, const std::vector<double>& y,
                     QPen appearance, bool smooth,
                     QString name = QString(), double region = -1, double peak = -1);

  void plotBackgroundEdge(DAQuiri::SUM4Edge edge,
                          const std::vector<double> &x,
                          double region,
                          QString button_name);

  void plotExtraButtons();
  void plotEnergyLabel(double peak_id, double peak_energy, QCPItemTracer *crs);

  void plotRegion(double region_id, const DAQuiri::RegionManager &region, QCPGraph *data_graph);
  void plotPeak(double region_id, double peak_id, const DAQuiri::Peak &peak);

  void make_SUM4_range(double region, double peak);
  void make_background_range(double region, bool left);

  void what_is_in_range(std::set<double>& peaks,
                        std::set<double>& rois,
                        std::set<double>& best_rois);

  void toggle_items(const std::set<double> &good_peak, const std::set<double> &prime_roi);
  void toggle_label(QCPAbstractItem *item, const std::set<double> &good_peaks);

  void toggle_sum4_line(QCPItemLine* line, bool only_one_region);
  void toggle_background_line(QCPItemLine* line, bool only_one_region,
                              bool adjusting_background, bool adjusting_sum4);

  void toggle_region_button(QPlot::Button* button, const std::set<double> &prime_roi);
  void toggle_peak_button(QPlot::Button* button);
  void toggle_sum4_button(QPlot::Button* button, bool only_one_region);
  void toggle_background_button(QPlot::Button* button, bool only_one_region);


  void toggle_graphs(const std::set<double>& good_peak, const std::set<double>& good_roi);

  bool region_is_selected(double region);
  bool peaks_in_region_selected(double region);
  bool only_one_region_selected();

  void make_range(double energy);

  void add_resids_rect();
  void add_resids_plots(const std::vector<double>& xx);
};
