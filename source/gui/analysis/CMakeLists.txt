set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/isotope.cpp
  ${dir}/widget_isotopes.cpp
  ${dir}/widget_plot_calib.cpp
  ${dir}/form_energy_calibration.cpp
  ${dir}/form_fwhm_calibration.cpp
  ${dir}/form_efficiency_calibration.cpp
  ${dir}/form_fit_results.cpp
  ${dir}/form_analysis_1d.cpp
  )

set(UI
  ${dir}/widget_isotopes.ui
  ${dir}/form_energy_calibration.ui
  ${dir}/form_fwhm_calibration.ui
  ${dir}/form_efficiency_calibration.ui
  ${dir}/form_fit_results.ui
  ${dir}/form_analysis_1d.ui
  )

set(HEADERS
  ${dir}/isotope.h
  ${dir}/widget_isotopes.h
  ${dir}/widget_plot_calib.h
  ${dir}/form_energy_calibration.h
  ${dir}/form_fwhm_calibration.h
  ${dir}/form_efficiency_calibration.h
  ${dir}/form_fit_results.h
  ${dir}/form_analysis_1d.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
set(${this_target}_ui ${${this_target}_ui} ${UI} PARENT_SCOPE)
