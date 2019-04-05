set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/FormFitter.cpp
  ${dir}/FormFitterSettings.cpp
  ${dir}/FormPeakInfo.cpp
  ${dir}/qp_fitter.cpp
  ${dir}/qp_range_selector.cpp
  ${dir}/rollback_dialog.cpp
  ${dir}/thread_fitter.cpp
  ${dir}/WidgetFitterAutoselect.cpp
  )

set(UI
  ${dir}/FormFitter.ui
  ${dir}/FormFitterSettings.ui
  ${dir}/FormPeakInfo.ui
  ${dir}/WidgetFitterAutoselect.ui
  )

set(HEADERS
  ${dir}/FormFitter.h
  ${dir}/FormFitterSettings.h
  ${dir}/FormPeakInfo.h
  ${dir}/qp_fitter.h
  ${dir}/qp_range_selector.h
  ${dir}/rollback_dialog.h
  ${dir}/thread_fitter.h
  ${dir}/WidgetFitterAutoselect.h
  )

add_subdirectory(widgets)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
set(${this_target}_ui ${${this_target}_ui} ${UI} PARENT_SCOPE)
