set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/uncertain_double_widget.cpp
  ${dir}/unbounded_parameter_widget.cpp
  ${dir}/positive_parameter_widget.cpp
  ${dir}/bounded_parameter_widget.cpp
  ${dir}/peak_skew_widget.cpp
  ${dir}/peak_step_widget.cpp
  )

set(HEADERS
  ${dir}/uncertain_double_widget.h
  ${dir}/unbounded_parameter_widget.h
  ${dir}/positive_parameter_widget.h
  ${dir}/bounded_parameter_widget.h
  ${dir}/peak_skew_widget.h
  ${dir}/peak_step_widget.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
