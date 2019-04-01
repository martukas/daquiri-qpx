set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/gaussian_precalc_values.cpp
  ${dir}/step.cpp
  ${dir}/skew.cpp
  ${dir}/poly_background.cpp
  ${dir}/peak.cpp
  ${dir}/gamma_region.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
