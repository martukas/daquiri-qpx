set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/gaussian_precalc_values.cpp
  ${dir}/skew.cpp
  ${dir}/step.cpp
  ${dir}/peak.cpp
  ${dir}/poly_background.cpp
  ${dir}/gamma_region.cpp
  )

set(HEADERS
  ${dir}/gaussian_precalc_values.h
  ${dir}/skew.h
  ${dir}/step.h
  ${dir}/peak.h
  ${dir}/poly_background.h
  ${dir}/gamma_region.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

