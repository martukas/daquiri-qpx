set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/value.cpp
  ${dir}/step.cpp
  ${dir}/poly_background.cpp
  ${dir}/peak.cpp
  ${dir}/region.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
