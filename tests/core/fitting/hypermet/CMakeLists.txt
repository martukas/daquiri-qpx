set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/value.cpp
  ${dir}/function_test.h
  ${dir}/step.cpp
  ${dir}/tail.cpp
  ${dir}/poly_background.cpp
  ${dir}/peak.cpp
  #${dir}/region.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
