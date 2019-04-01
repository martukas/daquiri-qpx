set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/clever_hist.h
  ${dir}/clever_hist.cpp
  ${dir}/function_test.h
  ${dir}/function_test.cpp
  ${dir}/weighted_data.cpp
  ${dir}/data_model.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
