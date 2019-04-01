set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/data_model.cpp
  ${dir}/weighted_data.cpp
  ${dir}/weight_strategies.cpp
  )

set(HEADERS
  ${dir}/data_model.h
  ${dir}/weighted_data.h
  ${dir}/weight_strategies.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
