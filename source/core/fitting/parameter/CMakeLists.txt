set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/abstract_param.cpp
  ${dir}/abstract_bounded_param.cpp
  ${dir}/unbounded_param.cpp
  ${dir}/positive_param.cpp
  ${dir}/sine_bounded_param.cpp
  ${dir}/atan_bounded_param.cpp
  )

set(HEADERS
  ${dir}/abstract_param.h
  ${dir}/abstract_bounded_param.h
  ${dir}/unbounded_param.h
  ${dir}/positive_param.h
  ${dir}/sine_bounded_param.h
  ${dir}/atan_bounded_param.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
