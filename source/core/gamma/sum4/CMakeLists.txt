set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/background_sample.cpp
  ${dir}/peak.cpp
  )

set(HEADERS
  ${dir}/background_sample.h
  ${dir}/peak.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

