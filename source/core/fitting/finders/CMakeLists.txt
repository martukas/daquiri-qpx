set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/finder_kon_naive.cpp
  ${dir}/finder_kon_calibrated.cpp
  ${dir}/finder.cpp
  )

set(HEADERS
  ${dir}/finder_kon_naive.h
  ${dir}/finder_kon_calibrated.h
  ${dir}/finder.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

