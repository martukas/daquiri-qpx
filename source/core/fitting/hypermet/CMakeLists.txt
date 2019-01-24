set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/EfficiencyCal.cpp
  ${dir}/NonlinearityCal.cpp
  ${dir}/Calibration.cpp
  ${dir}/Spectrum.cpp
  ${dir}/Value.cpp
  ${dir}/Tail.cpp
  ${dir}/Step.cpp
  ${dir}/Peak.cpp
  ${dir}/Region.cpp
  )

set(HEADERS
  ${dir}/EfficiencyCal.h
  ${dir}/NonlinearityCal.h
  ${dir}/Calibration.h
  ${dir}/Spectrum.h
  ${dir}/Value.h
  ${dir}/Tail.h
  ${dir}/Step.h
  ${dir}/Peak.h
  ${dir}/Region.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

