set(this_target ${PROJECT_NAME}_optimizerBFGS)

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
  ${dir}/BFGS.cpp
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
  ${dir}/Fittable.h
  ${dir}/BFGS.h
  ${dir}/more_math.h
  )

set(${this_target}_headers ${HEADERS})
set(${this_target}_sources ${SOURCES})

add_library(
  ${this_target} STATIC
  ${${this_target}_headers}
  ${${this_target}_sources}
)

include_directories(
  ${this_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/source
)

target_link_libraries(
  ${this_target}
  PUBLIC ${PROJECT_NAME}_core
)

enable_coverage(${this_target})

