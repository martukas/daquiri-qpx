set(this_target ${PROJECT_NAME}_optimizerBFGS)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/EfficiencyCal.cpp
  ${dir}/NonlinearityCal.cpp
  ${dir}/Calibration.cpp
  ${dir}/optimizer_bfgs.cpp
  ${dir}/Spectrum.cpp
  ${dir}/Peak.cpp
  ${dir}/Region.cpp
  )

set(HEADERS
  ${dir}/EfficiencyCal.h
  ${dir}/NonlinearityCal.h
  ${dir}/Calibration.h
  ${dir}/optimizer_bfgs.h
  ${dir}/Spectrum.h
  ${dir}/Peak.h
  ${dir}/Region.h
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

