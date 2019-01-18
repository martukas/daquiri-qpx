set(this_target ${PROJECT_NAME}_optimizerBFGS)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/optimizer_bfgs.cpp
  ${dir}/SpectrumIO.cpp
  )

set(HEADERS
  ${dir}/optimizer_bfgs.h
  ${dir}/SpectrumIO.h
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

