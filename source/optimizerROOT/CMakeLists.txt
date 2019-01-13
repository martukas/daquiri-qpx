include(${ROOT_USE_FILE})

set(this_target ${PROJECT_NAME}_optimizerROOT)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/optimizer_ROOT.cpp
  ${dir}/wrapper_ROOT.cpp
  )

set(HEADERS
  ${dir}/optimizer_ROOT.h
  ${dir}/wrapper_ROOT.h
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
  PUBLIC ${ROOT_INCLUDE_DIR}
)

target_link_libraries(
  ${this_target}
  PUBLIC ${PROJECT_NAME}_core
  PUBLIC ${ROOT_LIBRARIES}
)

enable_coverage(${this_target})

