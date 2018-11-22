set(this_target ${PROJECT_NAME}_importers)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/importers_autoreg.cpp
  ${dir}/ImporterAVA.cpp
  ${dir}/ImporterCNF.cpp
  ${dir}/ImporterMCA.cpp
  ${dir}/ImporterN42.cpp
  ${dir}/ImporterTKA.cpp
  ${dir}/string_to_chans.cpp
  ${dir}/xylib/canberra_cnf.cpp
  ${dir}/xylib/canberra_mca.cpp
  ${dir}/xylib/util.cpp
  ${dir}/xylib/xylib.cpp
  )

set(HEADERS
  ${dir}/importers_autoreg.h
  ${dir}/ImporterAVA.h
  ${dir}/ImporterCNF.h
  ${dir}/ImporterMCA.h
  ${dir}/ImporterN42.h
  ${dir}/ImporterTKA.h
  ${dir}/string_to_chans.h
  ${dir}/xylib/canberra_cnf.h
  ${dir}/xylib/canberra_mca.h
  ${dir}/xylib/util.h
  ${dir}/xylib/xylib.h
  )

set(${this_target}_headers ${HEADERS})
set(${this_target}_sources ${SOURCES})

add_library(
  ${this_target} STATIC
  ${${this_target}_headers}
  ${${this_target}_sources}
)

target_include_directories(
  ${this_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/source
  PRIVATE ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(
  ${this_target}
  PUBLIC ${PROJECT_NAME}_core
  PUBLIC ${PROJECT_NAME}_consumers
  PRIVATE pugixml
)

enable_coverage(${this_target})

