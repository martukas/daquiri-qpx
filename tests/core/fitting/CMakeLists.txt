set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/function_test.h
  ${dir}/weighted_data.cpp
  ${dir}/fittable_region.cpp
  ${dir}/fit_evaluation.cpp
  )

add_subdirectory(sum4)
add_subdirectory(optimizers)
add_subdirectory(hypermet)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
