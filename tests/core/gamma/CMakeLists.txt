set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_evaluation.cpp
  )


add_subdirectory(sum4)
add_subdirectory(hypermet)
add_subdirectory(nuclides)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
