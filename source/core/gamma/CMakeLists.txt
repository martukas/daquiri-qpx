set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_evaluation.cpp
  ${dir}/fit_settings.cpp
  ${dir}/region_manager.cpp
  ${dir}/fitter.cpp
  )

set(HEADERS
  ${dir}/fit_evaluation.h
  ${dir}/fit_settings.h
  ${dir}/region_manager.h
  ${dir}/fitter.h
  )

add_subdirectory(finders)
add_subdirectory(sum4)
add_subdirectory(hypermet)
add_subdirectory(nuclides)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
