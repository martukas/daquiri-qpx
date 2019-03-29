set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_settings.cpp
  ${dir}/weighted_data.cpp
  ${dir}/fittable_region.cpp
  ${dir}/fit_evaluation.cpp
  ${dir}/region_manager.cpp
  ${dir}/fitter.cpp
  ${dir}/uncertain.cpp
  )

set(HEADERS
  ${dir}/fit_settings.h
  ${dir}/weighted_data.h
  ${dir}/fittable_region.h
  ${dir}/fit_evaluation.h
  ${dir}/region_manager.h
  ${dir}/fitter.h
  ${dir}/uncertain.h
  )

add_subdirectory(sum4)
add_subdirectory(optimizers)
add_subdirectory(hypermet)
add_subdirectory(finders)
add_subdirectory(nuclides)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
