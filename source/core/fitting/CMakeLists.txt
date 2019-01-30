set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_settings.cpp
  ${dir}/finder.cpp
#  ${dir}/optimizer.cpp
  ${dir}/roi.cpp
  ${dir}/fitter.cpp
  ${dir}/uncertain.cpp
  )

set(HEADERS
  ${dir}/fit_settings.h
  ${dir}/finder.h
#  ${dir}/optimizer.h
  ${dir}/roi.h
  ${dir}/fitter.h
  ${dir}/uncertain.h
  )

add_subdirectory(sum4)
add_subdirectory(BFGS)
add_subdirectory(hypermet)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
