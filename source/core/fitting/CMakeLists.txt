set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/fit_settings.cpp
#  ${dir}/gaussian.cpp
#  ${dir}/hypermet.cpp
  ${dir}/finder.cpp
  ${dir}/sum4.cpp
  ${dir}/peak.cpp
#  ${dir}/optimizer.cpp
  ${dir}/roi.cpp
#  ${dir}/fitter.cpp
  )

set(HEADERS
  ${dir}/fit_settings.h
#  ${dir}/gaussian.h
#  ${dir}/hypermet.h
  ${dir}/finder.h
  ${dir}/sum4.h
  ${dir}/peak.h
#  ${dir}/optimizer.h
  ${dir}/roi.h
#  ${dir}/fitter.h
  )

add_subdirectory(BFGS)
add_subdirectory(hypermet)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
