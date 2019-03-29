set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/transition.cpp
  ${dir}/nuclide.cpp
  ${dir}/std_importer.cpp
  )

set(HEADERS
  ${dir}/transition.h
  ${dir}/nuclide.h
  ${dir}/std_importer.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

