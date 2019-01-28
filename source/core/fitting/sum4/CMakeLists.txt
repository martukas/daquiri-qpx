set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/sum4edge.cpp
  ${dir}/sum4.cpp
  )

set(HEADERS
  ${dir}/sum4edge.h
  ${dir}/sum4.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

