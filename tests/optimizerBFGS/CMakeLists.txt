set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/more_math.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
