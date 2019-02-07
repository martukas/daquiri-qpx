set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/weighted_data.cpp
  )

add_subdirectory(sum4)
add_subdirectory(hypermet)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
