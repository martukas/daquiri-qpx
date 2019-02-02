set(dir ${CMAKE_CURRENT_SOURCE_DIR})

#set(SOURCES
#  ${dir}/parameter.cpp
#  )

add_subdirectory(hypermet)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
