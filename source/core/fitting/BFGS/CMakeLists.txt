set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/BFGS.cpp
  ${dir}/dlib.cpp
  )

set(HEADERS
  ${dir}/BFGS.h
  ${dir}/dlib.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)

