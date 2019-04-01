set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/uncertain.cpp
  )

set(HEADERS
  ${dir}/uncertain.h
  )

add_subdirectory(optimizers)
add_subdirectory(data_model)

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
