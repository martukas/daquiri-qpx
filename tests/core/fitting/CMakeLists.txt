set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/function_test.cpp
  ${dir}/function_test.h
  ${dir}/simple_functions.h
  )

add_subdirectory(optimizers)
add_subdirectory(data_model)
add_subdirectory(parameter)

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
