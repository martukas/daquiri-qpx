set(dir ${CMAKE_CURRENT_SOURCE_DIR})

get_filename_component(TEST_DATA_PATH "test_data/prompt0.cnf" DIRECTORY)

set(SOURCES
  ${dir}/ImporterTKA.cpp
  ${dir}/string_to_chans.cpp
  )

#target_compile_definitions(this_target
#PRIVATE TEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DATA_PATH}/")

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
