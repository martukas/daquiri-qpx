set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/calib_function_factory.cpp
  ${dir}/polynomial.cpp
  ${dir}/sqrt_poly.cpp
#  ${dir}/calibration.cpp
  )

set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
