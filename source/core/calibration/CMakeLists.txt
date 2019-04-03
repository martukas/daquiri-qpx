set(dir ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCES
  ${dir}/calib_function_factory.cpp
  ${dir}/polynomial.cpp
  ${dir}/sqrt_poly.cpp
  ${dir}/calibration.cpp
  )

set(HEADERS
  ${dir}/calib_function.h
  ${dir}/calib_function_factory.h
  ${dir}/polynomial.h
  ${dir}/sqrt_poly.h
  ${dir}/calibration.h
  )

set(${this_target}_headers ${${this_target}_headers} ${HEADERS} PARENT_SCOPE)
set(${this_target}_sources ${${this_target}_sources} ${SOURCES} PARENT_SCOPE)
