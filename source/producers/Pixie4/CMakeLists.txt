set(this_target pixie4)

set(dir ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  message("64 bits compiler detected")
  set(EX_PLATFORM 64)
  set(EX_PLATFORM_NAME "x64")
else(CMAKE_SIZEOF_VOID_P EQUAL 8)
  message("32 bits compiler detected")
  set(EX_PLATFORM 32)
  set(EX_PLATFORM_NAME "x86")
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

if(EX_PLATFORM EQUAL 64)
  link_directories(${dir}/PLX/Library/64bit)
  set(PLX_LIB_DIRECTORY ${dir}/PLX/Library/64bit)
else(EX_PLATFORM EQUAL 64)
  link_directories(${dir}/PLX/Library/32bit)
  set(PLX_LIB_DIRECTORY ${dir}/PLX/Library/32bit)
endif(EX_PLATFORM EQUAL 64)

set(SOURCES
  ${dir}/pixie4.cpp
  ${dir}/pixie4_api_wrapper.cpp
  ${dir}/XIA/Boot.c
  ${dir}/XIA/eeprom.c
  ${dir}/XIA/globals.c
  ${dir}/XIA/pixie_c.c
  ${dir}/XIA/ua_par_io.c
  ${dir}/XIA/utilities.c
  )

set(HEADERS
  ${dir}/pixie4.h
  ${dir}/pixie4_api_wrapper.h
  ${dir}/XIA/boot.h
  ${dir}/XIA/defs.h
  ${dir}/XIA/sharedfiles.h
  ${dir}/XIA/utilities.h
  ${dir}/PLX/Include/PciRegs.h
  ${dir}/PLX/Include/PciTypes.h
  ${dir}/PLX/Include/PexApi.h
  ${dir}/PLX/Include/PlxApi.h
  ${dir}/PLX/Include/PlxDefCk.h
  ${dir}/PLX/Include/Plx.h
  ${dir}/PLX/Include/PlxIoctl.h
  ${dir}/PLX/Include/PlxNetIoctl.h
  ${dir}/PLX/Include/PlxNetTypes.h
  ${dir}/PLX/Include/PlxStat.h
  ${dir}/PLX/Include/Plx_sysdep.h
  ${dir}/PLX/Include/PlxTypes.h
  ${dir}/PLX/Include/Reg9054.h
  )

set_source_files_properties(${dir}/XIA/utilities.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(${dir}/XIA/Boot.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(${dir}/XIA/eeprom.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(${dir}/XIA/globals.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(${dir}/XIA/pixie_c.c PROPERTIES LANGUAGE CXX)
set_source_files_properties(${dir}/XIA/ua_par_io.c PROPERTIES LANGUAGE CXX)

add_library(
  ${this_target} STATIC
  ${SOURCES}
  ${HEADERS}
)

target_link_libraries(
  ${this_target}
  PUBLIC ${PROJECT_NAME}_core
  PRIVATE ${PLX_LIB_DIRECTORY}/PlxApi.a
  PRIVATE ${CMAKE_DL_LIBS}
)

target_include_directories(
  ${this_target}
  PUBLIC ${PROJECT_SOURCE_DIR}/source
  PRIVATE ${dir}/XIA
  PRIVATE ${dir}/PLX/Include
)

if(UNIX)
  message("  !! XIA and PLX using LINUX !! ")
  add_definitions(-DXIA_LINUX)
  add_definitions(-DPLX_LINUX)
  set(CMAKE_CXX_FLAGS "-fpermissive -w")
endif()
