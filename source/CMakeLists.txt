include(InstallConfig)
include(BuildInfo)

find_package(date REQUIRED)
find_package(fmt REQUIRED)
find_package(h5cpp 0.1.0 REQUIRED)
find_package(Threads REQUIRED)
add_definitions("-DGRAYLOG")

#find_package(GSL REQUIRED)
#find_package(dlib REQUIRED)
#find_package(MPFR REQUIRED)
#set(optimizer_LIBRARIES dlib::dlib GSL::gsl GSL::gslcblas)
#${MPFR_LIBRARIES})

add_subdirectory(core)
add_subdirectory(consumers)
add_subdirectory(producers)
add_subdirectory(importers)

add_subdirectory(cmd)
add_subdirectory(gui)
