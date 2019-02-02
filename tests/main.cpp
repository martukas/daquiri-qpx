#include <gtest/gtest.h>

#pragma push_macro("POSIX")
#undef POSIX
#include <h5cpp/hdf5.hpp>
#pragma pop_macro("POSIX")


#include <core/util/custom_logger.h>

int main(int argc, char **argv)
{
  hdf5::error::Singleton::instance().auto_print(false);
  CustomLogger::initLogger(Severity::Debug, nullptr, "unit_tests.log");

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
