#pragma once

#pragma GCC diagnostic push
#if defined(__GNUC__) && (__GNUC__ >= 7)
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Core>
#pragma GCC diagnostic pop
