#pragma once

#include <core/fitting/hypermet/Calibration.h>
#include <vector>

namespace Hypermet
{

// \todo these go somewhere else

double dead_time(double real_time, double live_time);
double rate(double live_time, double sum_counts);

}
