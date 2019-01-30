#pragma once

#include <core/fitting/uncertain.h>
#include <core/fitting/sum4/sum4.h>
#include <core/fitting/hypermet/Peak.h>

namespace DAQuiri
{

// \todo these go somewhere else
// \todo use uncert type
double dead_time(double real_time, double live_time);
double rate(double live_time, double sum_counts);

int value_quality(UncertainDouble ud, double error_threshold = 50);

int peak_good(const Peak& h, const SUM4& s);


}
