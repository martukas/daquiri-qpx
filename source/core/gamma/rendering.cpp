#include <core/gamma/rendering.h>

namespace DAQuiri
{

void PeakRendering::clear()
{
  peak.clear();
  full_fit.clear();
}

void PeakRendering::render(const Peak& /*h*/)
{

}

void RegionRendering::reserve(size_t count)
{
  channel.assign(count, 0.0);
  energy.assign(count, 0.0);
  background.assign(count, 0.0);
  back_steps.assign(count, 0.0);
  full_fit.assign(count, 0.0);
  sum4_background.assign(count, 0.0);
}

void RegionRendering::clear()
{
  channel.clear();
  energy.clear();
  background.clear();
  back_steps.clear();
  full_fit.clear();
  sum4_background.clear();
  peaks.clear();
}

void RegionRendering::render(const Region& r,
                             const Calibration& energy_calib)
{
  auto sum4back = SUM4Edge::sum4_background(r.LB_, r.RB_);

  clear();

  double start = r.left();
  double end = r.right();

  auto count = static_cast<size_t>(std::ceil((end - start) * subdivisions));
  double step = 1.0 / static_cast<double>(subdivisions);

  reserve(count);

  for (size_t i = 0; i < count; ++i)
  {
    double x = start + static_cast<double>(i) * step;
    channel[i] = x;
    energy[i] = energy_calib.transform(x);
    full_fit[i] = back_steps[i] = background[i] = r.background.eval(x);
    for (auto& p : r.peaks_)
    {
      auto vals = p.second.eval(x);
      back_steps[i] += vals.step_tail();
      full_fit[i] += vals.all();
    }
    sum4_background[i] = sum4back(x);
  }

  for (auto& p : r.peaks_)
  {
    auto& peak = peaks[p.first];
    const auto& hyp = p.second;
    peak.full_fit = back_steps;
    peak.peak.assign(count, 0.0);
    for (size_t i = 0; i < channel.size(); ++i)
    {
      auto vals = hyp.eval(channel[i]);
      peak.peak[i] += vals.peak_skews();
      peak.full_fit[i] += vals.peak_skews();
    }
  }
}


}
