#include <fmt/format.h>

template<typename T>
std::string visualize(const T& hist, size_t nstars, bool non_empty_only = false)
{
  if (hist.empty())
    return {};

  auto vmax{hist[0]};
  size_t start{0}, end{0};
  bool print{false};
  for (uint32_t i = 0; i < hist.size(); i++)
  {
    const auto& val = hist[i];
    vmax = std::max(vmax, val);
    if (val > 0)
    {
      end = i;
      if (!print)
      {
        start = i;
        print = true;
      }
    }
  }

  std::string largesti = fmt::format("{}", end);
  std::string pad = "{:<" + fmt::format("{}", largesti.size()) + "}";
  std::string pad2 = "{:<" + std::to_string(nstars + 2) + "}";

  std::stringstream ss;
  for (size_t i = start; i <= end; i++)
  {
    auto val = hist[i];
    if (!non_empty_only || (val > 0))
      ss << fmt::format(pad, i) << ": "
         << fmt::format(pad2, std::string((nstars * val) / vmax, '*'))
         << val << "\n";
  }
  return ss.str();
}

template<typename T>
std::string visualize(const T& x, const T& y, size_t nstars, bool non_empty_only = false)
{
  if (y.empty() || (x.size() != y.size()))
    return {};

  auto vmax{y[0]};
  size_t start{0}, end{0};
  size_t longest_x{0};
  bool print{false};
  for (uint32_t i = 0; i < y.size(); i++)
  {
    longest_x = std::max(longest_x, fmt::format("{}", x[i]).size());
    const auto& val = y[i];
    vmax = std::max(vmax, val);
    if (val > 0)
    {
      end = i;
      if (!print)
      {
        start = i;
        print = true;
      }
    }
  }

  std::string pad = "{:<" + std::to_string(longest_x) + "}";
  std::string pad2 = "{:<" + std::to_string(nstars + 2) + "}";

  std::stringstream ss;
  for (size_t i = start; i <= end; i++)
  {
    auto val = y[i];
    if (!non_empty_only || (val > 0))
      ss << fmt::format(pad, x[i]) << ": "
         << fmt::format(pad2, std::string((nstars * val) / vmax, '*'))
         << val << "\n";
  }
  return ss.str();
}
