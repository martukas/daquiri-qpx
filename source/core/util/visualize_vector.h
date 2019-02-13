#pragma once

#include <fmt/format.h>
#include <core/util/more_math.h>

template<typename T>
std::string visualize(const T& hist, size_t nstars, bool non_empty_only = false)
{
  if (hist.empty())
    return {};

  auto vmax{hist[0]};
  auto vmin{hist[0]};
  size_t start{0}, end{0};
  bool print{false};
  for (uint32_t i = 0; i < hist.size(); i++)
  {
    const auto& val = hist[i];
    vmax = std::max(vmax, val);
    vmin = std::min(vmin, val);
    if (val != 0.0)
    {
      end = i;
      if (!print)
      {
        start = i;
        print = true;
      }
    }
  }

  if (!print)
    return "all zeros";

  size_t largesti = fmt::format("{}", end).size();

  auto range = vmax - vmin;
  auto starsize = range / nstars;
  size_t neg_stars = (signum(vmin) < 0) ? (std::abs(vmin) / starsize) : 0;
  size_t pos_stars = std::abs(vmax) / starsize;

  std::string pad = "{:<" + std::to_string(largesti) + "}: ";
  std::string pad_neg = pad
      + (neg_stars ? "{:>" + std::to_string(neg_stars) + "}" : "")
      + std::string(pos_stars, ' ');
  std::string pad_pos = pad
      + std::string(neg_stars, ' ')
      + (pos_stars ? "{:<" + std::to_string(pos_stars) + "}" : "");

  std::stringstream ss;
  for (size_t i = start; i <= end; i++)
  {
    auto val = hist[i];
    if (non_empty_only && (val == 0.0))
      continue;
    bool is_neg = (signum(val) < 0);
    const std::string& def = is_neg ? pad_neg : pad_pos;
    size_t count = is_neg ? (neg_stars * val / vmin) : (pos_stars * val / vmax);
    char star = is_neg ? '-' : '+';
    ss << fmt::format(def, i, std::string(count, star)) << "  " << val << "\n";
  }
  return ss.str();
}

template<typename T>
std::string visualize(const T& x, const T& y, size_t nstars, bool non_empty_only = false)
{
  if (y.empty() || (x.size() != y.size()))
    return {};

  auto vmax{y[0]};
  auto vmin{y[0]};
  size_t start{0}, end{0};
  size_t longest_x{0};
  bool print{false};
  for (uint32_t i = 0; i < y.size(); i++)
  {
    longest_x = std::max(longest_x, fmt::format("{}", x[i]).size());
    const auto& val = y[i];
    vmax = std::max(vmax, val);
    vmin = std::min(vmin, val);
    if (val != 0.0)
    {
      end = i;
      if (!print)
      {
        start = i;
        print = true;
      }
    }
  }

  if (!print)
    return "all zeros";

  auto range = vmax - vmin;
  auto starsize = range / nstars;
  size_t neg_stars = (signum(vmin) < 0) ? (std::abs(vmin) / starsize) : 0;
  size_t pos_stars = std::abs(vmax) / starsize;

  std::string pad = "{:<" + std::to_string(longest_x) + "}: ";
  std::string pad_neg = pad
      + (neg_stars ? "{:>" + std::to_string(neg_stars) + "}" : "")
      + (pos_stars ? std::string(pos_stars, ' ') : " ");
  std::string pad_pos = pad
      + (neg_stars ? std::string(neg_stars, ' ') : " ")
      + (pos_stars ? "{:<" + std::to_string(pos_stars) + "}" : "");

  std::stringstream ss;
  for (size_t i = start; i <= end; i++)
  {
    auto val = y[i];
    if (non_empty_only && (val == 0.0))
      continue;
    bool is_neg = (signum(val) < 0);
    const std::string& def = is_neg ? pad_neg : pad_pos;
    size_t count = is_neg ? (neg_stars * val / vmin) : (pos_stars * val / vmax);
    char star = is_neg ? '-' : '+';
    auto stars = count ? std::string(count, star) : " ";
    ss << fmt::format(def, x[i], stars) << "  " << val << "\n";
  }
  return ss.str();
}

template<typename T>
std::string visualize_all(const T& x, const T& y, size_t nstars, bool non_empty_only = false)
{
  if (y.empty() || (x.size() != y.size()))
    return {};

  auto vmax{y[0]};
  auto vmin{y[0]};
  size_t longest_x{0};
  bool print{false};
  for (uint32_t i = 0; i < y.size(); i++)
  {
    longest_x = std::max(longest_x, fmt::format("{}", x[i]).size());
    const auto& val = y[i];
    vmax = std::max(vmax, val);
    vmin = std::min(vmin, val);
    if (val != 0.0)
      print = true;
  }

  if (!print)
    return "all zeros";

  auto range = vmax - vmin;
  auto starsize = range / nstars;
  size_t neg_stars = (signum(vmin) < 0) ? (std::abs(vmin) / starsize) : 0;
  size_t pos_stars = std::abs(vmax) / starsize;

  std::string pad = "{:<" + std::to_string(longest_x) + "}: ";
  std::string pad_neg = pad
      + (neg_stars ? "{:>" + std::to_string(neg_stars) + "}" : "")
      + (pos_stars ? std::string(pos_stars, ' ') : " ");
  std::string pad_pos = pad
      + (neg_stars ? std::string(neg_stars, ' ') : " ")
      + (pos_stars ? "{:<" + std::to_string(pos_stars) + "}" : "");

  std::stringstream ss;
  for (size_t i = 0; i < y.size(); i++)
  {
    auto val = y[i];
    if (non_empty_only && (val == 0.0))
      continue;
    bool is_neg = (signum(val) < 0);
    const std::string& def = is_neg ? pad_neg : pad_pos;
    size_t count = is_neg ? (neg_stars * val / vmin) : (pos_stars * val / vmax);
    char star = is_neg ? '-' : '+';
    auto stars = count ? std::string(count, star) : " ";
    ss << fmt::format(def, x[i], stars) << "  " << val << "\n";
  }
  return ss.str();
}
