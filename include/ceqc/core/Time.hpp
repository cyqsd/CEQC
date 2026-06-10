#pragma once
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

namespace ceqc::model {
using TimePoint = std::chrono::system_clock::time_point;

inline std::string formatUTC(const TimePoint& tp) {
  auto tt = std::chrono::system_clock::to_time_t(tp);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &tt);
#else
  gmtime_r(&tt, &tm);
#endif
  std::ostringstream os;
  os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
  return os.str();
}

inline std::tm toUTC(const TimePoint& tp) {
  auto tt = std::chrono::system_clock::to_time_t(tp);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &tt);
#else
  gmtime_r(&tt, &tm);
#endif
  return tm;
}

inline TimePoint makeUTC(int year, int month, int day, int hour, int minute, double second) {
  std::tm tm{};
  tm.tm_year = year - 1900;
  tm.tm_mon = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min = minute;
  tm.tm_sec = static_cast<int>(second);
#if defined(_WIN32)
  auto t = _mkgmtime(&tm);
#else
  auto t = timegm(&tm);
#endif
  auto tp = std::chrono::system_clock::from_time_t(t);
  auto frac = second - static_cast<int>(second);
  return tp + std::chrono::nanoseconds(static_cast<long long>(frac * 1e9 + 0.5));
}
}
